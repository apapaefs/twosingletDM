"""Isolated Odysseus regression pilot; never launch a production campaign."""
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time

ROOT = Path('/tmp/trsm-mg5-survey-arlac7oh')
RUNTIME = Path('/home/apapaefs/Projects/TwoSingletDM/runtime-v2-new')
REPO = Path('/home/apapaefs/Projects/TwoSingletDM/twosingletDM')
PROCESS = {'gg_heta0': 'gg_heta0', 'pp_eta0Z': 'pp_eta0Z',
           'hh': 'gg_hh_twoscalar', 'hhh': 'gg_hhh_twoscalar'}


def ignore(path, names):
    return [n for n in names if n in ('Events', 'HTML', 'ME5_debug', '.trsm_madevent.lock', '__pycache__')
            or n.endswith('_debug.log') or re.fullmatch(r'G\d+(?:_\d+)?', n)]


def run_case(case, key, *, legacy=False, reference=False):
    from tools.repair_mg5_survey import repair_survey
    directory = ROOT / case / PROCESS[key]
    assert not directory.exists(), directory
    shutil.copytree(RUNTIME / 'MG5_aMC_v3_5_15' / PROCESS[key], directory,
                    symlinks=True, ignore=ignore)
    for name in ('mpmodule.mod', 'libiregi.a', 'libcts.a'):
        path = directory / 'lib' / name
        if path.is_symlink():
            data = path.read_bytes()
            path.unlink()
            path.write_bytes(data)
    for name in ('Events', 'HTML'):
        (directory / name).mkdir(exist_ok=True)
    banner = (ROOT / 'original-banner.txt').read_text()
    slha = re.search(r'<slha>\s*(.*?)\s*</slha>', banner, re.S | re.I)[1] + '\n'
    (directory / 'Cards/param_card.dat').write_text(slha)
    card = directory / 'Cards/run_card.dat'
    text = (re.search(r'<mgruncard>\s*(.*?)\s*</mgruncard>', banner, re.S | re.I)[1] + '\n'
            if key == 'gg_heta0' else card.read_text())
    text, n = re.subn(r'(?m)^\s*\S+\s*=\s*iseed\b', ' 22 = iseed' if reference else ' 21 = iseed', text)
    assert n == 1
    if key != 'pp_eta0Z':
        text, n = re.subn(r'(?m)^\s*\S+\s*=\s*nhel\b', ' 1 = nhel', text)
        assert n == 1
    card.write_text(text)
    with (directory / 'Cards/me5_configuration.txt').open('a') as stream:
        stream.write('\nnotification_center = False\nautomatic_html_opening = False\n')
    receipt = repair_survey(directory) if legacy else None
    code = ROOT / ('code-legacy' if legacy else 'code-final')
    env = dict(os.environ, PYTHONPATH=str(code), TRSM_MG5_LOCATION=str(directory.parent),
               TRSM_MG5_CORES='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               PATH=str(RUNTIME / 'venv/bin') + ':' + os.environ['PATH'])
    for key_ in ('MAKEFLAGS', 'GNUMAKEFLAGS', 'MFLAGS', 'MAKELEVEL', 'MAKEOVERRIDES'):
        env.pop(key_, None)
    log = ROOT / (case + '-' + key + '.log')
    result_path = ROOT / (case + '-' + key + '.json')
    started = time.monotonic()
    with log.open('w') as stream:
        child = subprocess.Popen([str(RUNTIME / 'venv/bin/python'), str(Path(__file__)),
                                  '--worker', key, str(result_path), str(int(reference))],
                                 cwd=directory, env=env, stdout=stream, stderr=subprocess.STDOUT,
                                 start_new_session=True)
        try:
            code_ = child.wait(timeout=240)
        except subprocess.TimeoutExpired:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
            raise
    result = json.loads(result_path.read_text()) if result_path.exists() else {}
    result.update(case=case, process=key, exit_code=code_, elapsed_s=time.monotonic() - started,
                  log=str(log), maintenance_receipt=str(receipt) if receipt else None)
    result_path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)
    assert code_ == 0 and result['status'] == 'passed', result
    return result


def worker():
    import generate_mg5_trsm_xsecs as mg
    p = json.loads((ROOT / 'point.json').read_text())
    process, destination, reference = sys.argv[2:]
    if int(reference):
        mg.LOOP_INDUCED_SURVEY_ITERATIONS = 6
        original_run = mg.subprocess.run
        def reference_run(args, **kwargs):
            command = Path(args[1])
            command.write_text(command.read_text().replace('--points=300 ', '--points=3000 '))
            return original_run(args, **kwargs)
        mg.subprocess.run = reference_run
    value = mg.get_mg5_xsec(process, 'surveytest', p['lambdas'], p['k1'], p['k2'], p['k3'],
                            p['m2'], p['w2'], p['m3'], p['w3'], ecm=p['ecm'], w1=p['w1'], k233=p['k233'])
    directory = Path(mg.MGLocation) / mg.ProcLocation[process]
    files = list((directory / 'Events').glob('*/unweighted_events.lhe.gz'))
    assert len(files) == 1
    with gzip.open(files[0], 'rt') as stream:
        lhe = stream.read()
    init = re.search(r'<init>\s*(.*?)\s*</init>', lhe, re.S)[1].splitlines()
    process_line = init[1].split()
    error = float(process_line[1])
    assert lhe.count('<event>') == 1 and math.isfinite(value) and value > 0
    assert math.isfinite(error) and error >= 0
    weight = float(re.search(r'<event>\s*([^\n]+)', lhe)[1].split()[2])
    assert math.isfinite(weight) and weight > 0
    inputs = [path.read_text().splitlines() for path in directory.glob('SubProcesses/P*/input_app.txt')]
    assert inputs
    assert all(int(lines[0].split()[1]) >= int(lines[0].split()[2]) for lines in inputs)
    if process != 'pp_eta0Z':
        assert all(int(lines[4].split()[0]) == 0 for lines in inputs)
    result = dict(status='passed', cross_section_pb=value, error_pb=error, event_count=1,
                  weight=weight, lhe=str(files[0]),
                  wrapper_sha256=hashlib.sha256(Path(mg.__file__).read_bytes()).hexdigest(),
                  survey_inputs=inputs)
    Path(destination).write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    if '--worker' in sys.argv:
        worker()
    else:
        (ROOT / 'code-legacy').mkdir()
        shutil.copy2(REPO / 'generate_mg5_trsm_xsecs.py', ROOT / 'code-legacy/generate_mg5_trsm_xsecs.py')
        results = []
        for args, kwargs in [
            (('current', 'gg_heta0'), {}),
            (('legacy-repaired', 'gg_heta0'), {'legacy': True}),
            (('reference', 'gg_heta0'), {'reference': True}),
            (('current', 'pp_eta0Z'), {}),
            (('current', 'hh'), {}),
            (('current', 'hhh'), {}),
        ]:
            results.append(run_case(*args, **kwargs))
            (ROOT / 'validation.json').write_text(json.dumps({'cases': results}, indent=2) + '\n')
        current, legacy, reference = results[:3]
        assert current['cross_section_pb'] == legacy['cross_section_pb']
        sigma = math.hypot(current['error_pb'], reference['error_pb'])
        assert abs(current['cross_section_pb'] - reference['cross_section_pb']) <= 3 * sigma
        print('All native comparisons passed.', flush=True)
