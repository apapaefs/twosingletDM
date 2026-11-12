"""Apply the validated runtime-only repair, preserving both saved campaigns."""
from contextlib import ExitStack
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path('/tmp/trsm-mg5-survey-arlac7oh')
REPO = Path('/home/apapaefs/Projects/TwoSingletDM/twosingletDM')
RUNTIME = Path('/home/apapaefs/Projects/TwoSingletDM/runtime-v2-new')
sys.path.insert(0, str(REPO))
import run_trsm_seed_campaign as campaign
sys.path.append(str(ROOT / 'code-final'))
from tools.repair_mg5_survey import repair_survey, pending_updates

PROCESSES = [RUNTIME / 'MG5_aMC_v3_5_15' / name for name in
             ('gg_heta0', 'gg_hh_twoscalar', 'gg_hhh_twoscalar')]
CAMPAIGNS = [REPO / 'output' / name for name in
             ('odysseus-smoke-10x100-mg5', 'odysseus-100x1000-66661')]


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def protected_files():
    paths = {REPO / name for name in subprocess.check_output(
        ['git', 'ls-files', '-z'], cwd=REPO).decode().split('\0') if name}
    for process in PROCESSES + [RUNTIME / 'MG5_aMC_v3_5_15/pp_eta0Z']:
        paths.update(path for path in (process / 'Cards').glob('*') if path.is_file()
                     and path.name not in ('run_card.dat', 'run_card_default.dat'))
        paths.update(path for path in (process / 'Events').rglob('*') if path.is_file())
        paths.update((process / 'SubProcesses').rglob('*.f'))
        paths.update((process / 'SubProcesses').rglob('*.inc'))
        paths.update((process / 'bin/internal/ufomodel').rglob('*.py'))
    for directory in CAMPAIGNS:
        paths.update(path for path in directory.rglob('*') if path.is_file()
                     and path.suffix in ('.json', '.dat', '.tsv', '.csv', '.log'))
    return sorted(path for path in paths if path.is_file())


def preflight(state, label):
    directory = ROOT / label
    directory.mkdir()
    args = campaign.parse_args(['--campaign-dir', str(directory), '--resume'])
    campaign.hydrate_configuration(args, state)
    return campaign.run_preflight(args)


if __name__ == '__main__':
    validation = json.loads((ROOT / 'validation.json').read_text())
    assert len(validation['cases']) == 6 and all(case['status'] == 'passed' for case in validation['cases'])
    assert Path(campaign.__file__).resolve().parent == REPO
    assert Path(os.environ['TRSM_MG5_LOCATION']).resolve() == RUNTIME / 'MG5_aMC_v3_5_15'
    with ExitStack() as locks:
        states = []
        for directory in CAMPAIGNS:
            locks.enter_context(campaign.CampaignLock(directory / 'campaign_state.json', [__file__]))
            state = json.loads((directory / 'campaign_state.json').read_text())
            campaign.assert_no_live_workers(state)
            states.append(state)
        assert states[0]['runtime_receipt'] == states[1]['runtime_receipt']
        initial_receipt = preflight(states[0], 'preflight-before')
        assert initial_receipt == states[0]['runtime_receipt'], 'Existing runtime differs before repair'
        protected = {str(path): digest(path) for path in protected_files()}
        (ROOT / 'protected-before.json').write_text(json.dumps(protected, indent=2) + '\n')
        repairs = [str(repair_survey(process)) for process in PROCESSES]
        assert all(not pending_updates(process) for process in PROCESSES)
        final_receipt = preflight(states[0], 'preflight-after')
        changed = [name for name, value in protected.items() if digest(Path(name)) != value]
        assert not changed, changed
        assert final_receipt == initial_receipt
        fingerprint = campaign.configuration_fingerprint(final_receipt)
        assert all(state['runtime_fingerprint'] == fingerprint for state in states)
        result = dict(protected_files=len(protected), protected_files_unchanged=True,
                      runtime_fingerprint=fingerprint, campaign_receipts_unchanged=True,
                      repair_receipts=repairs, restarted_campaigns=False)
        (ROOT / 'repair-report.json').write_text(json.dumps(result, indent=2) + '\n')
        print(json.dumps(result, indent=2))
