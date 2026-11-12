#!/usr/bin/env python3
"""Run or resume CalcGW after a seed campaign has finished aggregating."""
import argparse
import json
import os
from pathlib import Path
import shlex
import socket
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from trsm_parallel import cancellation_signals, effective_jobs, job_limit, stop_processes
from trsm_scan_campaign import (CampaignLock, CampaignStateError, _pid_is_alive,
                                atomic_write_json, load_json, utc_now)


def completed_scan(directory, *, wait=False, poll_seconds=30):
    """Never read an aggregate while its scan supervisor can still replace it."""
    announced = False
    while True:
        state = load_json(directory / 'campaign_state.json', 'campaign state')
        if state.get('schema') != 'trsm_seed_campaign_v1':
            raise ValueError('Unsupported campaign state schema')
        if state.get('hostname') != socket.gethostname():
            raise ValueError('CalcGW follow-up must run on the campaign host')
        lock_path = directory / 'campaign_state.lock'
        active = False
        if lock_path.exists():
            try:
                lock = json.loads(lock_path.read_text())
            except FileNotFoundError:
                lock = None  # Supervisor completed between the two reads.
            if lock is not None:
                if not isinstance(lock, dict) or lock.get('hostname') != socket.gethostname() or not isinstance(lock.get('pid'), int):
                    raise ValueError('Cannot establish the campaign lock owner')
                active = _pid_is_alive(lock['pid'])
        if not active:
            # The supervisor may have published completion after our first read
            # and then released its lock. Recheck before declaring interruption.
            state = load_json(directory / 'campaign_state.json', 'campaign state')
        complete = (state.get('status') == 'complete' and bool(state.get('seeds'))
                    and all(s.get('status') == 'complete' for s in state['seeds']))
        if complete and not active:
            if not (directory / 'combined_points.tsv').is_file():
                raise ValueError('Completed campaign has no combined_points.tsv')
            return state
        if not wait or not active:
            raise ValueError('Campaign is not complete; finish or resume the scan before CalcGW')
        if not announced:
            print(f'Waiting for all seeds and aggregation in {directory}', flush=True)
            announced = True
        time.sleep(poll_seconds)


def calcgw_command(directory, state, jobs=None):
    config = state['configuration']
    requested = jobs if jobs is not None else state.get('execution', {}).get('jobs', 1)
    workers = effective_jobs(requested, 10**9)
    output = directory / 'combined_points_gw.tsv'
    command = [str(config.get('python_executable') or sys.executable),
               str(ROOT / 'reprocess_trsm_gw.py'), str(directory / 'combined_points.tsv'),
               '--output', str(output), '--thigh', str(config['ewpt_thigh']),
               '--jobs', str(workers)]
    # Match a campaign's explicitly selected BSMPT installation when present.
    if config.get('ewpt_executable'):
        command += ['--executable', str(Path(config['ewpt_executable']).with_name('CalcGW'))]
    if Path(str(output) + '.partial').exists() or output.exists() or output.with_suffix('.metadata.json').exists():
        command.append('--resume')
    return command


def run_command(command):
    """Forward cancellation to this follow-up's own process group only."""
    with cancellation_signals() as cancellation:
        process = subprocess.Popen(command, start_new_session=True)
        try:
            while process.poll() is None:
                if cancellation['signal'] is not None:
                    stop_processes([process], 30)
                    return 128 + cancellation['signal']
                try:
                    process.wait(timeout=1)
                except subprocess.TimeoutExpired:
                    pass
            return process.returncode
        except BaseException:
            stop_processes([process], 30)
            raise


def run_followup(campaign_dir, *, jobs=None, wait=False):
    directory = Path(campaign_dir).expanduser().resolve()
    if not (directory / 'campaign_state.json').is_file():
        raise ValueError('Campaign state does not exist')
    progress_path = directory / 'calcgw_followup.json'
    with CampaignLock(progress_path, sys.argv):
        progress = dict(schema='trsm_calcgw_followup_v1', pid=os.getpid(),
                        hostname=socket.gethostname(), started_utc=utc_now())

        def record(status, **fields):
            progress.update(status=status, updated_utc=utc_now(), **fields)
            atomic_write_json(progress_path, progress)

        try:
            record('waiting_for_scan')
            state = completed_scan(directory, wait=wait)
            command = calcgw_command(directory, state, jobs)
            print(shlex.join(command), flush=True)
            record('running_calcgw', command=command)
            code = run_command(command)
            if code:
                record('failed', returncode=code)
                return code
            output = directory / 'combined_points_gw.tsv'
            metadata = load_json(output.with_suffix('.metadata.json'), 'CalcGW metadata')
            counts = metadata.get('gw_enrichment', {}).get('counts', {})
            if not output.is_file() or counts.get('complete') is not True:
                raise ValueError('CalcGW exited without a complete output')
            record('complete', output=str(output), counts=counts, returncode=0)
            print(f'CalcGW follow-up complete: {output}', flush=True)
            return 0
        except (OSError, ValueError, CampaignStateError) as error:
            record('failed', error=str(error))
            raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign-dir', type=Path, required=True)
    parser.add_argument('--jobs', type=job_limit, help='Default: saved scan concurrency; capped by available CPUs')
    parser.add_argument('--wait', action='store_true', help='Wait for an active scan and its aggregation to finish')
    args = parser.parse_args(argv)
    try:
        return run_followup(args.campaign_dir, jobs=args.jobs, wait=args.wait)
    except (OSError, ValueError, CampaignStateError) as error:
        print(f'CalcGW follow-up error: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
