#!/usr/bin/env python3
"""Print or launch a configurable, resumable v2 scan campaign."""
import argparse
import json
import os
from pathlib import Path
import shlex
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from trsm_inputs import PHYSICS_VERSION
from trsm_parallel import job_limit
from generate_mg5_trsm_xsecs import ProcLocation


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT / 'config/next-scan-v2.json')
    parser.add_argument('--campaign-dir', type=Path, required=True)
    parser.add_argument('--seed-start', type=int)
    parser.add_argument('--nseeds', type=int)
    parser.add_argument('--nrandom', type=int)
    parser.add_argument('--nrandom-count-evo-thc', action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument('--jobs', type=job_limit)
    parser.add_argument('--checkpoint-every', type=int)
    parser.add_argument('--heartbeat-seconds', type=float)
    parser.add_argument('--shutdown-grace-seconds', type=float)
    parser.add_argument('--python-executable', type=Path)
    parser.add_argument('--ewpt-executable', type=Path)
    parser.add_argument('--ewpt-minima-executable', type=Path)
    parser.add_argument('--run-mg5', action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument('--mg5-without-dm', action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument('--mg5-process', dest='mg5_processes', action='append', choices=tuple(ProcLocation))
    parser.add_argument('--pilot', action='store_true', help='Default to two scans and two targets each; explicit values override these.')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--resume', action='store_true')
    mode.add_argument('--aggregate-only', action='store_true')
    parser.add_argument('--run', action='store_true', help='Execute the command; otherwise only print it.')
    args = parser.parse_args(argv)
    if args.pilot and (args.resume or args.aggregate_only):
        parser.error('--pilot applies only to fresh campaigns')
    return args


def build_command(args):
    python = str(args.python_executable or sys.executable)
    command = [python, str(ROOT / 'run_trsm_seed_campaign.py'),
               '--campaign-dir', str(args.campaign_dir.expanduser().resolve())]
    if args.resume or args.aggregate_only:
        command.append('--resume' if args.resume else '--aggregate-only')
        # Existing campaigns own their physics configuration; do not reread defaults.
        for name in ('seed_start', 'nseeds', 'nrandom', 'jobs', 'checkpoint_every',
                     'heartbeat_seconds', 'shutdown_grace_seconds', 'python_executable',
                     'ewpt_executable', 'ewpt_minima_executable'):
            value = getattr(args, name)
            if value is not None:
                command += ['--' + name.replace('_', '-'), str(value)]
        if args.nrandom_count_evo_thc is not None:
            command.append('--nrandom-count-evo-thc' if args.nrandom_count_evo_thc else '--no-nrandom-count-evo-thc')
        append_mg5_options(command, args)
        return command
    config = json.loads(args.config.read_text())
    if config.get('constraint_version') != PHYSICS_VERSION:
        raise ValueError('Scan configuration does not use the active physics version')
    values = {
        'seed_start': config['seed_start'],
        'nseeds': 2 if args.pilot else config['nseeds'],
        'nrandom': 2 if args.pilot else config.get('points_per_seed', config.get('draws_per_seed')),
        'jobs': 2 if args.pilot else config.get('jobs', 4),
        'checkpoint_every': config.get('checkpoint_every'),
        'heartbeat_seconds': config.get('heartbeat_seconds'),
        'shutdown_grace_seconds': config.get('shutdown_grace_seconds'),
        'python_executable': python,
        'ewpt_executable': config.get('ewpt_executable'),
        'ewpt_minima_executable': config.get('ewpt_minima_executable'),
    }
    for name, default in values.items():
        value = getattr(args, name) if getattr(args, name) is not None else default
        if value is not None:
            command += ['--' + name.replace('_', '-'), str(value)]
    if values['nrandom'] is None and args.nrandom is None:
        raise ValueError('Set points_per_seed (or legacy draws_per_seed) in the configuration')
    count = args.nrandom_count_evo_thc
    if count is None:
        count = config.get('nrandom_count_evo_thc', '--nrandom-count-evo-thc' in config.get('generator_arguments', []))
    if not isinstance(count, bool):
        raise ValueError('nrandom_count_evo_thc must be a JSON boolean')
    command.append('--nrandom-count-evo-thc' if count else '--no-nrandom-count-evo-thc')
    command += ['--run-ewpt', '--ewpt-thigh', str(config['ewpt_thigh_GeV'])]
    append_mg5_options(command, args, config)
    command += ['--generator-extra-arg=' + str(arg) for arg in config.get('generator_arguments', [])]
    return command


def append_mg5_options(command, args, config=None):
    config = config or {}
    if args.run_mg5 is False:
        config = {key: value for key, value in config.items()
                  if key not in ('mg5_without_dm', 'mg5_processes')}
    for name in ('run_mg5', 'mg5_without_dm'):
        value = getattr(args, name)
        if value is None:
            value = config.get(name)
        if value is not None:
            if not isinstance(value, bool):
                raise ValueError(f'{name} must be a JSON boolean')
            command.append('--' + ('' if value else 'no-') + name.replace('_', '-'))
    processes = args.mg5_processes if args.mg5_processes is not None else config.get('mg5_processes')
    if processes is not None:
        if not isinstance(processes, list) or not processes or any(p not in ProcLocation for p in processes):
            raise ValueError('mg5_processes must be a nonempty list of supported process names')
        for process in processes:
            command += ['--mg5-process', process]


def main(argv=None):
    args = parse_args(argv)
    command = build_command(args)
    # Use the launcher's own validation even for a dry run.
    from run_trsm_seed_campaign import parse_args as validate
    validate(command[2:])
    print(shlex.join(command), flush=True)
    if args.run:
        # Replacing the wrapper preserves signal handling and exit status.
        os.execv(command[0], command)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
