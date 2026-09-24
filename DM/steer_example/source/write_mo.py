#!/usr/bin/env python3
"""Write full-precision v2 cards, including the common SM inputs."""
import argparse
import os
from pathlib import Path
import shutil
import sys
import tempfile
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_single_point import read_points, write_card


def main():
    base = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=base / 'run/oks.dat')
    parser.add_argument('--output-dir', type=Path, default=base / 'run/cards')
    parser.add_argument('--overwrite', action='store_true',
                        help='Accepted for compatibility; replaced cards are always backed up')
    args = parser.parse_args()
    points = read_points(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    # Regeneration was part of the original interface. Stage the complete new
    # set first and retain changed/stale cards instead of requiring a new flag.
    with tempfile.TemporaryDirectory(prefix='.cards-', dir=args.output_dir) as temporary:
        stage = Path(temporary)
        for point in points:
            write_card(point, stage / f'MO_inp{point.index}.dat')
        changed = [path for path in args.output_dir.glob('MO_inp*.dat')
                   if not (stage / path.name).exists() or path.read_bytes() != (stage / path.name).read_bytes()]
        if changed:
            backup = Path(tempfile.mkdtemp(prefix='.previous-cards-', dir=args.output_dir))
            for path in changed:
                shutil.move(str(path), backup / path.name)
            print(f'Previous cards retained in {backup}')
        for path in stage.glob('MO_inp*.dat'):
            if not (args.output_dir / path.name).exists():
                os.replace(path, args.output_dir / path.name)
    print(f'Wrote {len(points)} v2 cards to {args.output_dir}')


if __name__ == '__main__':
    main()
