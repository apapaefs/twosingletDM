#!/usr/bin/env python3
"""Write full-precision v2 cards, including the common SM inputs."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from run_single_point import read_points, write_card


def main():
    base = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=base / 'run/oks.dat')
    parser.add_argument('--output-dir', type=Path, default=base / 'run/cards')
    args = parser.parse_args()
    points = read_points(args.input)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        parser.error('Card directory is not empty; use a new directory to avoid stale cards')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for point in points:
        write_card(point, args.output_dir / f'MO_inp{point.index}.dat')
    print(f'Wrote {len(points)} v2 cards to {args.output_dir}')


if __name__ == '__main__':
    main()
