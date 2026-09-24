#!/usr/bin/env python3
"""Plot cumulative relic, direct, gamma-line and optional CMB selections."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cutflow import stages


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('outdir', type=Path)
    parser.add_argument('--xvar', default='MX')
    parser.add_argument('--yvar', default='dm_omega')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--show', action='store_true')
    args = parser.parse_args()
    import matplotlib
    if not args.show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = json.loads((args.outdir / 'results.json').read_text())
    for key in (args.xvar, args.yvar):
        if not any(key in row for row in rows):
            parser.error(f'Column not found: {key}')
    steps = list(stages(rows))
    fig, axes = plt.subplots(1, len(steps), figsize=(3.3 * len(steps), 3.8), sharex=True, sharey=True, squeeze=False)
    for ax, step in zip(axes[0], steps):
        indices = set(step['indices'])
        valid = [r for r in rows if isinstance(r.get(args.xvar), (int, float)) and
                 isinstance(r.get(args.yvar), (int, float))]
        selected = [r for r in valid if r['index'] in indices]
        ax.scatter([r[args.xvar] for r in valid], [r[args.yvar] for r in valid],
                   s=18, color='lightgray', label='All inputs')
        ax.scatter([r[args.xvar] for r in selected], [r[args.yvar] for r in selected],
                   s=22, color='tab:blue', label='Remaining')
        ax.set_title(f"{step['stage']}\n{step['passed']} remain; {step['unassessed']} unassessed", fontsize=9)
        ax.set_xlabel(args.xvar)
        ax.grid(alpha=.2)
    axes[0][0].set_ylabel(args.yvar)
    axes[0][0].legend(fontsize=8)
    fig.tight_layout()
    target = args.output or args.outdir / 'constraint_cutflow.png'
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=160)
    print(target.resolve())
    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
