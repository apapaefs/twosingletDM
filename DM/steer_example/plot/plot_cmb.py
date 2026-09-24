#!/usr/bin/env python3
"""Compare native full-abundance and assessed Planck CMB ratios."""
import argparse
import json
from pathlib import Path
import sys
import matplotlib
if '--show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('outdir', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--show', action='store_true')
    args = parser.parse_args()
    rows = json.loads((args.outdir / 'results.json').read_text())
    enabled = [r for r in rows if r.get('dm_cmb_enabled')]
    valid = [r for r in enabled if r.get('dm_cmb_available')]
    missing = len(enabled) - len(valid)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
    for ax, key, title in zip(axes, ('dm_cmb_ratio_raw', 'dm_cmb_ratio'),
                              ('Native ratio (full DM abundance)', 'Ratio used for the CMB verdict')):
        ax.axhspan(1, 1e12, alpha=.06, color='crimson')
        ax.axhline(1, color='crimson', ls='--', lw=1, label='Exclusion: ratio > 1')
        ax.scatter([r['MX'] for r in valid], [r[key] for r in valid], s=32, color='tab:blue')
        if len(valid) <= 20:
            for r in valid:
                ax.annotate(str(r['index']), (r['MX'], r[key]), xytext=(5, 5), textcoords='offset points')
        ax.set_xscale('log')
        ax.set_yscale('symlog', linthresh=1e-8)
        ax.set_xlabel(r'$M_X$ [GeV]')
        ax.set_title(title)
        ax.grid(alpha=.2)
    values = [r[k] for r in valid for k in ('dm_cmb_ratio_raw', 'dm_cmb_ratio')]
    positive = [v for v in values if v > 0]
    axes[0].set_ylim(0 if any(v == 0 for v in values) else min(positive + [1]) / 3,
                     max(values + [1]) * 5)
    axes[0].set_ylabel(r'$R_{\rm CMB}$ (95% limit normalized to one)')
    axes[1].legend(fontsize=8)
    fig.suptitle(f'Planck CMB: {len(valid)} available, {missing} unavailable, {len(rows)-len(enabled)} disabled')
    fig.tight_layout()
    target = args.output or args.outdir / 'cmb_ratios.png'
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, dpi=160)
    print(target.resolve())
    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
