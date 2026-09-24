#!/usr/bin/env python3
"""Plot the same unrescaled SI or gamma-line limits used by the v2 assessment."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from test_trsm_DM import FERMI_LAT_R16_LINE_LIMITS
from trsm_direct_detection import DEFAULT_LIMIT_TABLE, load_si_limit_table


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kind', choices=('direct', 'indirect'))
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--limit-table', type=Path, default=DEFAULT_LIMIT_TABLE)
    args = parser.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.kind == 'direct':
        table = load_si_limit_table(args.limit_table)
        x, y = table.masses_gev, [v * 1e-36 for v in table.limits_pb]
        name = 'a_dirdet_limit_curve'
        xlabel, ylabel, title = 'DM mass [GeV]', 'SI upper limit per nucleon [cm²]', table.label + ' (90% CL)'
    else:
        x, y = zip(*FERMI_LAT_R16_LINE_LIMITS)
        name = 'a_indirect_limit_curve'
        xlabel, ylabel, title = 'Photon energy [GeV]', 'Line flux upper limit [cm⁻² s⁻¹]', 'Fermi-LAT R16 gamma-line limit'
    (args.output_dir / (name + '.dat')).write_text(
        f'# {xlabel}\t{ylabel}\n' + ''.join(f'{a:.17g}\t{b:.17g}\n' for a, b in zip(x, y)))
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.loglog(x, y, '.-', lw=1)
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.grid(alpha=.2)
    fig.tight_layout()
    target = args.output_dir / (name + '.png')
    fig.savefig(target, dpi=160)
    plt.close(fig)
    print(target.resolve())
    return 0


if __name__ == '__main__':
    sys.exit(main())
