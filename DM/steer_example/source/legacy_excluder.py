"""Compatibility for the original positional relic/DD/gamma-line assessment.

This interface has no solver-status or CMB inputs. Its historical 0.121
abundance reference and LZ2025 fit are deliberately confined to this module.
"""
import math
import os
from pathlib import Path
import sys

from run_single_point import read_points
from test_trsm_DM import fermi_lat_r16_line_limit, lz2025_source_direct_detection_base_limit


def main(argv):
    if len(argv) < 4 or (len(argv) - 4) % 2:
        raise ValueError('Usage: mO_excluder.py index m_DM Omega DirDet [E_gamma Phi_R16]...')
    index = int(argv[0])
    mass, omega, direct = map(float, argv[1:4])
    if not math.isfinite(mass) or mass <= 0 or any(not math.isfinite(v) or v < 0 for v in (omega, direct)):
        raise ValueError('Require positive finite mass and nonnegative finite Omega and SI cross section')
    base = Path(__file__).resolve().parents[1]
    output = Path(os.environ.get('MO_EXCLUDER_OUTPUT_DIR', base / 'output'))
    points = read_points(Path(os.environ.get('MO_EXCLUDER_OKS_FILE', base / 'run/oks.dat')))
    point = next((point for point in points if point.index == index), None)
    if point is None:
        raise ValueError(f'Index {index} not found in oks.dat')
    fraction = min(1., omega / .121)
    limit_base = lz2025_source_direct_detection_base_limit(mass)
    limit = limit_base / fraction if fraction else math.inf
    available, max_ratio = False, 0.
    best = (math.nan, math.nan, math.nan)
    for energy, flux in zip(map(float, argv[4::2]), map(float, argv[5::2])):
        bound = fermi_lat_r16_line_limit(energy) / fraction**2 if fraction else math.inf
        if not math.isfinite(bound) or bound <= 0 or not math.isfinite(flux) or flux < 0:
            continue
        available = True
        if flux / bound > max_ratio:
            max_ratio, best = flux / bound, (energy, flux, bound)
    relic_pass, direct_pass, indirect_pass = omega <= .121, direct <= limit, max_ratio <= 1
    values = [*point.columns().values(), mass, omega, direct, limit, limit_base,
              int(available), *best, max_ratio]
    line = '\t'.join(map(str, values)) + '\n'
    selected = {
        'relic_pass': relic_pass, 'relic_strict': abs(omega - .12) <= .001,
        'omexcl': not relic_pass, 'luxexcl': not direct_pass, 'luxpass': direct_pass,
        'all_dirpass': relic_pass and direct_pass, 'omgpass_dirfail': relic_pass and not direct_pass,
        'indirexcl': not indirect_pass, 'indirpass': indirect_pass,
        'indir_caughtit': relic_pass and direct_pass and not indirect_pass,
        'dir_caughtit': relic_pass and not direct_pass and indirect_pass,
        'dir_indir_caughtit': relic_pass and not direct_pass and not indirect_pass,
        'allall': relic_pass and direct_pass and indirect_pass,
        'dmexcl': not (relic_pass and direct_pass and indirect_pass),
    }
    output.mkdir(parents=True, exist_ok=True)
    # Do not mix historical assessments into a v2 output directory.
    if (output / 'results.json').exists() or (output / 'metadata.json').exists():
        raise ValueError('Use a separate directory for legacy positional assessments and v2 runs')
    (output / 'DM_data').mkdir(exist_ok=True)
    for filename in ('DM_data', f'DM_data_{index}'):
        (output / 'DM_data' / filename).write_text(line)
    for name, include in selected.items():
        with (output / (name + '.dat')).open('a') as stream:
            if include:
                stream.write(line)
    print('Legacy positional assessment: relic/DD/gamma lines; solver status and CMB were not supplied.')
    print(f'Point {index}: relic={relic_pass}, direct={direct_pass}, gamma={indirect_pass}')
    return 0
