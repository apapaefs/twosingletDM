#!/usr/bin/env python3
"""
Exclude dark matter points based on experimental constraints.
Evaluates relic density, direct detection, and indirect detection limits.
"""

import os
import sys
import math
import subprocess
from pathlib import Path

# ============================================================================
# Direct Detection Limits (LZ/LUX 2025 update)
# ============================================================================

def dir_det_excl(x):
    """
    Calculate direct detection exclusion limit for dark matter mass x (in GeV).
    Based on LZ and LUX constraints (2025 update).
    """
    result = 0.0

    if x < 10:
        result = (math.exp(-23.2949 + 2.2493 * pow((-3.60228 + math.log(x)), 2)))
        result = result * (1 - 0.465188 + 1.64349 * pow((-2.3605 + math.log(x)), 2))
    elif x < 40:
        result = 2.45e-12 * pow(10, (1.65 / pow(1 - math.log10(40), 2)) * pow(math.log10(x) - math.log10(40), 2))
    elif x < 60:
        result = 2.32775e-12 + 2.28365e-15 * x
    elif x < 80:
        result = 1.68775e-12 + 1.28365e-14 * x
    elif x < 100:
        result = math.exp(0.671326 * math.log(x) - 29.5696)
    else:
        result = 3.49167e-14 * math.exp(0.977674 * math.log(x))

    return result

# ============================================================================
# Fermi-LAT Limits (R16 ROI)
# ============================================================================

class FermiLatLimits:
    """Fermi-LAT gamma-ray line limits with log-linear interpolation."""

    # Data from arXiv:1506.00013
    _LIMITS_1 = [
        (0.214, 45.5), (0.234, 38.0), (0.255, 34.0), (0.278, 32.4), (0.303, 30.9),
        (0.329, 28.9), (0.358, 26.3), (0.388, 21.5), (0.421, 18.8), (0.456, 16.7),
        (0.493, 14.7), (0.533, 13.5), (0.576, 11.8), (0.620, 10.4), (0.668, 9.84),
        (0.718, 8.78), (0.770, 7.80), (0.826, 6.74), (0.885, 6.36), (0.947, 7.32),
        (1.01, 8.68), (1.08, 9.29), (1.16, 8.68), (1.24, 7.55), (1.32, 6.39),
        (1.41, 5.12), (1.50, 4.03), (1.60, 3.17), (1.70, 2.92), (1.81, 2.20),
        (1.93, 1.84), (2.06, 1.93), (2.19, 1.87), (2.33, 1.75), (2.49, 1.35),
        (2.65, 1.07), (2.82, 0.844), (3.00, 0.738), (3.20, 0.780), (3.40, 0.915),
        (3.62, 1.02), (3.85, 0.925), (4.09, 0.764), (4.35, 0.715), (4.63, 0.561),
        (4.91, 0.445), (5.22, 0.263),
    ]

    _LIMITS_2 = [
        (5.54, 21.6), (5.87, 20.6), (6.23, 17.0), (6.60, 22.5), (6.99, 26.0),
        (7.40, 33.4), (7.83, 25.3), (8.28, 18.6), (8.76, 22.1), (9.26, 12.7),
        (9.79, 7.59), (10.4, 8.14), (10.9, 12.4), (11.6, 16.0), (12.2, 7.72),
        (12.9, 7.83), (13.6, 8.42), (14.4, 5.91), (15.2, 6.52), (16.1, 6.66),
        (17.0, 8.18), (17.9, 10.2), (18.9, 5.79), (20.0, 3.65), (21.1, 6.56),
        (22.3, 3.66), (23.6, 3.74), (24.9, 2.97), (26.4, 3.78), (27.9, 4.56),
        (29.5, 7.05), (31.2, 4.37), (33.0, 3.28), (34.9, 4.17), (36.9, 4.71),
        (39.0, 3.18), (41.3, 3.07), (43.8, 4.71), (46.4, 5.66), (49.1, 6.40),
        (52.1, 4.56), (55.2, 3.96), (58.6, 4.85), (62.2, 3.32), (66.0, 1.82),
        (70.1, 1.90), (74.5, 3.63), (79.2, 1.48), (84.2, 0.951), (89.6, 0.947),
        (95.4, 0.891), (102, 2.29), (108, 4.89), (115, 4.92), (123, 3.84),
        (131, 3.11), (140, 1.48), (150, 0.765), (160, 0.764), (171, 1.11),
        (183, 1.70), (196, 2.22), (210, 2.85), (225, 1.59), (241, 1.93),
        (259, 0.867), (276, 0.843), (294, 1.32), (321, 1.45), (345, 1.17),
        (367, 0.646), (396, 0.613), (427, 0.560), (462, 0.487),
    ]

    @classmethod
    def get_limits(cls):
        """Get all Fermi-LAT limits with unit conversion (1e-8 and 1e-10)."""
        limits = []
        unit_1 = 1.0e-8
        for energy, flux in cls._LIMITS_1:
            limits.append((energy, flux * unit_1))

        unit_2 = 1.0e-10
        for energy, flux in cls._LIMITS_2:
            limits.append((energy, flux * unit_2))

        return limits

    @classmethod
    def get_limit(cls, energy_gev):
        """Interpolate Fermi-LAT limit at given energy using log-linear interpolation."""
        limits = cls.get_limits()

        if not limits or energy_gev < limits[0][0] or energy_gev > limits[-1][0]:
            return float('inf')

        for i in range(len(limits) - 1):
            low_energy, low_flux = limits[i]
            high_energy, high_flux = limits[i + 1]

            if energy_gev == low_energy:
                return low_flux
            if energy_gev <= high_energy:
                # Log-linear interpolation
                log_energy = (math.log(energy_gev) - math.log(low_energy)) / \
                             (math.log(high_energy) - math.log(low_energy))
                return math.exp(math.log(low_flux) +
                               log_energy * (math.log(high_flux) - math.log(low_flux)))

        return limits[-1][1]

# ============================================================================
# Exclusion Assessment
# ============================================================================

def relic_density_fraction(omega, relic_upper_limit):
    if not math.isfinite(omega) or omega < 0.0:
        raise ValueError("Relic density Omega must be finite and non-negative")
    if not math.isfinite(relic_upper_limit) or relic_upper_limit <= 0.0:
        raise ValueError("Relic-density upper limit must be finite and positive")
    return min(1.0, omega / relic_upper_limit)


def assess_indirect_limit(indirect_channels, omega, relic_upper_limit, rescale = True):
    """
    Assess indirect detection limit for a list of gamma-ray channels.
    Each channel is a tuple (energy_gev, flux_cm2s).
    Returns a dict with assessment results.
    """
    result = {
        'available': False,
        'excluded': False,
        'channels_seen': 0,
        'channels_used': 0,
        'max_ratio': 0.0,
        'energy_gev': float('nan'),
        'flux_cm2s': float('nan'),
        'limit_cm2s': float('nan'),
    }
    abundance_fraction = relic_density_fraction(omega, relic_upper_limit)

    for energy, flux in indirect_channels:
        result['channels_seen'] += 1

        limit = FermiLatLimits.get_limit(energy)
        if rescale:
            if abundance_fraction > 0.0:
                limit /= abundance_fraction**2
            else:
                limit = float('inf')

        if not math.isfinite(limit) or limit <= 0.0 or not math.isfinite(flux) or flux < 0.0:
            continue

        result['channels_used'] += 1
        result['available'] = True
        ratio = flux / limit


        if ratio > result['max_ratio']:
            result['max_ratio'] = ratio
            result['energy_gev'] = energy
            result['flux_cm2s'] = flux
            result['limit_cm2s'] = limit

    result['excluded'] = result['available'] and result['max_ratio'] > 1.0
    return result

# ============================================================================
# File I/O Utilities
# ============================================================================

def ensure_directory(path):
    """Create directory if it doesn't exist."""
    try:
        os.makedirs(path, exist_ok=True)
        return True
    except Exception as e:
        print(f"Could not create output directory {path}", file=sys.stderr)
        return False

def output_path(directory, filename):
    """Build output file path."""
    if not directory or directory == ".":
        return filename
    return os.path.join(directory, filename)

def append_line(path, line):
    """Append a line to a file."""
    try:
        with open(path, 'a') as f:
            f.write(line + '\n')
    except IOError as e:
        print(f"Failed to write to {path}: {e}", file=sys.stderr)

def shell_quote(text):
    """Quote a string for shell use."""
    escaped = text.replace("'", "'\\''")
    return f"'{escaped}'"

# ============================================================================
# Plotting Functions
# ============================================================================

def write_dirdet_limit_plot(output_dir):
    """Generate direct detection limit plot with gnuplot."""
    if not ensure_directory(output_dir):
        return

    data_path = output_path(output_dir, "a_dirdet_limit_curve.dat")
    script_path = output_path(output_dir, "a_dirdet_limit_curve.gnuplot")
    plot_path = output_path(output_dir, "a_dirdet_limit_curve.png")

    try:
        with open(data_path, 'w') as data:
            data.write("# mDM_GeV\tlimit_pb\tlimit_cm2\n")
            min_mass = 1.0
            max_mass = 1.0e4
            n_points = 800

            for i in range(n_points):
                t = i / (n_points - 1)
                mass = math.exp(math.log(min_mass) + t * (math.log(max_mass) - math.log(min_mass)))
                limit_pb = dir_det_excl(mass)
                data.write(f"{mass}\t{limit_pb}\t{limit_pb * 1.0e-36}\n")
    except IOError as e:
        print(f"Failed to write {data_path}: {e}", file=sys.stderr)
        return

    try:
        with open(script_path, 'w') as script:
            script.write("set terminal pngcairo size 1000,750 enhanced font 'Arial,12'\n")
            script.write(f"set output {shell_quote(plot_path)}\n")
            script.write("set logscale xy\n")
            script.write("set grid\n")
            script.write("set xlabel 'dark matter mass m_{{DM}} [GeV]'\n")
            script.write("set ylabel 'direct detection limit [cm^2]'\n")
            script.write("set title 'DirDetexcl(m_{{DM}}) limit used by mO_excluder.py'\n")
            script.write(f"plot {shell_quote(data_path)} using 1:3 with lines linewidth 2 title 'DirDetexcl'\n")
    except IOError as e:
        print(f"Failed to write {script_path}: {e}", file=sys.stderr)
        return

    try:
        subprocess.run(['gnuplot', script_path], check=True, capture_output=True)
        print(f"Wrote {data_path} and {plot_path}")
    except FileNotFoundError:
        print(f"Wrote {data_path} and {script_path}; install gnuplot or run the script manually to create {plot_path}")
    except subprocess.CalledProcessError as e:
        print(f"Wrote {data_path} and {script_path}; install gnuplot or run the script manually to create {plot_path}")

def write_indirect_limit_plot(output_dir):
    """Generate Fermi-LAT indirect detection limit plot."""
    if not ensure_directory(output_dir):
        return

    data_path = output_path(output_dir, "a_fermi_lat_line_limit_r16.dat")
    script_path = output_path(output_dir, "a_fermi_lat_line_limit_r16.gnuplot")
    plot_path = output_path(output_dir, "a_fermi_lat_line_limit_r16.png")

    limits = FermiLatLimits.get_limits()

    try:
        with open(data_path, 'w') as data:
            data.write("# Fermi-LAT 95% CL observed gamma-gamma line flux upper limits, R16 ROI\n")
            data.write("# Source: arXiv:1506.00013, Flux column Phi_gammagamma.\n")
            data.write("# E_gamma_GeV\tPhi_limit_cm-2_s-1\n")
            for energy, flux in limits:
                data.write(f"{energy}\t{flux}\n")
    except IOError as e:
        print(f"Failed to write {data_path}: {e}", file=sys.stderr)
        return

    try:
        with open(script_path, 'w') as script:
            script.write("set terminal pngcairo size 1000,750 enhanced font 'Arial,12'\n")
            script.write(f"set output {shell_quote(plot_path)}\n")
            script.write("set logscale xy\n")
            script.write("set grid\n")
            script.write("set xlabel 'line photon energy E_{{gamma}} [GeV]'\n")
            script.write("set ylabel '95% CL line flux upper limit [cm^{{-2}} s^{{-1}}]'\n")
            script.write("set title 'Fermi-LAT gamma-ray line limit, R16 ROI (arXiv:1506.00013)'\n")
            script.write(f"plot {shell_quote(data_path)} using 1:2 with linespoints linewidth 2 pointtype 7 pointsize 0.45")
            script.write(" title 'Observed Phi_{{gamma gamma}} limit'\n")
    except IOError as e:
        print(f"Failed to write {script_path}: {e}", file=sys.stderr)
        return

    try:
        subprocess.run(['gnuplot', script_path], check=True, capture_output=True)
        print(f"Wrote {data_path} and {plot_path}")
    except FileNotFoundError:
        print(f"Wrote {data_path} and {script_path}; install gnuplot or run the script manually to create {plot_path}")
    except subprocess.CalledProcessError as e:
        print(f"Wrote {data_path} and {script_path}; install gnuplot or run the script manually to create {plot_path}")

# ============================================================================
# Main
# ============================================================================

def main():
    if len(sys.argv) == 2 and sys.argv[1] == "--plot-dirdet-limits":
        output_dir = os.environ.get('MO_EXCLUDER_OUTPUT_DIR', '../output')
        write_dirdet_limit_plot(output_dir)
        return 0

    if len(sys.argv) == 2 and sys.argv[1] == "--plot-indirect-limits":
        output_dir = os.environ.get('MO_EXCLUDER_OUTPUT_DIR', '../output')
        write_indirect_limit_plot(output_dir)
        return 0

    if len(sys.argv) < 5 or (len(sys.argv) - 5) % 2 != 0:
        print("Sth wrong with inp!!", file=sys.stderr)
        print(f"i have {len(sys.argv)} arguments", file=sys.stderr)
        print(f"usage: {sys.argv[0]} index m_DM Omega DirDet [E_gamma Phi_R16]...", file=sys.stderr)
        print(f"       {sys.argv[0]} --plot-dirdet-limits", file=sys.stderr)
        print(f"       {sys.argv[0]} --plot-indirect-limits", file=sys.stderr)
        return 1

    index = int(sys.argv[1])
    m_dm = float(sys.argv[2])
    omega = float(sys.argv[3])
    dirdet = float(sys.argv[4])

    if not math.isfinite(m_dm) or m_dm <= 0.0:
        print("Dark matter mass must be finite and positive", file=sys.stderr)
        return 1
    if not math.isfinite(omega) or omega < 0.0:
        print("Relic density Omega must be finite and non-negative", file=sys.stderr)
        return 1
    if not math.isfinite(dirdet) or dirdet < 0.0:
        print(
            "Direct-detection cross section must be finite and non-negative",
            file=sys.stderr,
        )
        return 1

    indirect_channels = []
    for i in range(5, len(sys.argv), 2):
        energy = float(sys.argv[i])
        flux = float(sys.argv[i + 1])
        indirect_channels.append((energy, flux))

    print(f"processing {index}\t{m_dm}\t{omega}\t{dirdet}")

    output_dir = os.environ.get('MO_EXCLUDER_OUTPUT_DIR', '../output')
    oks_file = os.environ.get('MO_EXCLUDER_OKS_FILE', '../run/oks.dat')

    # Read oks.dat
    lx = lhx = lsx = mx = vevs = sint = mh2 = 0.0
    found_point = False

    try:
        with open(oks_file, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) < 8:
                    continue
                ii = int(parts[0])
                if ii == index:
                    lx = float(parts[1])
                    lhx = float(parts[2])
                    lsx = float(parts[3])
                    mx = float(parts[4])
                    vevs = float(parts[5])
                    sint = float(parts[6])
                    mh2 = float(parts[7])
                    found_point = True
                    break
    except (FileNotFoundError, IOError) as e:
        print(f"Could not open {oks_file}", file=sys.stderr)
        return 1

    if not found_point:
        print(f"Index {index} not found in oks.dat", file=sys.stderr)
        return 1

    if not ensure_directory(output_dir):
        return 1

    # Initialize flags
    dmisok = True
    omgok = True
    dirok = True
    indok = True

    relic_upper_limit = 0.121
    relic_strict_central = 0.12
    relic_strict_width = 0.001
    abundance_fraction = relic_density_fraction(omega, relic_upper_limit)
    lux = dir_det_excl(m_dm)
    dir_det_limit = lux

    rescale = True
    if rescale:
        if abundance_fraction > 0.0:
            dir_det_limit = lux / abundance_fraction
        else:
            dir_det_limit = float('inf')

    base_line = f"{index}\t{lx}\t{lhx}\t{lsx}\t{mx}\t{vevs}\t{sint}\t{mh2}\t{m_dm}"

    indirect = assess_indirect_limit(indirect_channels, omega, relic_upper_limit)

    dm_line = f"{base_line}\t{omega}\t{dirdet}\t{dir_det_limit}\t{lux}"

    indirect_line = (f"{dm_line}\t{1 if indirect['available'] else 0}"
                    f"\t{indirect['energy_gev']}\t{indirect['flux_cm2s']}"
                    f"\t{indirect['limit_cm2s']}\t{indirect['max_ratio']}")

    # Create empty output files
    for fname in ["relic_pass.dat", "relic_strict.dat", "omexcl.dat", "luxexcl.dat",
                  "luxpass.dat", "all_dirpass.dat", "omgpass_dirfail.dat", "indirexcl.dat",
                  "indirpass.dat", "indir_caughtit.dat", "dir_caughtit.dat",
                  "dir_indir_caughtit.dat", "dmexcl.dat", "allall.dat"]:
        try:
            open(output_path(output_dir, fname), 'a').close()
        except IOError:
            pass

    # Write DM_data
    dm_data_dir = output_path(output_dir, "DM_data")
    ensure_directory(dm_data_dir)
    try:
        with open(output_path(dm_data_dir, "DM_data"), 'w') as f:
            f.write(indirect_line + '\n')
    except IOError as e:
        print(f"Failed to write DM_data: {e}", file=sys.stderr)

    # Apply constraints
    if omega <= relic_upper_limit:
        append_line(output_path(output_dir, "relic_pass.dat"), indirect_line)

    if abs(omega - relic_strict_central) <= relic_strict_width:
        append_line(output_path(output_dir, "relic_strict.dat"), indirect_line)

    if omega > relic_upper_limit:
        dmisok = False
        omgok = False
        append_line(output_path(output_dir, "omexcl.dat"), indirect_line)

    if dirdet > dir_det_limit:
        dmisok = False
        dirok = False
        append_line(output_path(output_dir, "luxexcl.dat"), indirect_line)
    else:
        append_line(output_path(output_dir, "luxpass.dat"), indirect_line)

    if omgok and dirok:
        append_line(output_path(output_dir, "all_dirpass.dat"), indirect_line)

    if omgok and not dirok:
        append_line(output_path(output_dir, "omgpass_dirfail.dat"), indirect_line)

    if indirect['excluded']:
        dmisok = False
        indok = False
        append_line(output_path(output_dir, "indirexcl.dat"), indirect_line)
    else:
        append_line(output_path(output_dir, "indirpass.dat"), indirect_line)

    if omgok and dirok and not indok:
        append_line(output_path(output_dir, "indir_caughtit.dat"), indirect_line)

    if omgok and not dirok and indok:
        append_line(output_path(output_dir, "dir_caughtit.dat"), indirect_line)

    if omgok and not dirok and not indok:
        append_line(output_path(output_dir, "dir_indir_caughtit.dat"), indirect_line)

    if not dmisok:
        print(f"excluded from dm {omgok}{dirok}{indok}")
        append_line(output_path(output_dir, "dmexcl.dat"), indirect_line)

        if not omgok:
            print("Too high relic dens")

        if not dirok:
            print("Should be visible in LUX!")

        if not indok:
            print("Excluded by Fermi-LAT gamma-line indirect detection")
    else:
        print("Point agrees with DM data")
        append_line(output_path(output_dir, "allall.dat"), indirect_line)

    return 0

if __name__ == '__main__':
    sys.exit(main())
