#!/usr/bin/env python3
"""
Generate MicrOMEGAs input cards from scan parameters.
Reads oks.dat and creates MO_inp*.dat files for each point.
"""

import os
import sys

def ensure_directory(path):
    """Create directory if it doesn't exist."""
    try:
        os.makedirs(path, exist_ok=True)
        return True
    except Exception as e:
        print(f"Failed to create directory {path}: {e}", file=sys.stderr)
        return False

def main():
    input_file = "../run/oks.dat"
    output_dir = "../run/cards"

    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Failed to open {input_file}", file=sys.stderr)
        return 1
    except IOError as e:
        print(f"Failed to open {input_file}: {e}", file=sys.stderr)
        return 1

    if not ensure_directory(output_dir):
        print(f"Failed to create {output_dir}", file=sys.stderr)
        return 1

    for line in lines:
        line = line.strip()
        if not line:
            continue

        parts = line.split()
        if len(parts) < 8:
            continue

        ii = parts[0]
        lx_in = parts[1]
        lhx_in = parts[2]
        lsx_in = parts[3]
        mx_in = parts[4]
        vevs_in = parts[5]
        sint_in = parts[6]
        mh2_in = parts[7]

        print(ii)

        filename = f"MO_inp{ii}.dat"
        filepath = os.path.join(output_dir, filename)

        try:
            with open(filepath, 'w') as mO_inp:
                mO_inp.write(f"LX\t\t{lx_in}\n")
                mO_inp.write(f"LHX\t\t{lhx_in}\n")
                mO_inp.write(f"LSX\t\t{lsx_in}\n")
                mO_inp.write(f"MX\t\t{mx_in}\n")
                mO_inp.write(f"vevs\t{vevs_in}\n")
                mO_inp.write(f"SinT\t{sint_in}\n")
                mO_inp.write(f"Mh2\t\t{mh2_in}\n")
        except IOError as e:
            print(f"Failed to open {filepath}: {e}", file=sys.stderr)
            return 1

    return 0

if __name__ == '__main__':
    sys.exit(main())
