import os
import glob
import sys


def concat_dat_files(output_file, *input_files):
    """
    Concatenate two or more .dat files into a single output file.
    Removes duplicate lines from the output and renumbers indices sequentially.
    
    Args:
        output_file (str): Path to the output concatenated .dat file
        *input_files: Variable number of input .dat file paths
    """
    if len(input_files) < 2:
        raise ValueError("At least two input files are required for concatenation")
    
    seen_lines = set()
    duplicate_count = 0
    line_number = 1
    
    def get_line_key(line):
        """Extract all but the first number from a line for comparison."""
        parts = line.split()
        if len(parts) > 1:
            return ' '.join(parts[1:])
        return line
    
    with open(output_file, 'w') as outfile:
        for input_file in input_files:
            if not os.path.exists(input_file):
                print(f"Warning: File {input_file} not found, skipping...")
                continue
            
            with open(input_file, 'r') as infile:
                for line in infile:
                    stripped = line.rstrip('\n')
                    if not stripped.strip():
                        continue
                    line_key = get_line_key(stripped)
                    if line_key not in seen_lines:
                        parts = stripped.split()
                        if len(parts) > 1:
                            outfile.write(f"{line_number} {' '.join(parts[1:])}\n")
                        else:
                            outfile.write(f"{line_number}\n")
                        seen_lines.add(line_key)
                        line_number += 1
                    else:
                        duplicate_count += 1
    
    if duplicate_count > 0:
        print(f"Deleted {duplicate_count} lines due to repetition")
    print(f"Successfully concatenated {len(input_files)} files into {output_file}")


def concat_dat_files_from_pattern(output_file, pattern):
    """
    Concatenate .dat files matching a glob pattern.
    
    Args:
        output_file (str): Path to the output concatenated .dat file
        pattern (str): Glob pattern to match input files (e.g., '*.dat')
    """
    input_files = sorted(glob.glob(pattern))
    
    if len(input_files) < 2:
        raise ValueError(f"Found fewer than 2 files matching pattern {pattern}")
    
    with open(output_file, 'w') as outfile:
        for input_file in input_files:
            with open(input_file, 'r') as infile:
                outfile.write(infile.read())
                outfile.write('\n')
    
    print(f"Successfully concatenated {len(input_files)} files into {output_file}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python concat_dat.py <output_file> <input_file1> <input_file2> [input_file3] ...")
        sys.exit(1)
    
    output_file = sys.argv[1]
    input_files = sys.argv[2:]
    
    concat_dat_files(output_file, *input_files)
