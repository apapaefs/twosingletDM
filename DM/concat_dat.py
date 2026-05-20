import os
import glob
import shutil
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
    line_length = 0
    
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
                    if len(seen_lines) != 0:
                        if len(line_key.split()) != line_length:
                            print(f"Warning: Line in {input_file} has different number of columns than previous lines.")
                    
                    if line_key not in seen_lines:
                        parts = stripped.split()
                        if len(parts) > 1:
                            outfile.write(f"{line_number} {' '.join(parts[1:])}\n")
                        else:
                            outfile.write(f"{line_number}\n")
                        seen_lines.add(line_key)
                        line_number += 1
                        line_length = len(line_key.split())
                    else:
                        duplicate_count += 1
    
    if duplicate_count > 0:
        print(f"Deleted {duplicate_count} lines due to repetition")
    print(f"Successfully concatenated {len(input_files)} files into {output_file}")

# def concat_dat_files_from_pattern(output_file, pattern):
#     """
#     Concatenate .dat files matching a glob pattern.
    
#     Args:
#         output_file (str): Path to the output concatenated .dat file
#         pattern (str): Glob pattern to match input files (e.g., '*.dat')
#     """
#     input_files = sorted(glob.glob(pattern))
    
#     if len(input_files) < 2:
#         raise ValueError(f"Found fewer than 2 files matching pattern {pattern}")
    
#     concat_dat_files(output_file, *input_files)
    
#     print(f"Successfully concatenated {len(input_files)} files into {output_file}")

def concat_out_folders(output_folder, *input_folders):
    """
    Concatenate two or more output folders with .dat files into a single output folder.

    Args:
        output_folder (str): Path to the output folder
        *input_folders: Variable number of input folder paths
    """
    if os.path.exists(output_folder):
        print(f"Output folder {output_folder} already exists.")
        return
    
    os.makedirs(output_folder, exist_ok=True)
    list_of_file_names = []
    for input_folder in input_folders:
        for file in glob.glob(os.path.join(input_folder, "*.dat")):
            list_of_file_names.append(os.path.basename(file))

    if len(list_of_file_names) == 0:
        print("No .dat files found in the specified input folders.")
        return

    for file in list_of_file_names:
        input_files = [os.path.join(input_folder, file) for input_folder in input_folders]
        if len(input_files) < 2:
            print(f"Warning: File {file} found in only one input folder, copying...")
            # copy file to output folder
            shutil.copy(input_files[0], output_folder)
        else:
            concat_dat_files(os.path.join(output_folder, file), *input_files)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python concat_dat.py <output_file> <input_file1> <input_file2> [input_file3] ... or [pattern*.dat]")
        print("or for folders: python concat_out.py <output_folder> <input_folder1> <input_folder2> [input_folder3] ...")
        sys.exit(1)
    
    if sys.argv[1].endswith(".dat"):
        output_file = sys.argv[1]
        input_files = sys.argv[2:]
        concat_dat_files(output_file, *input_files)
    else:
        output_folder = sys.argv[1]
        input_folders = sys.argv[2:]
        concat_out_folders(output_folder, *input_folders)