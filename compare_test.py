#!/usr/bin/env python3

import argparse
import pandas as pd
import filecmp
import os
import subprocess

debug = True

def runcmd(*args, shell=True, capture_output=False, check=True, debug=debug):
    """run a shell command

    This function runs a shell command (converts arguments to a space-separated string first).
    It will print the command before running it if debug is set to true.

    Args:
        args: unnamed arguments to use for the shell command
        shell (boolean, optional): Run as a shell command. Defaults to True
        capture_output (boolean, optional): Capture the output in result. Defaults to False
        check (boolean, optional): Check for non-zero exit code. Defaults to True
        debug (boolean, optional): Defaults to the global value at runtime

    Returns:
        completed_process: instance returned from subprocess.run
    """
    cmd = " ".join([str(item) for item in args])
    if debug:
        print("--- Running: " + cmd)
    return subprocess.run(cmd, shell=shell, capture_output=capture_output, check=check)

def check_files_identical(file1_path, file2_path):
    """
    Compares two files for identity using an efficient standard library function.

    The shallow comparison checks:
    1. Same size (os.stat().st_size)
    2. Same modification time (os.stat().st_mtime)

    If the shallow comparison passes, it performs a full content comparison.
    """

    # Optional: Basic check to ensure files exist before comparison
    if not os.path.exists(file1_path) or not os.path.exists(file2_path):
        print("Error: One or both files not found.")
        return False

    # filecomp.cmp(f1, f2, shallow=True) is the default
    # If shallow is True, files with the same os.stat() signatures are considered
    # identical without reading them. This is fast and what you want.
    are_identical = filecmp.cmp(file1_path, file2_path, shallow=False)

    return are_identical
##################################################

argparser = argparse.ArgumentParser(
        description="Compare amrfinder output from two different runs of AMRFinderPlus",
        formatter_class=argparse.RawTextHelpFormatter
        )
argparser.add_argument(
        'directory',
        type=str,
        help="The directory created by random_pathogen_asm.py"
        )
args = argparser.parse_args()

try:
    df = pd.read_csv(args.directory + '/asm_data.tab', sep='\t')
except FileNotFoundError:
    print(f"Error: File not found at {file_path}")
    # Handle the error gracefully or re-raise
    raise

for acc in df.asm_acc:
    old = f"{args.directory}/amrfinder_old/{acc}.amrfinder_old"
    new = f"{args.directory}/amrfinder_new/{acc}.amrfinder_new"
    # check if different
    if check_files_identical(old, new):
        print("--- files for acc are identical ---")
    else:
        runcmd(f"vimdiff {old} {new}")
    
