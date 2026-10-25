#!/usr/bin/env python3

import argparse
import pandas as pd
import filecmp
import os
import shutil
import subprocess
import logging

debug = True

# Configure logging
logging.basicConfig(
    level=logging.DEBUG if debug else logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

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
        logger.debug(f"Running: {cmd}")
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
        logger.error(f"One or both files not found: {file1_path}, {file2_path}")
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
argparser.add_argument(
        '--status-file',
        type=str,
        default=None,
        help="File to record comparison classifications (default: <directory>/comparison_status.txt)"
        )
args = argparser.parse_args()

try:
    df = pd.read_csv(args.directory + '/asm_data.tab', sep='\t')
except FileNotFoundError:
    logger.error(f"File not found at {args.directory}/asm_data.tab")
    # Handle the error gracefully or re-raise
    raise

status_file = args.status_file if args.status_file else os.path.join(args.directory, "comparison_status.txt")

# Read existing classifications to skip already classified comparisons
classified = {}
if os.path.exists(status_file):
    with open(status_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            acc = parts[0].strip()
            status = parts[1].strip() if len(parts) > 1 else "classified"
            classified[acc] = status
    logger.info(f"Loaded {len(classified)} previously classified comparison(s) from {status_file}")

use_cdiff = shutil.which("cdiff") is not None

try:
    for acc in df.asm_acc:
        acc = str(acc)
        if acc in classified:
            logger.info(f"Skipping {acc}: already classified as {classified[acc]}")
            continue

        old = f"{args.directory}/amrfinder_old/{acc}.amrfinder_old"
        new = f"{args.directory}/amrfinder_new/{acc}.amrfinder_new"
        # check if different
        if check_files_identical(old, new):
            logger.info(f"Files for {acc} are identical")
        else:
            if use_cdiff:
                runcmd(f"cdiff {old} {new} | less -SR", check=False)
            else:
                runcmd(f"vimdiff {old} {new}", check=False)

            response = input("Is this comparison OK? [s]crutinize further or any other key to continue: ")
            if response.strip().lower().startswith("s"):
                status = "scrutinize"
            else:
                status = "OK"

            with open(status_file, "a") as f:
                f.write(f"{acc}\t{status}\n")
            classified[acc] = status
            logger.info(f"Recorded {acc} as {status}")
except (KeyboardInterrupt, EOFError):
    logger.info("Exiting on user interrupt.")
    
