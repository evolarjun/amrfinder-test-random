#!/usr/bin/env python3
"""
Prepare test cases for posting to developers based on comparison_status.txt.

Generates a directory containing:
  - The assembly files (*_genomic.fna)
  - The two output versions (.amrfinder_old and .amrfinder_new)
  - A text file (commands.txt) containing the two command-lines to reproduce
"""

import argparse
import logging
import os
import re
import shutil
from glob import glob
from pathlib import Path
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def taxgroup2orgopt(taxgroup):
    """Convert taxgroup name to the organism option expected by AMRFinderPlus."""
    if not taxgroup or pd.isna(taxgroup):
        return None
    match taxgroup:
        case 'Salmonella enterica':
            return 'Salmonella'
        case 'Acinetobacter baumannii':
            return 'Acinetobacter'
        case 'E.coli and Shigella':
            return 'Escherichia_coli_Shigella'
        case 'Klebsiella pneumoniae':
            return 'Klebsiella'
        case 'Serratia marcescens':
            return 'Serratia'
        case 'Campylobacter jejuni':
            return 'Campylobacter'
        case _:
            return str(taxgroup).replace(' ', '_')


def detect_amrfinder_paths(directory):
    """Extract amrfinder_old and amrfinder_new from pipeline_run.py if available."""
    old_val, new_val = "4.2.7/amrfinder", "4.2.17/amrfinder"
    pipeline_file = os.path.join(directory, "pipeline_run.py")
    if os.path.exists(pipeline_file):
        with open(pipeline_file, "r") as f:
            for line in f:
                m_old = re.search(r'^\s*amrfinder_old\s*=\s*["\']([^"\']+)["\']', line)
                if m_old:
                    old_val = m_old.group(1)
                m_new = re.search(r'^\s*amrfinder_new\s*=\s*["\']([^"\']+)["\']', line)
                if m_new:
                    new_val = m_new.group(1)
    return old_val, new_val


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare test cases for posting to developers from comparison_status.txt",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        'directory',
        type=str,
        help="Directory created by random_pathogen_asm.py (e.g. test1)"
    )
    parser.add_argument(
        '--status-file',
        type=str,
        default=None,
        help="Path to comparison_status.txt (default: <directory>/comparison_status.txt)"
    )
    parser.add_argument(
        '--output-dir', '-o',
        type=str,
        default=None,
        help="Directory to place prepared test cases (default: <directory>/test_cases)"
    )
    parser.add_argument(
        '--status',
        type=str,
        default='scrutinize',
        help="Status to filter for from status-file (default: scrutinize). Use 'all' for all lines."
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help="Include all accessions in comparison_status.txt regardless of status"
    )
    parser.add_argument(
        '--acc',
        type=str,
        default=None,
        help="Limit to a specific assembly accession"
    )
    parser.add_argument(
        '--amrfinder-old',
        type=str,
        default=None,
        help="Command/path for old AMRFinderPlus (default: from pipeline_run.py or 4.2.7/amrfinder)"
    )
    parser.add_argument(
        '--amrfinder-new',
        type=str,
        default=None,
        help="Command/path for new AMRFinderPlus (default: from pipeline_run.py or 4.2.17/amrfinder)"
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=8,
        help="Number of threads for commands (default: 8)"
    )
    parser.add_argument(
        '--split',
        action='store_true',
        help="Create a subfolder per accession inside output directory"
    )
    parser.add_argument(
        '--symlink',
        action='store_true',
        help="Create symbolic links instead of copying files"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    status_file = args.status_file if args.status_file else os.path.join(args.directory, "comparison_status.txt")
    if not os.path.exists(status_file):
        logger.error(f"Status file not found at: {status_file}")
        raise FileNotFoundError(f"Status file not found at: {status_file}")

    output_dir = args.output_dir if args.output_dir else os.path.join(args.directory, "test_cases")
    os.makedirs(output_dir, exist_ok=True)

    # Detect AMRFinder paths
    detected_old, detected_new = detect_amrfinder_paths(args.directory)
    amrfinder_old = args.amrfinder_old if args.amrfinder_old else detected_old
    amrfinder_new = args.amrfinder_new if args.amrfinder_new else detected_new

    # Load organism / taxgroup metadata
    taxgroup_map = {}
    asm_data_file = os.path.join(args.directory, "asm_data.tab")
    if os.path.exists(asm_data_file):
        try:
            df = pd.read_csv(asm_data_file, sep="\t")
            if "asm_acc" in df.columns and "taxgroup_name" in df.columns:
                for _, row in df.iterrows():
                    taxgroup_map[str(row["asm_acc"])] = row["taxgroup_name"]
        except Exception as e:
            logger.warning(f"Could not read {asm_data_file}: {e}")

    # Read target accessions from status file
    target_accessions = []
    with open(status_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            acc = parts[0].strip()
            status = parts[1].strip() if len(parts) > 1 else ""

            if args.acc and acc != args.acc:
                continue

            if args.all or args.status.lower() == "all":
                target_accessions.append(acc)
            elif not status or status.lower() == args.status.lower():
                target_accessions.append(acc)

    if not target_accessions:
        logger.warning(f"No accessions matched criteria (status='{args.status}') in {status_file}")
        return

    logger.info(f"Processing {len(target_accessions)} accession(s) into {output_dir}")

    all_commands = []

    for acc in target_accessions:
        # 1. Locate assembly file
        assembly_matches = glob(f"{args.directory}/ncbi_dataset/data/{acc}/{acc}*.fna*")
        if not assembly_matches:
            assembly_matches = glob(f"{args.directory}/**/{acc}*.fna*", recursive=True)

        if not assembly_matches:
            logger.error(f"Assembly file not found for {acc}")
            continue
        assembly_path = assembly_matches[0]

        # 2. Locate output versions
        old_output = os.path.join(args.directory, "amrfinder_old", f"{acc}.amrfinder_old")
        new_output = os.path.join(args.directory, "amrfinder_new", f"{acc}.amrfinder_new")

        if not os.path.exists(old_output):
            logger.error(f"Old output not found for {acc}: {old_output}")
            continue
        if not os.path.exists(new_output):
            logger.error(f"New output not found for {acc}: {new_output}")
            continue

        # Target directory for this accession
        target_dir = os.path.join(output_dir, acc) if args.split else output_dir
        os.makedirs(target_dir, exist_ok=True)

        # Copy or symlink files
        assembly_dest = os.path.join(target_dir, os.path.basename(assembly_path))
        old_dest = os.path.join(target_dir, os.path.basename(old_output))
        new_dest = os.path.join(target_dir, os.path.basename(new_output))

        for src, dst in [(assembly_path, assembly_dest), (old_output, old_dest), (new_output, new_dest)]:
            if args.symlink:
                if os.path.exists(dst) or os.path.islink(dst):
                    os.remove(dst)
                os.symlink(os.path.abspath(src), dst)
            else:
                shutil.copy2(src, dst)

        # Build commands
        orgopt = taxgroup2orgopt(taxgroup_map.get(acc))
        org_arg = f"--organism '{orgopt}' " if orgopt else ""

        asm_name = os.path.basename(assembly_dest)
        old_name = os.path.basename(old_dest)
        new_name = os.path.basename(new_dest)

        cmd_old = f"{amrfinder_old} --nucleotide {asm_name} --output {old_name} --threads {args.threads} {org_arg}--gpipe_org"
        cmd_new = f"{amrfinder_new} --nucleotide {asm_name} --output {new_name} --threads {args.threads} {org_arg}--gpipe_org"

        if args.split:
            # Write individual commands.txt in the accession folder
            with open(os.path.join(target_dir, "commands.txt"), "w") as f_cmd:
                f_cmd.write(f"{cmd_old}\n")
                f_cmd.write(f"{cmd_new}\n")

        all_commands.append((acc, cmd_old, cmd_new))
        logger.info(f"Prepared test case for {acc}")

    # Write combined commands.txt in the main output directory
    commands_file = os.path.join(output_dir, "commands.txt")
    with open(commands_file, "w") as f_cmd:
        for acc, cmd_old, cmd_new in all_commands:
            if len(all_commands) > 1:
                f_cmd.write(f"# {acc}\n")
            f_cmd.write(f"{cmd_old}\n")
            f_cmd.write(f"{cmd_new}\n\n")

    logger.info(f"Wrote command lines to {commands_file}")
    logger.info(f"Successfully prepared test cases in: {output_dir}")


if __name__ == '__main__':
    main()
