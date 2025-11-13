#!/usr/bin/env python3
# run like uv run random_pathogen_asm.py
debug = True

from pathlib import Path
import sys
import shutil
import os
import subprocess
import pandas as pd
from glob import glob


query = """
select i.biosample_acc, i.asm_acc, i.taxgroup_name, i.asm_level, i.Platform, i.Run, i.assembly_method, i.species_taxid
from `ncbi-pathogen-detect.pdbrowser.isolates` i tablesample system(20 percent)
where i.asm_acc is not NULL
order by rand()
limit 1000
"""

prefix = "test6random"

amrfinder_old = "/netmnt/vast01/pd/usr/work/bin/amrfinder"
amrfinder_new = "/home/aprasad/amrfinder/src/amrfinder"

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

def taxgroup2orgopt(taxgroup):
    """Fix some irregularities in converting taxgroup name to the short_name used by AMRFinderPlus

    Args:
        taxgroup: taxgroup returned from BigQuery / browsers

    Returns:
        organism_option: organism option to pass to amrfinder (with --gpipe_org enabled)
    """
    match taxgroup:
        case 'Salmonella enterica':
            return('Salmonella')
        case 'Acinetobacter baumannii':
            return('Acinetobacter')
        case 'E.coli and Shigella':
            return('Escherichia_coli_Shigella')
        case 'Klebsiella pneumoniae':
            return('Klebsiella')
        case 'Serratia marcescens':
            return('Serratia')
        case 'Campylobacter jejuni':
            return('Campylobacter')
        case _:
            return(taxgroup.replace(' ', '_'))
    

# Download a random list of assemblies unless it already exists
if not os.path.exists(prefix):
    # make directory for all analysis outputs
    Path(prefix).mkdir()
    shutil.copy(sys.argv[0], f"{prefix}/pipeline_run.py")

if not os.path.exists(prefix + '/asm_data.tab'):
    # 1. Connect to and use BigQuery to get a list of N assemblies, their accessions and their taxgroups
    from google.cloud import bigquery

    client = bigquery.Client()
    print("--- Running Query to get list of isolates ---")
    query_job = client.query(query)  # API request - starts the query
    rows = query_job.result()  # Wait for the query to finish and get the results
    column_names = [field.name for field in rows.schema]

    with open(prefix + "/asm_data.tab", 'w') as fh:
        # print header
        print("\t".join(column_names), file=fh)
        # print data
        for row in rows:
            values = [str(value) for value in row]
            print("\t".join(values), file=fh)

else:
    print("--- Found directory " + prefix + " not querying for assemblies again ---")


# 2. Download those assemblies using Datasets (or using FTP?)

if not os.path.exists(prefix + '/ncbi_dataset/data'):
    # Path(prefix + '/data').mkdir()
    result = runcmd(f"head -1 {prefix}/asm_data.tab | tr '\t' '\n' | grep -n asm_acc | cut -d: -f1", capture_output=True)
    colnum = int(result.stdout.strip())
    if colnum < 1:
        raise ExceptionType("Error: did not get a column for assembly accession from " + prefix + "/asm_data.tab")
    # get the accessions in a file
    runcmd(f"cut -f {colnum} {prefix}/asm_data.tab | awk " + "'NR>1{print}'" + f"> {prefix}/asm_acc.txt")
    # pass the accessions to datasets download to download the dehydrated data package (only assemblies for now)
    runcmd(f"cd {prefix}; datasets download genome accession --dehydrated --include genome --inputfile asm_acc.txt")
    # Now rehydrate
    runcmd(f"cd {prefix}; pwd; unzip ncbi_dataset.zip")
    runcmd(f"datasets rehydrate --directory {prefix}")
else:
    print("--- Found " + prefix + '/ncbi_dataset/data so skipping downloads ---')


# 3. Run amrfinder old and new
dir_old = prefix + "/amrfinder_old"
if not os.path.exists(dir_old):
    Path(dir_old).mkdir()
dir_new = prefix + "/amrfinder_new"
if not os.path.exists(dir_new):
    Path(dir_new).mkdir()

# This time read data using pandas so I can learn
data_file_path = prefix + '/asm_data.tab'
try:
    df = pd.read_csv(data_file_path, sep='\t')
except FileNotFoundError:
    print(f"Error: Could not find {data_file_path} to open.")
    exit()

for row in df.itertuples():
    # run_old
    assembly = glob(f"{prefix}/ncbi_dataset/data/{row.asm_acc}/{row.asm_acc}*_genomic.fna")
    # breakpoint()
    if not assembly:
        print(f"Error: Could not find assembly for {row.asm_acc}")
    else:
        assembly = assembly[0]
        orgopt = taxgroup2orgopt(row.taxgroup_name)

        output = f"{dir_old}/{row.asm_acc}.amrfinder_old"
        cmd = " ".join([
            f"{amrfinder_old} --nucleotide {assembly} --output {output}",
            f"--threads 4 --organism '{orgopt}' --gpipe_org"])
        if os.path.exists(output) and os.path.getsize(output) > 10:
            print(f"--- Found {output} so skipping running {amrfinder_old} on {row.asm_acc} ---")
        else:
            runcmd(cmd)

        output = f"{dir_new}/{row.asm_acc}.amrfinder_new"
        cmd = " ".join([
            f"{amrfinder_new} --nucleotide {assembly} --output {output}",
            f"--threads 4 --organism '{orgopt}' --gpipe_org"])
        if os.path.exists(output) and os.path.getsize(output) > 10:
            print(f"--- Found {output} so skipping running {amrfinder_new} on {row.asm_acc} ---")
        else:
            runcmd(cmd)


# asm_list = []
# with open(prefix + "/asm_acc.txt") as f:
#     asm_list = f.read().splitlines()
# for asm_acc in asm_list:
#     if not os.path.exists(f"{prefix}/amrfinder_old/





