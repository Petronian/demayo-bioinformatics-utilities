#!/bin/env python

# A simple Python script to screen SAM files, find the proportion of mapped and
# unmapped reads, and find the aligned species using BLAST.
#
# Information about usage limits:
# https://www.nlm.nih.gov/ncbi/workshops/2022-10_Basic-Web-BLAST/how-blast-works.html
#
# Only takes in FASTQ files and searches them up on BLAST.

# Builtin libraries.
import argparse
import re
import os
import gzip
from subprocess import run
from datetime import datetime
from sys import exit
from pathlib import Path
from tempfile import mkstemp

# Installed libraries.
try:
    import numpy as np
    import pandas as pd
    from joblib import Parallel, delayed
except ImportError:
    print("You need numpy, pandas, and joblib installed for this program to work:\n"
          "pip install numpy pandas joblib")
    exit(1)

# Constant definitions.
MAX_NTIDES = 1000000 # Not sure if this is needed anymore
UNMAPPED_READS = 4

# Define and parse arguments.
parser = argparse.ArgumentParser(prog = "NCBI SAMScreen",
                                 description = "Runs some QC metrics on alrready-aligned SAM files..")
parser.add_argument('csvfile', type = str)
parser.add_argument('-db', '--database', type = str, required = True)
parser.add_argument('-o', '--outdir', type = str, default = '.')
parser.add_argument('-e', '--blastn-outfile-extension', type = str, default = '.txt')
parser.add_argument('-n', '--nreads', type = int, default = 100)
parser.add_argument('-s', '--seed', type = int, default = None)
parser.add_argument('-p', '--num-cores', type = int, default = 1)
args, blast_args = parser.parse_known_args()


# Process the arguments.
db = args.database
n_reads_req = args.nreads
num_cores = args.num_cores
outfile_ext = args.blastn_outfile_extension

# Add . if extension is present (implying not folder) but no '.'
if outfile_ext and not outfile_ext.startswith('.'):
    outfile_ext = '.' + outfile_ext

# Other useful variables.
rng = np.random.default_rng(seed = args.seed)
sam_pattern = re.compile(pattern = r'^(?:\S+\s){1}(\S+\s)(?:\S+\s){7}(\w+)')
outdir = Path(args.outdir)

def check_sam_flags(samfile_line, flags, invert = False):
    # Will allow passthrough if 0
    result = (int(re.match(pattern = sam_pattern, string = samfile_line).group(1)) & flags) == flags
    return result if not invert else not result

def get_datetime():
    return datetime.now().ctime()

# Make a helper function to construct queries from SAM file line numbers.
# Assumes lines is a sorted structure. Could be a stack/queue.
def query_constructor(sam_fd, lines, read_type, debug):

    # For ease of logging later.
    read_type = read_type.lower()
    assert read_type in ('mapped', 'unmapped')

    # Scrub through the file and build the query.
    sam_fd.seek(0)
    n_reads_act = 0
    n_ntides_act = 0
    query = ""
    i = 0
    line = sam_fd.readline().decode()
    while line:
        if i == lines[0]:
            line = line.rstrip()
            # Clean up line with SAM file, only extract the read part.
            line = re.match(pattern = sam_pattern, string = line).group(2)
            if n_ntides_act + len(line) > MAX_NTIDES:
                debug(f"For {read_type} reads: Truncating number of reads to avoid exceeding "
                      f"BLAST limits (requested {n_reads_req}, actual {n_reads_act}). ")
                break
            query = query + line + "\n"
            n_ntides_act = n_ntides_act + len(line)
            n_reads_act = n_reads_act + 1
            del lines[0]
        if len(lines) == 0:
            break
        i = i + 1
        line = sam_fd.readline().decode()
    if len(lines) != 0:
        debug(f"For {read_type} reads: Selected less than the requested number of reads due "
              f"to an unspecified complication. (Requested {n_reads_req}, actual {n_reads_act}.)")
        
    return query.rstrip()

def query_runner(query, outfn, *args):
    fd, path = mkstemp(text = True)
    try:
        with os.fdopen(fd, mode = 'w') as f:
            f.write(query)
        args = ['blastn', '-query', path, '-out', outfn, '-db', db] + \
               [arg.replace('--', '-') for arg in blast_args]
        print(f'[{get_datetime()}] Running external program: ' + ' '.join(args))
        run(args)
    finally:
        os.remove(path)
    return os.path.exists(outfn)

# Define the function that will be used for SAM work.
def process_samfile(fn):

    # Create the return object.
    data_out = {
        'status': '',
        'percent_no_hit': None,
        'blast_query_nreads_mapped': None,
        'blast_query_nreads_unmapped': None,
        'blast_results_mapped': None,
        'blast_results_unmapped': None,
    }

    # Just to give the user some comfort.
    base_fn = Path(fn).stem
    print(f'[{get_datetime()} - {base_fn}] Now processing sample \'{base_fn}\'... ')

    # Screen to see if the file is a valid SAMfile.
    is_sam = re.match(pattern = r'.*\.sam\.?.*', string = fn) is not None
    if not is_sam:
        data_out['status'] = 'Error: not a valid SAMfile.'
        return

    # Screen to see if gzipped.
    is_gzipped = re.search(pattern = r'.*\.gz(ip)?$', string = fn) is not None

    # Read through the SAM file and extract relevant statistics.
    open_fn = open if not is_gzipped else gzip.open
    mapped_lines = []
    unmapped_lines = []
    with open_fn(fn, mode = 'rb') as fd:

        # Scan through the SAM file.
        num_lines = 0
        line = fd.readline().decode()
        while line:
            if not line.startswith('@'):
                if check_sam_flags(line, flags = UNMAPPED_READS):
                    unmapped_lines.append(num_lines)
                else:
                    mapped_lines.append(num_lines)
            num_lines = num_lines + 1
            line = fd.readline().decode()
        n_reads_file = len(mapped_lines) + len(unmapped_lines)
        n_unmapped_reads_file = len(unmapped_lines)
        print(f'[{get_datetime()} - {base_fn}] Found {n_reads_file} reads ({n_unmapped_reads_file} unmapped) to choose from.')
        if len(mapped_lines) < n_reads_req or len(unmapped_lines) < n_reads_req:
            data_out['status'] = data_out['status'] + 'Warning: fewer than requested number of reads in file. See CSV output. '

        # Report the percent hit and miss.
        data_out['percent_no_hit'] = n_unmapped_reads_file / n_reads_file * 100

        # Select random lines.
        sel_mapped_lines = sorted(rng.choice(mapped_lines, size = n_reads_req, replace = False) \
                                  if len(mapped_lines) > n_reads_req else mapped_lines)
        sel_unmapped_lines = sorted(rng.choice(unmapped_lines, size = n_reads_req, replace = False) \
                                    if len(unmapped_lines) > n_reads_req else unmapped_lines)

        # Report sizes of queries.
        data_out['blast_query_nreads_mapped'] = len(sel_mapped_lines)
        data_out['blast_query_nreads_unmapped'] = len(sel_unmapped_lines)

        # Build queries for each.
        def debug_fn(str): data_out['status'] = data_out['status'] + str + ' '
        query_mapped = query_constructor(fd, lines = sel_mapped_lines, read_type = 'mapped', debug = debug_fn)
        query_unmapped = query_constructor(fd, lines = sel_unmapped_lines, read_type = 'unmapped', debug = debug_fn)
        
        print(f'[{get_datetime()} - {base_fn}] Constructed queries. Submitting {data_out["blast_query_nreads_mapped"]} '
            f'mapped reads and {data_out["blast_query_nreads_unmapped"]} unmapped reads to BLAST...')

        # Submit to BLAST offline tool.
        blast_mapped_outfn = base_fn + "_blast_mapped_results" + outfile_ext
        blast_unmapped_outfn = base_fn + "_blast_unmapped_results" + outfile_ext
        mapped_ok = query_runner(query_mapped, outfn = blast_mapped_outfn)
        unmapped_ok = query_runner(query_unmapped, outfn = blast_unmapped_outfn)
        data_out['blast_results_mapped'] = blast_mapped_outfn if mapped_ok else "Error with blastn"
        data_out['blast_results_unmapped'] = blast_unmapped_outfn if unmapped_ok else "Error with blastn"

        print(f'[{get_datetime()} - {base_fn}] Handled BLAST responses. Done.')

        return data_out

if __name__ == "__main__":
    os.makedirs(outdir, exist_ok = True)
    meta_df = pd.read_csv(args.csvfile)
    fns = meta_df.iloc[:,0]
    print(f"[{get_datetime()}] Successfully read input CSV and found {len(fns)} SAM file(s).")

    outs = Parallel(n_jobs = num_cores)(delayed(process_samfile)(fn) for fn in fns)
    out_df = pd.DataFrame({fn: out for fn, out in zip(fns, outs)}).T.reset_index(drop = True)
    meta_df = pd.concat([meta_df, out_df], axis = 1)
    meta_df.to_csv(outdir / "results.csv", index = False)








