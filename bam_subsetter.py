#!/usr/bin/env python3

import argparse
import os
import glob
import subprocess
from multiprocessing import Pool
import time


def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Subset bam files according to segs defined by Numbat.')

    # Add the optional arguments
    parser.add_argument('--r_libs_location', default='/g/data/pq84/andrew/libs/r', help='Path to R_LIBS_USER directory. Default: /g/data/pq84/andrew/libs/r')
    parser.add_argument('--numbat_folder', required=True, help='Path to Numbat output folder')
    parser.add_argument('--working_dir', default='haplotypeCaller', help='Name for folder that will be made for output files. Default: haplotypecaller')
    parser.add_argument('--walltime', default='24:00:00', help='Walltime for the qsub job. Default: 24:00:00')
    parser.add_argument('--mem', default='8GB', help='Memory for the qsub job. Default: 8Gb')
    parser.add_argument('--cpus', default='1', help='Number of CPUs for the qsub job. Default: 1')
    parser.add_argument('--input_bam', required=True, help='Path to input gatk best practices corrected BAM file')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Create the working directory
    working_dir = os.path.join(args.numbat_folder, args.working_dir)
    os.makedirs(working_dir, exist_ok=True)

    # Set the R_LIBS_USER environment variable
    os.environ['R_LIBS_USER'] = args.r_libs_location

    # Find the highest number in the segs_consensus files
    segs_files = glob.glob(os.path.join(args.numbat_folder, 'segs_consensus_*.tsv'))
    if segs_files:
        max_iteration = max([int(os.path.splitext(file)[0].split('_')[-1]) for file in segs_files])
    else:
        max_iteration = 0

    # Construct the R script
    r_script = f'''
    library(glue)
    library(dplyr)
    library(numbat)

    working_dir <- "{os.path.join(args.numbat_folder, args.working_dir)}"
    setwd(working_dir)
    nb <- Numbat$new(glue('{args.numbat_folder}'), i = {max_iteration})
    values <- unique(nb$clone_post %>% select(GT_opt))$GT_opt
    
    write(values, file="values.txt")

    for (value in values) {{
        filtered_data <- nb$clone_post %>% filter(GT_opt == value) %>% select(cell)
        if (value == "") {{
            value <- "null"
        }}
        file_name <- glue(value, "_cell_barcodes.txt", sep = "_")
        write.table(filtered_data, file_name, col.names = FALSE, quote = FALSE, row.names = FALSE)
    }}
    '''

    # Save the R script to a file
    r_script_file = os.path.join(working_dir, f'{args.working_dir}.R')
    with open(r_script_file, 'w') as file:
        file.write(r_script)

    # Create the qsub file
    qsub_file = os.path.join(working_dir, f'{args.working_dir}.qsub')
    with open(qsub_file, 'w') as file:
        file.write(f'''#!/bin/bash
#PBS -P pq84
#PBS -q normalbw
#PBS -l walltime={args.walltime},mem={args.mem},jobfs=100GB,ncpus={args.cpus}
#PBS -l storage=gdata/pq84+gdata/u86+scratch/u86+gdata/xx92

module load R/4.2.1

R --vanilla -f {r_script_file}
''')

    # Submit the qsub job
    subprocess.run(f'qsub {qsub_file}', shell=True)

    # Wait until the job has completed
    job_completed = False
    while not job_completed:
        time.sleep(10)  # Wait for 10 seconds
        qstat_output = subprocess.run('qstat', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
        if args.working_dir.encode() not in qstat_output:
            job_completed = True

    # Read the values from the file
    with open(os.path.join(working_dir, 'values.txt'), 'r') as file:
        values = file.read().splitlines()

    # Run subset-bam_linux for each value in parallel
    with Pool() as pool:
        pool.starmap(run_subset_bam, zip(values, [args] * len(values), [working_dir] * len(values)))

# Define a function to execute the subset-bam_linux command
def run_subset_bam(value, args, working_dir):
    if value == "":
        value = "null"
    output_bam = f'{value}.bam'

    # Create the qsub file. The replace command is there because you can't use commas in qsub file names
    qsub_file = os.path.join(working_dir, f'{value.replace(",", "-")}.qsub')
    with open(qsub_file, 'w') as file:
        file.write(f'''#!/bin/bash
#PBS -P pq84
#PBS -q normalbw
#PBS -l walltime={args.walltime},mem={args.mem},jobfs=100GB,ncpus={args.cpus}
#PBS -l storage=gdata/pq84+gdata/u86+scratch/u86+gdata/xx92

module load samtools

subset-bam_linux --bam {args.input_bam} --cell-barcodes {os.path.join(working_dir, value + "_cell_barcodes.txt")} --out-bam {os.path.join(working_dir, output_bam)}

samtools index -@{args.cpus} {os.path.join(working_dir, value + ".bam")}
''')

    # Submit the qsub job
    subprocess.run(f'qsub {qsub_file}', shell=True)

#    # Wait until the job has completed
#    job_completed = False
#    while not job_completed:
#        time.sleep(10)  # Wait for 10 seconds
#        qstat_output = subprocess.run('qstat', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
#        if args.working_dir.encode() not in qstat_output:
#            job_completed = True
#
#    command = f'subset-bam_linux --bam {args.input_bam} --cell-barcodes {os.path.join(working_dir, value + "_cell_barcodes.txt")} --out-bam {os.path.join(working_dir, output_bam)}'
#    subprocess.run(command, shell=True)

if __name__ == '__main__':
    main()
