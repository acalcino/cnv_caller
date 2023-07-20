#!/usr/bin/env python3

import os
import argparse
import subprocess

def submit_qsub_job(args):
    # Get the directory of the script
    script_directory = os.path.dirname(os.path.realpath(__file__))

    # Set the location of the resources folder
    resources_directory = os.path.join(script_directory, "resources")

    # Set the location of the pileup_and_phase.R script
    pileup_and_phase_script = os.path.join(resources_directory, "pileup_and_phase.R")

    # Set the location of the genetic map file
    gmap_location = os.path.join(resources_directory, "genetic_map_hg38_withX.txt.gz")

    # Set the location of the SNP VCF file
    snpvcf_location = os.path.join(resources_directory, "genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz")

    # Check if the SNP VCF file exists
    if not os.path.exists(snpvcf_location):
        # Download the SNP VCF file
        vcf_url = "http://ufpr.dl.sourceforge.net/project/cellsnp/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
        vcf_file = os.path.join(resources_directory, "genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz")
        os.system(f"wget {vcf_url} -O {vcf_file}")

    # Set the location of the panel directory
    paneldir_location = os.path.join(resources_directory, "1000G_hg38")

    # Check if the panel directory exists
    if not os.path.exists(paneldir_location):
        # Download the 1000G_hg38.zip file
        zip_url = "http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip"
        zip_file = os.path.join(resources_directory, "1000G_hg38.zip")
        os.system(f"wget {zip_url} -O {zip_file}")

        # Unzip the downloaded file
        os.system(f"unzip {zip_file} -d {resources_directory}")

        # Remove the downloaded zip file
        os.remove(zip_file)

# Prepare the qsub command
    qsub_cmd = f"""\
#!/bin/bash
#PBS -P pq84
#PBS -q normalbw
#PBS -l walltime={args.walltime},mem={args.mem},jobfs=100GB,ncpus={args.cpus}
#PBS -l storage=gdata/pq84+gdata/u86+scratch/u86+gdata/xx92
"""

    # Export R_LIBS_USER if provided by the user
    if args.r_libs_user:
        qsub_cmd += f"export R_LIBS_USER={args.r_libs_user}\n"

    # Append the pileup_and_phase.R command
    qsub_cmd += f"""\
module load R/4.2.1

Rscript {pileup_and_phase_script} \
--label {args.project_name} \
--samples {args.project_name} \
--bams {args.bam} \
--barcodes {args.barcodes} \
--outdir {args.outdir} \
--gmap {gmap_location} \
--snpvcf {snpvcf_location} \
--paneldir {paneldir_location} \
--ncores {args.cpus}"""

    # Submit the qsub job
    qsub_process = subprocess.Popen("qsub", stdin=subprocess.PIPE)
    qsub_process.communicate(input=qsub_cmd.encode())

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Submit qsub job for pileup_and_phase.R script.')
    
    # Add arguments
    parser.add_argument('--project_name', required=True, help='Name of the project')
    parser.add_argument('--bam', required=True, help='Location of the input BAM file (processed according to GATK best practices and indexed)')
    parser.add_argument('--barcodes', required=True, help='Location of the barcodes.tsv.gz file')
    parser.add_argument('--outdir', required=True, help='Output directory - will be created. Use the full path.')
    parser.add_argument('--cpus', type=int, default=1, help='Number of cores')
    parser.add_argument('--walltime', default='01:00:00', help='Walltime (default: 01:00:00)')
    parser.add_argument('--mem', default='2GB', help='Memory (default: 2GB)')
    parser.add_argument('--r_libs_user', default=None, help='Location of R libraries (default: /g/data/pq84/andrew/libs/r)')
    
    # Parse the command-line arguments
    args = parser.parse_args()

    # Submit the qsub job
    submit_qsub_job(args)

if __name__ == '__main__':
    main()

