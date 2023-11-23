#!/usr/bin/env python3

import argparse
import subprocess
import os

def build_qsub_content(output_location, r_libs_location, walltime, memory, cpus):
    qsub_content = f'''#!/bin/bash
#PBS -P pq84
#PBS -q normalbw
#PBS -l walltime={walltime},mem={memory},jobfs=100GB,ncpus={cpus}
#PBS -l storage=gdata/pq84+gdata/u86+scratch/u86+gdata/xx92

#Make sure R is accessible. If not, uncomment the following line and load the R module
module load R/4.3.1

'''

    if r_libs_location:
        qsub_content += f'export R_LIBS_USER={r_libs_location}\n\n'

    qsub_content += f'R --vanilla < {output_location}/numbat_{project_name}.r\n'

    return qsub_content

def build_numbat_r_content(filtered_feature_matrix_location, project_name, allele_file_location, output_location, max_entropy=0.7, init_k=4, min_cells=10, tau=0.2, min_LLR=25, max_iter=3, ncores=12):
    r_content = f'''library(Seurat)
library(numbat)
setwd("{output_location}")
raw_counts <- Read10X(data.dir = "{filtered_feature_matrix_location}")
raw_seurat <- CreateSeuratObject(counts = raw_counts, project = "{project_name}", min.cells = 0, min.features = 0)
exp_rawdata <- GetAssayData(object = raw_seurat, assay = "RNA", layer = "counts")
allele <- read.table(gzfile("{allele_file_location}"), header = TRUE, sep = "\t")

#Make new expression reference from the blood cell types from the hca expression reference
blood_ref <- ref_hca[, c(1, 3:5, 8:11)]
out <- run_numbat(exp_rawdata, blood_ref, allele, genome = "hg38", t = 1e-05, ncores = {cpus}, plot = TRUE, out_dir = "{output_location}", max_entropy = {max_entropy}, init_k = {init_k}, min_cells = {min_cells}, tau = {tau}, min_LLR = {min_LLR}, max_iter = {max_iter})
'''

    return r_content

# Create the argument parser
parser = argparse.ArgumentParser(description="Script for running cnv_caller", formatter_class=argparse.RawTextHelpFormatter)

# Add additional help information about required libraries
help_text = "Prior to running the script, make sure the following libraries are installed:\n"
help_text += "- Seurat\n"
help_text += "- numbat\n"
parser.epilog = help_text

# Add the command line arguments
parser.add_argument("output_location", help="Output location. Use output folder created by pileup_and_phaser.py")
parser.add_argument("filtered_feature_matrix_location", help="Filtered feature matrix folder location")
parser.add_argument("allele_file_location", help="Allele file location (output of pileup_and_phaser.py). File should end in: _allele_counts.tsv.gz")
parser.add_argument("--r_libs_location", help="R libs location")
parser.add_argument("--walltime", default="02:00:00", help="Walltime (default: 02:00:00)")
parser.add_argument("--memory", default="20GB", help="Memory (default: 20GB)")
parser.add_argument("--cpus", type=int, default=2, help="CPUs (default: 2)")
parser.add_argument("--max_entropy", type=float, default=0.7, help="Max entropy (default: 0.7)")
parser.add_argument("--init_k", type=int, default=4, help="Init k (default: 4)")
parser.add_argument("--min_cells", type=int, default=10, help="Min cells (default: 10)")
parser.add_argument("--tau", type=float, default=0.2, help="Tau (default: 0.2)")
parser.add_argument("--min_LLR", type=int, default=25, help="Min LLR (default: 25)")
parser.add_argument("--max_iter", type=int, default=3, help="Max iter (default: 3)")
parser.add_argument("--project_name", default="project", help="Project name (default: project)")

# Parse the command line arguments
args = parser.parse_args()

# Extract the values from the parsed arguments
output_location = args.output_location
r_libs_location = args.r_libs_location
walltime = args.walltime
memory = args.memory
cpus = args.cpus
filtered_feature_matrix_location = args.filtered_feature_matrix_location
project_name = args.project_name
allele_file_location = args.allele_file_location
max_entropy = args.max_entropy
init_k = args.init_k
min_cells = args.min_cells
tau = args.tau
min_LLR = args.min_LLR
max_iter = args.max_iter

# Get the full path of the required files
script_directory = os.path.dirname(os.path.realpath(__file__))
gmap_location = os.path.join(script_directory, "genetic_map_hg38_withX.txt.gz")
snpvcf_location = os.path.join(script_directory, "genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf")
paneldir_location = os.path.join(script_directory, "1000G_hg38")

# Build the qsub file content
qsub_content = build_qsub_content(output_location, r_libs_location, walltime, memory, cpus)

# Build the numbat_{project_name}.r file content
numbat_r_content = build_numbat_r_content(filtered_feature_matrix_location, project_name, allele_file_location, output_location, max_entropy, init_k, min_cells, tau, min_LLR, max_iter)

# Write the qsub content to a new file named "numbat_{project_name}.qsub" in the output location
qsub_file_path = f'{output_location}/numbat_{project_name}.qsub'
with open(qsub_file_path, 'w') as f:
    f.write(qsub_content)

# Write the numbat_{project_name}.r content to a new file named "numbat_{project_name}.r" in the output location
numbat_r_file_path = f'{output_location}/numbat_{project_name}.r'
with open(numbat_r_file_path, 'w') as f:
    f.write(numbat_r_content)

# Submit the qsub file
subprocess.run(['qsub', qsub_file_path])

print(f"Qsub file created successfully: {qsub_file_path}")

