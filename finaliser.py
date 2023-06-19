#!/usr/bin/env python3

import argparse
import os
import glob
import subprocess
import time

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Get some final numbers from your analysis. Make sure you have loaded the bcftool module prior to running. Also use python v3.7 or greater.')

    # Add the optional arguments
    parser.add_argument('--subset_folder', required=True, help='Path to folder containing bam subsets')
    parser.add_argument('--genome', required=True, help='Path to genome fasta')
    parser.add_argument('--walltime', default='24:00:00', help='Walltime for the qsub job. Default: 24:00:00')
    parser.add_argument('--mem', default='8GB', help='Memory for the qsub job. Default: 8GB')
    parser.add_argument('--cpus', default='1', help='Number of CPUs for the qsub job. Default: 1')


    # Parse the command-line arguments
    args = parser.parse_args()

    # Get the parent folder of the subset folder
    parent_folder = os.path.dirname(args.subset_folder)

    # Find the highest number in the segs_consensus files
    segs_files = glob.glob(os.path.join(parent_folder, 'segs_consensus_*.tsv'))
    if segs_files:
        max_iteration = max([int(os.path.splitext(file)[0].split('_')[-1]) for file in segs_files])
    else:
        max_iteration = 0

    # Iterate over the segs_consensus files
    segs_file = os.path.join(parent_folder, f'segs_consensus_{max_iteration}.tsv')
    with open(segs_file, 'r') as file:
        lines = file.readlines()

    # Write the intervals.bed file
    chromosome_numbers = []
    intervals_file = os.path.join(args.subset_folder, 'intervals.bed')
    with open(intervals_file, 'w') as file:
        for line in lines[1:]:
            if 'neu' not in line:
                chromosome_number = line.split('\t')[1]
                chromosome_number = 'chr' + chromosome_number   #Adds chr to the chromosome number. Comment if unnecessary
                chromosome_numbers.append(chromosome_number)  # Append the chromosome number to the list
                start = int(line.split('\t')[5])
                stop = int(line.split('\t')[6])
                seg = line.split('\t')[2]
                length = stop - start
                file.write('\t'.join([chromosome_number, str(start), str(stop), seg, str(length), '+']) + '\n')

    # Create qsub files to run haplotypecaller and designate column names for final table output
    bams = []
    for file_name in os.listdir(args.subset_folder):
        if file_name.endswith(".bam"):
            bams.append(os.path.splitext(file_name)[0])

    column_names = [bam_name for bam_name in bams]

    table_data = {
        "Chromosome": chromosome_numbers,
    }

    for bam_name in bams:
        table_data[bam_name] = {
            "Het Line Counts": [],
            "Hom Line Counts": [],
            "Het SNP Line Counts": [],
            "Hom SNP Line Counts": []
        }

    for bam_name in bams:
        qsub_file_name = bam_name.replace(",", "-")
        qsub_file = os.path.join(args.subset_folder, f"{qsub_file_name}.qsub")
        with open(qsub_file, 'w') as file:
            file.write(f'''#!/bin/bash
#PBS -P pq84
#PBS -q normalbw
#PBS -l walltime={args.walltime},mem={args.mem},jobfs=100GB,ncpus={args.cpus}
#PBS -l storage=gdata/pq84+gdata/u86+scratch/u86+gdata/xx92

module load gatk
module load bcftools

gatk HaplotypeCaller -R {args.genome} -I {os.path.join(args.subset_folder, f"{bam_name}.bam")} -L {os.path.join(args.subset_folder, "intervals.bed")} --native-pair-hmm-threads {args.cpus} -O {os.path.join(args.subset_folder, f"{bam_name}.vcf")}

bcftools view -i 'GT="het"' {os.path.join(args.subset_folder, f"{bam_name}.vcf")} >{os.path.join(args.subset_folder, f"{bam_name}" + "_het.vcf")}
bcftools view -i 'GT="hom"' {os.path.join(args.subset_folder, f"{bam_name}.vcf")} >{os.path.join(args.subset_folder, f"{bam_name}" + "_hom.vcf")}
''')

        # Submit the qsub job and capture the output
        submit_process = subprocess.run(['qsub', qsub_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        submit_output = submit_process.stdout.decode('utf-8')
        #subprocess.run(['qsub', qsub_file])

        # Extract the job ID from the output
        job_id = submit_output.strip().split('.')[0]
        job_ids.append(job_id)  # Add the job ID to the list

    # Wait for the qsub jobs to finish
    while True:
        time.sleep(10)  # Sleep for 10 seconds before checking the status again
        process = subprocess.run(['qstat'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = process.stdout.decode('utf-8')
        if all(job_id not in output for job_id in job_ids):
        #if 'qw' not in output and 'r' not in output:  # Check if there are no pending or running jobs
            break

    print("All qsub jobs have finished.")

    #Creates a list called formatted_line which is used to produce the header for the final table
    intervals_file = os.path.join(args.subset_folder, 'intervals.bed')
    formatted_lines = []  # List to store the formatted lines

    with open(intervals_file, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            column1 = columns[0]
            column2 = columns[1]
            column3 = columns[2]
            formatted_line = f"{column1}:{column2}-{column3}"
            formatted_lines.append(formatted_line)

    #For the final table
    for bam_name in bams:
        het_vcf_file = os.path.join(args.subset_folder, f"{bam_name}_het.vcf")
        hom_vcf_file = os.path.join(args.subset_folder, f"{bam_name}_hom.vcf")
        het_line_counts = get_matching_line_counts(het_vcf_file, chromosome_numbers)
        hom_line_counts = get_matching_line_counts(hom_vcf_file, chromosome_numbers)
        het_snp_line_counts = get_snp_line_counts(het_vcf_file, chromosome_numbers)
        hom_snp_line_counts = get_snp_line_counts(hom_vcf_file, chromosome_numbers)

        table_data[bam_name]["Het Line Counts"] = het_line_counts
        table_data[bam_name]["Hom Line Counts"] = hom_line_counts
        table_data[bam_name]["Het SNP Line Counts"] = het_snp_line_counts
        table_data[bam_name]["Hom SNP Line Counts"] = hom_snp_line_counts

    # Create a dictionary to store the transposed data
    transposed_data = {}

    # Transpose the data
    for bam_name in bams:
        for chromosome_number in chromosome_numbers:
            het_count = table_data[bam_name]["Het Line Counts"][chromosome_number]
            hom_count = table_data[bam_name]["Hom Line Counts"][chromosome_number]
            het_snp_count = table_data[bam_name]["Het SNP Line Counts"][chromosome_number]
            hom_snp_count = table_data[bam_name]["Hom SNP Line Counts"][chromosome_number]
            row_values = [f"Het: {het_count}", f"Hom: {hom_count}", f"Het SNP: {het_snp_count}", f"Hom SNP: {hom_snp_count}"]
            if bam_name in transposed_data:
                transposed_data[bam_name].append(row_values)
            else:
                transposed_data[bam_name] = [row_values]

    # Print the transposed table
    header = "\t".join(["seg"] + formatted_lines)
    print(header)
    for bam_name in bams:
        rows = []
        for i in range(len(chromosome_numbers)):
            row_values = [bam_name] + [transposed_data[bam_name][j][i] for j in range(len(transposed_data[bam_name]))]
            rows.append("\t".join(row_values))
        print("\n".join(rows))

    #Save the table to a file
    output_file = os.path.join(args.subset_folder, "final_numbers.tsv")

    with open(output_file, "w") as file:
        header = "\t".join(["seg"] + formatted_lines)
        file.write(header + "\n")

        for bam_name in bams:
            rows = []
            for i in range(len(chromosome_numbers)):
                row_values = [bam_name] + [transposed_data[bam_name][j][i] for j in range(len(transposed_data[bam_name]))]
                rows.append("\t".join(row_values))
            file.write("\n".join(rows) + "\n")

    print(f"Transposed table saved to {output_file}.")

def get_matching_line_counts(file_path, chromosome_numbers):
    line_counts = {chromosome_number: 0 for chromosome_number in chromosome_numbers}
    with open(file_path, 'r') as file:
        for line in file:
            for chromosome_number in chromosome_numbers:
                if line.startswith(chromosome_number + '\t'):
                    line_counts[chromosome_number] += 1
                    break
    return line_counts

def get_snp_line_counts(file_path, chromosome_numbers):
    line_counts = {chromosome_number: 0 for chromosome_number in chromosome_numbers}
    for chromosome_number in chromosome_numbers:
        command = f"bcftools view -i 'TYPE=\"snp\"' {file_path} | grep -P '^{chromosome_number}\t' | wc -l"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            count = int(result.stdout.strip())
            line_counts[chromosome_number] = count
    return line_counts


if __name__ == '__main__':
    main()

