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

    # Create the 'intervals' directory if it doesn't exist
    intervals_directory = 'intervals'
    os.makedirs(intervals_directory, exist_ok=True)

    # Create the 'vcfs' directory if it doesn't exist
    vcfs_directory = os.path.join(parent_folder, 'vcfs')
    os.makedirs(vcfs_directory, exist_ok=True)

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
#                chromosome_number = 'chr' + chromosome_number   #Adds chr to the chromosome number. Comment if unnecessary
                chromosome_numbers.append(chromosome_number)  # Append the chromosome number to the list
                start = int(line.split('\t')[5])
                stop = int(line.split('\t')[6])
                seg_cons = line.split('\t')[26]
                length = stop - start
                file.write('\t'.join([chromosome_number, str(start), str(stop), seg_cons, str(length), '+']) + '\n')

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
        job_ids = []  # Initialize the list here
        job_id = submit_output.strip().split('.')[0]
        job_ids.append(job_id)  # Add the job ID to the list

    # Wait for the qsub jobs to finish
#    while True:
#        time.sleep(10)  # Sleep for 10 seconds before checking the status again
#        process = subprocess.run(['qstat'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#        output = process.stdout.decode('utf-8')
#        if all(job_id not in output for job_id in job_ids):
#        #if 'qw' not in output and 'r' not in output:  # Check if there are no pending or running jobs
#            break

    # Wait for the qsub jobs to finish
    # Hopefully this works better
    for job_id in job_ids:
        while True:
            time.sleep(10)
            process = subprocess.run(['qstat', job_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output = process.stdout.decode('utf-8')
            if job_id not in output:
                break

    print("All qsub jobs have finished.")

    # Read the 'intervals.bed' file and create separate '.bed' files
    intervals_bed_path = os.path.join(args.subset_folder, 'intervals.bed')
    with open(intervals_bed_path, 'r') as intervals_bed_file:
        for line in intervals_bed_file:
            fields = line.strip().split('\t')
            # Extract the fourth column value as the filename
            filename = fields[3] + '.bed'
            # Create the full path for the new file in the 'intervals' directory
            output_path = os.path.join(intervals_directory, filename)
            # Write the line to the new '.bed' file
            with open(output_path, 'w') as output_file:
                output_file.write(line)

     # Iterate through VCF files to make sure any vcfs are indexed
    for vcf_filename in os.listdir(parent_folder):
        if vcf_filename.endswith("_het.vcf") or vcf_filename.endswith("_hom.vcf"):
            # Extract the value from the VCF filename
            value = os.path.splitext(vcf_filename)[0]

            # Zip the VCF file using bgzip
            subprocess.run(['bgzip', '-f', os.path.join(parent_folder, vcf_filename)])

            # Index the VCF file using bcftools
            subprocess.run(['bcftools', 'index', os.path.join(parent_folder, vcf_filename + '.gz')])

    # Iterate through indexed VCF files to create vcf files for each cnv for each clone
    for vcf_gz_filename in os.listdir(parent_folder):
        if vcf_gz_filename.endswith("_het.vcf.gz") or vcf_gz_filename.endswith("_hom.vcf.gz"):
            # Extract the value from the VCF filename
            base_filename = os.path.splitext(vcf_gz_filename)[0]  # Remove the .gz extension
            value_gz = base_filename.split('_')[0]  # Split on '_' and take the first part

            # Construct the command to extract values from 'bulk_clones_final.tsv.gz'
            bulk_clones_final_path = os.path.join(os.path.dirname(parent_folder), 'bulk_clones_final.tsv.gz')
            cmd = (
                f"zcat {bulk_clones_final_path} | "
                f"awk -F '\t' '$29 == \"{value_gz}\"' | "
                f"awk -F '\t' '$77 != \"neu\" {{print $47}}' | uniq"
            )
            output = subprocess.check_output(cmd, shell=True, text=True)

            # Split the 'output' into lines and iterate through each line
            for output_line in output.split('\n'):
                if output_line:
                    # Construct the tabix command for this output value
                    tabix_cmd = (
                        f"tabix -R {intervals_directory}/{output_line}.bed {vcf_gz_filename} > {vcfs_directory}/{value_gz}_{output_line}_h{base_filename.split('_h', 1)[-1]}"
                    )
                    # Run the tabix command
                    subprocess.run(tabix_cmd, shell=True)
 
                    print(f"Processed {vcf_gz_filename} for {output_line}")
 
    # Create vcfs for every location in intervals.bed for null_het.vcf.gz and null_hom.vcf.gz
    # Define paths to null VCF files
    null_het_vcf = os.path.join(parent_folder, 'null_het.vcf.gz')
    null_hom_vcf = os.path.join(parent_folder, 'null_hom.vcf.gz')
 
    # Iterate through bed files in intervals directory
    for bed_filename in os.listdir(intervals_directory):
        if bed_filename.endswith(".bed"):
            bed_path = os.path.join(intervals_directory, bed_filename)
            bed_name = os.path.splitext(bed_filename)[0]
 
            # Construct the tabix commands for this output value
            tabix_cmd_het = (
                f"tabix -R {intervals_directory}/{bed_name}.bed {null_het_vcf} > {vcfs_directory}/null_{bed_name}_het.vcf"
            )
            tabix_cmd_hom = (
                f"tabix -R {intervals_directory}/{bed_name}.bed {null_hom_vcf} > {vcfs_directory}/null_{bed_name}_hom.vcf"
            )
            # Run the tabix command
            subprocess.run(tabix_cmd_het, shell=True)
            subprocess.run(tabix_cmd_hom, shell=True)
 
            print(f"Generated null VCFs for {bed_name}")    

    # Define a dictionary to store counts
    counts = {}
 
    # Iterate through indexed VCF files
    for het_hom_vcf_gz_filename in os.listdir(vcfs_directory):
        if het_hom_vcf_gz_filename.endswith(".vcf"):
            # Extract the value from the VCF filename
            value_gz, output_line = os.path.splitext(het_hom_vcf_gz_filename)[0].split('_')[:2]

            # Determine if it's 'het' or 'hom' based on the filename
            file_type = 'het' if 'het' in het_hom_vcf_gz_filename else 'hom'

            # Read the VCF file and count lines
            with open(os.path.join(vcfs_directory, het_hom_vcf_gz_filename), 'r') as vcf_file:
                line_count = sum(1 for line in vcf_file)

            # Update the counts dictionary
            if value_gz not in counts:
                counts[value_gz] = {'het': {}, 'hom': {}}

            if output_line not in counts[value_gz][file_type]:
                counts[value_gz][file_type][output_line] = 0

            # Update the count
            counts[value_gz][file_type][output_line] += line_count

    # Define a dictionary to store null ratios
    null_ratios = {}

    # Calculate null ratios
    for bed_filename in os.listdir(intervals_directory):
        if bed_filename.endswith(".bed"):
            bed_name = os.path.splitext(bed_filename)[0]

            # Calculate null ratios
            null_het_file = f"{vcfs_directory}/null_{bed_name}_het.vcf"
            null_hom_file = f"{vcfs_directory}/null_{bed_name}_hom.vcf"

            with open(null_het_file, 'r') as null_het_vcf, open(null_hom_file, 'r') as null_hom_vcf:
                het_count = sum(1 for line in null_het_vcf)
                hom_count = sum(1 for line in null_hom_vcf)
                ratio = het_count / hom_count if hom_count != 0 else 'N/A'

                null_ratios[bed_name] = ratio

    # Define the output file path
    output_file_path = os.path.join(parent_folder, 'output_counts.tsv')

    # Write the counts to the output file
    with open(output_file_path, 'w') as output_file:
        output_file.write("clone\tcnv\thet number\thom number\tratio\tnull ratio\n")
        for value_gz, output_lines in counts.items():
            for output_line, het_count in output_lines['het'].items():
                if not value_gz.startswith('null'):  # Exclude rows starting with "null"
                    hom_count = output_lines['hom'].get(output_line, 0)
                    ratio = het_count / hom_count if hom_count != 0 else 'N/A'

                    # Fetch the null ratio
                    null_ratio = null_ratios.get(output_line, 'N/A')

                    output_file.write(f"{value_gz}\t{output_line}\t{het_count}\t{hom_count}\t{ratio}\t{null_ratio}\n")

    print(f"Counts written to {output_file_path}")


if __name__ == '__main__':
    main()

