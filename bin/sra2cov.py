import subprocess
import os
import sys
import argparse

def run_command(command, description):
    print(f"Running: {description}")
    try:
        result = subprocess.run(command, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.stdout.decode())
        print(result.stderr.decode())
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while {description}: {e}\n{e.stderr.decode()}")
        sys.exit(1)

def check_executable(executable):
    command = f"command -v {executable}"
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print(f"Error: {executable} is not installed or not in your PATH.")
        return False
    print(f"{executable} is available and executable.")
    return True

def download_sra(accession, outdir, threads):
    read1_files = []
    read2_files = []
    single_end_files = []

    accessions = accession.split(';')
    for acc in accessions:
        read1 = os.path.join(outdir, f"{acc}_1.fastq")
        read2 = os.path.join(outdir, f"{acc}_2.fastq")
        single_end_read = os.path.join(outdir, f"{acc}.fastq")

        # Check if files already exist
        if os.path.exists(read1) and os.path.exists(read2):
            print(f"Files {read1} and {read2} already exist. Skipping download for {acc}.")
            read1_files.append(read1)
            read2_files.append(read2)
        elif os.path.exists(single_end_read):
            print(f"File {single_end_read} already exists. Assuming single-end read for {acc}. Skipping download.")
            single_end_files.append(single_end_read)
        else:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            command_download = f"fasterq-dump {acc} --split-files --threads {threads} --outdir {outdir}"
            run_command(command_download, f"downloading SRA data for {acc}")

            if os.path.exists(single_end_read):
                print("Only single-end reads available.")
                single_end_files.append(single_end_read)
            else:
                print("Paired-end reads available.")
                read1_files.append(read1)
                read2_files.append(read2)
    
    return read1_files, read2_files, single_end_files

def run_fastp_single(input1, output1, threads):
    html_report = os.path.join(os.path.dirname(output1), "fastp_report.html")
    json_report = os.path.join(os.path.dirname(output1), "fastp_report.json")
    
    command_fastp = f"fastp -i {input1} -o {output1} -w {threads} -h {html_report} -j {json_report}"
    run_command(command_fastp, "trimming adapters and low-quality regions with fastp (single-end)")

def run_fastp_paired(input1, input2, output1, output2, threads):
    html_report = os.path.join(os.path.dirname(output1), "fastp_report.html")
    json_report = os.path.join(os.path.dirname(output1), "fastp_report.json")
    
    command_fastp = f"fastp -i {input1} -I {input2} -o {output1} -O {output2} -w {threads} -h {html_report} -j {json_report}"
    run_command(command_fastp, "trimming adapters and low-quality regions with fastp (paired-end)")

def map_reads_single(trimmed_read1, bowtie_index, threads, outdir, accession):
    sam_output = os.path.join(outdir, f"{accession}.sam")
    command_bowtie = f"bowtie2 -x {bowtie_index} -U {trimmed_read1} -S {sam_output} --no-unal -p {threads}"
    run_command(command_bowtie, "mapping reads with Bowtie2 (single-end)")
    return sam_output

def map_reads_paired(trimmed_read1, trimmed_read2, bowtie_index, threads, outdir, accession):
    sam_output = os.path.join(outdir, f"{accession}.sam")
    command_bowtie = f"bowtie2 -x {bowtie_index} -1 {trimmed_read1} -2 {trimmed_read2} -S {sam_output} --no-unal -p {threads}"
    run_command(command_bowtie, "mapping reads with Bowtie2 (paired-end)")
    return sam_output

def sam_to_sorted_bam(sam_file, threads, outdir, accession):
    bam_output = os.path.join(outdir, f"{accession}.bam")
    sorted_bam_output = os.path.join(outdir, f"{accession}.sorted.bam")
    
    command_sam_to_bam = f"samtools view -bS -q 25 {sam_file} -@ {threads} -o {bam_output}"
    run_command(command_sam_to_bam, "converting SAM to BAM with minimum quality 25")
    
    command_sort_bam = f"samtools sort {bam_output} -@ {threads} -o {sorted_bam_output}"
    run_command(command_sort_bam, "sorting BAM file")
    
    command_index_bam = f"samtools index {sorted_bam_output} -@ {threads}"
    run_command(command_index_bam, "indexing sorted BAM file")
    
    return sorted_bam_output

def merge_bam_files(bam_files, threads, outdir, merged_bam_name):
    merged_bam = os.path.join(outdir, merged_bam_name)
    command_merge_bam = f"samtools merge -@ {threads} {merged_bam} " + " ".join(bam_files)
    run_command(command_merge_bam, "merging BAM files")
    
    command_index_bam = f"samtools index {merged_bam} -@ {threads}"
    run_command(command_index_bam, "indexing merged BAM file")
    
    return merged_bam

def calculate_breadth_coverage(sorted_bam, outdir, accession):
    breadth_output = os.path.join(outdir, f"{accession}_breadth.txt")
    command_breadth = f"breadth_depth.py {sorted_bam} > {breadth_output}"
    run_command(command_breadth, "calculating breadth of coverage with breadth_depth.py")
    return breadth_output

def calculate_reads_per_contig(sorted_bam, outdir, accession):
    idxstats_output = os.path.join(outdir, f"{accession}_idxstats.txt")
    command_idxstats = f"samtools idxstats {sorted_bam} > {idxstats_output}"
    run_command(command_idxstats, "calculating number of reads per contig with samtools idxstats")
    return idxstats_output

def remove_tmp_files(outdir, accessions):
    for accession in accessions.split(';'):
        fastq_files = [
            os.path.join(outdir, f"{accession}_1.fastq"),
            os.path.join(outdir, f"{accession}_2.fastq"),
            os.path.join(outdir, f"{accession}.fastq"),
            os.path.join(outdir, f"{accession}_1.trimmed.fastq"),
            os.path.join(outdir, f"{accession}_2.trimmed.fastq"),
            os.path.join(outdir, f"{accession}.trimmed.fastq")
        ]
        sam_file = os.path.join(outdir, f"{accession}.sam")
        
        for file in fastq_files:
            if os.path.exists(file):
                os.remove(file)
                print(f"Removed {file}")
        
        if os.path.exists(sam_file):
            os.remove(sam_file)
            print(f"Removed {sam_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check for tools, download SRA reads, trim with fastp, map with Bowtie2, and process BAM files.")
    parser.add_argument("--accession", "-c", required=True, help="Semicolon-separated list of SRA accession numbers")
    parser.add_argument("--bowtie-index", "-i", required=True, help="Bowtie2 indexed reference")
    parser.add_argument("--threads", "-t", type=int, default=1, help="Number of threads to use for downloading, trimming, mapping, and BAM processing")
    parser.add_argument("--output", "-o", required=True, help="Output directory to save the reads")
    parser.add_argument("--sample-name", "-s", required=True, help="Name for the final BAM file and output files")
    parser.add_argument("--remove-tmp", action="store_true", help="Remove temporary fastq and sam files after processing")

    args = parser.parse_args()
    
    # Check if required executables are available
    tools = ["fasterq-dump", "fastp", "bowtie2", "samtools", "breadth_depth.py"]
    for tool in tools:
        if not check_executable(tool):
            print(f"Please install {tool} and make sure it is in your PATH.")
            sys.exit(1)
    
    # Define file paths
    outdir = args.output
    accessions = args.accession
    threads = args.threads
    bowtie_index = args.bowtie_index
    sample_name = args.sample_name
    
    # Create output directory if it does not exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Download SRA data
    read1_files, read2_files, single_end_files = download_sra(accessions, outdir, threads)
    
    bam_files = []
    
    for i, accession in enumerate(accessions.split(';')):
        if i < len(read1_files):
            read1 = read1_files[i]
            read2 = read2_files[i]
            
            # Run fastp on paired-end reads
            trimmed_read1 = read1.replace(".fastq", ".trimmed.fastq")
            trimmed_read2 = read2.replace(".fastq", ".trimmed.fastq")
            run_fastp_paired(read1, read2, trimmed_read1, trimmed_read2, threads)
            
            # Map paired-end reads with Bowtie2
            sam_file = map_reads_paired(trimmed_read1, trimmed_read2, bowtie_index, threads, outdir, accession)
        
        elif i < len(single_end_files):
            single_end_read = single_end_files[i]
            
            # Run fastp on single-end reads
            trimmed_read = single_end_read.replace(".fastq", ".trimmed.fastq")
            run_fastp_single(single_end_read, trimmed_read, threads)
            
            # Map single-end reads with Bowtie2
            sam_file = map_reads_single(trimmed_read, bowtie_index, threads, outdir, accession)
        
        # Convert SAM to sorted BAM
        sorted_bam = sam_to_sorted_bam(sam_file, threads, outdir, accession)
        bam_files.append(sorted_bam)
    
    # Merge BAM files if there are multiple accessions
    if len(bam_files) > 1:
        merged_bam = merge_bam_files(bam_files, threads, outdir, f"{sample_name}.bam")
    else:
        merged_bam = bam_files[0]
        # Rename BAM file to sample_name
        renamed_bam = os.path.join(outdir, f"{sample_name}.bam")
        os.rename(merged_bam, renamed_bam)
        merged_bam = renamed_bam
        
        # Index the renamed BAM file
        command_index_bam = f"samtools index {merged_bam} -@ {threads}"
        run_command(command_index_bam, "indexing renamed BAM file")
    
    # Calculate breadth of coverage
    calculate_breadth_coverage(merged_bam, outdir, sample_name)
    
    # Calculate number of reads per contig
    calculate_reads_per_contig(merged_bam, outdir, sample_name)
    
    # Remove temporary files if flag is set
    if args.remove_tmp:
        remove_tmp_files(outdir, accessions)
    
    print("Processing complete.")
