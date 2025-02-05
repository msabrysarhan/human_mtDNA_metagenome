#!/usr/bin/env python3
"""
Script Name: extract_matches.py
Description: This script blast fasta sequences against pre-built BLAST database.
Author: Mohamed S. Sarhan
Affiliation: Institute for Biomedicine, Eurac Research, Bolzano 39100, Italy
Contact: mohamed.sarhan@eurac.edu; m.sabrysarhan@gmail.com
Date Created: Feb 16, 2024
Version: 1.0
"""

import argparse
import subprocess
import os
from Bio import SeqIO
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def check_file_exists(file_path):
    """Check if the file exists."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

def check_blast_database(database_path):
    """Check if the BLAST database exists and is usable."""
    try:
        subprocess.run(['blastdbcmd', '-db', database_path, '-info'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except subprocess.CalledProcessError:
        raise FileNotFoundError(f"BLAST database not found or inaccessible: {database_path}")

def count_contigs(contigs_file):
    """Count the number of contigs in the input file."""
    return sum(1 for _ in SeqIO.parse(contigs_file, "fasta"))

def run_blast(contigs_file, blast_database, blast_output_file, threads):
    """Run BLAST search."""
    logger.info("Running BLAST search...")
    subprocess.run(['blastn', '-db', blast_database, '-query', contigs_file, '-out', blast_output_file, '-outfmt', '6', '-num_threads', str(threads)], check=True)

def filter_blast_results(blast_output_file, filtered_output_file, similarity, matches, bit_score):
    """Filter BLAST results based on provided criteria and write to a new file."""
    logger.info("Filtering BLAST results...")
    with open(blast_output_file, 'r') as blast_handle, open(filtered_output_file, 'w') as filtered_handle:
        for line in blast_handle:
            columns = line.strip().split('\t')
            contig_id, hit_id, percent_identity, alignment_length, mismatches, gap_openings, q_start, q_end, s_start, s_end, evalue, hit_bit_score = columns
            if float(percent_identity) >= similarity and float(alignment_length) >= matches and float(hit_bit_score) >= bit_score:
                filtered_handle.write(line)

def write_output(filtered_output_file, contigs_file, output_file):
    """Write filtered records to output FASTA file."""
    logger.info("Writing output...")
    contigs = []
    contig_sequences = {record.id: (str(record.seq), record.description) for record in SeqIO.parse(contigs_file, "fasta")}
    with open(filtered_output_file, "r") as filtered_handle, open(output_file, "w") as output_handle:
        for line in filtered_handle:
            columns = line.strip().split('\t')
            contig_id = columns[0]
            if contig_id not in contigs:
                contigs.append(contig_id)
                contig_seq, contig_desc = contig_sequences[contig_id]
                output_handle.write(f">{contig_id} {contig_desc}\n{contig_seq}\n")

def main(args):
    # Check if input files exist
    check_file_exists(args.contigs)
    check_blast_database(args.database)

    # Count the number of contigs
    contig_count = count_contigs(args.contigs)
    logger.info(f"Number of contigs in {args.contigs}: {contig_count}")

    # Run BLAST search
    blast_output_file = args.output + ".tab"
    run_blast(args.contigs, args.database, blast_output_file, args.threads)

    # Filter BLAST results
    filtered_output_file = args.output + ".filtered.tab"
    filter_blast_results(blast_output_file, filtered_output_file, args.similarity, args.matches, args.bit_score)
    logger.info(f"Filtered BLAST results written to: {filtered_output_file}")

    # Write output
    write_output(filtered_output_file, args.contigs, args.output)

    # Remove temporary files if requested
    if args.remove_temp:
        os.remove(blast_output_file)
        os.remove(filtered_output_file)
        logger.info("Removing temporary files")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process genome data.')
    parser.add_argument('--contigs', '-c', type=str, help='Path to contigs file', required=True)
    parser.add_argument('--database', '-d', type=str, help='Path to BLAST database', required=True)
    parser.add_argument('--output', '-o', type=str, help='Prefix for output files', required=True)
    parser.add_argument('--threads', '-t', type=int, default=1, help='Number of threads')
    parser.add_argument('--similarity', '-s', type=float, default=50, help='Minimum percentage identity (default: 50)')
    parser.add_argument('--matches', '-m', type=float, default=20, help='Minimum match length (default: 20) as a float')
    parser.add_argument('--bit_score', '-b', type=float, default=50, help='Minimum bit score (default: 50) as a float')
    parser.add_argument('--remove_temp', action='store_true', help='Remove temporary files after processing')
    args = parser.parse_args()
    main(args)
