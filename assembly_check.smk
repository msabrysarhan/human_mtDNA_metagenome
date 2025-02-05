
import pandas as pd
import os


sample_list = "NCBIaccession.txt"
sample_df = pd.read_csv(sample_list, sep='\t')
sample_df[['study', 'sample']] = sample_df['sampleID'].str.split('__', expand=True)


studies = sample_df['study'].tolist()
studies_list = sample_df['study'].unique().tolist()
samples = sample_df['sample'].tolist()


rule all:
    input:
        stats = ["summary/all_samples_rCRS_hg.tsv", "summary/all_samples_rCRS_stats.tsv", "summary/all__rCRS_fasta.depth"],

rule download_assemblies:
    params: 
        study_name = "{study}"
    input: 
        download_script = "download_metagenomic_assemblies.sh"
    output: 
        bz_assembly = temp("assemblies/{study}.tar.bz2")
    threads: 1
    resources: mem_mb=10000
    shell:
        """
        mkdir -p assemblies
        grep {params.study_name} {input.download_script} > assemblies/download_{wildcards.study}.sh
        bash assemblies/download_{wildcards.study}.sh
        rm assemblies/download_{wildcards.study}.sh
        """


rule extract_assemblies:
    input: 
        bz_assembly = "assemblies/{study}.tar.bz2",
    output: 
        assemblies = directory("contigs/{study}"),
    resources: mem_mb=20000
    shell:
        """
        mkdir -p {output.assemblies}
        tar -xvjf {input.bz_assembly} -C contigs/
        """



rule blast_assemblies:
    input:
        assembly = "contigs/{study}",
        mtDNA_ref = "reference/rCRS.nhr"
    params:
        assembly = "contigs/{study}/{study}__{sample}.fa"
    output:
        blast_fasta = "results/{study}/{sample}_rCRS.fasta",
    threads: 1
    resources: mem_mb=10000
    shell:
        """
        python bin/get_mito.py -c {params.assembly} \
        -d reference/rCRS -o {output.blast_fasta} -b 500 -s 99 --remove_temp
        """



rule rCRS_download:
    params:
        mtDNA_ref = "https://www.phylotree.org/resources/rCRS.fasta"
    output:
        mtDNA_ref = "reference/rCRS.fasta"
    threads: 1
    resources: mem_mb=10000
    shell:
        """
        wget -O {output.mtDNA_ref} {params.mtDNA_ref}
        """

rule rCRS_index:
    input:
        mtDNA_ref = "reference/rCRS.fasta"
    params:
        mtDNA_ref = "reference/rCRS"
    output:
        mtDNA_ref = "reference/rCRS.nhr"
    threads: 1
    resources: mem_mb=10000
    shell:
        """
        makeblastdb -in {input.mtDNA_ref} -out {params.mtDNA_ref} -dbtype nucl
        bwa index {input.mtDNA_ref} -p {params.mtDNA_ref}
        """


rule mito2vcf:
    input:
        fasta_files= "results/{study}/{sample}_rCRS.fasta",
        mtDNA_ref = "reference/rCRS.nhr",
    params:
        ref_dir = "reference/rCRS",
        quality_threshold = 30,
        max_depth = 1000
    output:
        bam = "results/{study}/{sample}_rCRS.bam",
        vcf = "results/{study}/{sample}_rCRS.vcf",
        fasta = "results/{study}/{sample}_rCRS.consensus.fasta"
    shell:
        """
        bwa mem {params.ref_dir} {input.fasta_files} | \
        samtools view -Sbq {params.quality_threshold} | \
        samtools sort > {output.bam}
        
        samtools index {output.bam}
        samtools consensus {output.bam} > {output.fasta}
        
        bcftools mpileup -B -Ou -d {params.max_depth} -q {params.quality_threshold} -f {params.ref_dir}.fasta {output.bam} | \
        bcftools call -mv --ploidy 1 -Ou -o {output.vcf}
        """

rule haplogrep3:
    input:
        vcf = "results/{study}/{sample}_rCRS.vcf",
    output:
        haplogrep = "results/{study}/{sample}_rCRS.hg",
    shell:
        """
        haplogrep3 classify --in {input.vcf} --out {output.haplogrep} --tree phylotree-rcrs@17.2 --chip
        """

rule combined_coverage:
    input:
        fasta = [f"results/{study}/{sample}_rCRS.fasta" for study, sample in zip(studies, samples)],
        mtDNA_ref = "reference/rCRS.nhr",
    params:
        ref_dir = "reference/rCRS",
    output:
        all_fasta = temp("summary/all__rCRS_fasta.fasta"),
        bam = temp("summary/all__rCRS_fasta.bam"),
        depth = "summary/all__rCRS_fasta.depth"
    shell:
        """
        cat {input.fasta} > {output.all_fasta}

        bwa mem {params.ref_dir} {output.all_fasta} | samtools view -Sbq 30 | samtools sort > {output.bam}
        samtools depth -a {output.bam} > {output.depth}
        """

rule summary:
    input:
        fasta = [f"results/{study}/{sample}_rCRS.fasta" for study, sample in zip(studies, samples)],
        hg = [f"results/{study}/{sample}_rCRS.hg" for study, sample in zip(studies, samples)],
    output:
        stats = "summary/all_samples_rCRS_stats.tsv",
        hg = "summary/all_samples_rCRS_hg.tsv",
    shell:
        """
        seqkit stats -T -a  {input.fasta} > {output.stats}
        cat {input.hg} | grep -v "Haplogroup" > {output.hg}
        """
