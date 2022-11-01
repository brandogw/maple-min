
rule get_fastq:
    input:
        copy_from = get_fastq
    output:
        copy_to = "data/processed/reads/{sample}/{sample}_R{direction}.fastq.gz"
    shell:
        """
        cp {input.copy_from} {output.copy_to} 
        """


rule merge_paired_end:
    input:
        fwd = "data/processed/reads/{sample}/{sample}_R1.fastq.gz",
        rvs = "data/processed/reads/{sample}/{sample}_R2.fastq.gz" 
    output:
        merged = "data/processed/merged/{sample}.merged.fastq.gz",
    params:
        flags = '-m 10',
    shell:
        """
        NGmerge -1 {input.fwd} -2 {input.rvs} -o {output.merged} -z {params.flags}

        """


rule minimap2:
    input:
        sequence = "data/processed/merged/{sample}.merged.fastq.gz",
        ref = config['ref']
    output:
        aligned = "data/processed/aligned/{sample}.sam"
    params:
        flags = '-a -A2 -B4 -O4 -E2 --secondary=no'
    threads: 3
    resources:
        threads = 3,
        mem_mb = 10000,
        time_min = 60
    shell:
        """
        minimap2 -t {threads} {params.flags} {input.ref} {input.sequence} > {output.aligned} 

        """

# sam to bam conversion
rule sam2bam:
    input:
        sam = "data/processed/aligned/{sample}.sam"
    output:
        bam = "data/processed/bam/{sample}.bam",
        bai = "data/processed/bam/{sample}.bam.bai"
    shadow: "minimal"
    threads: 1
    resources:
        threads = 3,
        mem_mb = 10000
    shell:
        """
        samtools view -b {input.sam} | samtools sort -m 4G > {output.bam}
        samtools index {output.bam}
        """
