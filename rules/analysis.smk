rule mutation_analysis:
    input:
        bam = "data/processed/bam/{sample}.bam",
        ref = "data/ref/ref.fasta",
    output:
        files= directory("data/processed/mutation_stats/{sample}")
    params:
        target = config["target"]
    shell:
        """
        python3 utils/mutation_analysis.py -r {input.ref} -b {input.bam} -t {params.target} -s {output.files} 
        """