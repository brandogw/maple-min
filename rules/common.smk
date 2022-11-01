import snakemake.common
import pandas as pd
import subprocess

###### Config file and sample sheets #####
configfile: "config/config.yaml"
samples = pd.read_csv(config["samples"]).set_index("Sample", drop=False)
project = config['project']
basemount = config['basemount']

def get_fastq(wildcards):
    sample = str(wildcards.sample)
    direction = str(wildcards.direction)
    subfolder_tag = project+"{:03d}".format(int(wildcards.sample))
    fastq_dir = basemount+'Basespace/Projects/'+ project +'/Samples/'+ subfolder_tag +'/Files/'+subfolder_tag + '_S'+ sample + '_L001_R'+direction+'_001.fastq.gz'
    print("Getting files from: " + fastq_dir)
    return fastq_dir