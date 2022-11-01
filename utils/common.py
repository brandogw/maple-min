import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
import pysam
from pathlib import Path

# Loads fasta file into pandas dataframe
def load_fasta(ref_dir):
    # Open FASTA file with all reference sequences
    with open(ref_dir) as fasta_file:  # Will close handle cleanly
        identifiers = []
        seqs = []
        rev_comp = [] 
        length = []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            identifiers.append(seq_record.id)
            seqs.append(str(seq_record.seq.upper()))
            rev_comp.append(str(seq_record.seq.reverse_complement().upper()))
            length.append(len(seq_record.seq))
    d = {'sequence': seqs,'reverse_compliment':rev_comp, 'length': length, 'identifiers': identifiers}
    references = pd.DataFrame(data=d)
    return references

# Loads BAM file into pandas dataframe
def load_bam(bam_dir):
    with pysam.AlignmentFile(bam_dir, 'rb') as bam_file:
        query_name = []
        reference_start = []
        reference_end = [] 
        reference_name = []
        cigartuples = []
        query_alignment_sequence = [] 
        query_alignment_qualities = []

        for entry in bam_file:
            query_name.append(entry.query_name)
            reference_start.append(entry.reference_start)
            reference_end.append(entry.reference_end)
            reference_name.append(entry.reference_name)
            cigartuples.append(entry.cigartuples)
            query_alignment_sequence.append(entry.query_alignment_sequence)
            query_alignment_qualities.append(entry.query_alignment_qualities)

    d = {'query_name': query_name,
        'reference_start':reference_start,
        'reference_end':reference_end,
        'reference_name':reference_name,
        'cigartuples':cigartuples,
        'query_alignment_sequence':query_alignment_sequence,
        'query_alignment_qualities':query_alignment_qualities
        }
    df_bam = pd.DataFrame(data=d)
    return df_bam

# Finds sequence in complete read and outputs protein indices
def find_target_indices(references):
    position =[]
    complete_read = str(references[(references.identifiers == 'completeRead')]['sequence'].to_numpy()[0])
    for n in references.index:
        target = str(references['sequence'][n])
        starting_position =complete_read.find(target)
        ending_position = starting_position + len(target)
        position.append([starting_position, ending_position])
    return position




if __name__ == '__main__':
    print('Common scripts for mini-maple sequencing pipeline. Run snakemake for full pipeline.')