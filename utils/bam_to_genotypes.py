import pysam, argparse
from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np
import pandas as pd

aa_list = "ACDEFGHIKLMNPQRSTVWY*"
nt_list = "ATGC"

# Analyze BAM entry and produce sequences for categorizing mutations
# INPUTS:
#     - bam_entry- panads entry with cigartuples, reference_start, query sequence, and quality fields
#     - target_ref - pandas entry with reference the BAM is compared against
# OUTPUTS:
#     - None if any of the variables are not present
#     OR
#     - refAln   CTGGTATCGCCAGGCGCCGGGCAAAGAAC (from complete read)
#     - alignStr ||||||||||||||||.|||||||||||| 
#     - queryAln CTGGTATCGCCAGGCGGCGGGCAAAGAAC (from bam entry)
#     - qScores - quality scores for query alignment 
#     - insertions - list of insertions found for bam entry
#     - deletions - list of deletions found for bam entry

def analyze_bam(bam_entry, target_ref, complete_read): # Pandas Entry

    completeRead_seq = str(complete_read['sequence'].to_numpy()[0])

    starting_pos_in_read = target_ref['target_pos_in_completeRead'].to_numpy()[0][0]
    ending_pos_in_read = target_ref['target_pos_in_completeRead'].to_numpy()[0][1]

    cigartuples = bam_entry['cigartuples'] #summary of results from alignment
    start = bam_entry['reference_start'] # where the BAM entry starts, this will serve as an index as more cigartuples get analyzed
    query_sequence = bam_entry['query_alignment_sequence'] #alignment sequence
    query_quality = bam_entry['query_alignment_qualities'] # quality of reads

    # Initialize varaibles to analyze cigartuples
    query_index = 0 # index for sequence position for query sequence (from BAM)
    ref_index = start # index for sequence position for reference (complete read) sequence (from fasta)
    
    # Initialization for strings for construction of alignment data
    refAln = ' ' * start  # variable for construction of reference sequence
    alignStr = ' ' * start # variable for construction of string showing relationship between strains
    queryAln = ' ' * start # variable for construction of query sequence
    queryQualities = [-1] * start 

    insertions = []
    deletions = [] 
     
    # Loop through each logged alignment difference
    for cigar_entry in cigartuples:
        operation = cigar_entry[0]
        length = cigar_entry[1]
        
        if operation == 0: #no indel
            from_reference = completeRead_seq[ref_index:ref_index + length]
            from_bam = query_sequence[query_index:query_index+length]

            aligned_segment = ''.join(['|' if r==q else '.' for r,q in zip(from_reference,from_bam)])
            
            refAln += from_reference
            queryAln += from_bam
            queryQualities += query_quality[query_index:query_index+length]

            alignStr += aligned_segment
            ref_index += length
            query_index += length

        if operation == 1: #insertion
            if (starting_pos_in_read <= ref_index) and (ending_pos_in_read > ref_index):
                insertions.append((ref_index - starting_pos_in_read, query_sequence[query_index:query_index+length]))
            query_index += length
        
        if operation == 2: #deletion, '-' added to sequence to maintain alignment to reference       
            refAln += completeRead_seq[ref_index:ref_index + length]
            queryAln += '-' * length
            alignStr += ' ' * length
            queryQualities += [0] * length

            if ( starting_pos_in_read <= ref_index + length ) and ( ref_index < ending_pos_in_read ): # record deletions as tuples of position and length
                deletions.append((ref_index-starting_pos_in_read, length))
            ref_index += length

    ref = refAln[starting_pos_in_read:ending_pos_in_read]
    alignStr = alignStr[starting_pos_in_read:ending_pos_in_read]
    seq = queryAln[starting_pos_in_read:ending_pos_in_read]
    qScores = queryQualities[starting_pos_in_read:ending_pos_in_read]

    alignment_data = [ref, alignStr, seq, qScores]
    insertions = insertions
    deletions = deletions

    return (alignment_data, insertions, deletions)

def bam_ready_check(target_ref, complete_read, bam_entry):
    bam_pass =True
    failure_reason = ""
    failure_index = None
    bam_reference_name = bam_entry.reference_name
    bam_start = bam_entry['reference_start']
    bam_end = bam_entry['reference_end']
    bam_length = bam_end-bam_start

    id = str(list(complete_read.identifiers)[0]) 
    starting_pos_in_read = target_ref['target_pos_in_completeRead'].to_numpy()[0][0]
    ending_pos_in_read = target_ref['target_pos_in_completeRead'].to_numpy()[0][1]
    template_length = ending_pos_in_read-starting_pos_in_read

    if bam_reference_name != id:
        failure_reason = 'alignment uses wrong reference sequence'
        bam_pass = False
    
    if bam_start > starting_pos_in_read:
        failure_reason = 'alignment starts past trimmed reference start'
        bam_pass = False

    if bam_end < ending_pos_in_read:
        failure_reason = 'alignment ends before trimmed reference end'
        bam_pass = False
    
    if bam_length < template_length:
        failure_reason = 'truncated read'
        bam_pass = False

    return failure_reason, bam_pass



def indel_to_genotype(target_ref, insertions, deletions):
    starting_pos_in_read = target_ref['target_pos_in_completeRead'].to_numpy()[0][0]
    ending_pos_in_read = target_ref['target_pos_in_completeRead'].to_numpy()[0][1]
    refProteinStart = 0
    refProteinEnd = ending_pos_in_read - starting_pos_in_read
    flagged_codons = []    # list of amino acid positions that are affected by indel (for insertion, insertion is within a codon; for deletion, at least one base of codon deleted)
    for index, _ in insertions:
        if refProteinStart <= index < refProteinEnd:
            protIndex = index-refProteinStart
            if protIndex%3 == 0: continue # ignore if insertion occurs between codons
            else: flagged_codons.append( int(protIndex/3) )
    for index, length in deletions:
        if (refProteinStart <= index < refProteinEnd) or (refProteinStart <= index+length < refProteinEnd):
            protIndexStart = index-refProteinStart
            protIndexEnd = (index+length)-refProteinStart
            firstCodon = int(protIndexStart/3)
            lastCodon = int(protIndexEnd/3)
            flagged_codons.extend([i for i in range(firstCodon,lastCodon+1)])

    in_output = ', '.join([str(index)+'ins'+nt_list for index,nt_list in insertions]) 
    out_output = ', '.join([str(index)+'del'+str(length) for index,length in deletions])     # string of all deletions for genotype output

    return flagged_codons, in_output, out_output


# Analyze alignment data and look for nucleotide and amino acid mutations
# INPUTS:
#     - alignment_data = [refAln, alignStr, queryAln, qScores]
#                - refAln   CTGGTATCGCCAGGCGCCGGGCAAAGAAC (from complete read)
#                - alignStr ||||||||||||||||.|||||||||||| 
#                - queryAln CTGGTATCGCCAGGCGGCGGGCAAAGAAC (from bam entry)
#     - protein_indicies: indicies in the truncated sequence (not complete read) that contains the protein. In most cases, should be [0 protein_length].
#     - flagged_codons: codons with indels to make sure they get handled properly
#     - QSminimum: minimium quality score to log as a genotype
# OUTPUTS:

def id_mutations(target_ref, alignment_data, flagged_codons, QSminimum=25):
    # Initialize variables
    ## Protein indicies
    starting_pos_in_read = target_ref['target_pos_in_completeRead'].to_numpy()[0][0]
    ending_pos_in_read = target_ref['target_pos_in_completeRead'].to_numpy()[0][1]

    refProteinStart = 0
    refProteinEnd = ending_pos_in_read - starting_pos_in_read

    ## Alignment Data
    ref = alignment_data[0]
    alignStr = alignment_data[1]
    seq = alignment_data[2]
    qScores = alignment_data[3]

    ## Mutation Logs
    NTmutArray = np.zeros((int(len(alignStr)), len(nt_list)), dtype=int)
    AAmutArray = np.zeros((int(len(alignStr)/3), len(aa_list)), dtype=int)
    codonsChecked = []
    NTsubstitutions = []
    AAnonsynonymous = []
    AAsynonymous = []

    # Find mismatches in alignment data
    mismatches = [i for i,a in enumerate(alignStr) if a=='.']
    for i in mismatches:

        if qScores[i] < QSminimum: 
            continue
        else:
            # Analyze NT mutations
            wtNT = ref[i] # original base pair
            mutNT = seq[i] # mutated base pair
            NTmutArray[i,nt_list.find(mutNT)] += 1 # Log mutation in array
            NTsubstitutions.append(wtNT+str(i+1)+mutNT) # genotype output 1-index

            # Analyze AA mutations
            if refProteinStart <= i < refProteinEnd:
                protIndex = i-refProteinStart
                codon = int(protIndex/3) # 0-index amino acid position
                if codon in codonsChecked: continue
                codonsChecked.append(codon)
                codonPosi = protIndex%3 
                codonIndices = list(range(i-codonPosi, i+(3-codonPosi)))

                # check that codon doesn't contain any bases influenced by an indel
                if codon in flagged_codons: continue

                #check that all three quality scores in codon are above threshold
                if qScores:
                    QStooLow = False
                    codonQS = qScores[codonIndices[0]:codonIndices[2]]
                    for qs in codonQS:
                        if qs < QSminimum:
                            QStooLow = True
                    if QStooLow: continue

                # nucleotide to AA
                wtAA = str(Seq(ref[codonIndices[0]:codonIndices[2]+1]).translate())
                mutAA = str(Seq(seq[codonIndices[0]:codonIndices[2]+1]).translate())

                if wtAA!=mutAA:
                    AAmutArray[codon, aa_list.find(mutAA)] += 1
                    AAnonsynonymous.append(wtAA+str(codon+1)+mutAA) # genotype output 1-index
                else:
                    AAsynonymous.append(wtAA+str(codon+1))
    
    return (NTsubstitutions, AAnonsynonymous, AAsynonymous, NTmutArray, AAmutArray)

# def genotypes_to_pandas():

#     genotypesDF = pd.DataFrame(genotypesList, columns=genotypesColumns)


if __name__ == '__main__':
    print("mini-MAPLE BAM analysis tools. Please run mutation_analysis.py instead.")