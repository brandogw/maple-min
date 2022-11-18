import common
import bam_to_genotypes as bg
import numpy as np
import pandas as pd
import csv_outputs, os, argparse, shutil

aa_list = "ACDEFGHIKLMNPQRSTVWY*"
nt_list = "ATGC"

def main(save_dir, ref_dir,bam_dir,target_name, qsminimum):
    # Make directory for output files
    if os.path.isdir(save_dir):
        shutil.rmtree(save_dir)
    os.mkdir(save_dir)

    # Load references from FASTA file 
    references = common.load_fasta(ref_dir)
    references['target_pos_in_completeRead'] = common.find_target_indices(references) # find position of target in complete read
    
    # Load BAM file
    df_bam = common.load_bam(bam_dir)

    # Identify target sequence from FASTA file
    complete_read = references[(references.identifiers == 'completeRead')] # name of fasta entry with complete read
    target_ref = references[(references.identifiers == target_name)] # name of fasta entry with target read, this can be flexible. 

    # Initialize array for general stats about mutations
    target_sequence = str(list(target_ref.sequence)[0]) # target seuence
    target_length = int(len(target_sequence)) # target nucleotide length
    protein_length = int( len(target_sequence) / 3 ) # length of protein
    nt_mutations_by_position = np.zeros((target_length, len(nt_list)), dtype=int) # array that counts mutations in each nt pos
    nt_total_mutations = np.zeros(target_length, dtype=int) # tallies number of mutations
    aa_mutations_by_position = np.zeros((protein_length, len(aa_list)), dtype=int) # array that counts mutations in each aa pos
    aa_total_mutations = np.zeros(protein_length, dtype=int) # tallies number of nonsyn aa mutations

    # Initialize arrays for final csv outputs          
    bam_query_names = []
    average_q_scores = []
    genotypes =[] # genotype list
    mutations_in_genotype =[] # number of mutations in genotype
    in_genotypes =[] # insertion genotypes found
    del_genotypes =[] # deletion genotypes found
    aa_nonsynonymous_muts =[] # nonsynonymous amino acid mutations
    aa_synonymous_muts =[] # synonymous amino acid mutations
    aa_nonsynonymous_in_genotype=[] # number of nonsynonymous mutations in genotype
    failure_list = []

    # Go through each BAM entry
    for bam_index in range(len(df_bam.index)):
        bam_entry = df_bam.loc[bam_index] #load based on location in df

        # Check for common alignment failure reasons before analysis
        failure_reason, pass_bam_check = bg.bam_ready_check(target_ref, complete_read, bam_entry) # check ref matches, start of protein, end of protein is in read
        alignment_data, insertions, deletions = bg.analyze_bam(bam_entry, target_ref, complete_read)

        if pass_bam_check and bam_entry['cigartuples']:
            # Analyze BAM, create alignment files, and label genotypes
            codons_with_indels, in_genotype, del_genotype = bg.indel_to_genotype(target_ref, insertions, deletions) 

            # id mutations flagged by BAM analysis
            (genotype, 
            aa_nonsynonymous_mut, 
            aa_synonymous_mut, 
            nt_mutation_by_position, 
            aa_mutation_by_position) = bg.id_mutations(target_ref, alignment_data, codons_with_indels, qsminimum)
            
            # Add general mutations stats from bam entry
            total_mutations_in_entry = sum(sum(nt_mutation_by_position)) 
            nt_mutations_by_position += nt_mutation_by_position
            nt_total_mutations[total_mutations_in_entry] += 1

            total_aa_in_entry = sum(sum(aa_mutation_by_position))
            aa_mutations_by_position += aa_mutation_by_position
            aa_total_mutations[total_aa_in_entry] += 1

            # Add genotypes from bam entry into array for entry into pandas dataframe
            bam_query_names.append(bam_entry['query_name'])
            average_q_scores.append(np.average(np.array(bam_entry['query_alignment_qualities'])))
            genotypes.append(', '.join(genotype))
            mutations_in_genotype.append(len(genotype))
            in_genotypes.append(in_genotype)
            del_genotypes.append(del_genotype)
            aa_nonsynonymous_muts.append(', '.join(aa_nonsynonymous_mut))
            aa_synonymous_muts.append(', '.join(aa_synonymous_mut))
            aa_nonsynonymous_in_genotype.append(len(aa_nonsynonymous_mut))
        else:
            failure_list.append([bam_entry['query_name'], failure_reason])
            csv_outputs.export_failure_alignments(save_dir, alignment_data, bam_entry['query_name'], failure_reason)


    # Pandas dataframe headers
    d = {'seq_ID':bam_query_names, 
        'avg_quality_score':average_q_scores, 
        'NT_substitutions':genotypes, 
        'NT_substitutions_count':mutations_in_genotype, 
        'NT_insertions':in_genotypes, 
        'NT_deletions':del_genotypes,
        'AA_substitutions_nonsynonymous':aa_nonsynonymous_muts,
        'AA_substitutions_synonymous':aa_synonymous_muts,
        'AA_substitutions_nonsynonymous_count':aa_nonsynonymous_in_genotype}

    # Create pandas dataframe
    df_genotypes = pd.DataFrame(data=d)
    df_failures = pd.DataFrame(failure_list, columns=['seq_ID', 'failure_reason'])
    df_genotypes, df_consolidated, df_alignments = csv_outputs.mods_for_export(df_genotypes)
    
    csv_outputs.export_alignment(save_dir, bam_dir, df_alignments, alignment_data)
    total_sequences = csv_outputs.export_genotypes(save_dir, target_sequence, nt_mutations_by_position, nt_total_mutations, df_genotypes, df_consolidated)
    csv_outputs.export_failures(save_dir,df_failures)
    csv_outputs.export_aa_stats(save_dir, target_ref, aa_mutations_by_position, aa_total_mutations, total_sequences)


        
if __name__ == '__main__':

    # Set command line interface variables
    parser = argparse.ArgumentParser(description='Extract mutational data from BAM files.')
    parser.add_argument('-r','--reference', help='Path to FASTA file', required=True)
    parser.add_argument('-b','--bam', help='Path to BAM file', required=True)
    parser.add_argument('-q','--qualityscore', help='Minimum quality score of the reads', default=25, required=False)
    parser.add_argument('-s','--savedir', help='Directory to output CSV files', required=True)
    parser.add_argument('-t','--targetseq', help='Target sequence to compare completeRead with', required=True)
    args = vars(parser.parse_args())

    ref_dir = args['reference']
    bam_dir = args['bam']
    qsminimum = args['qualityscore']
    save_dir = args['savedir']
    target_name = args['targetseq']

    main(save_dir, ref_dir,bam_dir,target_name, qsminimum)