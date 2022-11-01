import numpy as np
import pandas as pd
import pysam, os
from Bio.Seq import Seq
from Bio import SeqIO

highestAbundanceGenotypes = 5
aa_list = "ACDEFGHIKLMNPQRSTVWY*"
nt_list = "ATGC"

def mods_for_export(df_genotypes_raw):
    # Count colomns for each genotype
    genotype_columns = list(df_genotypes_raw.columns)
    df_genotypes = df_genotypes_raw
    df_genotypes['genotype->seq'] = df_genotypes.groupby(by=genotype_columns[2:]).ngroup()
    df_genotypes['count'] = df_genotypes.groupby(by='genotype->seq')['seq_ID'].transform('count') # count duplicate genotypes, add to new column

     # remove duplicate genotype rows, keep the seq_ID with the highest average quality score
    df_consolidated = df_genotypes.sort_values('avg_quality_score', ascending=False).drop_duplicates('genotype->seq', keep='first')[ ['genotype->seq','count'] + genotype_columns]
    sort = ['count','NT_substitutions_count']
    ascendBool = [False,True]
    df_consolidated = df_consolidated.sort_values(sort, ascending=ascendBool)

    # move wildtype row(s) to the beginning, if they exist. rename as wild type only if there aren't any barcodes in the genotype, as this would result in many different 'wildtype' rows
    wildtypeDF = df_consolidated.loc[(df_consolidated['NT_substitutions']=='')&(df_consolidated['NT_insertions']=='')&(df_consolidated['NT_deletions']=='')]
    if len(wildtypeDF) > 0:
        wildtype_in_df = True
    else: wildtype_in_df = False
    if wildtype_in_df:
        df_consolidated = df_consolidated.drop(index=wildtypeDF.index)
        df_consolidated = pd.concat([wildtypeDF, df_consolidated]).reset_index(drop=True)
        df_consolidated.rename(index={0:'wildtype'}, inplace=True)
    else:
        df_consolidated.reset_index(drop=True, inplace=True)
    if df_consolidated.index[0] == 0: # make barcode IDs 1-indexed if necessary
        df_consolidated.index += 1

    # now that genotype IDs are established, add column that correlates every sequence ID with a genotype ID from the condensed genotypes DF
    df_consolidated = df_consolidated.reset_index().rename(columns={'index':'genotype_ID'})
    seq_to_genotype_dict = dict(zip(df_consolidated['genotype->seq'], df_consolidated['genotype_ID']))
    df_genotypes['genotype_ID'] = df_genotypes['genotype->seq'].map(seq_to_genotype_dict)

    # iterate through x genotypes with highest counts and genotypes of specific ID # (both defined in config file) , get a representative sequence for each (that w highest avg_quality_score, or essentially random if there are no quality scores), and write alignments to file
    df_alignments = df_consolidated.iloc[0:highestAbundanceGenotypes+1]

    return df_genotypes, df_consolidated, df_alignments

def export_alignment(export_dir, bam_dir,df_alignments, alignment_data):
    alignments_dir = os.path.join(export_dir,'alignments.txt')

    ref = alignment_data[0]
    alignStr = alignment_data[1]
    seq = alignment_data[2]

    with open(alignments_dir, 'w') as txtOut:
        with pysam.AlignmentFile(bam_dir, 'rb') as bam_file:
            nameIndexedBAM = pysam.IndexedReads(bam_file)
            nameIndexedBAM.build()
            for row in df_alignments.itertuples():
                if row.genotype_ID=='wildtype':
                    continue
                seqID = row.seq_ID
                iterator = nameIndexedBAM.find(seqID)
                for BAMentry in iterator:
                    break
                txtOut.write(f'Genotype {row.genotype_ID} representative sequence. Sequence ID: {seqID}\n')
                for string in [ref, alignStr, seq]:
                    txtOut.write(string+'\n')
                txtOut.write('\n')
            txtOut.write('')

def export_genotypes(export_dir, target_sequence, NTmutArray, NTmutDist, genotypesDF, genotypesDFcondensed):
    genotypes_path = os.path.join(export_dir,'genotypes.csv')
    seqid_path = os.path.join(export_dir,'seq-IDs.csv')
    nt_muts_freq_dir = os.path.join(export_dir,'NT-muts-frequencies.csv')
    nt_muts_distribution_dir = os.path.join(export_dir,'NT-muts-distribution.csv')

    ntIDs = list(target_sequence)
    ntPositions = [f'{str(i)}' for i in range(0, len(target_sequence))]
    WTnts = [ID+ntPosi for ID,ntPosi in zip(ntIDs,ntPositions)]
    NTmutDF = pd.DataFrame(NTmutArray, columns=list(nt_list))
    NTmutDF['wt_nucleotides'] = pd.Series(WTnts)
    NTmutDF.set_index('wt_nucleotides', inplace=True)
    NTmutDF = NTmutDF.transpose()

    totalSeqs = int(NTmutDist.sum())

    NTdistDF = pd.DataFrame(NTmutDist, columns=['seqs_with_n_NTsubstitutions'])
    NTdistDF.index.name = 'n'

    genotypesDFcondensed.drop(columns=['genotype->seq', 'seq_ID', 'avg_quality_score']).to_csv(genotypes_path, index=False)
    genotypesDF.drop(columns=genotypesDF.columns.difference(['seq_ID', 'genotype_ID'])).to_csv(seqid_path, index=False)

    NTmutDF.index.name = 'NT_mutation_count'
    if totalSeqs>0:
        NTmutDF = NTmutDF.divide(totalSeqs)
    NTmutDF.to_csv(nt_muts_freq_dir)
    NTdistDF.to_csv(nt_muts_distribution_dir)

    return totalSeqs

def export_failures(export_dir, failuresDF):
    path = os.path.join(export_dir,'failures.csv')
    failuresDF.to_csv(path, index=False)

def export_aa_stats(export_dir, target_ref, AAmutArray, AAmutDist, totalSeqs):

    target_sequence = str(list(target_ref.sequence)[0])
    target_length = int(len(target_sequence))

    aa_muts_freq_dir = os.path.join(export_dir,'AA-muts-frequencies.csv')
    aa_muts_distribution_dir = os.path.join(export_dir,'AA-muts-distribution.csv')

    resiIDs = list(str(Seq(target_sequence).translate()))
    resiPositions = [str(i) for i in range(1, int((target_length/3)+1) )]

    WTresis = [ID+posi for ID,posi in zip(resiIDs,resiPositions)]
    AAmutDF = pd.DataFrame(AAmutArray, columns=list(aa_list))
    AAmutDF['wt_residues'] = pd.Series(WTresis)
    AAmutDF.set_index('wt_residues', inplace=True)
    AAmutDF = AAmutDF.transpose()
    AAmutDF.index.name = 'AA_mutation_count'
    if totalSeqs > 0:
        AAmutDF = AAmutDF.divide(totalSeqs)
    AAmutDF.to_csv(aa_muts_freq_dir)

    AAdistDF = pd.DataFrame(AAmutDist, columns=['seqs_with_n_AAsubstitutions'])
    AAdistDF.index.name = 'n'
    AAdistDF.to_csv(aa_muts_distribution_dir)


if __name__ == '__main__':
    print('Output scripts for mini-maple sequencing pipeline. Run snakemake for full pipeline.')