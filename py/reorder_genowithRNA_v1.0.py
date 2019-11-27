# 2017-5-12
# This script is used to change the order of sampleID in genotype file (chunk.1.geon) to the same with RNA sample (chr1_score)
# 2017-6-20
#    v0.2 big change match method by using RNA-seq individual ID to map DNA individual ID
#    eg. GTEX-12ZZX-0005_junction.bed : ind_ID : 12ZZX
#    the sample order of new_file.geno is based on the RNA_sample order.
#    If the RNA sample don't have genotype, its genotype will be skipping.
#    In order to match RNA with DNA,
#    other reorder_RNAwithgeno.py should run to remove unmatched RNA individual Sample before or after this script
#
# 2018-8-30
#    v0.3
#    use for GTEX v7
#    change sample GT column begin in  11th column!
# 2018-8-31
#    using pandas instand of numpy

import numpy as np
import sys, time
import pandas as pd

# def sample_map_geno(map_in_addr):
#
#     """input sample.map.txt
#     map the RNAseq sample name -> geno sample name
#     return a dict[RNAseq] = geno"""
#     f_map = open(map_in_addr)
#     map_dict = {}
#     for each_line in f_map:
#         each_line_l = each_line.strip().split('\t')
#         RNA_id = each_line_l[1].split('.')[0]
#         geno_id = each_line_l[2]
#         map_dict[RNA_id] = geno_id
#     return map_dict

def reorder_genotype_id(in_junct_merge_addr, geno_in_addr, geno_out_addr):
    f_rna = open(in_junct_merge_addr)
    rna_header = f_rna.next()
    f_rna.close()
    indiv_id_l = [rna_id.split('-')[1] for rna_id in rna_header.split('\t')[4:]]  # RNA order
    # dna_sample_id_l = [map_dict[rna_id] for rna_id in rna_sample_id_l]

    f_geno = open(geno_in_addr)
    # geno_lines = f_geno.readlines()
    # geno_header = geno_lines[0].strip().split('\t')
    # geno_array = np.array([line.strip().split() for line in f_geno])
    geno_array = np.genfromtxt(f_geno, dtype=str)
    geno_dict = {}
    # i = 13
    i = 10
    n = geno_array.shape[1]
    n_snp = geno_array.shape[0] - 1
    while i < n:
        geno_col = geno_array[:, i]
        geno_indiv_name = geno_col[0].split('-')[1]
        geno_dict[geno_indiv_name] = geno_col
        i += 1
    new_geno_array = geno_array[:, :13]
    for indiv_id in indiv_id_l:
        if indiv_id in geno_dict:
            new_geno_array = np.hstack((new_geno_array, geno_dict[indiv_id].reshape(-1,1)))
    np.savetxt(geno_out_addr, new_geno_array, fmt='%s')

def reoder_genotype_id_v2(in_junct_merge_addr, geno_in_addr, geno_out_addr):
    # pandas version
    f_rna = open(in_junct_merge_addr)
    rna_header = f_rna.next()
    f_rna.close()
    indiv_id_l = [rna_id.split('-')[1] for rna_id in rna_header.split('\t')[4:]]  # RNA order

    # df = pd.read_table(geno_in_addr, sep=' ')
    df = pd.read_table(geno_in_addr, sep='\s+')
    i = 10
    # new_df = df.iloc[:, :i]
    new_col_l = list(df.columns)[:i]
    geno_sample_id_l = list(df.columns)[i:]

    n = 0 # total match number

    for indiv_id in indiv_id_l:
        new_indiv_id = 'GTEX-%s'%indiv_id
        if new_indiv_id in geno_sample_id_l:
            new_col_l.append(new_indiv_id)
            # new_df[new_indiv_id] = df[new_indiv_id]
            n += 1
    new_df = df[new_col_l]
    new_df.to_csv(geno_out_addr, sep=' ', index=False)
    print 'Match N RNA sample: %d' % n

def main():
    print ' '.join(sys.argv)
    print 'start:', time.ctime()

    in_junct_merge_addr = sys.argv[1]
    # map_in_addr = sys.argv[2]
    geno_in_addr = sys.argv[2]
    geno_out_addr = sys.argv[3]

    # map_dict = sample_map_geno(map_in_addr)
    # reorder_genotype_id(in_junct_merge_addr, geno_in_addr, geno_out_addr)
    reoder_genotype_id_v2(in_junct_merge_addr, geno_in_addr, geno_out_addr)
    print 'end:', time.ctime()




main()