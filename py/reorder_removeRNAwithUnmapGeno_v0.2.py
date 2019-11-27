# 2017-6-20
# This script remove RNA-sample without mapping to the Geno file
# v0.2 2018-8-31
# change column index to 10 for GTEX v7
# using pandas instead of numpy
import numpy as np
import sys, time
import pandas as pd


def reorder_rna_sample(in_merge_junct_addr, geno_in_addr,out_merge_junct_addr):
    f_geno = open(geno_in_addr)
    geno_header = f_geno.next()
    f_geno.close()
    # j = 13
    j = 10   # genotype column start from 11th in table
    indiv_id_set = set([geno_id.split('-')[1] for geno_id in geno_header.split()[j:]])  # RNA order

    junct_array = np.genfromtxt(in_merge_junct_addr, dtype=str)

    i = 4   # junction count column start from 5th in table
    n = junct_array.shape[1]
    new_junct_array = junct_array[:, :4]
    unmap_id = []

    while i < n:
        junct_col = junct_array[:, i]
        junct_indiv_id = junct_col[0].split('-')[1]
        if junct_indiv_id in indiv_id_set:
            new_junct_array = np.hstack((new_junct_array, junct_col.reshape(-1,1)))
        else:
            unmap_id.append(junct_indiv_id)
        i += 1

    np.savetxt(out_merge_junct_addr, new_junct_array, fmt='%s', delimiter='\t')
    n_unmap = len(unmap_id)
    print 'No genotype RNA individual %d'%n_unmap
    print ' '.join(unmap_id)


def reorder_rna_sample_v2(in_merge_junct_addr, geno_in_addr, out_merge_junct_addr):
    f_geno = open(geno_in_addr)
    geno_header = f_geno.next()
    f_geno.close()
    # j = 13
    j = 10
    indiv_id_set = set([geno_id.split('-')[1] for geno_id in geno_header.split()[j:]])  # RNA order
    rna_df = pd.read_table(in_merge_junct_addr, sep='\t')
    i = 4
    # new_rna_df = rna_df.iloc[:, :i]
    new_rna_column = list(rna_df.columns[:i])
    junct_cols = rna_df.columns[i:]
    unmap_id = []

    n = 0  # matched sample in RNA-seq
    for junct_col in junct_cols:
        junct_indiv_id = junct_col.split('-')[1]
        if junct_indiv_id in indiv_id_set:
            new_rna_column.append(junct_col)
            # new_rna_df[junct_col] = rna_df[junct_col]
            n += 1
        else:
            unmap_id.append(junct_indiv_id)
    new_rna_df = rna_df[new_rna_column]
    new_rna_df.to_csv(out_merge_junct_addr, sep='\t', index=False)
    n_unmap = len(unmap_id)
    print 'With genotype RNA individual %d'%n
    print 'No genotype RNA individual %d'%n_unmap
    print ' '.join(unmap_id)


def main():
    print ' '.join(sys.argv)
    print 'start:',time.ctime()

    in_merge_junct_addr = sys.argv[1]
    geno_in_addr = sys.argv[2]
    out_merge_junct_addr = sys.argv[3]

    reorder_rna_sample_v2(in_merge_junct_addr, geno_in_addr, out_merge_junct_addr)
    print 'end:', time.ctime()


main()
