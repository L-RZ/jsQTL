# this scritp is use to find the 95% SNP in each exon-skipping event
# using p-value -> -log(likehood) -> -log(detla_likelihood)
# p-valvue = chi(likehood, df =1)
# v0.1.1 2017-10-13
# v0.1.2 2018-4-3
# add paramater p_cutoff as input
#
# 2018-07-26 find a but in proportion_p function line 14:
# Error      When p-value > 1, the p-value is not append to the junct_p list. This will make the length of junct_p is
#            less than junct_list, which leads to line 43 get error result.
# v0.1.4 fixed v0.1.2 bug.
#
# v0.2.0 big changes 2019-5-23
#        add R^2 as filter to fine-mapping in credible set
#        add Index SNP & Index SNP p-value & Index Pos


import os, sys, time
import numpy as np
import scipy.stats as stats
import argparse


def make_r2_dict(in_addr):
    f_in_r2 = open(in_addr)
    r2_dict = {}
    header = f_in_r2.next()
    for line in f_in_r2:
        line_l = line.strip().split()
        snpA_id = line_l[2]
        snpB_id = line_l[6]
        r2 = float(line_l[9])
        if snpA_id in r2_dict:
            r2_dict[snpA_id][snpB_id] = r2
        else:
            r2_dict[snpA_id] = {snpB_id: r2}
    return r2_dict

def proportion_p(junct_list, min_line_l, r2_dict, p_cutoff, r2_cutoff=0.4):
    # junct_p = np.array([l.strip().split('\t')[4] for l in junct_list if float(l.strip().split('\t')[4]) <= 1], dtype=float) v0.1.2
    # junct_p = np.array([l.strip().split('\t')[5] for l in junct_list], dtype=float)

    # r2_cutoff = 0.4 # set

    min_p_snpid = min_line_l[1]
    min_p = min_line_l[4]

    if min_p_snpid not in r2_dict:
        print 'ERROR:'
        print min_line_l
        return
    r2_dict_sub = r2_dict[min_p_snpid]
    junct_p = []
    junct_list_sub = []

    for l in junct_list:
        l_l = l.strip().split('\t')
        l_snpid = l_l[1]
        p = float(l_l[4])
        if l_snpid in r2_dict_sub:
            r2 = r2_dict_sub[l_snpid]
            if r2 > r2_cutoff:
                junct_list_sub.append(l.strip()+'\t'+str(r2))
                if float(p) > 1:
                    p = 1
                junct_p.append(p)

    junct_p = np.array(junct_p, dtype=float)
    sub_n = len(junct_p)
    if np.nanmin(junct_p) < p_cutoff:
        base_line_p = np.nanmax(junct_p)
        base_line_chi2 = stats.chi2.isf(base_line_p, 1)
        chi2_array = stats.chi2.isf(junct_p, 1)  # chi2 = 2log(M1/M0)
        ld_rate = np.exp((chi2_array - base_line_chi2) / 2.0)
        ld_rate_sum = np.nansum(ld_rate)

        ld_proportion = ld_rate / ld_rate_sum
        ld_proportion_index = ld_proportion
        ld_proportion_index[ld_proportion_index == np.nan] = 0
        ld_array = np.vstack([ld_proportion, ld_proportion_index])
        ld_array_arg = np.argsort(-ld_array[1, :])   # sort max -> min
        ld_array = ld_array[:, ld_array_arg]
        ld_array_cum = np.cumsum(ld_array[1, :])
        cum_s = ld_array_cum / ld_array_cum[-1]

        # if junct_list[0].split()[0] == '1_111731834_111731408_111730992_-|DENND2D':
        #     print(junct_list[0])

        ld_prop_cutoff = np.nanmax(ld_array[0, cum_s > 0.95])

        is_ld_prop_95 = ld_proportion >= ld_prop_cutoff
        n_is_ld_prop_95 = np.sum(is_ld_prop_95)

        index_snpid = min_p_snpid
        index_p = min_p
        index_prob = np.exp((stats.chi2.isf(float(min_p), 1) - base_line_chi2) / 2.0) / ld_rate_sum

        res = ['\t'.join([i[0].strip(), str(n_is_ld_prop_95), str(i[1]), index_snpid, str(index_p), str(index_prob)])
               for i in zip(junct_list_sub, ld_proportion, is_ld_prop_95)
               if i[2]]
        res_line = '\n'.join(res) + '\n'
        return res_line


def min_p_test(new_p, new_line, min_p, min_line):
    if new_p < min_p:
        return new_p, new_line
    else:
        return min_p, min_line


def junct_region(p_addr, out_addr, r2_dict, p_cutoff=1e-9, r2_cutoff=0.4):
    f_p = open(p_addr)
    fpath, fname = os.path.split(p_addr)
    name, ext = os.path.splitext(fname)

    f_out = open(out_addr, 'w')

    junct_id_old = ''
    junct_tmp_list = []
    for line in f_p:
        line_l = line.split('\t')
        if line_l[0] == 'GeneID':
            # header = line.strip() + '\tN_credible\tP_credible\n'
            header = line.strip() + '\tR2\tN_credible\tP_credible\tSNP_index\tP_index\tP_cred_index\n'
            f_out.write(header)
            break

    f_p = open(p_addr)
    min_line_l = []
    min_p = 1
    for line in f_p:
        line_l = line.strip().split('\t')
        junct_id = line_l[0]
        if line_l[0] != 'GeneID':
            line_p = float(line_l[4])
            if junct_id_old == '':
                junct_id_old = junct_id
                junct_tmp_list = [line]

                min_p, min_line_l = min_p_test(line_p, line_l, min_p, min_line_l)

            elif junct_id == junct_id_old:
                junct_tmp_list.append(line)
                min_p, min_line_l = min_p_test(line_p, line_l, min_p, min_line_l)
            else:
                if min_p < p_cutoff:  # run FM if minp of exon-skipping is < p_cutoff
                    out_line = proportion_p(junct_tmp_list, min_line_l, r2_dict, p_cutoff, r2_cutoff)
                    f_out.write(out_line)
                junct_tmp_list = [line]
                junct_id_old = junct_id
                min_p = line_p
                min_line_l = line_l

        # else:# for test
        #     print junct_tmp_list[-1] #for test
    else:
        if len(junct_tmp_list) != 0 and min_p < p_cutoff:
            out_line = proportion_p(junct_tmp_list, min_line_l, r2_dict, p_cutoff, r2_cutoff)
            if out_line:
                f_out.write(out_line)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('in_p_addr', help='input association result')
    parser.add_argument('in_ld_addr', help='input R^2 ld file')
    parser.add_argument('out_addr', help='output fine-mapping output')
    parser.add_argument('--R2', help='R2 cutoff (FM > R2). defaut: 0.4', default=0.4, type=float)
    parser.add_argument('--p_cutoff', help='p cutoff for index SNP. defaut: 1e-9', default=1e-9, type=float)
    args = parser.parse_args()

    # print ' '.join(sys.argv)
    # in_addr = sys.argv[1]
    # out_addr = sys.argv[2]
    # p_cutoff = float((sys.argv[3]))

    in_addr = args.in_p_addr
    in_ld_addr = args.in_ld_addr
    out_addr = args.out_addr
    r2_cutoff = args.R2
    p_cutoff = args.p_cutoff

    r2_dict = make_r2_dict(in_ld_addr)

    junct_region(in_addr, out_addr, r2_dict, p_cutoff, r2_cutoff)
    print 'DONE', time.ctime()
    # print 'DONE', time.ctime(), sys.argv[2]
main()




