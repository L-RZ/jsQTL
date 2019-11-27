# to find the skipping exon junction read and sum them up
#  exon1---exon2---exon3---exon4---
# 1st count the junction reads(exon1_exon2)  as non-skipping-exon-junction-reads : NSEJ
# 2st sum and count the jucnton reads ( exon1_exon3 + exon1_exon4 + ....exon1_exonN) : sum-skipping-exon-junction-reads: SSEJ
# junction count file should sorted by start pos and then end pos for "+" strand gene
# v0.1 without gtf annotation
#      with same start recursion junction test
# v1.1 skipping finding strage is changing. The group2 junction skipping is changing
#      add filter: output first sum of all sample the junction read > 10 as the ns_jun
#
#      NO add the gtf as reference to make sure the first junction is the annotated
# 2017-9-7
# v1.2 find the exon junction reads with the gtf annotation
# should be
# v1.2.1 junction should be in the same gene

# v1.3.0 to find the "-" strand junction which is different from the "+" strand

# >>--|exon1|----|exon2|----|exon3|----|exon4|---->>>
#           <--1->   #group1       NS
#           <-----2--------->      SS
#           <--------3----------------->   SS
#
#                      <--1->  #group2     NS
#                      <---2----------->   SS
# chr1 pos  1     3    5    7     9    10    11
#
#junct_id_1: 1_1_3_7_10_+|gene_id
#junct_id_2: 1_5_7_10_+|gene_id
#
#

# <<--|exon4|----|exon3|----|exon2|----|exon1|----<<<
#                                 <--1->   #group1  NS
#                      <-----2--------->            SS
#           <--------3----------------->            SS
#
#                      <--1->  #group2              NS
#           <---2----------->                       SS
# chr1 pos  1     3    5    7     9    10    11
#
#junct_id_1: 1_10_9_5_1_-|gene_id
#junct_id_2: 1_7_5_1_-|gene_id

# junction_extract_merge
#strand +
# chr     start   end     strand  GTEX-111YS-0006_junction        GTEX-1122O-0005_junction        GTEX-1128S-0005_junction
# 1       12227   12497   +       0       0       0
# 1       12227   12595   +       2       1       1
# 1       12227   12613   +       9       2       6
# 1       12227   12646   +       0       0       0
# 1       12456   12613   +       0       0       0
# 1       12697   13221   +       0       0       0
# 1       12697   13225   +       0       0       0
# 1       12697   13403   +       1       2       0
#strand -
# Chr     start   end     strand  GTEX-111YS-0006_junction        GTEX-1122O-0005_junction        GTEX-1128S-0005_junction
# 1       249224225       249227187       -       0       0       0
# 1       249153897       249157157       -       0       0       1
# 1       249152520       249157157       -       0       0       0
# 1       249153453       249153715       -       0       0       0
# 1       249152520       249153715       -       0       0       0




import numpy as np
import sys
def make_gtf_exon_dict(gtf_addr):
    f_gtf = open(gtf_addr)
    gtf_exon_dict = {}
    gtf_exon_s_dict = {}
    gtf_exon_e_dict = {}
    for line in f_gtf:
        if not line.startswith('#'):
            line_l = line.split('\t')
            if line_l[2] == 'exon':
                chr_id = line_l[0]
                exon_s = line_l[3]
                exon_e = line_l[4]
                exon_strand = line_l[6]
                exon_param = line_l[8].strip()
                gene_id = exon_param.split(';')[-1].split()[1][1:-1]  # changed gene_name
                if chr_id not in gtf_exon_dict:
                    gtf_exon_dict[chr_id] = {}
                    gtf_exon_s_dict[chr_id] = {}
                    gtf_exon_e_dict[chr_id] = {}
                key = '_'.join([exon_s, exon_e, exon_strand])
                gtf_exon_dict[chr_id][key] = [gene_id]
                exon_s_strand = '_'.join([exon_s, exon_strand])
                if exon_s_strand not in gtf_exon_s_dict[chr_id]:
                    gtf_exon_s_dict[chr_id][exon_s_strand] = []
                if gene_id not in gtf_exon_s_dict[chr_id][exon_s_strand]:
                    gtf_exon_s_dict[chr_id][exon_s_strand].append(gene_id)

                exon_e_strand = '_'.join([exon_e, exon_strand])
                if exon_e_strand not in gtf_exon_e_dict[chr_id]:
                    gtf_exon_e_dict[chr_id][exon_e_strand] = []
                if gene_id not in gtf_exon_e_dict[chr_id][exon_e_strand]:
                    gtf_exon_e_dict[chr_id][exon_e_strand].append(gene_id)

    return gtf_exon_dict, gtf_exon_s_dict, gtf_exon_e_dict

def format_output(jun_list_tmp):
    ns_jun_l = jun_list_tmp[0][:4]
    ns_jun_chr = ns_jun_l[0]
    ns_jun_s = ns_jun_l[1]
    ns_jun_e = ns_jun_l[2]
    ns_jun_str = '_'.join(ns_jun_l)
    ns_exp = jun_list_tmp[0][4:]
    ss_jun_list_exp = []
    jun_e_list = [ns_jun_e]
    for ss_each_line_l in jun_list_tmp[1:]:
        ss_jun_l = ss_each_line_l[:4]
        ss_jun_chr = ss_each_line_l[0]
        ss_jun_s = ss_each_line_l[1]
        ss_jun_e = ss_each_line_l[2]
        ss_jun_str = '_'.join(ss_jun_l)
        if ns_jun_chr == ss_jun_chr and ns_jun_s == ss_jun_s:
            jun_e_list.append(ss_jun_e)
            ss_jun_list_exp.append(ss_each_line_l[4:])
    id_out = '%s_%s_%s_%s' % (ns_jun_chr, ns_jun_s, '_'.join(jun_e_list), ns_jun_l[3])
    ss_array_sum = np.sum(np.array(ss_jun_list_exp, dtype=int), axis=0)
    ss_jun_exp_out_str = '\t'.join(map(str, ss_array_sum.tolist()))
    ns_exp_str = '\t'.join(ns_exp)
    ss_line_out = '\t'.join([id_out, 'S', ss_jun_exp_out_str])
    ns_line_out = '\t'.join([id_out, 'N', ns_exp_str])

    return ss_line_out, ns_line_out

def format_output_v2(jun_list_tmp, gtf_exon_s_dict, gtf_exon_e_dict):
    # output first sum of all sample the junction read > 10 as the ns_jun

    found_ns = False
    ss_jun_list_exp = []

    for each_line_l in jun_list_tmp:
        if not found_ns:                # find the first junction  as ns non-skipping
            ns_exp = each_line_l[4:]    # sum all this all sample junc read counts
            ns_exp_sum = np.array(ns_exp, dtype=int).sum()
            if ns_exp_sum >= 10:
                found_ns = True
                # ns_jun_l = jun_list_tmp[0][:4]    #wrong
                ns_jun_l = each_line_l[:4]
                ns_jun_chr = ns_jun_l[0]
                ns_jun_s = ns_jun_l[1]
                ns_jun_e = ns_jun_l[2]
                ns_jun_strand = ns_jun_l[3]

                jun_e_list = [ns_jun_e]
                ns_jun_s_key = '_'.join([ns_jun_s, ns_jun_strand])
                ns_jun_e_key = '_'.join([ns_jun_e, ns_jun_strand])

                ns_jun_s_gene_id_list = gtf_exon_s_dict[ns_jun_chr][ns_jun_e_key]
                ns_jun_e_gene_id_list = gtf_exon_e_dict[ns_jun_chr][ns_jun_s_key]
                ns_gene_id_set = set(ns_jun_s_gene_id_list) & set(ns_jun_e_gene_id_list)
                if not ns_gene_id_set:
                    break
                ns_gene_id = '|'.join(list(ns_gene_id_set))
        else:
            ss_each_line_l = each_line_l  # find  remained and sum up as ss sum-skipping
            ss_jun_l = ss_each_line_l[:4]
            ss_jun_chr = ss_each_line_l[0]
            ss_jun_s = ss_each_line_l[1]
            ss_jun_e = ss_each_line_l[2]
            ss_jun_str = '_'.join(ss_jun_l)
            ss_jun_strand = ss_each_line_l[3]

            ss_jun_e_key = '_'.join([ss_jun_e, ss_jun_strand])
            ss_jun_e_gene_id_list = gtf_exon_s_dict[ns_jun_chr][ss_jun_e_key]
            ss_gene_id_set = set(ss_jun_e_gene_id_list)

            if not ss_gene_id_set & ns_gene_id_set:
                continue

            if ns_jun_chr == ss_jun_chr and ns_jun_s == ss_jun_s:
                jun_e_list.append(ss_jun_e)
                ss_jun_list_exp.append(ss_each_line_l[4:])
    if len(ss_jun_list_exp) == 0:
        ss_line_out = ''
        ns_line_out = ''
        return ss_line_out, ns_line_out

    try:
        id_out = '%s_%s_%s_%s|%s' % (ns_jun_chr, ns_jun_s, '_'.join(jun_e_list), ns_jun_l[3], ns_gene_id)
        ss_jun_list_exp_array = np.array(ss_jun_list_exp, dtype=int)     # if there is no ns, it will be NameError.
        if len(ss_jun_list_exp_array.shape) == 1:                        # only one ss junction
            ss_array_sum = ss_jun_list_exp_array
        else:
            try:
                ss_array_sum = np.sum(ss_jun_list_exp_array, dtype=int, axis=0)
            except TypeError:   # for test
                print ss_jun_list_exp_array
                sys.exit()

        ss_jun_exp_out_str = '\t'.join(map(str, ss_array_sum.tolist()))  # if there is no ss, it will be NameError.
        ns_exp_str = '\t'.join(ns_exp)
        ss_line_out = '\t'.join([id_out, 'S', ss_jun_exp_out_str])
        ns_line_out = '\t'.join([id_out, 'N', ns_exp_str])
    except NameError:
        ss_line_out = ''
        ns_line_out = ''

    return ss_line_out, ns_line_out

def format_output_v3_posi(jun_list_tmp, gtf_exon_s_dict, gtf_exon_e_dict):
    # output first sum of all sample the junction read > 10 as the ns_jun
    found_ns = False
    ss_jun_list_exp = []
    for each_line_l in jun_list_tmp:
        if not found_ns:                # find the first junction  as ns non-skipping
            ns_exp = each_line_l[4:]    # sum all this all sample junc read counts
            ns_exp_sum = np.array(ns_exp, dtype=int).sum()
            if ns_exp_sum >= 10:
                found_ns = True
                # ns_jun_l = jun_list_tmp[0][:4]    #wrong
                ns_jun_l = each_line_l[:4]
                ns_jun_chr = ns_jun_l[0]
                ns_jun_s = ns_jun_l[1]
                ns_jun_e = ns_jun_l[2]
                ns_jun_strand = ns_jun_l[3]

                jun_e_list = [ns_jun_e]
                ns_jun_s_key = '_'.join([ns_jun_s, ns_jun_strand])
                ns_jun_e_key = '_'.join([ns_jun_e, ns_jun_strand])

                ns_jun_s_gene_id_list = gtf_exon_s_dict[ns_jun_chr][ns_jun_e_key]
                ns_jun_e_gene_id_list = gtf_exon_e_dict[ns_jun_chr][ns_jun_s_key]
                ns_gene_id_set = set(ns_jun_s_gene_id_list) & set(ns_jun_e_gene_id_list)
                if not ns_gene_id_set:
                    break
                ns_gene_id = '|'.join(list(ns_gene_id_set))
        else:
            ss_each_line_l = each_line_l  # find  remained and sum up as ss sum-skipping
            ss_jun_l = ss_each_line_l[:4]
            ss_jun_chr = ss_each_line_l[0]
            ss_jun_s = ss_each_line_l[1]
            ss_jun_e = ss_each_line_l[2]
            ss_jun_str = '_'.join(ss_jun_l)
            ss_jun_strand = ss_each_line_l[3]

            ss_jun_e_key = '_'.join([ss_jun_e, ss_jun_strand])
            ss_jun_e_gene_id_list = gtf_exon_s_dict[ns_jun_chr][ss_jun_e_key]
            ss_gene_id_set = set(ss_jun_e_gene_id_list)

            if not ss_gene_id_set & ns_gene_id_set:
                continue

            if ns_jun_chr == ss_jun_chr and ns_jun_s == ss_jun_s:
                jun_e_list.append(ss_jun_e)
                ss_jun_list_exp.append(ss_each_line_l[4:])
    if len(ss_jun_list_exp) == 0:
        ss_line_out = ''
        ns_line_out = ''
        return ss_line_out, ns_line_out

    try:
        id_out = '%s_%s_%s_%s|%s' % (ns_jun_chr, ns_jun_s, '_'.join(jun_e_list), ns_jun_l[3], ns_gene_id)
        ss_jun_list_exp_array = np.array(ss_jun_list_exp, dtype=int)     # if there is no ns, it will be NameError.
        if len(ss_jun_list_exp_array.shape) == 1:                        # only one ss junction
            ss_array_sum = ss_jun_list_exp_array
        else:
            try:
                ss_array_sum = np.sum(ss_jun_list_exp_array, dtype=int, axis=0)
            except TypeError:   # for test
                print ss_jun_list_exp_array
                sys.exit()

        ss_jun_exp_out_str = '\t'.join(map(str, ss_array_sum.tolist()))  # if there is no ss, it will be NameError.
        ns_exp_str = '\t'.join(ns_exp)
        ss_line_out = '\t'.join([id_out, 'S', ss_jun_exp_out_str])
        ns_line_out = '\t'.join([id_out, 'N', ns_exp_str])
    except NameError:
        ss_line_out = ''
        ns_line_out = ''

    return ss_line_out, ns_line_out

def format_output_v3_negi(jun_list_tmp, gtf_exon_s_dict, gtf_exon_e_dict):
    # output first sum of all sample the junction read > 10 as the ns_jun
    found_ns = False
    ss_jun_list_exp = []
    for each_line_l in jun_list_tmp:
        if not found_ns:                # find the first junction  as ns non-skipping
            ns_exp = each_line_l[4:]    # sum all this all sample junc read counts
            ns_exp_sum = np.array(ns_exp, dtype=int).sum()
            if ns_exp_sum >= 10:
                found_ns = True
                # ns_jun_l = jun_list_tmp[0][:4]    #wrong
                ns_jun_l = each_line_l[:4]
                ns_jun_chr = ns_jun_l[0]    # junct             [1]       [2]
                ns_jun_s = ns_jun_l[2]      # for strand - :     E|<-----<|S              change direction
                ns_jun_e = ns_jun_l[1]      # for exon:  S--exon2_E       S_exon1---E     in gtf unchange direction
                ns_jun_strand = ns_jun_l[3]

                jun_e_list = [ns_jun_e]
                ns_jun_s_key = '_'.join([ns_jun_s, ns_jun_strand])
                ns_jun_e_key = '_'.join([ns_jun_e, ns_jun_strand])

                ns_jun_s_gene_id_list = gtf_exon_s_dict[ns_jun_chr][ns_jun_s_key]
                ns_jun_e_gene_id_list = gtf_exon_e_dict[ns_jun_chr][ns_jun_e_key]
                ns_gene_id_set = set(ns_jun_s_gene_id_list) & set(ns_jun_e_gene_id_list)
                if not ns_gene_id_set:
                    break
                ns_gene_id = '|'.join(list(ns_gene_id_set))
        else:
            ss_each_line_l = each_line_l  # find  remained and sum up as ss sum-skipping
            ss_jun_l = ss_each_line_l[:4]
            ss_jun_chr = ss_each_line_l[0]
            ss_jun_s = ss_each_line_l[2]   #  change
            ss_jun_e = ss_each_line_l[1]   #  change
            ss_jun_str = '_'.join(ss_jun_l)
            ss_jun_strand = ss_each_line_l[3]

            ss_jun_e_key = '_'.join([ss_jun_e, ss_jun_strand])
            ss_jun_e_gene_id_list = gtf_exon_e_dict[ss_jun_chr][ss_jun_e_key]
            ss_gene_id_set = set(ss_jun_e_gene_id_list)

            if not ss_gene_id_set & ns_gene_id_set:
                continue

            if ns_jun_chr == ss_jun_chr and ns_jun_s == ss_jun_s:
                jun_e_list.append(ss_jun_e)
                ss_jun_list_exp.append(ss_each_line_l[4:])
    if len(ss_jun_list_exp) == 0:
        ss_line_out = ''
        ns_line_out = ''
        return ss_line_out, ns_line_out

    try:
        id_out = '%s_%s_%s_%s|%s' % (ns_jun_chr, ns_jun_s, '_'.join(jun_e_list), ns_jun_l[3], ns_gene_id)
        ss_jun_list_exp_array = np.array(ss_jun_list_exp, dtype=int)     # if there is no ns, it will be NameError.
        if len(ss_jun_list_exp_array.shape) == 1:                        # only one ss junction
            ss_array_sum = ss_jun_list_exp_array
        else:
            try:
                ss_array_sum = np.sum(ss_jun_list_exp_array, dtype=int, axis=0)
            except TypeError:   # for test
                print ss_jun_list_exp_array
                sys.exit()

        ss_jun_exp_out_str = '\t'.join(map(str, ss_array_sum.tolist()))  # if there is no ss, it will be NameError.
        ns_exp_str = '\t'.join(ns_exp)
        ss_line_out = '\t'.join([id_out, 'S', ss_jun_exp_out_str])
        ns_line_out = '\t'.join([id_out, 'N', ns_exp_str])
    except NameError:
        ss_line_out = ''
        ns_line_out = ''

    return ss_line_out, ns_line_out

def find_junction(in_addr, out_addr, gtf_exon_s_dict, gtf_exon_e_dict):
    f_in = open(in_addr)
    dat_in = f_in.readlines()
    head_line = dat_in[0].strip().split('\t')
    sample_id = '\t'.join(head_line[4:])
    ns_jun_index = 1  # header is 0
    result = ['id\tisoform_type\t%s' % sample_id]
    # raw_result = ['id\tSKIP_type\tisoform_type\tjunction_type\t%s' % sample_id]


    jun_s_tmp = 0
    jun_chr_tmp_posi = 0
    jun_list_tmp_posi = []

    jun_e_tmp = 0
    jun_chr_tmp_negi = 0
    jun_list_tmp_negi = []

    for each_line in dat_in[1:]:
        each_line_l = each_line.strip().split('\t')
        jun_chr = each_line_l[0]
        jun_s = each_line_l[1]
        jun_e = each_line_l[2]
        jun_strand = each_line_l[3]

        jun_s_key = '_'.join([jun_s, jun_strand])
        jun_e_key = '_'.join([jun_e, jun_strand])

        if not jun_s_key in gtf_exon_e_dict[jun_chr] or not jun_e_key in gtf_exon_s_dict[jun_chr]:
            continue

        if jun_strand in '+':
            if jun_s_tmp == 0 and jun_chr_tmp_posi == 0:
                jun_list_tmp_posi.append(each_line_l)
                jun_s_tmp = jun_s
                jun_chr_tmp_posi = jun_chr
            elif jun_s_tmp == jun_s and jun_chr_tmp_posi == jun_chr:
                jun_list_tmp_posi.append(each_line_l)
            else:
                if len(jun_list_tmp_posi) != 1:
                    # ss_out, ns_out = format_output(jun_list_tmp)
                    ss_out, ns_out = format_output_v3_posi(jun_list_tmp_posi, gtf_exon_s_dict, gtf_exon_e_dict)
                    if ss_out != '':
                        result.append(ss_out)
                        result.append(ns_out)

                jun_s_tmp = jun_s
                jun_chr_tmp_posi = jun_chr
                jun_list_tmp_posi = [each_line_l]
        else:
            if jun_e_tmp == 0 and jun_chr_tmp_negi == 0:
                jun_list_tmp_negi.append(each_line_l)
                jun_e_tmp = jun_e
                jun_chr_tmp_negi = jun_chr
            elif jun_e_tmp == jun_e and jun_chr_tmp_negi == jun_chr:
                jun_list_tmp_negi.append(each_line_l)
            else:
                if len(jun_list_tmp_negi) != 1:
                    # ss_out, ns_out = format_output(jun_list_tmp)
                    ss_out, ns_out = format_output_v3_negi(jun_list_tmp_negi, gtf_exon_s_dict, gtf_exon_e_dict)
                    if ss_out != '':
                        result.append(ss_out)
                        result.append(ns_out)
                jun_e_tmp = jun_e
                jun_chr_tmp_negi = jun_chr
                jun_list_tmp_negi = [each_line_l]


    else:
        if len(jun_list_tmp_posi) != 1:
            ss_out, ns_out = format_output_v3_posi(jun_list_tmp_posi, gtf_exon_s_dict, gtf_exon_e_dict)
            if ss_out != '':
                result.append(ss_out)
                result.append(ns_out)

        if len(jun_list_tmp_negi) !=1:
            ss_out, ns_out = format_output_v3_negi(jun_list_tmp_negi, gtf_exon_s_dict, gtf_exon_e_dict)
            if ss_out != '':
                result.append(ss_out)
                result.append(ns_out)
  
    f_out_res = open(out_addr + '_sumLJ.txt','w')
    f_out_res.write('\n'.join(result))
    f_out_res.close()

'''
output in start with recursion 
id	isoform_type	GTEX-1122O-0003_junction	GTEX-11EM3-0001_junction	GTEX-11EMC-0002_junction
2_231109795_231110480_231110578_231112631_231112637_231113600_231120167_231134247_+	S	146	56	167
2_231109795_231110480_231110578_231112631_231112637_231113600_231120167_231134247_+	N	0	0	0
2_231110655_231111964_231112631_231112637_231113600_231120167_+	S	271	34	272
2_231110655_231111964_231112631_231112637_231113600_231120167_+	N	4	0	3
'''

def main():
    print ' '.join(sys.argv)
    in_jun_addr = sys.argv[1] # .rmUnMapID.txt
    gtf_addr = sys.argv[2]    # file.gtf
    output_addr = sys.argv[3] 

    gene_dict, exon_s_dict, exon_e_dict = make_gtf_exon_dict(gtf_addr)
    find_junction(in_jun_addr, output_addr, exon_s_dict, exon_e_dict)
main()

