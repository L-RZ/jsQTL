# This script is to classify the expression of different isoform in a gene by the genotype

# v3.9.0s_wm:
# Mann-Whitney U test (Wilcoxon rnak-sum test)
# merge test AA vs (Aa, aa)
# only test MAF > 0.05 ( edit in snp_pos_gp_dict() )
# add new format geno.gz file in geno

# v3.9.0s_wm_1:
# add p-value threshold parameter
# -p
# clean rm unused function

import os, time, sys, getopt, matplotlib, HTSeq, gzip
from scipy import stats
import numpy as np
from textwrap import wrap
matplotlib.use('Agg')
import matplotlib.pylab as plt
# from matplotlib.backends.backend_pdf import PdfPages
# import statsmodels.api as sm

def find_xm_gene_pos(gtf_in_addr, xm = 0.5):
    # use for gene_name is junction_boundary :chr_s_e_strand  1_1000_2000_+
    '''
    :param gtf_in_addr: file.gtf address str
    :param gene_name:  a set of gene names
    :param xm: expand half range * Mb
    :return: dict[gene_id] = [ChrID, start + half_range, end+range]
    '''
    f_gtf = open(gtf_in_addr)
    gene_pos_dict = {}
    for gtf_line in f_gtf:
        if gtf_line.startswith('#'):
            continue
        gtf_line_l = gtf_line.split('\t')
        if gtf_line_l[2] == 'gene':
            gene_id = gtf_line_l[8].split('"')[1] # gene_id  : ENSG00000079263.14
            chr_id = gtf_line_l[0]
            start = int(gtf_line_l[3]) - xm * 1e6
            if start < 1:
                start = 1
            end = int(gtf_line_l[4]) + xm * 1e6
            gene_pos_dict[gene_id] = [chr_id, start, end]  # [Chr(str), start(int), end(int)]

    return gene_pos_dict


def find_xm_snp_pos( chr, s, e, xm = 1):
    #use for junction_s_e is junction_boundary :chr_s_e_strand  1_1000_2000_+
    #need to fix
    s = int(s) - xm * 1e6
    if s < 0:
        s = 0
    e = int(e) + xm * 1e6
    return chr, s, e

def find_xm_snp_pos_s(chr, s, e, xm=0.5):
    # use for junction_s_e is junction_boundary :chr_s_e_strand  1_1000_2000_+
    # need to fix
    s = int(s)
    e = int(e)

    if s > e:  # strand -
        s, e = e, s
    s = s - xm * 1e6
    if s < 0:
        s = 0
    e = e + xm * 1e6
    return chr, s, e


def snp_pos_gp_dict(geno_in_addr):
    # output HTSeq genome with snp:pos and genotype
    if geno_in_addr.endswith('geno.gz'):
        f_geno = gzip.open(geno_in_addr)
    else:
        f_geno = open(geno_in_addr)
    # sample_id = f_geno.next().split()[13:]
    sample_id = f_geno.next().split()[8:]  # GTEX v7
    gas_snpid = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    geno_dict = {}
    geno_freq_dict = {}
    for each_line in f_geno:
        each_line_l = each_line.split()
        chr = each_line_l[0]
        pos = int(each_line_l[1])
        snp_id = each_line_l[2]
        # v3.9.0 for WM  MAF > 0.05
        af = float(each_line_l[8])
        maf = 0.5 - abs(0.5 - af)
        if maf < 0.05:
            continue

        gas_snpid[HTSeq.GenomicInterval(chr, pos, pos+1)] += snp_id   # HTSeq genome[pos] =  snp_id : set([snp_id])
        # geno_dict[snp_id] = np.array(each_line_l[13:], dtype=int)   # geno_dict[snp_id] = snp_genotype :np.array([0,1,2,...])
        # geno_freq_dict[snp_id] = each_line_l[7]
        # geno_dict[snp_id] = (each_line_l[7], np.array(each_line_l[13:], dtype=int))  # v6: version
        geno_dict[snp_id] = (each_line_l[8], np.array(each_line_l[10:], dtype=int))    # v7: 8th freq_A1; 10th is GT
    f_geno.close()
    # return gas_snpid, geno_dict, geno_freq_dict
    return gas_snpid, geno_dict


def genotype_zscore(in_count_sum_addr, geno_in_addr, output_dir_addr, gas_snpid, geno_dict, p_threshold):

    f_count_sum = open(in_count_sum_addr)
    count_header = f_count_sum.next()
    count_header_l = count_header.split()

    fpath, fname = os.path.split(geno_in_addr)
    fname = fname.strip('.gz')
    if not os.path.isdir(output_dir_addr):
        try:
            os.makedirs(output_dir_addr)       # creat more dirs
        except OSError:
            pass

    #for genotype zscore
    # new_fname_gp_count = '%s_gp_counts.txt' % fname
    # f_addr_gp_count = os.path.join(output_dir_addr, new_fname_gp_count)
    # f_count_out = open(f_addr_gp_count, 'w')

    #for  ttest_p
    # new_fname_ttestp = '%s_gp_glm_p.txt' % fname
    new_fname_ttestp = '%s_gp_p.txt' % fname       # 2018-9-3 edit
    f_addr_ttest_p = os.path.join(output_dir_addr, new_fname_ttestp)
    f_test_out = open(f_addr_ttest_p, 'w')

    header_test = ['GeneID', 'snpID', 'EXP_FREQ_A1', 'Nob', 'p']  #

    header_test_str = '\t'.join(header_test) + '\n'
    f_test_out.write(header_test_str)

    # for sign count and gt
    new_fname_count_gt = '%s_gp_count_sig.txt' % fname
    f_addr_count_gt_sig = os.path.join(output_dir_addr, new_fname_count_gt)
    f_count_gt_sig = open(f_addr_count_gt_sig, 'w')
    header_count = ['GeneID', 'snpID'] + count_header_l[2:]
    f_count_gt_sig.write('\t'.join(header_count) + '\n')


    # for plot mathat p
    # pdf version
    # new_fname_mathat_fig = '%s_manhat.pdf' % fname
    # fig_addr_matha =os.path.join(output_dir_addr,new_fname_mathat_fig)
    # # pdf_manhat = PdfPages(fig_addr_matha)

    fid_manhat_dir_name = '%s_manhat' % fname
    fig_manhat_dir_addr = os.path.join(output_dir_addr, fid_manhat_dir_name)
    if not os.path.isdir(fig_manhat_dir_addr):
        os.mkdir(fig_manhat_dir_addr)


    tmp_count_list = []
    for each_line in f_count_sum:
            each_line_l = each_line.strip().split()
            line_id = each_line_l[0]
            jun_id = line_id

            type_l = ['S', 'N']
            tmp_count_list.append(each_line_l[2:])  # count from 3rd column

            if len(tmp_count_list) == 2:
                tmp_count = np.array(tmp_count_list, dtype=int)
                tmp_count_list = []
                jun_id_l = line_id.split('_')
                jun_chr = jun_id_l[0]

                jun_s = jun_id_l[1]
                jun_e = jun_id_l[-2]

                test_p_result = []
                # zscore_list = []
                count_gt_list = []

                geno_chrid = gas_snpid.chrom_vectors.items()[0][0]

                range_chr, range_start, range_end = find_xm_snp_pos_s(jun_chr, jun_s, jun_e)
                if jun_chr == geno_chrid:

                    snp_iv = HTSeq.GenomicInterval(range_chr, range_start, range_end)

                    if_plot_manhat = False  # only plot significatied manhan

                    for iv, val in gas_snpid[snp_iv].steps():
                        if len(val) != 0:
                            snp_id = list(val)[0]

                            # if snp_id in '1_155033308_G_A_b37':  # test
                            #     pass    # test

                            # snp_genotype = geno_dict[snp_id]
                            # snp_freq = geno_freq_dict[snp_id]
                            snp_freq, snp_genotype = geno_dict[snp_id]

                            skip_rate = make_ratio(tmp_count[0, ], tmp_count[1,])

                            filter_ind_gt = snp_genotype != -1  # filter 1: genotype != -1
                            filter_ind_count = tmp_count[0, ] != -1  # filter 2: count    != -1
                            filter_ind_skip_rate_nonan = ~np.isnan(skip_rate)  # filter 3: N_junct / S_junct != nan
                            use_ind = filter_ind_gt * filter_ind_count * filter_ind_skip_rate_nonan  # intersection of filter1 and filter2 and filter3

                            Nob = np.sum(use_ind)
                            use_ind_snp_gt = snp_genotype[use_ind]
                            use_ind_count = tmp_count[:, use_ind]
                            use_ind_skip_rate = skip_rate[use_ind]

                            gt_ind_len_l = []
                            gt_ind_rate_l = []

                            if len(set(use_ind_skip_rate)) < 2:
                                continue

                            for i in range(3):  # i in [ 0, 1, 2]
                                gt_ind_rate_l.append(use_ind_skip_rate[use_ind_snp_gt == i])
                                gt_ind_len_l.append(np.sum(use_ind_snp_gt == i))  # test num zscore in each genotype

                            if np.sum(np.array(gt_ind_len_l) > 0) < 2:  # if there are only one genotype don't test
                                continue

                            p_mw_merge = genotype_mw_u_merge_test(gt_ind_rate_l)
                            # p_mw = genotype_mw_u_test(gt_ind_rate_l)


                            out_line = '\t'.join(map(str, [Nob, p_mw_merge]))

                            # out_line = '\t'.join(map(str, [Nob, beta, intercept, p, r_value]))
                            snp_test_p_line_str = '\t'.join([jun_id, snp_id, snp_freq, out_line])
                            test_p_result.append(snp_test_p_line_str)

                            # plot density

                            if p_mw_merge < p_threshold:

                                if_plot_manhat = True

                                gt_out_line = '\t'.join(map(str, snp_genotype))
                                out_line_str = '\t'.join([jun_id, snp_id, gt_out_line])
                                count_gt_list.append(out_line_str)

                                # f_count_out.write(
                                # '\n'.join(['\t'.join(map(str, i_line)) for i_line in gt_ind_rate_l]) + '\n')

                        else:
                            continue
                    else:
                        if len(test_p_result) == 0:
                            tmp_count = []
                            continue

                        # f_count_out.write('\n'.join(['\t'.join(map(str, i_line)) for i_line in gt_ind_rate_l]) + '\n')

                        f_test_out.write('\n'.join(test_p_result) + '\n')

                        # change to only plot manhat plot with signif
                        if if_plot_manhat:
                            # pdf_manhat = plot_p_value_v2(test_p_result, pdf_manhat)

                            # don't plot manhat
                            manhat_fig_name = '_'.join(jun_id_l[:3] + [jun_id_l[-1]]) + "_manhat.png"
                            manhat_fig_addr = os.path.join(fig_manhat_dir_addr, manhat_fig_name)
                            plot_p_value_v3_1(test_p_result, manhat_fig_addr, p_threshold)

                            f_count_gt_sig.write('\t'.join([jun_id, 'S'] + map(str, tmp_count[0])) + '\n')
                            f_count_gt_sig.write('\t'.join([jun_id, 'N'] + map(str, tmp_count[1])) + '\n')

                            f_count_gt_sig.write('\n'.join(count_gt_list) + '\n')

                else:
                    print 'No SNP in Chr%s of this chunk: %s' % (geno_chrid, fname)
                    # f_count_out.close()
                    f_test_out.close()
                    # pdf_manhat.close()
                    # pdf_denst.close()
                    # os.remove(f_addr_gp_count)
                    os.remove(f_addr_ttest_p)
                    os.remove(f_addr_count_gt_sig)
                    # os.rmdir(fig_denst_dir_addr)
                    # os.remove(fig_addr_matha)
                    # os.remove(fig_addr_denst)
                    return

    f_test_out.write(header_test_str)
    # f_count_out.close()
    f_test_out.close()
    # pdf_manhat.close()
    # pdf_denst.close()


def make_ratio(list1, list2):
    l1 = np.array(list1, dtype=float)  # skipping count
    l2 = np.array(list2, dtype=float)
    l_ratio = l1 / (l2+l1)  # empty [] in empty [] out
    return l_ratio


def make_ratio_log(list1, list2):
    l1 = np.array(list1, dtype=float)
    l1_log = np.log10(l1)
    l2 = np.array(list2, dtype=float)
    l2_log = np.log10(l2)
    l_ratio_log = l1_log - l2_log  # empty [] in empty [] out
    return l_ratio_log

def genotype_mw_u_test(genotype_ratio_list):
    """
    using Mann-Whitney U test to each 2 in 3 genotypes and return a list contains 3 elements(0vs1,1vs2,0vs2)
    :param genotype_ratio_list: [genotype0_list_ratio, genotype1_list_ratio, genotype2_list_ratio]
    :return: np.min([p_genotype0vs1, p_genotype1vs2, p_genotype1vs3])
    """
    test_list = [(0, 1), (1, 2), (0, 2)]
    # ans  = []
    # ans2 = []
    ans3 = []
    for i, j in test_list:
        x = genotype_ratio_list[i]
        y = genotype_ratio_list[j]
        if (len(x) * len(y) == 0) | (set(x) == set(y)):
            # ans.append(2)
            # ans2.append(2)
            ans3.append(2)
        else:
            # statistic, p = stats.ranksums(x, y)
            # x_y = list(x) + list(y)
            # x_y_r = stats.rankdata(x_y, method='ordinal')
            # x_r = x_y_r[:len(x)]
            # y_r = x_y_r[len(x):]
            # statistic2, p2 = stats.mannwhitneyu(x_r, y_r, alternative='two-sided')
            statistic3, p3 = stats.mannwhitneyu(x, y, alternative='two-sided')
            # ans.append(p)
            # ans2.append(p2)
            ans3.append(p3)
    # min_p = np.min(ans)
    # min_p2 = np.min(ans2)
    min_p3 = np.min(ans3)
    # return min_p, min_p2, min_p3
    return min_p3

def genotype_mw_u_merge_test(genotype_ratio_list):
    """
    using MW_test to test AA vs (Aa, aa)
    :param snp_ratio_list: [genotype0_list_ratio, genotype1_list_ratio, genotype2_list_ratio]
    :return: ks-pvalue
    """
    len_l = map(len, genotype_ratio_list)
    if max(len_l) == len_l[1]:
        if min(len_l) == len_l[2]:
            x = genotype_ratio_list[0]
            y = np.append(genotype_ratio_list[1], genotype_ratio_list[2])
        else:
            x = genotype_ratio_list[2]
            y = np.append(genotype_ratio_list[1], genotype_ratio_list[0])
    elif len_l[0] >= len_l[1]:
        x = genotype_ratio_list[0]
        y = np.append(genotype_ratio_list[1], genotype_ratio_list[2])
    else:
        x = np.append(genotype_ratio_list[0], genotype_ratio_list[1])
        y = genotype_ratio_list[2]
    mw, p_mw = stats.mannwhitneyu(x, y, alternative='two-sided')
    # if str(p) == 'nan':
    #     print len_l
    #     print genotype_ratio_list
    #     print x,y
    return p_mw


def plot_p_value_v3_1(p_result_list, manhat_addr, p_threshold):
    # plot Manhattan plot for splice junction
    # p_result: list [GeneID, snpID, Nob, beta, z, p, loglikeh, chi2, ci_0025_0975, ...]
    # p_result: list[GeneID, snpID, EXP_FREQ_A1, Nob, p]  # for np ks
    # plot the figure in png
    p_result_list = np.array([i.split('\t') for i in p_result_list])
    try:
        gene_id = p_result_list[1,0]
    except IndexError:
        print p_result_list
        sys.exit()
    snp_id_list = p_result_list[:,1]
    snp_pos = np.array([snp_id.split('_')[1] for snp_id in snp_id_list],dtype=int)
    p_array = []
    # for p in p_result_list[:, 5]:
    for p in p_result_list[:, 4]:  #  for np_ks
        if p == '': p = np.nan
        p_array.append(p)

    p_array = -np.log10(np.array(p_array, dtype=float))   # -log10(p)

    pos_max = snp_pos.max()
    pos_min = snp_pos.min()

    p_min = np.nanmin(p_array)
    p_max = np.nanmax(p_array)
    if p_max > -np.log10(p_threshold): # test   #  -log10(p) > 7
        # fig = plt.figure()

        # plt.subplot(311)
        plt.plot(snp_pos, p_array, 'k.')
        plt.title("\n".join(wrap(gene_id, 100)), fontsize=8)
        plt.ylabel('-log10(P)')
        plt.xlim(pos_min-1000, pos_max+1000)
        plt.ylim(0, p_max+1)



        plt.savefig(manhat_addr, dpi=200)
        plt.close()

    else:
        pass
        # print gene_id,'No signature'

    # return pdf_manhat


def usage():
    print """

    -i isof_exp_sum.tab     # sum each isof_type express tab
    -o output_dir           # /output_dir
    -n genotype.geno        # trunk.N.geno
    -p p-value (float)      # p-value threshold for plot
    -h                      # help
    """


def main():
    try:
        if len(sys.argv) < 4:
            usage()
            sys.exit()
        # opts, args = getopt.getopt(sys.argv[1:], 'hi:o:t:n:')
        opts, args = getopt.getopt(sys.argv[1:], 'hi:o:n:p:')
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    for opt, value in opts:
        if opt == "-i":
            isof_zscore_addr = value
        elif opt == "-o":
            output_dir_addr = value
        elif opt == "-p":
            p_threshold = float(value)
        elif opt == "-n":
            geno_in_addr = value
        elif opt == "-h":
            usage()
            sys.exit()
        else:
            usage()
            sys.exit()

    # print ' '.join(sys.argv)
    print 'START:\t%s' % time.ctime()
    # gene_pos_dict = find_xm_gene_pos(gtf_in_addr)
    # print 'gtf Done:\t%s' % time.ctime()

    snp_pos_dict, snp_gp_dict = snp_pos_gp_dict(geno_in_addr)
    print 'snp in gene Done:\t%s' % time.ctime()

    # genotype_zscore(isof_zscore_addr, geno_in_addr, output_dir_addr, gene_pos_dict, snp_pos_dict, snp_gp_dict)
    genotype_zscore(isof_zscore_addr, geno_in_addr, output_dir_addr, snp_pos_dict, snp_gp_dict, p_threshold)

    print 'All Done\t%s' % time.ctime()





if __name__ == '__main__':
    main()

