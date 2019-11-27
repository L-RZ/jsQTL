# calculate the z-score for each individual by two type isoform's junction reads
# 20170522 filter the min reads counts: 25% sample's reads counts >= 1
#          added postfix
# 20170626 change the filter method
# 20170705 change the filter method to left in out
#          individual sum count must > 5
#          remain at least 3 individual

#20170818  using the chi2 test when max nmber of test_array is > 1e7 instand of the fisher test
#          to reduce the memory and running time.
#20171007  add # at the end of the output file

import scipy.stats as stats
import numpy as np
import os,sys,time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
# from matplotlib.backends.backend_pdf import PdfPages


def filter(remain_expr_array):
    # using left in out method
    # input : remian_expr_array
    #  [[ 1 3 4 4 4 4 ],    S_junction skips
    #   [ 2 4 5 5 5 5 ]]    N_junction NonSkip
    # output
    #  min_p : float
    #  p_result_list_srt: list [ p_str,..]
    #  test_array_list  : list [ '[test_array]', ...]

    new_expr_array = np.array(remain_expr_array, dtype=float) # to avoid /0
    n_ind = new_expr_array.shape[1]
    p_result = []
    test_array_list = []
    for i in xrange(n_ind):
        # is_test = np.arange(n_ind) == i
        is_other= np.arange(n_ind) != i
        test_ind_array = new_expr_array[:,i]
        other_ind_array = np.sum(new_expr_array[:, is_other], axis=1)   # sum two rows

        test_array = np.array([test_ind_array, other_ind_array], dtype=int)

        if np.max(test_array) > 1e7:                              # change to chi2 test   20170818
            try:
                chi, p, df, exp = stats.chi2_contingency(test_array)  # 20181127 avoid error "frequencies has a zero element at %s." % (zeropos,))
            except ValueError:
                odd, p = stats.fisher_exact(test_array)
        else:
            odd, p = stats.fisher_exact(test_array)
        p_result.append(p)

        test_array_list.append(str(test_array.tolist()))
    # p_result = np.array(p_result) * n_ind  #  correct p-value
    # p_result[ p_result >1 ] = 1
    p_result = 1 - (1 - np.array(p_result)) ** n_ind #  correct p-value
    min_p = p_result.min()
    p_result_list_srt = map(str, p_result.tolist())
    return min_p, p_result_list_srt, test_array_list

def filter_fun(in_addr, postfix):

    '''
    input_file format: 
    LongJunt_ShortJunc1_ShortJunc2  isofrom_type(N:normal S:skipping) sample1_counts sample2_counts sample3_counts
    1_29206_30976_+_1_29206_30564_+_1_30667_30976_+_M	N	0.0	0.0	0.0
    1_29206_30976_+_1_29206_30564_+_1_30667_30976_+_M	S	0.0	0.0	0.0
    '''

    f_in = open(in_addr)
    path, ext = os.path.splitext(in_addr)
    dat = f_in.readlines()
    tmp = []
    # result_exp = [dat[0]]
    exp_header = dat[0]

    # result_parameter = ['\t'.join(['id', 'min_p_val', 'n_remian', 'p_val','[[test_S,test_N],[other_S,other_N]\n'])]
    parameter_header = '\t'.join(['id', 'min_p_val', 'n_remian', 'p_val', '[[test_S,test_N],[other_S,other_N]\n'])
    ### change to output on time

    out_addr_exp = '%s_%s.txt' % (path, postfix)
    out_addr_para = '%s_%s_para.txt' % (path, postfix)

    f_out = open(out_addr_exp, 'w')
    f_out_para = open(out_addr_para, 'w')

    f_out.write(exp_header)
    f_out_para.write(parameter_header)


    # header_l =  dat[0].split()
    # new_header = '\t'.join([header_l[0]] + header_l[2:])
    # result = [new_header]

    # n_sample = len(header_l[2:])
    # n_sample_filter = 0.25 * n_sample   # 25% sample


    for each_line in dat[1:]:
        each_line_l = each_line.strip().split()
        line_id = each_line_l[0]

        type_l = ['S', 'N']
        tmp.append(each_line_l[2:])
        if len(tmp) == 2:
            #  filtered low juntion reads counts events
            # n_expr = np.sum(np.array(tmp, dtype=float) > 0, axis=1)
            # if (n_expr > n_sample_filter).all():


            # if line_id=='22_18604468_18606772_18606928_18606938_18606965_18606970_18609145_+':    #outlier is NA
            # if line_id in '22_17585700_17586481_17586743_17589197_+_0.0620568992935_83_0.00468140442133':   # outlier is error
            # if line_id in '22_22930674_23006934_23101221_+':    # for test
            #     print  line_id


            tmp_array = np.array(tmp, dtype=int)  # expression array

            ind_expr_sum = np.sum(tmp_array, axis=0)
            ind_expr_sum_filter = ind_expr_sum > 5  # individual sum count must > 5
            ind_expr_sum_filter_out = ind_expr_sum <= 5
            if np.sum(ind_expr_sum_filter) < 3:  # remain at least 3 individual
                tmp = []
                continue

            remain_expr_array = tmp_array[:, ind_expr_sum_filter]

            n_remian = np.sum(ind_expr_sum_filter)

            # n_expr = np.sum(remain_expr_array > 0, axis=1)     # at last 3 individual count > 0 at two isform
            # if (n_expr < 3).any():
            #     tmp = []
            #     continue

            tmp_array[:, ind_expr_sum_filter_out] = -1           # change filtered value to -1

            min_p, p_result_list_srt, test_array_list = filter(remain_expr_array)

            para_line_l = map(str, [min_p, n_remian]) + p_result_list_srt
            para_line_str = '\t'.join([line_id, '\t'.join(para_line_l)]) + '\n'

            coutn_line_l = map(str, [min_p, n_remian]) + test_array_list
            count_line_str = '\t'.join([line_id, '\t'.join(coutn_line_l)]) + '\n'
            # result_test_count.append(count_line)

            # result_parameter.append(para_line_str)
            # result_parameter.append(count_line_str)

            f_out_para.write(para_line_str)
            f_out_para.write(count_line_str)

            i = 0
            if min_p < 0.05:
                for expr in tmp_array:
                    ans_str = '\t'.join(map(str, expr))
                    ans_line_str = '\t'.join([line_id, type_l[i], ans_str]) + '\n'
                    # result_exp.append(ans_line_str)
                    f_out.write(ans_line_str)
                    i += 1
            tmp = []

    # out_addr_exp = '%s_%s.txt' % (path, postfix)
    # out_addr_para = '%s_%s_para.txt' % (path, postfix)

    # out_addr_count = '%s_%s_count.txt' % (path, postfix)

    # f_out = open(out_addr_exp, 'w')
    # f_out.write(''.join(result_exp))
    f_out.write('#')
    f_out.close()

    # f_out_para = open(out_addr_para, 'w')
    # f_out_para.write(''.join(result_parameter))
    f_out_para.close()

    # f_out_count = open(out_addr_count, 'w')
    # f_out_count.write(''.join(result_test_count))
    # f_out_count.close()


def usage():

    print '''
    zscore_calc.py isfo_juncts.txt postfix
    
    out: isfo_juncts_postifx.txt  
         isfo_juncts_postifx_para.txt

    '''


def main():
    print ' '.join(sys.argv)
    print 'start:',time.ctime()
    in_addr = sys.argv[1]
    postfix = sys.argv[2]
    filter_fun(in_addr, postfix)
    print 'end:',time.ctime()
main()

