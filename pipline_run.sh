# jsQTL analysis pipline demo
# copy folder "junct_bed" and "py" to working folder.
# copy toy.bed toy.bim toy.fam toy.geno to working folder.
# run the pipline_run.sh
# sh pipline_run.sh

# require python 2.7, samtools, regtools
 

###########################
# public software or tools
###########################
## 1st extract and count the unique junction reads from the bam or cram 
## 
## extract unique junction reads
# samtools view -bS -h -q 10 ${INPUT_CRAM} --reference ${INPUT_REF} -o ${INPUT_CRAM}.uniq.bam
# samtools index ${INPUT_CRAM}.uniq.bam ${INPUT_CRAM}.uniq.bai
## count the junction reads and output bed file in junct_bed/
# regtools junction extract -a 4 -i 50 -o ${OUTPUT_BED} ${INPUT_CRAM}.uniq.bam 

##########################
# self-script
# 2nd
# merge_sortStrand_splitChr_dsub.sh
tissue=tissue1

if [ ! -d "junct_extract" ]; then
  mkdir junct_extract
fi

python py/fix_jun_read_bed_v0.3.py junct_extract/${tissue}_junct_extract_merge.txt junct_bed/*.bed

# strand + 
head -n1 junct_extract/${tissue}_junct_extract_merge.txt > junct_extract/${tissue}_junct_extract_merge.sorted+.txt
awk 'NR!=1 && $4=="+" {print $0}' junct_extract/${tissue}_junct_extract_merge.txt| \
sort -t $'\t' -k 1,1 -k 4,4 -k 2,2n -k 3,3n >>junct_extract/${tissue}_junct_extract_merge.sorted+.txt

# strand -  reverse
head -n1 junct_extract/${tissue}_junct_extract_merge.txt >junct_extract/${tissue}_junct_extract_merge.sorted-.txt 
awk 'NR!=1 && $4=="-" {print $0}' junct_extract/${tissue}_junct_extract_merge.txt| \
sort -t $'\t' -k 1,1 -k 4,4 -k 3,3nr -k 2,2nr >>junct_extract/${tissue}_junct_extract_merge.sorted-.txt

# split by Chr1-22
cd junct_extract
python ../py/split_junct_by_chr.py ${tissue}_junct_extract_merge.sorted-.txt
python ../py/split_junct_by_chr.py ${tissue}_junct_extract_merge.sorted+.txt
cd ..

# 3rd reorder geno 
# reorder genotype sample to match RNA sample and remove unmatched sample
dir=geno_${tissue}
if [ ! -d "${dir}" ]; then
  mkdir ${dir}
fi

python py/reorder_genowithRNA_v1.0.py junct_extract/${tissue}_junct_extract_merge.txt toy.geno ${dir}/toy.${tissue}.geno

# 4th reorder RNA sample to match genotype sample and remove unmatched sample
# reorder_rnaID_strand_dsub.sh
num=1 # chromosome number 1-22. In the toy example, there is only chromosom 1. 

python py/reorder_removeRNAwithUnmapGeno_v0.2.py ./junct_extract/${tissue}_junct_extract_merge.sorted+_${num}.txt \
	${dir}/toy.${tissue}.geno \
	./junct_extract/${tissue}_junct_extract_merge.sorted+_${num}.rmUnMapID.txt
python py/reorder_removeRNAwithUnmapGeno_v0.2.py ./junct_extract/${tissue}_junct_extract_merge.sorted-_${num}.txt \
	${dir}/toy.${tissue}.geno \
	./junct_extract/${tissue}_junct_extract_merge.sorted-_${num}.rmUnMapID.txt


# 5th 
# find Junction skipping event and filter with leave-one-out test
# findJun_gtf_filter_dsub.sh

output_dir=whole_genome_sumLongJun_gtf_strand
if [ ! -d "${output_dir}" ]; then
  mkdir ${output_dir}
fi

python py/find_skip_exon_sum_v1.3.1.py  \
junct_extract/${tissue}_junct_extract_merge.sorted+_${num}.rmUnMapID.txt \
toy.gtf \
${output_dir}/${tissue}_junct_exonskipping_+_${num}
python py/find_skip_exon_sum_v1.3.1.py  \
junct_extract/${tissue}_junct_extract_merge.sorted-_${num}.rmUnMapID.txt \
toy.gtf \
${output_dir}/${tissue}_junct_exonskipping_-_${num}

cd ${output_dir}
python ../py/filter_calc_1.3.py ${tissue}_junct_exonskipping_+_${num}_sumLJ.txt filterInOut
python ../py/filter_calc_1.3.py ${tissue}_junct_exonskipping_-_${num}_sumLJ.txt filterInOut

### 5.5th 
# filter count to n0.05(leave-one-out)

python ../py/filter_calc_ckeck_Ninout_v0.1.py ${tissue}_junct_exonskipping_-_${num}_sumLJ_filterInOut_para.txt \
${tissue}_junct_exonskipping_-_${num}_sumLJ_filterInOut_para_n0.05.txt
python ../py/filter_calc_outlier_v0.1.py ${tissue}_junct_exonskipping_-_${num}_sumLJ_filterInOut.txt \
${tissue}_junct_exonskipping_-_${num}_sumLJ_filterInOut_para_n0.05.txt \
${tissue}_junct_exonskipping_-_${num}_sumLJ_filterInOut_outlier5.txt

python ../py/filter_calc_ckeck_Ninout_v0.1.py ${tissue}_junct_exonskipping_+_${num}_sumLJ_filterInOut_para.txt \
${tissue}_junct_exonskipping_+_${num}_sumLJ_filterInOut_para_n0.05.txt
python ../py/filter_calc_outlier_v0.1.py ${tissue}_junct_exonskipping_+_${num}_sumLJ_filterInOut.txt \
${tissue}_junct_exonskipping_+_${num}_sumLJ_filterInOut_para_n0.05.txt \
${tissue}_junct_exonskipping_+_${num}_sumLJ_filterInOut_outlier5.txt
cd ..

# 6th MW test jsQTL
# asj_chr_SKrate_np_ks_test_dsub.sh

output_dir=chrAll_SK_np_mw
if [ ! -d "${dir}" ]; then
  mkdir ${output_dir}
fi
python py/heter_skipExon_expr_test_v3.9.0s_mw_1.py \
-i whole_genome_sumLongJun_gtf_strand/${tissue}_junct_exonskipping_-_${num}_sumLJ_filterInOut_outlier5.txt \
-n geno_${tissue}/toy.${tissue}.geno \
-o ${output_dir}/${tissue}_asj_-_chr${num}_sumLJ_filterInOut_outlier5_SK_np_mw_test \
-p 0.01
python py/heter_skipExon_expr_test_v3.9.0s_mw_1.py \
-i whole_genome_sumLongJun_gtf_strand/${tissue}_junct_exonskipping_+_${num}_sumLJ_filterInOut_outlier5.txt \
-n geno_${tissue}/toy.${tissue}.geno \
-o ${output_dir}/${tissue}_asj_+_chr${num}_sumLJ_filterInOut_outlier5_SK_np_mw_test \
-p 0.01

# 7th 
# summary the output file 
# loop_cat_sort_uniq_dsub.sh
output_dir=${tissue}_asj_chrAll_n0.05
if [ ! -d "${output_dir}" ]; then
  mkdir ${output_dir}
fi

cat chrAll_SK_np_mw/${tissue}_asj_+_chr${num}_sumLJ_filterInOut_outlier5_SK_np_mw_test/*.geno_gp_p.txt \
	chrAll_SK_np_mw/${tissue}_asj_-_chr${num}_sumLJ_filterInOut_outlier5_SK_np_mw_test/*.geno_gp_p.txt \
 	> ${output_dir}/chr${num}_${tissue}.geno_gp_p.txt

sort -k1,1 -k2,2 -u ${output_dir}/chr${num}_${tissue}.geno_gp_p.txt > ${output_dir}/chr${num}_${tissue}.geno_gp_p_sorted_uniq.txt
tail -n 1 ${output_dir}/chr${num}_${tissue}.geno_gp_p_sorted_uniq.txt > ${output_dir}/chr${num}_${tissue}.geno_gp_p_sorted_uniq_new.txt
awk '$1!="GeneID" {print $0}' ${output_dir}/chr${num}_${tissue}.geno_gp_p_sorted_uniq.txt >> ${output_dir}/chr${num}_${tissue}.geno_gp_p_sorted_uniq_new.txt
# head -n -1 ${output_dir}/chr${num}_${tissue}.geno_gp_p_sorted_uniq.txt >> ${output_dir}/chr${num}_${tissue}.geno_gp_p_sorted_uniq_new.txt
mv ${output_dir}/chr${num}_${tissue}.geno_gp_p_sorted_uniq_new.txt ${output_dir}/chr${num}_${tissue}.geno_gp_p_sorted_uniq.txt
rm ${output_dir}/chr${num}_${tissue}.geno_gp_p.txt

cat chrAll_SK_np_mw/${tissue}_asj_+_chr${num}_sumLJ_filterInOut_outlier5_SK_np_mw_test/*.geno_gp_count_sig.txt \
 	chrAll_SK_np_mw/${tissue}_asj_-_chr${num}_sumLJ_filterInOut_outlier5_SK_np_mw_test/*.geno_gp_count_sig.txt  \
 	|sort -u -k 1,1 -k 2,2 -r > ${output_dir}/chr${num}_${tissue}.geno_gp_count_sig_uniq.txt

python py/extract_sig_p_v0.1.py ${output_dir}/chr${num}_${tissue}.geno_gp_count_sig_uniq.txt ${output_dir}/chr${num}_${tissue}.geno_gp_p_sorted_uniq.txt

python py/find_minp_snp_v0.2.py ${output_dir}/chr${num}_${tissue}.geno_gp_p_sorted_uniq.txt \
${output_dir}/chr${num}_${tissue}.geno_gp_p_sorted_uniq_np_minp.txt

# 
# 8th Fine-mapping with LD

output_dir=${tissue}_asj_chrAll_p0.01_n0.05_r2
if [ ! -d "${output_dir}" ]; then
  mkdir ${output_dir}
fi

cd ${output_dir}
# create a file with index SNPid list
awk 'NR !=1 && $5 < 0.01 {print $2}' \
../${tissue}_asj_chrAll_n0.05/chr${num}_${tissue}.geno_gp_p_sorted_uniq_np_minp.txt | sort -u \
> chr${num}_${tissue}.geno_gp_p_sorted_uniq_np_top95_minp_snpID.txt
# make index SNP LD with other SNP
plink --bfile  ../toy \
--ld-snp-list chr${num}_${tissue}.geno_gp_p_sorted_uniq_np_top95_minp_snpID.txt \
--ld-window 9999999 \
--ld-window-kb 500 \
--ld-window-r2 0 \
--out chr${num}_${tissue}.geno_gp_p_sorted_uniq_np_top95_minp_snpID \
--r2 in-phase dprime with-freqs
# Fine-mapping
python ../py/find_sig_snp_v0.2.0.py \
../${tissue}_asj_chrAll_n0.05/chr${num}_${tissue}.geno_gp_p_sorted_uniq.txt \
chr${num}_${tissue}.geno_gp_p_sorted_uniq_np_top95_minp_snpID.ld \
chr${num}_${tissue}.geno_gp_p_sorted_uniq_np_top95.txt \
--p_cutoff 0.01
cd ..

