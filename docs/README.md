# piRNA项目

## data preprocessing

```bash
# raw data quality check
fastqc -o $PWD/data/data/dm_test_raw_date/fastqc_reports $PWD/data/dm_test_raw_date/*.fastq.gz -t 48
bash $PWD/scripts/process_fastq_with_fastp.sh $PWD/data/dm_test_raw_date $PWD/data/data/dm_test_clean_date
bash $PWD/scripts/rename.sh $PWD/data/dm_smRNA_clean_data $PWD/data/sample_map.txt
fastqc -o $PWD/data/data/dm_test_clean_date/fastqc_reports $PWD/data/dm_smRNA_clean_data/*.fq.gz -t 48
for dir in $PWD/data/dm_smRNA_clean_data/fastp_reports $PWD/data/dm_smRNA_clean_data/fastqc_reports ;do
    $HOME/miniforge3/envs/trimgalore/bin/multiqc -o $dir/multiqc_report $dir
done
# calculate the frequency of reads in each file and generate length distribution plot
bash $PWD/scripts/count_length.sh $PWD/data/dm_test_clean_date
python3 $PWD/scripts/plot_length_distribution.py $PWD/data/dm_test_clean_date/length_distribution.txt -o $PWD/data/dm_test_clean_date/plots

# remove reads about dm rRNA
bash $PWD/scripts/remove_rRNA.sh $PWD/data/dm_smRNA_clean_data $HOME/software/piPipes/common/dm6/BowtieIndex $PWD/data/dm_smRNA_clean_data_x_rRNA
bash $PWD/scripts/count_length.sh $PWD/data/dm_smRNA_clean_data_x_rRNA_merge
python3 $PWD/scripts/plot_length_distribution.py $PWD/data/dm_smRNA_clean_data_x_rRNA_merge/length_distribution.txt -o $PWD/data/dm_smRNA_clean_data_x_rRNA_merge/plots

```

<!-- ## piRNA分析--piPipes

```bash
# Batch processing of multiple samples with piPipes
bash $PWD/scripts/piPipes_run.sh $PWD/data/dm_smRNA_clean_data $PWD/results/dm_smRNA_piPipes_out $PWD/logs/dm_smRNA_piPipes_logs 20 dm6

# test
piPipes_debug small \
        -i data/dm_smRNA_clean_data/DXH15-Rep1.fq.gz \
        -g dm6 \
        -o results/BDSC41992-Rep2.piPipes_out \
        -c 70 \
        -N uniqueXmiRNA \
        1> logs/DXH15-Rep1.piPipes.stdout \
        2> logs/DXH15-Rep1.piPipes.stderr

# Tailor pipeline test
bash $PWD/scripts/tailor_run.sh $PWD/data/dm_smRNA_clean_data $PWD/results/dm_smRNA_tailor_out 20 dm3

run_tailing_pipeline.sh \
	-i $PWD/data/dm_smRNA_clean_data/DXH15-RNAi-Rep1.fq.gz \
	-g $HOME/software/Tailor/annotation/dm3/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa \
	-t $HOME/software/Tailor/annotation/dm3.genomic_features \
	-o $PWD/results/dm_smRNA_tailer_out/DXH15-RNAi-Rep1 \
	-c 24
``` -->

## percentage of piRNA 5‘ end shift and 3‘ end position

```bash
# bmori
./scripts/process_pirna.sh data/bm_smRNA_clean_data data/bm_smRNA_pirna_out data/bmo.v3.0.fa --skip-existing

# dm
./scripts/process_pirna.sh data/dm_smRNA_clean_data_x_rRNA data/dm_smRNA_clean_data_x_rRNA_out data/dme.v3.0.fa --skip-existing

# dm_merge
./scripts/process_pirna.sh data/dm_smRNA_clean_data_x_rRNA_merge data/dm_smRNA_clean_data_x_rRNA_merge_out data/dme.v3.0.fa --skip-existing

## test
fastx_collapser -i data/bm_smRNA_pirna_out/Control.fa -o data/bm_smRNA_pirna_out/Control.fa.rmdup
seqmap 0 data/bmo.v3.0.fa data/bm_smRNA_pirna_out/Control.fa.collapsed data/bm_smRNA_pirna_out/Control.fa.collapsed.test.map /do_not_output_probe_without_match /forward_strand /cut:3,23 /output_all_matches

## piRNA difference analysis
for i in data/bm_smRNA_pirna_out/*.fa.collapsed
do
    # awk 'NR>1 {print $3}' $i.map | uniq > $i.map.uniq
    seqtk subseq $i $i.map.uniq > $i.map.uniq.fa
done

# Control vs DHX15-KD3
Rscript scripts/pirna_edger_pairwise.R data/bm_smRNA_pirna_out/Control.fa.collapsed.map.uniq.fa data/bm_smRNA_pirna_out/DHX15-KD3.fa.collapsed.map.uniq.fa
csvtk join -f 1 results/pirna_edger_pairwise_out/Control.fa.collapsed.map_vs_DHX15-KD3.fa.collapsed.map_count_matrix.csv results/pirna_edger_pairwise_out/Control.fa.collapsed.map_vs_DHX15-KD3.fa.collapsed.map_edger_results.csv | \
csvtk filter -f "FDR<0.05" > results/pirna_edger_pairwise_out/Control_vs_DHX15-KD3_edger_results_with_counts_FDR_0.05.csv

cat results/pirna_edger_pairwise_out/Control_vs_DHX15-KD3_edger_results_with_counts_FDR_0.05.csv | \
csvtk filter2 -f '$control + $treatment >= 100' | \
csvtk sort -k logFC:nr | \
csvtk mutate2 -n Length -e 'len($1)' > results/pirna_edger_pairwise_out/Control_vs_DHX15-KD3_edger_results_with_counts_FDR_0.05_length.csv

cat results/pirna_edger_pairwise_out/Control_vs_DHX15-KD3_edger_results_with_counts_FDR_0.05_length.csv | \
csvtk filter -f "logFC>0" | \
csvtk plot hist -f Length --bins 34 --title "Upregulated piRNA length distribution" -o results/pirna_edger_pairwise_out/Control_vs_DHX15-KD3_upregulated_pirna_length_distribution.pdf

cat results/pirna_edger_pairwise_out/Control_vs_DHX15-KD3_edger_results_with_counts_FDR_0.05_length.csv | \
csvtk filter -f "logFC<0" | \
csvtk plot hist -f Length --bins 9 --title "Downregulated piRNA length distribution" -o results/pirna_edger_pairwise_out/Control_vs_DHX15-KD3_downregulated_pirna_length_distribution.pdf


# Control vs SUGP1-KD2
Rscript scripts/pirna_edger_pairwise.R data/bm_smRNA_pirna_out/Control.fa.collapsed.map.uniq.fa data/bm_smRNA_pirna_out/SUGP1-KD2.fa.collapsed.map.uniq.fa
csvtk join -f 1 results/pirna_edger_pairwise_out/Control.fa.collapsed.map_vs_SUGP1-KD2.fa.collapsed.map_count_matrix.csv results/pirna_edger_pairwise_out/Control.fa.collapsed.map_vs_SUGP1-KD2.fa.collapsed.map_edger_results.csv | \
csvtk filter -f "FDR<0.05" > results/pirna_edger_pairwise_out/Control_vs_SUGP1-KD2_edger_results_with_counts_FDR_0.05.csv

cat results/pirna_edger_pairwise_out/Control_vs_SUGP1-KD2_edger_results_with_counts_FDR_0.05.csv | \
csvtk filter2 -f '$control + $treatment >= 100' | \
csvtk sort -k logFC:nr | \
csvtk mutate2 -n Length -e 'len($1)' > results/pirna_edger_pairwise_out/Control_vs_SUGP1-KD2_edger_results_with_counts_FDR_0.05_length.csv

cat results/pirna_edger_pairwise_out/Control_vs_SUGP1-KD2_edger_results_with_counts_FDR_0.05_length.csv | \
csvtk filter -f "logFC>0" | \
csvtk plot hist -f Length --bins 34 --title "Upregulated piRNA length distribution" -o results/pirna_edger_pairwise_out/Control_vs_SUGP1-KD2_upregulated_pirna_length_distribution.pdf

cat results/pirna_edger_pairwise_out/Control_vs_SUGP1-KD2_edger_results_with_counts_FDR_0.05_length.csv | \
csvtk filter -f "logFC<0" | \
csvtk plot hist -f Length --bins 9 --title "Downregulated piRNA length distribution" -o results/pirna_edger_pairwise_out/Control_vs_SUGP1-KD2_downregulated_pirna_length_distribution.pdf

# 提取差异piRNA
python scripts/csv2fa.py results/pirna_edger_pairwise_out/Control_vs_SUGP1-KD2_edger_results_with_counts_FDR_0.05_length.csv

python3 scripts/csv2fa.py results/pirna_edger_pairwise_out/Control_vs_DHX15-KD3_edger_results_with_counts_FDR_0.05_length.csv


# 使用seqmap进行piRNA映射
for file in data/bm_smRNA_pirna_diffrencial/*.fa.collapsed
do
    FILENAME=$(basename "$file" | sed 's/\.fa\.collapsed//')
    FILE_LOG="logs/process_pirna_20250424_230656/${FILENAME}.log"
    OUTPUT_DIR="data/bm_smRNA_pirna_diffrencial"
    seqmap 0 data/bmo.v3.0.fa "$OUTPUT_DIR/${FILENAME}.fa.collapsed" "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map" \
    /do_not_output_probe_without_match /forward_strand /cut:3,23 /output_all_matches >> "$FILE_LOG" 2>&1
    # 使用awk过滤出map文件的第1，2，4列，3减去第二列，然后过滤掉第二列小于-2的行
    awk '{print $4"\t"3-$2"\t"$1}' "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map" | awk '$2 > -3' > "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map.filtered"
    mv "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map.filtered" "$OUTPUT_DIR/${FILENAME}.fa.collapsed.map"
done

export PYTHONPATH="$PWD/scripts:$PYTHONPATH"
python3 -m pirna_analysis.main -i data/bm_smRNA_pirna_diffrencial/*.fa.collapsed.map -o data/bm_smRNA_pirna_diffrencial_results -r data/bmo.v3.0.fa -v

python3 /home/gyk/project/ld_pirna/scripts/create_dual_sequence_logos.py -i /home/gyk/project/ld_pirna/results/pirna_edger_pairwise_out -o /home/gyk/project/ld_pirna/results/pirna_edger_pairwise_out/logos_weighted_down --control_name Control --treatment_name DHX15-KD3 --logfc_filter "<0"

python3 /home/gyk/project/ld_pirna/scripts/create_dual_sequence_logos.py -i /home/gyk/project/ld_pirna/results/pirna_edger_pairwise_out -o /home/gyk/project/ld_pirna/results/pirna_edger_pairwise_out/logos_weighted_up --logfc_filter ">0" --auto_detect --sequence_filter TCTTCGGTAGTATAGTGGTCAGTATCCC

# piRNA source analysis
for i in data/bm_smRNA_pirna_diffrencial/*.fa.collapsed; do
    sample_name=$(basename "$i" .fa.collapsed)
    # bowtie -v 0 -a --best --strata -p 20 -x $HOME/reference/bm_ncbi/GCF_030269925.1/bowtie_index/bt  -f $i -S "data/bm_smRNA_pirna_diffrencial_out/${sample_name}.fa.collapsed.sam"
    # samtools view -Sb "data/bm_smRNA_pirna_diffrencial_out/${sample_name}.fa.collapsed.sam" > "data/bm_smRNA_pirna_diffrencial_out/${sample_name}.fa.collapsed.bam"
    # bamToBed -i "data/bm_smRNA_pirna_diffrencial_out/${sample_name}.fa.collapsed.bam" > "data/bm_smRNA_pirna_diffrencial_out/${sample_name}.fa.collapsed.bed"
    # rm data/bm_smRNA_pirna_diffrencial_out/${sample_name}.fa.collapsed.bam data/bm_smRNA_pirna_diffrencial_out/${sample_name}.fa.collapsed.sam
done

gff2bed < reference/genomic_chr_fixed_introns_level2_transposon.gff > reference/feature_level2.bed
gff2bed < reference/genomic_chr_fixed_introns_level3.gff > reference/feature_level3.bed
for i in data/bm_smRNA_pirna_diffrencial/*.fa.collapsed; do
    sample_name=$(basename "$i" .fa.collapsed)
    intersectBed -a "data/bm_smRNA_pirna_diffrencial_out/${sample_name}.fa.collapsed.bed" \
                -b reference/feature_level3.bed -wa -wb -f 1 -s | csvtk cut -H -t -f 14 | sort | uniq -c | csvtk space2tab | csvtk tab2csv | csvtk cut -H -f 2,1 > "data/bm_smRNA_pirna_diffrencial_out/${sample_name}.fa.collapsed.feature.level3.csv"
done

csvtk join -H data/bm_smRNA_pirna_diffrencial_out/*.level2.csv -p -e | csvtk rename2 -F -f '*' -p '\.fa.collapsed.feature.level2-c2' > data/bm_smRNA_pirna_diffrencial_out/feature_level2.csv

csvtk join -H data/bm_smRNA_pirna_diffrencial_out/*.level3.csv -p -e | csvtk rename2 -F -f '*' -p '\.fa.collapsed.feature.level3-c2' > data/bm_smRNA_pirna_diffrencial_out/feature_level3.csv

Rscript scripts/plot_featrue_distribution.R -i data/bm_smRNA_pirna_diffrencial_out/feature_level2.csv -o results/diff_piRNA_analysis/feature_level2.pdf

Rscript scripts/plot_featrue_distribution.R -i data/bm_smRNA_pirna_diffrencial_out/feature_level3.csv -o results/diff_piRNA_analysis/feature_level3.pdf

intersectBed -a data/bm_smRNA_pirna_diffrencial_out/Control_vs_DHX15-KD3-control.fa.collapsed.bed \
-b reference/feature_level2.bed -wa -wb -f 1 -s | csvtk cut -H -t -f 4,14 | sort | uniq -c | csvtk space2tab | csvtk cut -H -t -f 2,3,1 | csvtk spread -H -t -k 2 -v 3 > "data/bm_smRNA_pirna_diffrencial_out/${sample_name}.fa.collapsed.feature.tsv"

intersectBed -a data/bm_smRNA_pirna_diffrencial_out/Control_vs_DHX15-KD3-control.fa.collapsed.bed \
-b reference/feature_level2.bed -wa -wb -f 1 -s | csvtk cut -H -t -f 14 | sort | uniq -c | csvtk space2tab | csvtk tab2csv | csvtk cut -H -f 2,1 | head
```

## B. mori transposon analysis
**分析思路**

1. First, perform quality control on raw sequencing data to check whether QC has already been applied (software: FastQC).

```bash
# Quality control for raw data
cd data/bm_mRNA_raw_data
mkdir -p fastqc_output multiqc_report
fastqc -o fastqc_output -t 70 *.fq.gz

# Summarize QC reports with MultiQC
multiqc fastqc_output -o multiqc_report
```

3. Based on QC results, some adapters were not removed in part of the data; adapters were trimmed using fastp.

```bash
scripts/process_fastq_with_fastp.sh -p data/bm_mRNA_raw_data data/bm_mRNA_clean_data
cd data/bm_mRNA_clean_data
mkdir -p multiqc_report
# 使用MultiQC汇总质控报告
multiqc fastp_reports -o multiqc_report
```

2. Align reads to the transposon genome using STAR, convert to BAM with samtools, and sort.

```bash
# Create STAR index
cd ~/reference/bm_ncbi/GCF_030269925.1
STAR --runThreadN 70 --runMode genomeGenerate --genomeDir STARIndex \
--genomeFastaFiles GCF_030269925.1_ASM3026992v2_genomic_chr.fa \
--sjdbGTFfile genomic_chr.gtf --sjdbOverhang 100 --genomeSAindexNbases 13
ln -sf $(pwd)/STARIndex $HOME$/project/ld_pirna/reference/STARIndex
# maping
cd $HOME$/project/ld_pirna
ls data/bm_mRNA_clean_data/*_R1.clean.fastq.gz | \
parallel -j 5 'r1={}; r2=${r1/_R1/_R2}; sample_name=$(basename $r1 _R1.clean.fastq.gz); \
output_dir=work/trasposon_analysis/STAR/${sample_name}; mkdir -p $output_dir; \
STAR --runThreadN 14 --genomeDir reference/STARIndex \
     --readFilesIn $r1 $r2 --readFilesCommand zcat \
     --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \
     --outFileNamePrefix $output_dir/ --outSAMtype BAM Unsorted'
```

3. 使用 TEcounts 进行转座子计数

```bash
## python /home/gyk/project/ld_pirna/scripts/transposon_analysis/hite_to_tetranscripts_gtf.py <输入gff文件> <输出gtf文件>
python /home/gyk/project/ld_pirna/scripts/transposon_analysis/hite_to_tetranscripts_gtf.py \
/home/gyk/project/ld_pirna/reference/HiTE.full_length.gff \
/home/gyk/project/ld_pirna/reference/HiTE.full_length.tetranscripts.gtf

mkdir -p work/trasposon_analysis/TEcount

# TEcount
find work/trasposon_analysis/STAR -name "*Aligned.out.bam" |
parallel -j 5 'sample_dir=$(dirname {});
sample_name=$(basename $sample_dir);
apptainer exec sifs/tetranscripts.sif TEcount \
  --format BAM --mode multi -b {} \
  --GTF reference/genomic_chr.gtf \
  --TE reference/HiTE.full_length.tetranscripts.gtf \
  --outdir work/trasposon_analysis/TEcount --project $sample_name'

# 合并所有样本的转座子计数结果
# 先进入TEcount目录
cd /home/gyk/project/ld_pirna/work/trasposon_analysis/TEcount

# 合并文件并重命名列
# 修正后的命令
csvtk join -t -f "1" \
	DHX15.cntTable RNPS1.cntTable SUGP1-Rep1.cntTable SUGP1-Rep2.cntTable control.cntTable \
	--na "0" | \
	csvtk rename -t -f "2,3,4,5,6" -n "DHX15,RNPS1,SUGP1-Rep1,SUGP1-Rep2,control" \
	--out-file merged_TEcount.tsv
```

4. 使用edger进行差异分析

```r
# 运行已编写的差异分析脚本
# 参数：计数矩阵文件、样本配置文件、输出目录
Rscript scripts/transposon_analysis/te_differential_analysis.R \
  work/trasposon_analysis/TEcount/merged_TEcount.tsv \
  data/sample_config.txt \
  work/trasposon_analysis/edgeR_results

# 或不带参数运行，使用脚本中的默认路径
# Rscript scripts/transposon_analysis/te_differential_analysis.R
```