# TEcount 工作流

[English](README.md) | 简体中文

用 STAR 比对家蚕双端 mRNA 测序数据，再用 TEcount 对基因和转座子计数。已有 BAM 可以直接计数；同一份样本表可以混用两种入口。工作流不包含 fastp、STAR 建索引或 edgeR。

## 运行

从仓库根目录运行。默认只有 TEcount 在 Apptainer 中运行，使用 `sifs/tetranscripts.sif`，可用 `--tecount_image` 指定其他本地镜像文件。宿主机的 `PATH` 需包含 Nextflow；用 FASTQ 时还需要 `STAR` 和 `zcat`，合并时需要 `csvtk`，转换 HiTE 注释时需要 `python3`。若宿主机已安装兼容的 `TEcount`，可用 `-profile native` 代替容器，不需要镜像文件。

```bash
nextflow run workflow/TEcount_workflow \
  --samplesheet workflow/TEcount_workflow/samplesheet.example.csv \
  --merge -resume
```

示例包含已有 RNPS1 BAM 和 control 的 clean FASTQ。若只想验证 RNPS1，不需要 `--merge` 或 STAR 索引，可从仓库根目录运行：

```bash
mkdir -p work/tecount_workflow
printf 'sample,fastq_1,fastq_2,bam\nRNPS1,,,%s/work/trasposon_analysis/STAR/RNPS1/Aligned.out.bam\n' "$PWD" \
  > work/tecount_workflow/rnps1.samples.csv
nextflow run workflow/TEcount_workflow \
  --samplesheet work/tecount_workflow/rnps1.samples.csv \
  --outdir work/tecount_workflow/rnps1_check -resume
```

样本表的输入路径相对样本表所在目录，也可用绝对路径；如果把示例表移到别处，请相应调整其中的路径。每行只能填写 `bam` 或同时填写 `fastq_1`、`fastq_2`（双端 gzip FASTQ）。BAM 必须让同一 read 的比对记录相邻，例如此流程生成的 STAR `Aligned.out.bam`；不要直接传入按坐标排序的 BAM。

参考文件默认取自仓库的 `reference/`，可用 `--star_index`、`--gene_gtf`、`--te_gtf` 覆盖。若原 `reference/STARIndex` 链接已失效，传 `--star_index /path/to/star_index` 指向可用索引。`--te_gff` 会用仓库现有脚本转换 HiTE GFF，并优先于 `--te_gtf`。命令行传入的参考路径相对启动目录解析。

## 输出和验证

- 默认输出到 `work/tecount_workflow/`，可用 `--outdir` 指定其他目录；任务缓存默认在独立的 `work/.tecount_nextflow/`，可用 Nextflow 的 `-work-dir` 指定。
- `STAR/<sample>/`：从 FASTQ 生成的未排序 BAM 和 STAR 日志。比对保留原流程的两个多重比对参数 `100`，线程数默认 `14`，可用 `--star_cpus` 调整。
- `TEcount/<sample>.cntTable` 和 `TEcount/<sample>.log`：单样本计数及日志；计数使用 `--mode multi`。
- `--merge` 额外输出 `TEcount/merged_TEcount.tsv`，按 `gene/TE` 内连接（与原 `csvtk join` 默认行为一致，只保留各样本都有的条目），列按样本表中的顺序排列。单样本时仅将计数列改为样本名。

验证 RNPS1 时，按 `gene/TE` 比较新旧计数；TEcount 的第二列表头会因为临时 BAM 文件名不同而改变。重跑会更新本次生成的同名发布文件，但不会删除本次未生成的旧文件。若样本或输出选项有变化，建议使用新的 `--outdir`，也不要让并行运行共用同一个输出目录。默认不会写入历史结果目录 `work/trasposon_analysis/`。
