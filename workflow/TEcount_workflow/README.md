# TEcount Workflow

English | [简体中文](README.zh-CN.md)

Use STAR to align paired-end silkworm mRNA sequencing reads, then use TEcount to count reads assigned to genes and transposable elements. Existing BAM files can skip STAR. A single samplesheet can mix both input types. The workflow does not run fastp, build the STAR index, or run edgeR.

## Run the workflow

Run from the repository root. By default, only TEcount runs in Apptainer, using `sifs/tetranscripts.sif`. Pass `--tecount_image` to use a different local image. Nextflow must be on the host's `PATH`. FASTQ input also requires `STAR` and `zcat`; merging requires `csvtk`; converting HiTE annotations requires `python3`. If a compatible `TEcount` is installed on the host, use `-profile native` to run without an image.

```bash
nextflow run workflow/TEcount_workflow \
  --samplesheet workflow/TEcount_workflow/samplesheet.example.csv \
  --merge -resume
```

The example samplesheet uses an existing BAM for RNPS1 and cleaned FASTQ files for the control. To check RNPS1 alone, run the following from the repository root. Neither `--merge` nor a STAR index is needed.

```bash
mkdir -p work/tecount_workflow
printf 'sample,fastq_1,fastq_2,bam\nRNPS1,,,%s/work/trasposon_analysis/STAR/RNPS1/Aligned.out.bam\n' "$PWD" \
  > work/tecount_workflow/rnps1.samples.csv
nextflow run workflow/TEcount_workflow \
  --samplesheet work/tecount_workflow/rnps1.samples.csv \
  --outdir work/tecount_workflow/rnps1_check -resume
```

Paths in the samplesheet are relative to the directory containing the samplesheet, unless they are absolute. If you move the example samplesheet, adjust the paths accordingly. Each row must specify either `bam` or both `fastq_1` and `fastq_2` (paired-end gzip-compressed FASTQ files). BAM records for the same read must be adjacent, as in the `Aligned.out.bam` produced by this workflow. Do not pass a coordinate-sorted BAM directly.

Reference files default to the repository's `reference/` directory. Override them with `--star_index`, `--gene_gtf`, or `--te_gtf`. If `reference/STARIndex` is a broken symlink, pass `--star_index /path/to/star_index` to point to a working index. `--te_gff` converts a HiTE GFF with the repository's existing script and takes precedence over `--te_gtf`. Relative reference paths passed on the command line are resolved from the launch directory.

## Outputs and verification

- Published outputs default to `work/tecount_workflow/`; override this with `--outdir`. The separate Nextflow task cache defaults to `work/.tecount_nextflow/`; override it with `-work-dir`.
- `STAR/<sample>/`: unsorted BAM and STAR logs for FASTQ inputs. The two multimapping parameters remain at `100`. STAR uses `14` threads by default; adjust this with `--star_cpus`.
- `TEcount/<sample>.cntTable` and `TEcount/<sample>.log`: per-sample counts and logs. TEcount runs with `--mode multi`.
- `--merge` also writes `TEcount/merged_TEcount.tsv`. It inner-joins by `gene/TE` (matching the original `csvtk join` default), keeping only entries present in all samples. Columns follow samplesheet order. For a single sample, only the count-column header is renamed to the sample name.

To verify RNPS1, compare the new and historical counts by `gene/TE`. The second header field in the raw TEcount table differs because the staged BAM has a different name. Reruns replace published files produced with the same names but do not remove files that the new run does not produce. Use a new `--outdir` if the sample set or output options change, and do not share an output directory across concurrent runs. By default, the workflow leaves the historical `work/trasposon_analysis/` results untouched.
