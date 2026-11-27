nextflow.enable.dsl = 2

process CONVERT_TE_GFF {
    input:
    path 'converter.py'
    path 'te.gff'

    output:
    path 'HiTE.full_length.tetranscripts.gtf'

    script:
    """
    python3 converter.py te.gff HiTE.full_length.tetranscripts.gtf
    """
}

process STAR_ALIGN {
    tag "${sample}"
    cpus params.star_cpus
    publishDir { "${params.outdir}/STAR/${sample}" }, mode: 'copy', overwrite: true, failOnError: true

    input:
    tuple val(sample), val(rowIndex), path('R1.fastq.gz'), path('R2.fastq.gz'), path('star_index')

    output:
    tuple val(sample), val(rowIndex), path('Aligned.out.bam'), emit: bam
    path 'Log*.out', emit: logs

    script:
    """
    STAR --runThreadN ${task.cpus} --genomeDir star_index \\
        --readFilesIn R1.fastq.gz R2.fastq.gz --readFilesCommand zcat \\
        --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 \\
        --outSAMtype BAM Unsorted
    """
}

process TE_COUNT {
    tag "${sample}"
    container { image ?: null }
    publishDir "${params.outdir}/TEcount", mode: 'copy', overwrite: true, failOnError: true

    input:
    tuple val(sample), val(rowIndex), path('sample.bam')
    path 'genes.gtf'
    path 'transposons.gtf'
    val image

    output:
    tuple val(sample), val(rowIndex), path("${sample}.cntTable"), emit: table
    path "${sample}.log", emit: log

    script:
    """
    TEcount --format BAM --mode multi -b sample.bam \\
        --GTF genes.gtf --TE transposons.gtf \\
        --outdir . --project '${sample}' 2>&1 | tee '${sample}.log'
    """
}

process MERGE_COUNTS {
    publishDir "${params.outdir}/TEcount", mode: 'copy', overwrite: true, failOnError: true

    input:
    tuple val(samples), path(tables)

    output:
    path 'merged_TEcount.tsv'

    script:
    def names = samples.join(',')
    def columns = (2..(samples.size() + 1)).join(',')
    def quotedTables = tables.collect { "'${it}'" }.join(' ')
    def rename = "csvtk rename -t -f '${columns}' -n '${names}'"
    def command = samples.size() == 1
        ? "${rename} ${quotedTables}"
        : "csvtk join -t -f 1 --na 0 ${quotedTables} | ${rename}"
    """
    ${command} > merged_TEcount.tsv
    """
}

workflow {
    def checkedPath = { base, value, label, directory = false, hint = '' ->
        def path = base.resolve(value.toString()).normalize()
        if (!(directory ? path.isDirectory() : path.isFile()) || !java.nio.file.Files.isReadable(path)) {
            error "${label} not found or unreadable: ${path}${hint}"
        }
        path.toRealPath()
    }

    if (!params.samplesheet) {
        error 'Specify --samplesheet (see samplesheet.example.csv).'
    }

    def sheetPath = launchDir.resolve(params.samplesheet.toString()).normalize()
    def sheet = checkedPath.call(launchDir, sheetPath, 'Samplesheet')
    def sheetBase = sheetPath.parent
    def geneGtf = checkedPath.call(launchDir, params.gene_gtf, 'Gene GTF')
    def image = workflow.containerEngine == 'apptainer'
        ? checkedPath.call(launchDir, params.tecount_image, 'TEcount container image').toString()
        : ''

    if (params.te_gff) {
        def gff = checkedPath.call(launchDir, params.te_gff, 'TE GFF')
        def converter = checkedPath.call(projectDir, '../../scripts/transposon_analysis/hite_to_tetranscripts_gtf.py', 'HiTE converter')
        // With only value inputs, conversion returns a reusable value channel.
        teGtf = CONVERT_TE_GFF(converter, gff)
    } else {
        teGtf = channel.value(checkedPath.call(launchDir, params.te_gtf, 'TE GTF'))
    }

    samples = channel.value(sheet)
        .splitCsv(header: true, strip: true, quote: '"')
        .toList()
        .flatMap { rows ->
            if (!rows) {
                error 'Samplesheet has no samples.'
            }
            def seen = [] as Set
            def starIndex
            rows.withIndex().collect { row, rowIndex ->
                def sample = row.sample
                if (!(sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) || !seen.add(sample)) {
                    error "Invalid or duplicate sample name: ${sample}"
                }
                if (row.bam && !row.fastq_1 && !row.fastq_2) {
                    return [sample, rowIndex, checkedPath.call(sheetBase, row.bam, "${sample} BAM")]
                }
                if (!row.bam && row.fastq_1 && row.fastq_2) {
                    starIndex = starIndex ?: checkedPath.call(launchDir, params.star_index, 'STAR index', true,
                        '. Use --star_index to specify a valid index directory.')
                    def read1 = checkedPath.call(sheetBase, row.fastq_1, "${sample} FASTQ 1")
                    def read2 = checkedPath.call(sheetBase, row.fastq_2, "${sample} FASTQ 2")
                    if (read1 == read2) {
                        error "${sample}: fastq_1 and fastq_2 must reference different files."
                    }
                    return [
                        sample,
                        rowIndex,
                        read1,
                        read2,
                        starIndex
                    ]
                }
                error "${sample}: specify either bam or both fastq_1 and fastq_2."
            }
        }
        .branch { row ->
            bam: row.size() == 3
            fastq: row.size() == 5
        }

    aligned = STAR_ALIGN(samples.fastq)
    counted = TE_COUNT(samples.bam.mix(aligned.bam), geneGtf, teGtf, image)

    if (params.merge) {
        ordered = counted.table.toList().map { entries ->
            entries.sort { it[1] }
            tuple(entries.collect { it[0] }, entries.collect { it[2] })
        }
        MERGE_COUNTS(ordered)
    }
}
