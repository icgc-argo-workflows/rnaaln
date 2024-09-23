process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' :
        'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' }"

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)

    output:
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    path  "versions.yml"                      , emit: versions

    tuple val(meta), path('*d.out.bam')              , optional:true, emit: bam
    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')               , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab
    tuple val(meta), path('*.SJ.out.tab')            , optional:true, emit: spl_junc_tab
    tuple val(meta), path('*.ReadsPerGene.out.tab')  , optional:true, emit: read_per_gene_tab
    tuple val(meta), path('*.out.junction')          , optional:true, emit: junction
    tuple val(meta), path('*.out.sam')               , optional:true, emit: sam
    tuple val(meta), path('*.wig')                   , optional:true, emit: wig
    tuple val(meta), path('*.bg')                    , optional:true, emit: bedgraph

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.study_id}.${meta.patient}.${meta.sample}.${meta.id}"
    def reads1 = [], reads2 = [] // first half in reads1, second half in reads2
    meta.single_end ? [reads].flatten().each{reads1 << it} : reads.eachWithIndex{ v, ix -> (ix < reads.size() / 2 ? reads1 : reads2) << v }
    def attrRG          = '' // read group information
    if (meta.read_group.contains(",")) {
        readGroups = meta.read_group.split(",")
        readGroups.each { readGroup ->
            attrRG += "ID:${readGroup.trim()} SM:$meta.sample CN:$meta.sequencing_center PL:$meta.platform, "
            }
        attrRG = attrRG[0..-3] // Remove the trailing comma and space
    } else {
        attrRG = "ID:$meta.read_group SM:$prefix CN:$meta.sequencing_center PL:$meta.platform"
    }
    def uncompressionCommand = '' // file uncompression command
    if (reads.toList()[0].toString().endsWith('.gz')) {
        uncompressionCommand = '--readFilesCommand zcat'
    } else if (reads.toList()[0].toString().endsWith('.bz2')) {
        uncompressionCommand = '--readFilesCommand bzcat'
    }

    def readgroupinfo = meta.id.split(',').collect { "ID:$it" }.join(' , ')

    """
    echo ${attrRG}

    STAR \\
        --genomeDir $index \\
        --sjdbGTFfile $gtf \\
        --readFilesIn ${reads1.join(",")} ${reads2.join(",")} \\
        $uncompressionCommand \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        --twopassMode Basic \\
        --outFilterMultimapScoreRange 1 \\
        --outFilterMultimapNmax 20 \\
        --outFilterMismatchNmax 10 \\
        --alignIntronMax 500000 \\
        --alignMatesGapMax 1000000 \\
        --sjdbScore 2 \\
        --alignSJDBoverhangMin 1 \\
        --genomeLoad NoSharedMemory \\
        --outFilterMatchNminOverLread 0.33 \\
        --outFilterScoreMinOverLread 0.33 \\
        --outSAMstrandField intronMotif \\
        --outSAMmode Full \\
        --outSAMattributes NH HI NM MD AS XS \\
        --outSAMunmapped Within \\
        --limitSjdbInsertNsj 2000000 \\
        --outSAMtype BAM Unsorted SortedByCoordinate \\
        --outSAMheaderHD @HD VN:1.4 \\
        --quantMode TranscriptomeSAM \\
        --outSAMattrRGline $readgroupinfo \\
        $args

    mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam

    samtools reheader -P -c \'sed -e "s/^@RG.*/${meta.read_group.replaceAll(/'/,'')}/g"\' ${prefix}.Aligned.sortedByCoord.out.bam > ${prefix}.star_Aligned.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}Xd.out.bam
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    touch ${prefix}.sortedByCoord.out.bam
    touch ${prefix}.toTranscriptome.out.bam
    touch ${prefix}.Aligned.unsort.out.bam
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.unmapped_1.fastq.gz
    touch ${prefix}.unmapped_2.fastq.gz
    touch ${prefix}.tab
    touch ${prefix}.SJ.out.tab
    touch ${prefix}.ReadsPerGene.out.tab
    touch ${prefix}.Chimeric.out.junction
    touch ${prefix}.out.sam
    touch ${prefix}.Signal.UniqueMultiple.str1.out.wig
    touch ${prefix}.Signal.UniqueMultiple.str1.out.bg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
}

//  // samtools reheader -P -c \'sed -e "s/^@RG.*/${meta.read_group.replaceAll(/'/,'')}/g"\' ${prefix}.Aligned.sortedByCoord.out.bam | awk '!seen[$0]++' > ${prefix}.star_Aligned.bam