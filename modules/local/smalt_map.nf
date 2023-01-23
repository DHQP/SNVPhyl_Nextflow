/* Map reads to reference genome & create BAM file */
process SMALT_MAP {
    tag "$meta.id"
    label 'process_medium'
    container = "staphb/smalt:0.7.6"

    input:
    tuple val(meta), file(reads)
    path(ref_fai)
    path(ref_sma)
    path(ref_smi)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bams
    path("versions.yml"),                    emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    REFNAME=\$(basename ${ref_sma} .sma)
    smalt map -f bam -n 4 -l pe -i 1000 -j 20 -r 1 -y 0.5 -o ${prefix}.bam \${REFNAME} ${reads[0]} ${reads[1]}
    smalt version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smalt: \$( smalt version | grep "Version:" | sed 's/Version://' )
    END_VERSIONS
    """
}