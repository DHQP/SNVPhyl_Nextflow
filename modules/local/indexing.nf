process INDEXING {
    tag "${refgenome}"
    label 'process_low'
    container = "staphb/smalt:0.7.6"

    input:
    path(refgenome)

    output:
    path('*.fai'),        emit: ref_fai
    path('*.sma'),        emit: ref_sma
    path('*.smi'),        emit: ref_smi
    path("versions.yml"), emit: versions

    script:
    """
    REF_BASENAME=\$(basename ${refgenome} .fasta)
    smalt index -k 13 -s 6 \${REF_BASENAME} ${refgenome}
    samtools faidx ${refgenome}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smalt: \$( smalt version | grep "Version:" | sed 's/Version://' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}