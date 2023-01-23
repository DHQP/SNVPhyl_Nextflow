/* Making SNV matrix */
process MAKE_SNV {
    label 'process_low'
    container = "staphb/snvphyl-tools:1.8.2"

    input:
    path(snvAlignment_phy)

    output:
    path('snvMatrix.tsv'), emit: snvMatrix
    path("versions.yml"),  emit: versions

    script:
    """
    snv_matrix.pl ${snvAlignment_phy} -o snvMatrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snvphyl-tools:
        perl: \$(perl --version | grep "This is perl" | sed 's/.*(\\(.*\\))/\\1/' | cut -d " " -f1)
    END_VERSIONS
    """
}