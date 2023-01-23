/* Verifying Mapping Quality */
process VERIFYING_MAP_Q {
    tag ""
    label 'process_medium'
    container = "staphb/snvphyl-tools:1.8.2"

    input:
    path(sorted_bams)
    val(bam_line)

    output:
    path("mappingQuality.txt"), emit: mapping_quality
    path("versions.yml"),       emit: versions

    script:
    """
    verify_mapping_quality.pl -c 4 --min-depth 10 --min-map 80 --output mappingQuality.txt ${bam_line} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snvphyl-tools:
        perl: \$(perl --version | grep "This is perl" | sed 's/.*(\\(.*\\))/\\1/' | cut -d " " -f1)
    END_VERSIONS
    """
}