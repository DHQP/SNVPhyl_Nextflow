/* Filter Stats */
process FILTER_STATS {
    label 'process_low'
    container = "staphb/snvphyl-tools:1.8.2"

    input:
    path(snvTable)

    output:
    path('filterStats.txt'), emit: filter_stats
    path("versions.yml"),    emit: versions

    script:
    """
    filter-stats.pl -i ${snvTable} -a > filterStats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snvphyl-tools:
        perl: \$(perl --version | grep "This is perl" | sed 's/.*(\\(.*\\))/\\1/' | cut -d " " -f1)
    END_VERSIONS
    """
}