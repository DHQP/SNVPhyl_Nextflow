/* CONSOLIDATED_ALL */
process CONSOLIDATE_FILTERED_DENSITY {
    label 'process_low'
    container = 'quay.io/jvhagey/phoenix:base_v1.0.0'

    input:
    path(filtered_densities)
    path(invalid_positions)

    output:
    path('filtered_density_all.txt'),   emit: filtered_densities
    path('new_invalid_positions.bed'),  emit: new_invalid_positions
    path("versions.yml"),               emit: versions

    script:
    """
    find ./ -name '*_filtered_density.txt' -exec cat {} + > filtered_density_all.txt
    catWrapper.py new_invalid_positions.bed filtered_density_all.txt ${invalid_positions}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}