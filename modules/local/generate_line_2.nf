/* Generate a line for the next process */
process GENERATE_LINE_2 {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.0.0'

    input:
    path(consolidated_bcf)

    output:
    path("consolidation_line.txt"), emit: consolidation_line
    //path("versions.yml"),           emit: versions

    shell:
    '''
    #! /bin/bash -l
    for f in *_consolidated.bcf
    do
    fname=\$(basename $f _consolidated.bcf)
    echo "--consolidate_vcf \$fname=$f " | tr -d "\n" >> consolidation_line.txt
    done
    '''
}