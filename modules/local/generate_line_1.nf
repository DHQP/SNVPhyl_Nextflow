/* Generate a line for the next process */
process GENERATE_LINE_1 {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.0.0'

    input:
    path(sorted_bams)

    output:
    path("bam_line.txt"), emit: bam_lines_file
    //path("versions.yml"), emit: versions

    shell:
    '''
    #! /bin/bash -l
    count=0 
    for f in *_sorted.bam
    do
    ((count++))
    echo "--bam bam\$count=./$f " | tr -d "\n" >> bam_line.txt
    done
    '''
}