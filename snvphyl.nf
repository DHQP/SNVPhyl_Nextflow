//snvphyl.nf

/*
========================================================================================
   SNVPhyl Nextflow Workflow
========================================================================================
   Github   : https://git.biotech.cdc.gov/mmb/snvphyl/
   Contact  : Jill Hagey, qpk9@cdc.gov

This script was based on SNVPhyl https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5628696/
and it has a github: https://github.com/phac-nml/snvphyl-galaxy

----------------------------------------------------------------------------------------

*/

nextflow.enable.dsl=2

// Color of text
ANSI_RESET = "\u001B[0m"
ANSI_RED = "\u001B[31m"
ANSI_GREEN = "\u001B[32m"
ANSI_YELLOW = "\u001B[33m";

/*
========================================================================================
   Create Parameters
========================================================================================
*/

// Initialize required parameters
params.outdir = "./results"
params.refgenome = "reference.fasta"
params.input_reads = "./FASTQs/"
params.window_size = "11"
params.density_threshold = "2"

/*
========================================================================================
   Printing OutPipeline Details
========================================================================================
*/

println(ANSI_GREEN + """\
         \n
         =====================================================================================================================================
                                                     S N V P H Y L - N F   P I P E L I N E
         =====================================================================================================================================
         Currently the SNVPhyl workflow is for use with Illumina paired end samples
         Author: Jill Hagey
         Email: qpk9@cdc.gov
         Version: v.1.0.0\n
         Based on original SNVPhyl Pipeline by Aaron Petkau et al. https://github.com/phac-nml/snvphyl-galaxy\n
         """.stripIndent()
         + ANSI_RESET)
println(ANSI_GREEN + """\
         =====================================================================================================================================
                                                       I N P U T   P A R A M E T E R S
         =====================================================================================================================================
         --refgenome            Reference Genome [Default "reference.fasta"]                                    : ${params.refgenome}
         --input_reads          Input folder with reads [Default "./FASTQs/"]                                   : ${params.input_reads}
         --outdir               Output Directory [Default "./results"]                                          : ${params.outdir}
         --window_size          Window size for identifying high-density SNV regions [Default "11"]             : ${params.window_size}
         --density_threshold    SNV threshold for identifying high-density SNV regions [Default "2"]            : ${params.density_threshold}
         """.stripIndent()
         + ANSI_RESET)

/*
========================================================================================
   Workflow
========================================================================================
*/

workflow {
    
    // Create channel from path for Reference Genome
    ref_ch = Channel.fromPath(params.refgenome , checkIfExists: true )  
    // Create channel from file-pairs for Input Fastq files
    reads_ch = Channel.fromFilePairs("${params.input_reads}/*_R{1,2}*.{fastq,fastq.gz,fq,fq.gz}", checkIfExists: true )

    //1. index process takes 1 input channel as a argument
    INDEXING(ref_ch)
    //2. find repeats process takes 1 input channel as a argument
    FIND_REPEATS(ref_ch)

    //3. smalt map process takes 2 input channels as arguments
    SMALT_MAP(INDEXING.out.indexes.combine(reads_ch))

    //4. sorting and indexing bam files from smalt process takes 1 input channel as an arguments
    SORT_INDEX_BAMS(SMALT_MAP.out)

    //5. Generating mapping_quality.txt file
    GENERATE_LINE_1(SORT_INDEX_BAMS.out.sorted_bams.collect())
    VERIFYING_MAP_Q(SORT_INDEX_BAMS.out.sorted_bams.collect(), GENERATE_LINE_1.out.splitText())

    //6. freebays variant calling process takes 2 input channels as arguments
    FREEBAYES(SORT_INDEX_BAMS.out.sorted_bams_and_sampleID.combine(ref_ch))
    //7. filter freebays variant file process takes 1 input channel as an argument
    FILTER_FREEBAYES(FREEBAYES.out.vcf_files)
    // Zip up the freebayes vcf
    BGZIP_FREEBAYES_VCF(FILTER_FREEBAYES.out)
    //8. Convert vcf freebays variant file to bcf process takes 1 input channel as an argument
    FREEBAYES_VCF_TO_BCF(BGZIP_FREEBAYES_VCF.out)

    //9. mplileup process takes 1 input channel as argument
    MPILEUP(SORT_INDEX_BAMS.out.sorted_bams_and_sampleID.combine(ref_ch))
    // Zip up the mpileup vcf
    BGZIP_MPILEUP_VCF(MPILEUP.out)
    //10. mplileup variant calls takes 1 input channel as an argument
    BCFTOOLS_CALL(BGZIP_MPILEUP_VCF.out)

    //Joining channels of multiple outputs
    combined_ch = BCFTOOLS_CALL.out.mpileup_bcf.join(FREEBAYES_VCF_TO_BCF.out)
    //11. consolidate variant calling files process takes 2 input channels as arguments
    CONSOLIDATE_BCFS(combined_ch)

    // Concat filtered densities to make new invalid_postions
    catwrap_ch = Channel.fromPath("$PWD/catWrapper.py" , checkIfExists: true ) 
    CONSOLIDATE_FILTERED_DENSITY(CONSOLIDATE_BCFS.out.filtered_densities.collect(), FIND_REPEATS.out, catwrap_ch)

    // Making string that looks like... this is needed for the next process
    //--consolidate_vcf 2021JQ-00457-WAPHL-M5130-211029=2021JQ-00457-WAPHL-M5130-211029_consolidated.bcf --consolidate_vcf 2021JQ-00459-WAPHL-M5130-211029=2021JQ-00459-WAPHL-M5130-211029_consolidated.bcf --consolidate_vcf 2021JQ-00460-WAPHL-M5130-211029=2021JQ-00460-WAPHL-M5130-211029_consolidated.bcf
    GENERATE_LINE_2(CONSOLIDATE_BCFS.out.consolidated_bcfs.collect())

    // Get line out of file we just made that has the --consolidate_vcf line...
    //13. consolidate variant calling files process takes 2 input channels as arguments
    VCF2SNV_ALIGNMENT(GENERATE_LINE_2.out.splitText(), CONSOLIDATE_BCFS.out.consolidated_bcfs.collect(), CONSOLIDATE_FILTERED_DENSITY.out.new_invalid_positions, ref_ch, CONSOLIDATE_BCFS.out.consolidated_bcf_index.collect())
    //14. Filter Stats
    FILTER_STATS(VCF2SNV_ALIGNMENT.out.snvTable)
    //15. Using phyml to build tree process takes 1 input channel as an argument
    PHYML(VCF2SNV_ALIGNMENT.out.snvAlignment)
    //16. Make SNVMatix.tsv
    MAKE_SNV(VCF2SNV_ALIGNMENT.out.snvAlignment)
}

/*
========================================================================================
   Program info
========================================================================================
*/

// This is where the results will be
println(ANSI_YELLOW + "\nPipeline Starting! \nThe files and directory for results are in a folder called " + params.outdir + "\n" + ANSI_RESET)

/*
========================================================================================
   Processes
========================================================================================
*/

process INDEXING {
    tag {"INDEXING ${refgenome}"}

    //publishDir "${params.outdir}", mode: 'copy'

    input:
    path(refgenome)

    output:
    path '*.{sma,smi,fai}', emit: indexes

    script:
    """
    REF_BASENAME=\$(basename ${refgenome} .fasta)
    smalt index -k 13 -s 6 \${REF_BASENAME} ${refgenome}
    samtools faidx ${refgenome}
    """
}
/* Find repeats and create invaild_positions.bed */
process FIND_REPEATS {
    tag {"FIND_REPEATS ${refgenome}"}
    label 'process_snvphyl_low'

    //publishDir "${params.outdir}", mode: 'copy'
    
    input:
     path(refgenome)

    output:
     path( "invalid_positions.bed" )

    script:
    """
    find-repeats.pl ${refgenome} --min-length 150 --min-pid 90 > invalid_positions.bed
    """
}
/* Map reads to reference genome & create BAM file */
process SMALT_MAP {
    tag {"SMALT_MAP ${sample_id}"}

    //publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple path(reffai), path(refsma), path(refsmi), val(sample_id), file(reads)

    output:
    tuple val(sample_id), path( "${sample_id}.bam" )

    script:
    """
    REFNAME=\$(basename ${refsma} .sma)
    smalt map -f bam -n 4 -l pe -i 1000 -j 20 -r 1 -y 0.5 -o ${sample_id}.bam \${REFNAME} ${reads.get(0)} ${reads.get(1)}
    """
}
/* Create BAMs and Sort */
process SORT_INDEX_BAMS {
    tag {"SORT_INDEX_BAMS ${sample_id}"}

    //publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bams)

    output:
    tuple val(sample_id), path( "${sample_id}_sorted.bam" ), emit: sorted_bams_and_sampleID
    path( "${sample_id}_sorted.bam" ), emit: sorted_bams

    script:
    """
    samtools sort -O bam -o ${sample_id}_sorted.bam ${bams}
    samtools index ${sample_id}_sorted.bam
    """
}
/* Generate a line for the next process */
process GENERATE_LINE_1 {
    label 'process_gen_line'

    input:
     path(sorted_bams)

    output:
     file("bam_line.txt")

    shell:
    '''
    #! /bin/bash -l
    count=0 
    for f in *_sorted.bam
    do
    ((count++))
    echo "--bam bam$count=./$f " | tr -d "\n" >> bam_line.txt
    done
    '''

}
/* Verifying Mapping Quality */
process VERIFYING_MAP_Q {
    tag {"SORT_INDEX_BAMS"}

    publishDir "${params.outdir}", mode: 'copy'

    input:
     path(sorted_bams)
     val(bam_line)

    output:
     file("mappingQuality.txt")

    script:
    """
    verify_mapping_quality.pl -c 4 --min-depth 10 --min-map 80 --output mappingQuality.txt ${bam_line} 
    """
}
/* Freebayes variant calling */
process FREEBAYES {
    tag {"FREEBAYES ${sample_id}"}

    //publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bams), path(refgenome)

    output:
    tuple val(sample_id), path( "${sample_id}_freebayes.vcf" ), emit: vcf_files

    script:
    """
    freebayes --bam ${sorted_bams} --ploidy 1 --fasta-reference ${refgenome} --vcf ${sample_id}_freebayes.vcf 
    """
}
/* Mpileup */
process MPILEUP {
    tag {"MPILEUP ${sample_id}"}

    //publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bams), path(refgenome)

    output:
    tuple val(sample_id), path( "${sample_id}_mpileup.vcf" )

    script:
    """
    bcftools mpileup --threads 4 --fasta-ref ${refgenome} -A -B -C 0 -d 1024 -q 0 -Q 0 --output-type v -I --output ${sample_id}_mpileup.vcf ${sorted_bams}
    """
}
/* Zip mpileup vcf*/
process BGZIP_MPILEUP_VCF {
    tag {"BGZIP_MPILEUP_VCF ${sample_id}"}
    label 'process_bgzip'

    input:
     tuple val(sample_id), path(mpileup_vcf)

    output:
     tuple val(sample_id), path("${sample_id}_mpileup.vcf.gz")

    script:
    """
    bgzip -f ${mpileup_vcf}
    """
}
/* Bcftools call  */
process BCFTOOLS_CALL {
    tag {"BCFTOOLS_CALL ${sample_id}"}
    label 'process_bcf'

    //publishDir "${params.outdir}", mode: 'copy'

    input:
     tuple val(sample_id), path(mpileup_vcf_gz)

    output:
     tuple val(sample_id), path("${sample_id}_mpileup.bcf"), emit: mpileup_bcf

    script:
    """
    bcftools index -f ${mpileup_vcf_gz}
    bcftools call --ploidy 1 --threads 4 --output ${sample_id}_mpileup.bcf --output-type b --consensus-caller ${mpileup_vcf_gz}
    """
}
/* Filter freebayes vcf */
process FILTER_FREEBAYES {
    tag {"FILTER_FREEBAYES ${sample_id}"}
    label 'process_snvphyl_low'

    //publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(freebayes_vcf)

    output:
    tuple val(sample_id), path( "${sample_id}_freebayes_filtered.vcf" )

    script:
    """
    filterVcf.pl --noindels ${freebayes_vcf} -o ${sample_id}_freebayes_filtered.vcf
    """
}
/* Zip freebayes vcf*/
process BGZIP_FREEBAYES_VCF {
    tag {"BGZIP_FREEBAYES_VCF ${sample_id}"}
    label 'process_bgzip'

    input:
     tuple val(sample_id), path(freebayes_filtered_vcf)

    output:
     tuple val(sample_id), path("${sample_id}_freebayes_filtered.vcf.gz")

    script:
    """
    bgzip -f ${freebayes_filtered_vcf}
    """
}
/* Filtered freebayes vcf to bcf */
process FREEBAYES_VCF_TO_BCF {
    tag {"FREEBAYES_VCF_TO_BCF ${sample_id}"}
    label 'process_bcf'

    input:
     tuple val(sample_id), path(freebayes_filtered_vcf_gz)

    output:
     tuple val(sample_id), path( "${sample_id}_freebayes_filtered.bcf" ), path( "${sample_id}_freebayes_filtered.bcf.csi" )

    script:
    """
    bcftools index -f ${freebayes_filtered_vcf_gz}
    bcftools view --output-type b --output-file ${sample_id}_freebayes_filtered.bcf ${freebayes_filtered_vcf_gz}
    bcftools index -f ${sample_id}_freebayes_filtered.bcf
    """
}
/* Consolidate bcfs */
process CONSOLIDATE_BCFS {
    tag {"CONSOLIDATE_BCF ${sample_id}"}

    publishDir "${params.outdir}", mode: 'copy', pattern: "Log_File.txt"

    input:
     tuple val(sample_id), path(mpileup_bcf), path(freebayes_filtered_bcf), path(freebayes_filtered_csi)

    output:
     val(sample_id), emit: sample_id
     path( "${sample_id}_consolidated.bcf" ), emit: consolidated_bcfs
     path( "${sample_id}_consolidated.vcf" ), emit: consolidated_vcfs
     path( "${sample_id}_consolidated.bcf.csi" ), emit: consolidated_bcf_index
     path( "${sample_id}_filtered_density.txt" ), emit: filtered_densities

    script:
    """
    consolidate_vcfs.pl --coverage-cutoff 10 --min-mean-mapping 30 --snv-abundance-ratio 0.75 --vcfsplit ${freebayes_filtered_bcf} --mpileup ${mpileup_bcf} --filtered-density-out ${sample_id}_filtered_density.txt --window-size ${params.window_size} --density-threshold ${params.density_threshold} -o ${sample_id}_consolidated.bcf > ${sample_id}_consolidated.vcf
    bcftools index -f ${sample_id}_consolidated.bcf
    """
}
/* CONSOLIDATED_ALL */
process CONSOLIDATE_FILTERED_DENSITY {
    label 'process_gen_line'

    input:
     path(filtered_densities)
     path(invalid_positions)
     path(catWrapper)

    output:
     path 'filtered_density_all.txt', emit: filtered_densities
     path 'new_invalid_positions.bed', emit: new_invalid_positions

    script:
    """
    find ./ -name '*_filtered_density.txt' -exec cat {} + > filtered_density_all.txt
    python ${catWrapper} new_invalid_positions.bed filtered_density_all.txt ${invalid_positions}
    """
}
/* Generate a line for the next process */
process GENERATE_LINE_2 {
    label 'process_gen_line'

    input:
     path(consolidated_bcf)

    output:
     file("consolidation_line.txt")

    shell:
    '''
    #! /bin/bash -l
    for f in *_consolidated.bcf
    do
    fname=$(basename $f _consolidated.bcf)
    echo "--consolidate_vcf $fname=$f " | tr -d "\n" >> consolidation_line.txt
    done
    '''
}
/* VCF2SNV_ALIGNMENT Call variants */
process VCF2SNV_ALIGNMENT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
     val(consolidate_bcfs)
     path(bcf)
     path(new_invalid_positions)
     path(refgenome)
     path(consolidated_bcf_index)

    output:
     path 'snvAlignment.phy', emit: snvAlignment
     path 'vcf2core.tsv', emit: vcf2core
     path 'snvTable.tsv', emit: snvTable

    script:
    """
    vcf2snv_alignment.pl --reference reference --invalid-pos ${new_invalid_positions} --format fasta --format phylip --numcpus 4 --output-base snvalign --fasta ${refgenome} ${consolidate_bcfs} 
    mv snvalign-positions.tsv snvTable.tsv
    mv snvalign-stats.csv vcf2core.tsv
    if [[ -f snvalign.phy ]]; then
        mv snvalign.phy snvAlignment.phy
        sed -i "s/'//" snvAlignment.phy
        sed -i "s/'//" snvAlignment.phy
    else
        touch snvAlignment.phy
    fi    
    """
}
/* Filter Stats */
process FILTER_STATS {
    label 'process_snvphyl_low'

    publishDir "${params.outdir}", mode: 'copy'

    input:
     path(snvTable)

    output:
     path('filterStats.txt')

    script:
    """
    filter-stats.pl -i ${snvTable} -a > filterStats.txt
    """
}
/* PHYML to make tree */
process PHYML {
    publishDir "${params.outdir}", mode: 'copy'

    input:
     path(snvAlignment_phy)

    output:
     path 'phylogeneticTree.newick', emit: phylogeneticTree
     path 'phylogeneticTreeStats.txt', emit: phylogeneticTreeStats

    script:
    """
    phyml -i ${snvAlignment_phy} --datatype nt --model GTR -v 0.0 -s BEST --ts/tv e --nclasses 4 --alpha e --bootstrap -4 --quiet
    mv snvAlignment.phy_phyml_stats.txt phylogeneticTreeStats.txt
    mv snvAlignment.phy_phyml_tree.txt phylogeneticTree.newick
    """
}
/* Making SNV matrix */
process MAKE_SNV {
    label 'process_snvphyl_low'

    publishDir "${params.outdir}", mode: 'copy'

    input:
     path(snvAlignment_phy)

    output:
     path 'snvMatrix.tsv', emit: snvMatrix

    script:
    """
    snv_matrix.pl ${snvAlignment_phy} -o snvMatrix.tsv
    """
}

/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
       Pipeline Execution Summary
       ---------------------------
       Success     : ${workflow.success}
       workDir     : ${workflow.workDir}
       exit status : ${workflow.exitStatus}
       """ : """
       Failed: ${workflow.errorReport}
       exit status : ${workflow.exitStatus}
       """
   )
       //Completed at: ${workflow.complete}
       //Duration    : ${workflow.duration}
}

/*
========================================================================================
   THE END
========================================================================================
*/
