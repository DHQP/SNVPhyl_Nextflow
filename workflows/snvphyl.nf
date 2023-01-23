/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSnvphyl.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INDEXING                     } from '../modules/local/indexing'
include { FIND_REPEATS                 } from '../modules/local/find_repeats'
include { SMALT_MAP                    } from '../modules/local/smalt_map'
include { SORT_INDEX_BAMS              } from '../modules/local/sort_index_bams'
include { GENERATE_LINE_1              } from '../modules/local/generate_line_1'
include { VERIFYING_MAP_Q              } from '../modules/local/verifying_map_q'
include { FREEBAYES                    } from '../modules/local/freebayes'
include { FILTER_FREEBAYES             } from '../modules/local/filter_freebayes'
include { BGZIP_FREEBAYES_VCF          } from '../modules/local/bgzip_freebayes_vcf'
include { FREEBAYES_VCF_TO_BCF         } from '../modules/local/freebayes_vcf_to_bcf'
include { MPILEUP                      } from '../modules/local/mpileup'
include { BGZIP_MPILEUP_VCF            } from '../modules/local/bgzip_mpileup_vcf'
include { BCFTOOLS_CALL                } from '../modules/local/bcftools_call'
include { CONSOLIDATE_BCFS             } from '../modules/local/consolidate_bcfs'
include { CONSOLIDATE_FILTERED_DENSITY } from '../modules/local/consolidate_filtered_density'
include { GENERATE_LINE_2              } from '../modules/local/generate_line_2'
include { FILTER_STATS                 } from '../modules/local/filter_stats'
include { VCF2SNV_ALIGNMENT            } from '../modules/local/vcf2snv_alignment'
include { PHYML                        } from '../modules/local/phyml'
include { MAKE_SNV                     } from '../modules/local/make_snv'


/*
========================================================================================
    IMPORT NF-CORE SUBWORKFLOWS
========================================================================================
*/

include { INPUT_CHECK                    } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SNVPHYL {

    ch_versions = Channel.empty() // Used to collect the software versions

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //1. index process takes 1 input channel as a argument
    INDEXING(
        params.refgenome
    )

    //2. find repeats process takes 1 input channel as a argument
    FIND_REPEATS(
        params.refgenome
    )
    ch_versions = ch_versions.mix(FIND_REPEATS.out.versions)

    //3. smalt map process takes 2 input channels as arguments
    SMALT_MAP(
        INPUT_CHECK.out.reads, INDEXING.out.ref_fai, INDEXING.out.ref_sma, INDEXING.out.ref_smi
    )
    ch_versions = ch_versions.mix(SMALT_MAP.out.versions)

    //4. sorting and indexing bam files from smalt process takes 1 input channel as an arguments
    SORT_INDEX_BAMS(
        SMALT_MAP.out.bams
    )
    ch_versions = ch_versions.mix(SORT_INDEX_BAMS.out.versions)

    //5. Generating mapping_quality.txt file
    GENERATE_LINE_1(
        SORT_INDEX_BAMS.out.sorted_bams.collect()
    )
    //ch_versions = ch_versions.mix(GENERATE_LINE_1.out.versions)

    VERIFYING_MAP_Q(
        SORT_INDEX_BAMS.out.sorted_bams.collect(), GENERATE_LINE_1.out.bam_lines_file.splitText()
    )
    ch_versions = ch_versions.mix(VERIFYING_MAP_Q.out.versions)

    //6. freebays variant calling process takes 2 input channels as arguments
    FREEBAYES(
        SORT_INDEX_BAMS.out.sorted_bams_and_sampleID, params.refgenome
    )
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)

    //7. filter freebays variant file process takes 1 input channel as an argument
    FILTER_FREEBAYES(
        FREEBAYES.out.vcf_files
    )
    ch_versions = ch_versions.mix(FILTER_FREEBAYES.out.versions)

    // Zip up the freebayes vcf
    BGZIP_FREEBAYES_VCF(
        FILTER_FREEBAYES.out.filtered_vcf
    )
    ch_versions = ch_versions.mix(BGZIP_FREEBAYES_VCF.out.versions)

    //8. Convert vcf freebays variant file to bcf process takes 1 input channel as an argument
    FREEBAYES_VCF_TO_BCF(
        BGZIP_FREEBAYES_VCF.out.filtered_zipped_vcf
    )
    ch_versions = ch_versions.mix(FREEBAYES_VCF_TO_BCF.out.versions)

    //9. mplileup process takes 1 input channel as argument
    MPILEUP(
        SORT_INDEX_BAMS.out.sorted_bams_and_sampleID, params.refgenome
    )
    ch_versions = ch_versions.mix(MPILEUP.out.versions)

    // Zip up the mpileup vcf
    BGZIP_MPILEUP_VCF(
        MPILEUP.out.mpileup
    )
    ch_versions = ch_versions.mix(BGZIP_MPILEUP_VCF.out.versions)

    //10. mplileup variant calls takes 1 input channel as an argument
    BCFTOOLS_CALL(
        BGZIP_MPILEUP_VCF.out.mpileup_zipped
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CALL.out.versions)

    //Joining channels of multiple outputs
    combined_ch = BCFTOOLS_CALL.out.mpileup_bcf.join(FREEBAYES_VCF_TO_BCF.out.filtered_bcf)
    //11. consolidate variant calling files process takes 2 input channels as arguments
    CONSOLIDATE_BCFS(
        combined_ch
    )
    ch_versions = ch_versions.mix(CONSOLIDATE_BCFS.out.versions)

    // Concat filtered densities to make new invalid_postions
    CONSOLIDATE_FILTERED_DENSITY(
        CONSOLIDATE_BCFS.out.filtered_densities.collect(), FIND_REPEATS.out.repeats_bed_file
    )
    ch_versions = ch_versions.mix(CONSOLIDATE_FILTERED_DENSITY.out.versions)

    // Making string that looks like... this is needed for the next process
    //--consolidate_vcf 2021JQ-00457-WAPHL-M5130-211029=2021JQ-00457-WAPHL-M5130-211029_consolidated.bcf --consolidate_vcf 2021JQ-00459-WAPHL-M5130-211029=2021JQ-00459-WAPHL-M5130-211029_consolidated.bcf --consolidate_vcf 2021JQ-00460-WAPHL-M5130-211029=2021JQ-00460-WAPHL-M5130-211029_consolidated.bcf
    GENERATE_LINE_2(
        CONSOLIDATE_BCFS.out.consolidated_bcfs.collect()
    )
    //ch_versions = ch_versions.mix(GENERATE_LINE_2.out.versions)

    // Get line out of file we just made that has the --consolidate_vcf line...
    //13. consolidate variant calling files process takes 2 input channels as arguments
    VCF2SNV_ALIGNMENT(
        GENERATE_LINE_2.out.consolidation_line.splitText(), CONSOLIDATE_BCFS.out.consolidated_bcfs.collect(), CONSOLIDATE_FILTERED_DENSITY.out.new_invalid_positions, params.refgenome, CONSOLIDATE_BCFS.out.consolidated_bcf_index.collect()
    )
    ch_versions = ch_versions.mix(VCF2SNV_ALIGNMENT.out.versions)

    //14. Filter Stats
    FILTER_STATS(
        VCF2SNV_ALIGNMENT.out.snvTable
    )
    ch_versions = ch_versions.mix(FILTER_STATS.out.versions)

    //15. Using phyml to build tree process takes 1 input channel as an argument
    PHYML(
        VCF2SNV_ALIGNMENT.out.snvAlignment
    )
    ch_versions = ch_versions.mix(PHYML.out.versions)

    //16. Make SNVMatix.tsv
    MAKE_SNV(
        VCF2SNV_ALIGNMENT.out.snvAlignment
    )
    ch_versions = ch_versions.mix(MAKE_SNV.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowSnvphyl.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowSnvphyl.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
