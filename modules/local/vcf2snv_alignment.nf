/* VCF2SNV_ALIGNMENT Call variants */
def VERSION = '1.8.2' // Version information not provided by tool on CLI

process VCF2SNV_ALIGNMENT {
    label 'process_medium'
    container = "staphb/snvphyl-tools:1.8.2"

    input:
    val(consolidate_bcfs)
     path(bcf)
     path(new_invalid_positions)
     path(refgenome)
     path(consolidated_bcf_index)

    output:
    path('snvAlignment.phy'), emit: snvAlignment
    path('vcf2core.tsv'),     emit: vcf2core
    path('snvTable.tsv'),     emit: snvTable
    path("versions.yml"),     emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snvphyl-tools:
        perl: \$(perl --version | grep "This is perl" | sed 's/.*(\\(.*\\))/\\1/' | cut -d " " -f1)
    END_VERSIONS
    """
}