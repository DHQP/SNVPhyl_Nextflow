/* PHYML to make tree */
process PHYML {
    label 'process_low'
    //container = "staphb/phyml:3.3.20220408"
    container = "https://depot.galaxyproject.org/singularity/phyml:3.3.20211231--hee9e358_0"

    input:
    path(snvAlignment_phy)

    output:
    path('phylogeneticTree.newick'),   emit: phylogeneticTree
    path('phylogeneticTreeStats.txt'), emit: phylogeneticTreeStats
    path("versions.yml"),              emit: versions

    script:
    """
    phyml -i ${snvAlignment_phy} --datatype nt --model GTR -v 0.0 -s BEST --ts/tv e --nclasses 4 --alpha e --bootstrap -4 --quiet
    mv snvAlignment.phy_phyml_stats.txt phylogeneticTreeStats.txt
    mv snvAlignment.phy_phyml_tree.txt phylogeneticTree.newick

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phyml: \$( phyml --version | grep -P -o '[0-9]+.[0-9]+.[0-9]+' )
    END_VERSIONS
    """
}