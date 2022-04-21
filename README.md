# SNVPhyl_Nextflow

This Nextflow version of SNVPhyl is based on the original SNVPhyl pipeline written by Aaron Petkau. For documentation on the pipeline and how it works see the following links. 

- [Galaxy Script](https://github.com/phac-nml/snvphyl-galaxy/blob/development/docs/workflows/SNVPhyl/1.0.1/snvphyl-workflow-1.0.1.ga)
- [Workflow Image](https://snvphyl.readthedocs.io/en/latest/images/snvphyl-overview-galaxy.png)
- [SNVPhyl Paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5628696/)
- [SNVPhyl Documenation](https://phac-nml.github.io/irida-documentation/administrator/galaxy/pipelines/phylogenomics/)

## Workflow

![SNVPhyl DAG](https://github.com/DHQP/SNVPhyl_Nextflow/blob/main/SNVPhyl_DAG.png)

## Install and Dependencies

This version of the software was run with:
- `Nextflow version 21.04.3 build 5560 created 21-07-2021 15:09 UTC (11:09 EDT)`  
- `singularity-ce version 3.8.0`  

Once you have nextflow and singularity installed then get this repo with either:

`git clone https://github.com/DHQP/SNVPhyl_Nextflow.git`

or 

```
wget https://github.com/DHQP/SNVPhyl_Nextflow/releases/download/1.0.0/SNVPhyl_Nextflow.tar.gz
tar -xvzf SNVPhyl_Nextflow.tar.gz
```

## Running Pipeline

To run the pipeline do the following:

```
nextflow run snvphyl.nf --outdir ./results -c snvphyl.config --refgenome reference.fasta
```

**Make sure to pick an appropriate reference fasta file for your samples!** The one included in this repository is just for example purposes. 

### Inputs  

```
==============================================================================================================================
                                              I N P U T   P A R A M E T E R S
==============================================================================================================================
--refgenome            Reference Genome [Default "reference.fasta"]                                    : reference.fasta
--input_reads          Input folder with reads [Default "./FASTQs/"]                                   : ./FASTQs/
--outdir               Output Directory [Default "./results"]                                          : ./test
--window_size          Window size for identifying high-density SNV regions [Default "11"]             : 11
--density_threshold    SNV threshold for identifying high-density SNV regions [Default "2"]            : 2

```  

The config file is currently run locally with containers using Singularity. Alternatively you can adjust the config file to run with a [different executors](https://www.nextflow.io/docs/latest/executor.html) or [Docker](https://www.nextflow.io/docs/latest/docker.html). 

Prior to running you will need to change the `singularity.cacheDir = '$PATH/Singularity_Containers'` in [line 25](https://github.com/DHQP/SNVPhyl_Nextflow/blob/d400b20b4c147f11b3c1f456fcef83215bf16b56/snvphyl.config#L25) of the config file to the location of your singularity cache directory. 

### Output  

The following files will be output by SNVPhyl and are explained in the [original SNVPhyl Documentation](https://snvphyl.readthedocs.io/en/latest/user/output/).
1. filterStats.txt   
2. phylogeneticTree.newick    
3. phylogeneticTreeStats.txt
5. snvAlignment.phy  
6. snvTable.tsv
7. snvMatrix.tsv    
8. vcf2core.tsv

Additionally, two other files are created by this workflow. 

10. Log_File.txt -- Contains information on the software versions used for this run.
11. report.html  -- This is the nextflow report that will tell you the CPU usage, time to run each job and RAM requirements. 

