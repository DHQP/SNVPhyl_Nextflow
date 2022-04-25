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


---
## Notices and Disclaimers

### Public Domain
This repository constitutes a work of the United States Government and is not subject to domestic copyright protection under 17 USC § 105. This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/). All contributions to this repository will be released under the CC0 dedication. By submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

### License

Unless otherwise specified, the repository utilizes code licensed under the terms of the Apache Software License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under the terms of the Apache Software License version 2, or (at your option) any later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

Any source code forked from other open source projects will inherit its license.

### Privacy

This repository contains only non-sensitive, publicly available data and information. All material and community participation is covered by the [Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md) and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md). For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

### Contributing

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo) and submitting a pull request. (If you are new to GitHub, you might start with a [basic tutorial](https://help.github.com/articles/set-up-git)). By contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users under the terms of the [Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or later.

All comments, messages, pull requests, and other submissions received through CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

### Records

This repository is not a source of government records, but is a copy to increase collaboration and collaborative potential. All government records will be published through the [CDC web site](http://www.cdc.gov).

