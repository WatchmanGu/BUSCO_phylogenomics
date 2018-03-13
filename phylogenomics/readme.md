This folder contains scripts to extract BUSCO genes from several BUSCO run folders (! having used the same BUSCO dataset) and produce a species phylogeny.

- `Dockerfile`: to build a Docker container able to run the script. (see https://www.docker.com/what-docker) (BUSCO is included if you need it. But not required to run the script. Use python3 within the docker to run BUSCO)
- `run.sh`: the main script to proceed to the analysis. The final output is final.nwk
- `extract_buscos_phylo.py`: the main script to identify and extract BUSCO sequences
- `fetch_best_sequence.py`: an additional script required by the previous one
- `superalignment.py`: it fuses all single gene alignement together
- `config.ini`: BUSCO configuration to be used within the container

You will need to adapt several scripts for your own use, each files are commented. For feedbacks and help, please use the main BUSCO gitlab project: https://gitlab.com/ezlab/busco

BUSCO main gitlab repository: https://gitlab.com/ezlab/busco

**How to cite BUSCO**

*BUSCO applications from quality assessments to gene prediction and phylogenomics.*
Robert M. Waterhouse, Mathieu Seppey, Felipe A. Simão, Mose Manni, Panagiotis Ioannidis, Guennadi Klioutchnikov, Evgenia V. Kriventseva, and Evgeny M. Zdobnov
*Mol Biol Evol*, published online Dec 6, 2017 
doi: 10.1093/molbev/msx319 

*BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs.*
Felipe A. Simão, Robert M. Waterhouse, Panagiotis Ioannidis, Evgenia V. Kriventseva, and Evgeny M. Zdobnov
*Bioinformatics*, published online June 9, 2015 
doi: 10.1093/bioinformatics/btv351

Copyright (c) 2016-2018, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license. See LICENSE.md file.
