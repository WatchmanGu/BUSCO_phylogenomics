This folder contains scripts to extract BUSCO genes from several BUSCO run folders (! having used the same BUSCO dataset) and produce a species phylogeny.

- `Dockerfile`: to build a Docker container able to run the script. (see https://www.docker.com/what-docker) (BUSCO is included if you need it. But not required to run the script. Use python3 within the docker)
- `run.sh`: the main script to proceed to the analysis. The final output is final.nwk
- `extract_buscos_phylo.py`: the main script to identify and extract BUSCO sequences
- `fetch_best_sequence.py`: an additional script required by the previous one
- `superalignment.py`: it fuses all single gene alignement together
- `config.ini`: BUSCO configuration to be used within the container

You will need to adapt several scripts for your own use, each files are commented. For feedbacks and help, please use the main BUSCO gitlab project: https://gitlab.com/ezlab/busco

Copyright (c) 2016-2018, Evgeny Zdobnov (ez@ezlab.org)
Licensed under the MIT license

