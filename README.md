# Global Atlas of the Pan-Urban Virome(GPUV)
The Global Atlas of the Pan-Urban Virome from over 12,000 metagenomes from 280 cities across 86 countries or regions. Our analysis cataloged 23,800 viral species, over 94% of which were previously uncharacterized. 
### Workflow
Here is the Workflow of GPUV
<center><img src="img/data_workflow.png" width="90%"></center>
The following is a list of the primary software applications:
- [fastp](https://github.com/OpenGene/fastp), and [bowtie2](https://github.com/BenLangmead/bowtie2) were used to removing low-quality sequences and trimming adapters sequences from host such as human reference genome.
- [metaSPAdes](https://github.com/ablab/spades), and [megahit](https://github.com/voutcn/megahit) were used reads assembleing.
- [geNomad](https://portal.nersc.gov/genomad/),  and [CheckV](https://bitbucket.org/berkeleylab/checkv) were used to remove sequences from cellular organisms and plasmids, as necessary.
