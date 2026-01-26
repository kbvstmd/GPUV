Viromics Bioinformatics Workflow
This repository contains the bioinformatics pipeline for processing metagenomic sequencing data, focusing on the identification, characterization, and analysis of viral sequences.

Workflow Overview
Preprocessing: Quality control and host sequence depletion.

Assembly: Metagenome assembly of filtered reads.

Identification: Viral contig discovery and quality assessment.

Quantification: vOTU clustering and abundance estimation.

Annotation: Gene prediction and host/lifestyle inference.

Phylogenetics: Evolutionary analysis of marker genes (e.g., TerL).

1. Quality Control and Host Removal
Raw paired-end reads are quality-controlled using fastp. Host-derived reads are removed by mapping to the host reference genome using bowtie2; only unmapped reads are retained.

Bash
# Quality control
fastp -i R1.fq.gz -I R2.fq.gz \
      -o clean_R1.fq.gz -O clean_R2.fq.gz

# Host removal (e.g., human, mouse, or plant genome)
bowtie2 -x host_db \
        -1 clean_R1.fq.gz -2 clean_R2.fq.gz \
        --very-sensitive --un-conc-gz dehost.fq.gz \
        -S host.sam
2. Metagenome Assembly
Host-filtered reads are assembled into contigs using MEGAHIT or metaSPAdes.

Bash
# Option A: MEGAHIT (Memory efficient)
megahit -1 dehost.1.fq.gz -2 dehost.2.fq.gz -o megahit_out

# Option B: metaSPAdes (Optimized for complex communities)
metaspades.py -1 dehost.1.fq.gz -2 dehost.2.fq.gz -o metaspades_out
3. Viral Contig Identification
Viral sequences are identified from the assembly using a combination of machine learning and marker-based tools.

Bash
# geNomad identification
genomad end-to-end contigs.fna genomad_out genomad_db

# DeepVirFinder identification
python DeepVirFinder.py \
       -i contigs.fna \
       -o dvf_out \
       -m DVF_model
4. Viral Genome Quality Assessment
The quality and completeness of identified viral contigs are assessed using CheckV.

Bash
checkv end_to_end viral_contigs.fna checkv_out -d checkv_db
5. vOTU Clustering
Viral contigs are clustered into Viral Operational Taxonomic Units (vOTUs) based on Average Nucleotide Identity (ANI) thresholds.

Bash
python cluster.py \
  --fna all_contigs.fna \
  --ani ani.tsv \
  --out clusters.tsv \
  --min_ani 95 \
  --min_qcov 0 \
  --min_tcov 85
6. Gene Prediction
Protein-coding genes are predicted using Prodigal in metagenomic mode (-p meta).

Bash
prodigal -p meta \
         -i vOTU.fna \
         -a vOTU_prot.faa \
         -d vOTU_gene.fna \
         -f gff \
         -o genes.gff
7. Relative Abundance Estimation
Quality-controlled reads are mapped back to vOTU sequences using BWA-MEM to estimate abundance based on coverage.

Bash
# Indexing and mapping
bwa index vOTU.fna
bwa mem vOTU.fna dehost.1.fq.gz dehost.2.fq.gz > vOTU.sam

# Post-processing
samtools view -bS vOTU.sam | samtools sort -o vOTU.sorted.bam
samtools index vOTU.sorted.bam

# Coverage calculation
msamtools coverage vOTU.sorted.bam > vOTU_abundance.tsv
8. Host and Lifestyle Prediction
Inference of putative bacterial/archaeal hosts and viral lifestyles.

Bash
# Host Prediction (iPHoP)
iphop predict --fa vOTU.fna --out iphop_out --db iphop_db --threads 16

# Lifestyle Prediction (BACPHLIP)
bacphlip -i vOTU.fna -o bacphlip_out.tsv -t 16
9. Phylogenetic Analysis
Evolutionary analysis using the Large Terminase Subunit (TerL) as a phylogenetic marker.

Bash
# Multiple sequence alignment (MUSCLE)
muscle -align terl.faa -output terl.aln.faa

# Alignment trimming (TrimAL)
trimal -in terl.aln.faa -out terl.aln.trim.faa -automated1

# Phylogenetic tree inference (IQ-TREE)
iqtree2 -s terl.aln.trim.faa \
        -m MFP \
        -B 1000 \
        -T AUTO \
        -pre terl_iqtree
The resulting *.treefile (Newick format) can be visualized using iTOL or FigTree.
