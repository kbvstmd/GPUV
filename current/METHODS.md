# Methods for the current GPUV figures

## Catalogue, diversity and geographic comparisons

The dataset comprised 12,082 sequencing runs with 372 city labels across 85 countries or regions. The curated catalogue contained 22,494 vOTUs. The study workflow covers read quality control and host filtering, metagenomic assembly, contig annotation, quality assessment, species-level clustering and evidence-based curation. The catalogue and individual figures use different analysis cohorts, including cohorts with complete environmental or geographic labels.

- **Figure 1** contrasts genome-length distributions among the urban and reference collections. The long-genome inset shows counts and median lengths. Accumulation curves use sequencing runs as the sampling unit. The set-intersection panel uses operational genus clusters. Reference overlap uses each reference collection as its denominator.
- **Figures 2 and S2** show Shannon and Simpson diversity or sample-level taxonomic composition by environment. The Hellinger ordination for Fig. 2b–c uses 9,520 runs. ANOSIM and PERMANOVA used 999 permutations. The city and environment sharing networks retain edges with BH-adjusted *P* below 0.05 under a fixed–fixed null. The city classifier measures internal performance after feature selection on the available data.
- **Figure S1** describes the study workflow, CheckV groups and continent-level accumulation and richness. The quality and richness panels use Kruskal–Wallis comparisons.
- **Figure S3** compares Shannon diversity and composition across continents and climate groups. Its continental Hellinger ordination uses 10,106 runs. The latitude panel shows city-level residual-adjusted summaries and mixed-model predictions accounting for environmental category, sequencing depth, project and city. Continent-level occurrence in the set-intersection panel requires detection in at least 10% and at least five samples. The classifier panels describe the internal evaluation cohort.

## Host assignments and phage phylogeny

- **Figures 3 and S4** summarize predicted host associations, assignment resolution and enrichment by environment. The Fig. 3 matrix uses BH-adjusted probabilities by predicted host group. Multiple retained predictions describe assignment uncertainty. Fig. S4 uses its plotted table of 22,846 entries, including unassigned entries.
- **Figures 4, S5 and S6** summarize TerL phylogenies and jumbo phage characteristics. Genomes over 200 kb are used in the jumbo phage length comparison. Figure S5 partitions unrooted branch length by source and allocates each edge contribution across descendant tips for annotation-level summaries. Its branch-length fractions describe representation on the tree. Figure S6 uses Hellinger ordination on branch-positive samples and BH-adjusted comparisons for the selected genomic traits.

## Candidate annotations and detection summaries

- **Figure 5** uses 446 nonredundant candidate references for class summaries, based on 452 loci on 333 contigs. The d inset includes 316 assignments among 333 candidate-bearing contigs and 7,913 assignments among 22,161 contigs with no candidate detected. The environmental panel shows run-weighted frequencies with Wilson intervals alongside equally weighted project estimates. The geographic panel shows category-level sample percentages at the displayed locations. Gene-reference matches are evidence of a detected annotation; attribution to a specific carrier requires additional context.
- **Figure 6** summarizes 149 retained vOTUs from five families, including 41 associated with human-virus references. Sample detection uses a 0.01% abundance threshold for a represented member. Figure 6b uses 4,999 project bootstrap resamples. Figure 6c uses Wilson intervals and positive-only abundance summaries. Figure 6d uses Cochran–Mantel–Haenszel comparisons stratified by project, country and city, with BH correction across the planned test family. These labels describe catalogue and reference assignments. The trees in 6f–g refer to their assembly-level comparison sets.

## Virus–bacterium shared loci

**Figure S7a–f** reports the discovery-set filtering, functional categories, paired protein comparisons and composition distances. The paired GC and tetranucleotide comparisons share the displayed 307-locus cohort. The panels summarize association and similarity evidence in the discovery assemblies.

**Figure S7g–h** counts the known shared-locus carrier vOTUs represented in the filtered abundance matrix. Among 6,895 eligible runs, each represented carrier vOTU contributed at most once per run if its abundance was greater than zero. The displayed all-run analyses retain zeros and use 6,301 runs with environmental labels and 6,815 with continental labels. The boxplots summarize raw counts and the overall test annotations use 9,999 permutations of the group-explained sum of squares for log(1 + count). The included tables distinguish the all-run and positive-run cohorts and the unadjusted and stratified analyses. The count endpoint is the carrier vOTU. Shared-locus boundaries were established in the discovery assemblies.

## Plot reproduction

The repository tables are plotting inputs and statistical summaries. The R and Python redraws use them directly and leave the archived model estimates fixed. The publication layout is assembled separately from the editable PDF masters. See the [redraw and assembly instructions](README.md).
