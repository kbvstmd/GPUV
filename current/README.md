# Current GPUV figures and plotting inputs

This snapshot records the six main and seven supplementary figures assembled on 26 November 2026. The PNGs are previews of the publication PDFs. The [figure manifest](figure_manifest.tsv) records each preview's dimensions and the SHA-256 hash of its source PDF.

## Figures

| Figure | Result |
| --- | --- |
| [1](figures/Fig_1.png) | Catalogue geography, genome lengths, accumulation and reference overlap |
| [2](figures/Fig_2.png) | Diversity and sharing across environments and cities |
| [3](figures/Fig_3.png) | Predicted phage–host associations |
| [4](figures/Fig_4.png) | TerL phylogeny and jumbo phage characteristics |
| [5](figures/Fig_5.png) | Candidate resistance annotations and environmental distribution |
| [6](figures/Fig_6.png) | Eukaryotic DNA virus assignments and environmental detection |
| [S1](figures/Fig_S1.png) | Study workflow and catalogue characteristics |
| [S2](figures/Fig_S2.png) | Environmental diversity and community composition |
| [S3](figures/Fig_S3.png) | Continental and climate-associated variation |
| [S4](figures/Fig_S4.png) | Host assignment and resolution |
| [S5](figures/Fig_S5.png) | TerL branch-length partition and allocation |
| [S6](figures/Fig_S6.png) | Environmental distribution and genomic traits of four jumbo phage clades |
| [S7](figures/Fig_S7.png) | Shared loci, composition and carrier counts per run |

## Redraw from the included tables

Run these commands from the repository root with Python 3.11 or later, R with Cairo support, and Arial installed.

```sh
python -m pip install -r current/requirements.txt
python current/code/plot.py
Rscript current/code/redraw/Fig_S5/plot.R
```

The scripts write PDFs under `current/pdf/redrawn/`. The included data support the redraws listed below. Some redraws use a simplified layout or show a labeled summary of one panel. The previews above show the publication layout.

| Panels | Plotting inputs |
| --- | --- |
| 2b–c and S3b | [Hellinger coordinates, distance-rank quantiles and omnibus tests](data/diversity/) |
| 5a, 5d inset, 5e–f | [Annotation counts, host assignment counts, prevalence and location heatmap](data/Fig_5/) |
| 6a–d | [Family counts, confidence intervals and environmental summaries](data/Fig_6/) |
| S4a–b | [Assignment counts and resolution percentages](data/Fig_S4/) |
| S5a–c | [TerL branch length tables](data/Fig_S5/) |
| S7g–h | [Grouped carrier count distributions and statistical summaries](data/Fig_S7/carrier_counts/) |

Coordinate tables retain group and plot coordinates. The S7 count tables contain group frequencies, which reproduce the boxplot statistics and outliers. See [figure methods](METHODS.md) for the analysis definitions, denominators and limitations.

## Assemble publication PDFs

The [assembly script](code/build_pdfs.py) uses the editable PDF masters from the full figure package. Set `GPUV_FIGURE_MASTERS` to its `data/masters` directory. The repository contains PNG previews and selected plotting inputs. The separate figure package contains the editorial masters and full panel PDFs.

```powershell
$env:GPUV_FIGURE_MASTERS = 'C:/path/to/figure-package/data/masters'
python current/code/build_pdfs.py --out current/pdf
```

Set `ARIAL_BOLD_FILE` to the installed Arial Bold font when rebuilding S3g on a system where the Windows font path is unavailable. Use `--fig S4 S5` for a subset and `--reuse-panels` to assemble edited panel PDFs at their original dimensions.
