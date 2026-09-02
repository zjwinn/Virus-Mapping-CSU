# Identification of Novel Loci in Winter Wheat Elucidates a Potential Shift in Viral Disease Complexes in the U.S. High Plains

This repository contains the analysis scripts, phenotypic data, and supplementary materials supporting the manuscript:

> **Identification of Novel Loci in Winter Wheat Elucidates a Potential Shift in Viral Disease Complexes in the U.S. High Plains**  
> Zachary J. Winn, Mik Hammers, Emily Hudson-Arns, Matthew West, Nichole Barber, Robyn Roberts, Mary Guttieri, Guorong Zhang, Scott D. Haley, and R. Esten Mason.  
> *Submitted / Under Review (2026).*

---

## Overview

Wheat (*Triticum aestivum* L.) production in the U.S. High Plains faces major challenges from viral diseases vectored by the wheat curl mite (*Aceria tosichella*). This repository documents a decade-long study (2017-2025) capturing an epidemiological shift from wheat streak mosaic virus (WSMV) to Triticum mosaic virus (TriMV), the breakdown of historical resistance loci (*Wsm2* and *Cmc4*), and the discovery of novel genetic resistance on chromosome 4AL.

---

## Repository Structure

* `main_analysis.Rmd` / `main_analysis.R`: Complete quantitative genetic pipeline including spatial two-stage mixed linear modeling, inverse normal transformations, and Q+K genome-wide association studies (GWAS) for continuous and binomial traits.
* `predictive_haplotyping.R`: Linkage disequilibrium (LD) decay, Silhouette/K-means clustering, haplotype group BLUE estimations for the 4AL locus, and cultivar classifications.
* `visualizing_PCA_for_publication.R`: Population structure analysis (PCA) evaluating market classes and breeding methodologies across PreDCS and PostDCS panels.
* `virus_frequency_images.R`: Spatial and temporal tracking of diagnostic virus frequencies across Colorado and adjacent High Plains regions.
* `manuscript/`: Final manuscript figures (Figures 1-9), summary tables (Tables 1-5), and supplementary information files (SI 1-7).
* `multiple_environment_adjustment/`: Across-environment best linear unbiased estimates (BLUEs) for PreDCS and PostDCS mapping panels.
* `PreDCS_GWAS/` & `PostDCS_GWAS/`: Summary statistics, Manhattan plots, and representative QTN associations.

---

## Citation & Contact

If you use code or data from this repository, please cite the corresponding paper upon publication.

For questions or inquiries, please contact:
* **Zachary J. Winn** - [zwinn@outlook.com](mailto:zwinn@outlook.com)
* **R. Esten Mason** - [esten.mason@colostate.edu](mailto:esten.mason@colostate.edu)
