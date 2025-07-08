# Dynamic Trans-Omics Mechanisms Underpinning Retinol Tolerance


![Experimental Design](Experimental%20design.png)


### Overview

Retinol remains an essential component in anti-aging skincare; however, a subset of users develop intolerance, characterized by compromised barrier integrity and inflammation. This project presents a comprehensive multi-omics investigation of retinol-induced skin adaptation using a prospective 28-day longitudinal study design, aiming to elucidate the molecular mechanisms of host-microbiome metabolic interactions during retinol tolerance establishment.

### Key Findings

#### 1. Biphasic Skin Phenotype Response
- **Acute Phase:** Increased stratum corneum hydration, reduced sebum secretion, lower skin pH, improved wrinkle metrics
- **Re-equilibration Phase:** Sustained barrier restoration and functional stabilization

#### 2. Functional Microbial Community Remodeling
- 969 microbial species detected with stable overall α-diversity
- **Intolerant Group:** Overrepresentation of *Cutibacterium acnes*
- **Tolerant Group:** Enrichment of protective species (*Sphingomonas hankookensis*, *Acinetobacter johnsonii*)

#### 3. Temporal Metabolomic Fluctuations
- 4,457 metabolites identified showing initial metabolic turbulence followed by partial recovery
- **Early Adverse Phase:** Upregulated lipid/fatty acid metabolism, suppressed TCA cycle and sphingolipid activity
- **Tolerance Establishment:** Enhanced TCA cycle, sphingolipid, ascorbate, and pentose metabolism

### Project Structure

```plaintext
retinol-tolerance-study/
│
├── Experimental design.png # Experimental workflow
├── data/ # Research data
│ ├── intolerant_group/ # Intolerant group data
│ │ ├── kegg_gene.csv # KEGG gene functional annotations
│ │ ├── metabolites.csv # Metabolite abundance data
│ │ ├── mOTUs3_species.txt # Microbial species abundance
│ │ └── skin_trait.csv # Skin phenotype data
│ └── tolerant_group/ # Tolerant group data
│ ├── kegg_gene.csv
│ ├── metabolites.csv
│ ├── mOTUs3_species.txt
│ └── skin_trait.csv
│
├── R_code/ # Analysis scripts
│ ├── part1_Skin_phenotype.R # Skin phenotype analysis (Section 3.1)
│ ├── part2_Microbial_species.R # Microbial species analysis (Section 3.2)
│ ├── part3_Microbial_gene.R # Microbial gene analysis (Section 3.3)
│ ├── part4_Metabolites.R # Metabolite analysis (Section 3.4)
│ └── part5_Multi-omics.R # Multi-omics integration (Section 3.5)
│
├── requirements.txt # R environment dependencies
├── README.md # Project documentation
└── LICENSE # License
```

### Environment Setup

#### R Version
- R version 4.3.2

#### Package Installation

# Install required R packages
```r
install.packages(c(
  "ggpubr", "ggplot2", "vegan", "dplyr", "ape", "cluster",
  "pairwiseAdonis", "reshape2", "tidyr", "randomcoloR",
  "ReporterScore", "KEGGREST", "ggpmisc", "pheatmap",
  "readr", "patchwork", "grid", "ggrepel", "forcats",
  "ggsci", "readxl", "ggcor"
))
```

### Analysis Workflow

#### Step 1: Skin Phenotype Analysis
```r
source("R_code/part1_Skin_phenotype.R")
```
Analyze dynamic changes in skin phenotypes during retinol use, identifying biphasic response patterns.

#### Step 2: Microbial Species Analysis
```r
source("R_code/part2_Microbial_species.R")
```
Assess microbial community α/β diversity changes and identify significantly different key species.

#### Step 3: Microbial Gene Function Analysis
```r
source("R_code/part3_Microbial_gene.R")
```
Perform functional pathway enrichment analysis based on KEGG database to reveal microbial functional remodeling.

#### Step 4: Metabolite Analysis
```r
source("R_code/part4_Metabolites.R")
```
Untargeted metabolomics analysis to identify temporal regulatory patterns of key metabolic pathways.

#### Step 5: Multi-omics Integration
```r
source("R_code/part5_Multi-omics.R")
```
Integrate phenotype, microbial, and metabolomic data to construct host-microbiome interaction networks.


### Citation
If you use the data or code from this project, please cite our research:
```r
[Author Information]. Dynamic Trans-Omics Mechanisms Underpinning Retinol Tolerance: Stage-specific Reconstruction of Skin Barrier Function and Host--Microbiome Metabolic Interactions. [Journal Information], [Year].
```

### Contact
For questions or collaboration inquiries, please contact:
Email: [huangyx23@mails.tsinghua.edu.cn]

### License
This project is licensed under the MIT License - see the LICENSE file for details.
