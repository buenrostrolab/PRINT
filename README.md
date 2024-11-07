# PRINT

![README_overview](https://github.com/user-attachments/assets/14c4619b-9d05-451f-bde0-b71b1b467fce)


### 1. Overview

In this Github repo we present code for reproducing key analyses results in our Hu, Horlbeck, Zhang et al. paper. Our framework can be divided into two major parts including 

<img src="https://user-images.githubusercontent.com/44768711/193936026-b49715d8-7ec9-4e23-8aa9-330c1f93f2e7.png" width="350" align="left">

(1) Multi-scale footprinting: takes ATAC-Seq (bulk or single cells) data as input, and try to detect DNA-protein interaction across spatial scales. 

The multi-scale footprint pattern at any genomic location delineates local chromatin structure and can be used to infer TF and nucleosome binding. Additionally, we have implemented the infrastructure for generating pseudo-bulks using single cell data, as well as running multi-scale footprinting using the pseudo-bulked data. This provides us with the unique opportunity to track chromatin structure dynamics across pseudo-time.

(2) seq2PRINT: builds a deep learning model that predict footprint patterns from genomic DNA sequence. This model is then used to extract sequence attribution scores and predict TF binding. seq2PRINT can be scaled up to learn footprints and infer TF binding across hundreds of cell states with the help of LoRA.

![README_seq2PRINT](https://github.com/user-attachments/assets/001977bc-0aed-4d95-93e2-43700c2b96f8)

If you want to try our tool on your own data, we recommend using our newest Python package scPrinter, which includes both of the above components. scPrinter can be found at https://github.com/buenrostrolab/scPrinter.

### 2. Key Components

* Correction of Tn5 insertion bias and obtain single-base pair resolution chromatin accessibility

* Calling footprint signal across spatial scales to resolve local chromatin structure

* Tracking nucleosome dynamics across pseudo-time

* Predicting multi-scale footprints using DNA sequence as input

* Decoding sequence syntax of footprints

* Infer TF binding within cis-regulatory elements (CREs)

* Detecting de novo TF motifs

### 4. References

Hu, Horlbeck, Zhang et al., Multi-scale footprints reveal the organization of cis-regulatory elements

### 5. Installation

Currently the framework can be installed by cloning the github repo. Relevant input data can be found at https://zenodo.org/records/13963610. Accession number of public datasets can be found in the Data Availability section and Supplementary Table 8 in our paper.

### 6. Interative data browsers

Interactive visualization using Shinyapps can be found at https://buenrostrolab.shinyapps.io/ACAMShiny/ (human bone marrow) and https://buenrostrolab.shinyapps.io/aging/ (mouse HSC aging).

### 7. Support

If you have any questions, please feel free to open an issue. We appreciate everyone's contribution!
