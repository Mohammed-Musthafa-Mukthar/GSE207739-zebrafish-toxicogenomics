# GSE207739 zebrafish toxicogenomics

Reanalysis of the zebrafish toxicogenomics dataset GSE207739 from Morash et al. (2023).

The aim of this project is to reproduce selected parts of the published analysis and then work back from the processed data to the raw RNA-seq data.

## Data

The first part of the analysis uses the processed expression matrices available from GEO for Raloxifene and Resorcinol.

Each dataset contains 30 samples across six dose groups, with five replicates per group.

The raw sequencing data will be analysed separately later.

## Current analysis

I started by checking the processed expression files, extracting the sample information, and running PCA.

Raloxifene:
- PC1: 47.3%
- PC2: 13.6%

Resorcinol:
- PC1: 37.0%
- PC2: 11.1%

The PCA showed a strong pattern associated with the replicate labels A-E in both datasets. Samples with the same replicate label tended to group together across dose levels.

At this point I am treating this as replicate-associated variation rather than calling it a batch effect, since the available metadata do not clearly explain what the A-E labels represent experimentally.

## Scripts

`R/bmd_data.R`

Reads the processed GEO files, extracts sample information, checks the data, runs PCA, and saves the PCA coordinates and figures.

`R/pca_replicates.R`

Uses the same PCA results and plots dose group by colour and replicate label by shape.

## Next steps

- reproduce the BMDExpress analysis
- compare benchmark dose estimates with the published values
- download the raw RNA-seq data from SRA
- run read QC and STAR alignment against GRCz11
- reproduce the differential expression analysis
- compare pathway results with the paper