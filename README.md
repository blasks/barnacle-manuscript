# Barnacle Manuscript

This repository organizes the scripts, analyses, figures, and latex files associated with the publication entitled _Simultaneous acclimation to nitrogen and iron scarcity in open ocean cyanobacteria revealed by sparse tensor decomposition of metatranscriptomes_.

## Barnacle

We developed a Python package to implement the sparse tensor decomposition model used in the manuscript to analyze cyanobacterial gene expression data. The code for this package is published in a separate [GitHub repository](https://github.com/blasks/barnacle). Details on library usage can be found on the [Barnacle documentation website](https://barnacle-py.readthedocs.io). The documentation includes:

- [An overview of tensors and tensor decomposition.](https://barnacle-py.readthedocs.io/en/latest/model.html)
- [A gallery of examples that demonstrate tensor analysis with Barnacle.](https://barnacle-py.readthedocs.io/en/latest/examples.html)
- [An API reference of Barnacle modules.](https://barnacle-py.readthedocs.io/en/latest/autodoc/barnacle.html)

For a more technical discussion of the sparse tensor decomposition model implemented in Barnacle, please see the Methods section of the paper.

## Data

All data needed to evaluate the conclusions in the paper can be found in this repository and the Supplementary Materials. 

Supplementary Data files are available in an associated [Zenodo repository](https://doi.org/10.5281/zenodo.12210994). 

The sequences reported have been deposited in the NCBI Sequence Read Archive under the following BioProject accession numbers: 
- [Gradients 1 (2016)](https://simonscmap.com/catalog/cruises/KOK1606)
   - [PRJNA492143](https://www.ncbi.nlm.nih.gov/bioproject/492143) (transect samples)
- [Gradients 2 (2017)](https://simonscmap.com/catalog/cruises/MGL1704)
   - [PRJNA816919](https://www.ncbi.nlm.nih.gov/bioproject/816919) (transect samples)
   - [PRJNA1090086](https://www.ncbi.nlm.nih.gov/bioproject/1090086) (depth samples)
   - [PRJNA1090042](https://www.ncbi.nlm.nih.gov/bioproject/1090042) (20L nutrient incubations)
   - [PRJNA1088528](https://www.ncbi.nlm.nih.gov/bioproject/1088528) (4L sulfur incubations)
- [Gradients 3 (2019)](https://simonscmap.com/catalog/cruises/KM1906)
   - [PRJNA1091352](https://www.ncbi.nlm.nih.gov/bioproject/1091352) (transect samples)
   - [PRJNA1090467](https://www.ncbi.nlm.nih.gov/bioproject/1090467) (depth samples)
   - [PRJNA1090899](https://www.ncbi.nlm.nih.gov/bioproject/1090899) (diel samples)




