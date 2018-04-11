# FlowSorted.Blood.EPIC
[![Travis-CI Build Status](https://travis-ci.org/immunomethylomics/FlowSorted.Blood.EPIC.svg?branch=master)](https://travis-ci.org/immunomethylomics/FlowSorted.Blood.EPIC)
Illumina Human Methylation data from EPIC on magnetic sorted adult blood cell populations

The FlowSorted.Blood.EPIC package contains Illumina HumanMethylationEPIC DNA methylation microarray data from the immunomethylomics group (manuscript in preparation), consisting of 37 blood cell references and 12 samples, formatted as an RGChannelSet object for integration and normalization using most of the existing Bioconductor packages.

This package contains data similar to the FlowSorted.Blood.450k package consisting of data from peripheral blood samples generated from adult men. However, when using the newer EPIC microarray minfi estimates of cell type composition using the FlowSorted.Blood.450k package are less precise compared to actual cell counts. Hence, this package consists of appropriate data for deconvolution of adult blood samples used in for example EWAS relying in the newer EPIC technology.

Researchers may find this package useful as these samples represent different cellular populations ( T lymphocytes (CD4+ and CD8+), B cells (CD19+), monocytes (CD14+), NK cells (CD56+) and Neutrophils of cell sorted blood generated with high purity estimates. As a test of accuracy 12 experimental mixtures were reconstructed using fixed amounts of DNA from purified cells. These data can be integrated with the minfi Bioconductor package to estimate cellular composition in users' whole blood Illumina EPIC samples using a modified version of the algorithm constrained projection/quadratic programming described in Houseman et al. 2012. For a slightly more accurate estimations we also offered an IDOL optimized CpG selection for cell deconvolution (see the references) and a modified version of estimateCellCounts named estimateCellCounts2 which allows using customized sets of probes from IDOL.

References: D Koestler et al. (2016). Improving cell mixture deconvolution by identifying optimal DNA methylation libraries (IDOL). BMC bioinformatics. 17, 120.

EA Houseman et al. (2012) DNA methylation arrays as surrogate measures of cell mixture distribution. BMC Bioinformatics 13, 86.


