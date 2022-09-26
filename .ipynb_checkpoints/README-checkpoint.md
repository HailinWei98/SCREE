# SCREE
<strong>SCREE</strong>(<strong>S</strong>ingle-cell <strong>C</strong>RISPR sc<strong>RE</strong>en data analyses and p<strong>E</strong>rturbation modeling) is a workflow for single-cell CRISPR screens data processing and analysis. SCREE has two separate section: pre-process and analysis. The pre-process section is to perform reads alignment and quantification based on <a href="https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome">cellranger</a> and <a href="https://support.10xgenomics.com/single-cell-atac/software/overview/welcome">cellranger-atac</a>. The analysis section is an R package including functions of sgRNA assignment, sgRNA information visualization, single-cell quality control, perturbation enrichment ratio calculation, perturbation efficiency evaluation, regulatory score estimation and downstream functional analysis. All the outputs of SCREE can be integrated into an html file. <br/>

<img src="image/workflow.png" width="700" height="500"></img>

## Document
We are hosting SCREEN documentation, instruction and tutorials at <a href="https://hailinwei98.github.io/SCREEN.html">SCREEN Website</a>.

## Dependency
	R >= 4.0.3
	Seurat >= 4.0.3
	cellranger >= 6.0
	cellranger-atac
	setuptools
	devtools

## Installation of SCREE
Before installation and operation of SCREE, you should install cellranger (for scRNA-seq based data) or cellranger-atac (scATAC-seq based data) first to make sure that SCREE can perform alignment and quantification accurately.

## Installation of R package
	$ R
	> library(devtools)
	> install_github("HailinWei98/SCREE")
