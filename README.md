# RCPA

The R package for Consensus Pathway Analysis (RCPA) implements a complete analysis pipeline including: i) download and process data from NCBI Gene Expression Omnibus, ii) perform differential analysis using techniques developed for both microarray and sequencing data, iii) perform systems-level analysis using different methods for enrichment analysis and topology-based (TB) analysis, iv) perform meta-analysis and consensus analysis, and v) visualize analysis results and explore significantly impacted pathways across multiple analyses. The package supports the analysis of more than 1,000 species, two pathway databases, three differential analysis techniques, eight pathway analysis tools, six meta-analysis methods, and two consensus analysis techniques.

<img align="center" src="./img/pipeline.png?raw=true">

# How to install
- The package can be installed from CRAN or this repository.
- Using CRAN: `install.packages('RCPA')`
- Using devtools:
    - Install the package using: `devtools::install_github('tinnlab/RCPA')`  

# Tutorial
- A detailed tutorial on how to use scDHA package is available in Jupyter Notebook format, which can be accessed at this [Example](https://github.com/tinnlab/RCPA/blob/main/examples/RCPA-Protocol.ipynb).