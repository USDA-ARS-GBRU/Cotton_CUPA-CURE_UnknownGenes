# Cotton_CUPA-CURE_UnknownGenes

Cotton suffers from a mostly functionally unannotated proteome. This problem is unsolvable by our current automated annotation tools. We systematically identified all unannotated proteins and their corresponding genes and filtered out those that were not appropriate for downstream analysis. We devised one mutually beneficial solution, which involves manual annotation relying on undergraduate students presented in a course-based format (an effort entitled the Cotton Unknown Protein Analysis Course Based Undergraduate Research Experience or CUPA-CURE). This approach ultimately provides valuable annotations, in the form of a microPublication, that would otherwise remain synthesized. 

This repository contains code for analyzing the unknown proteins and workbooks for delivering the CUPA-CURE.


## Citation
[![DOI](https://zenodo.org/badge/217529920.svg)](https://zenodo.org/badge/latestdoi/217529920)
Add a DOI link when it's available.

* How to cite - update prior to publication and release. 

*This is a project supported by the U.S. Department of Agriculture - Agricultural Research Service (USDA-ARS) - Genomics and Bioinformatics Research Unit (GBRU).* 

## Description

These scripts are here for extracting unknown genes from the annotation .txt file. A seperate .R file is provided for generating the graphs displayed in the publication.

## Getting Started

### Dependencies

* Windows 10 Powershell

### Executing program
* Enter your input and output files into the script provided for each step to generate the list of unknown genes.
* Each step should be run individually in order.

#### Basic use-case

* Extracting unnannoted or unknown genes from an annotation file.
* One of the analyzed genome annotation files is included as a sample data set. 
* The first command with the output and input variables replaced is provided below.
```
Get-Content "Ghirsutum_527_v2.1.annotation_info.txt" | Where-Object {$_ -match 'unknown|DUF|putative|uncharacterised' } | Out-File -FilePath "Output.txt"
```

#### Detailed Usage

* It is important to note that a genomes annotation file may differ from the cotton genomes analyzed. Where that is the case the column or key words being searched must be adjusted.
* This can be done by editing the [##] in order to adjust the columns being searched.
* An additional command is provided at the bottom of the code document detailing a use case where only the first 12 columns are being searched.

```
Get-Content "INPUT" |
Where-Object {($_.Split("`t")[0..11]) -match 'unknown|DUF|putative|uncharacterised'} |
Out-File -FilePath "OUTPUT"
```

## Help

* Keep in mind that powershell works on a zero index. In other words the first column is called by [0] the second column by [1] etc.


### Corresponding Contact

Contact info for current maintaining author

Jonathan Zirkel jwzirkel@ncsu.edu
Amanda Hulse-Kemp amanda.hulse-kemp@usda.gov

## License

This software is a work of the United States Department of Agriculture, Agricultural Research Service and is released under a Creative Commons CC0 public domain attribution - see the LICENSE.txt file for details

## Acknowledgments

Inspiration, code snippets, etc.
* [readme_template](https://gist.github.com/DomPizzie/7a5ff55ffa9081f2de27c315f5018afc)

### Funding Support
This is a project supported by the U.S. Department of Agriculture - Agricultural Research Service (USDA-ARS) - Genomics and Bioinformatics Research Unit (GBRU) through CRIS Project No. 6066-21310-006-000-D.
Additional project support was through <and any additional agreements or grants>.
