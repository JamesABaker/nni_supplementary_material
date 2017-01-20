

# Supplementary Code
## Charged residues next to transmembrane regions revisited: “Positive-inside rule” is complemented by the “negative inside depletion/outside enrichment rule”

## Abstract from manuscript
Transmembrane helices (TMHs) are frequently occurring among protein architectures as means for protein attachment to or embedding into biological membranes. Physical constraints such as the membrane’s hydrophobicity and electrostatic potential apply uniform requirements to TMHs and their flanking regions; consequently, they are mirrored in their sequence patterns (besides TMHs being a span of generally hydrophobic residues) on top of variations enforced by the specific protein’s biological functions. With statistics derived from a large body of protein sequences, we demonstrate that, in addition to the positive charge preference at the cytoplasmic inside (“positive-inside rule”), negatively charged residues preferentially occur or are even enriched at the non-cytoplasmic flank (negative outside rule) or, at least, they are suppressed at the cytoplasmic flank (a “negative-notinside rule”). As negative residues are generally rare within or near TMHs, the statistical significance is sensitive with regard to details of TMH alignment and residue frequency normalization besides to dataset size and, therefore, this trend was obscured in previous work. We observe variations among taxa as well as for organelles along the secretory pathway. The effect is most pronounced for TMHs from single-pass transmembrane (bitopic) proteins compared to those with multiple TMH (polytopic proteins) and especially for the class of simple TMHs that evolved for the sole role as membrane anchors. The charged residue flank bias is only one of the TMH sequence features with a role in the anchorage mechanisms, others apparently being the leucine intra-helix propensity skew towards the cytoplasmic side, tryptophan flanking as well as the cysteine and tyrosine inside preference. These observations will stimulate new prediction methods for TMHs and protein topology from sequence as well as new engineering designs for artificial membrane proteins.

## About these files
Provided here are the files used to perform the studies outlined in the abstract. These scripts can be used to mine database files from uniprot into transmembrane sequence distribution data that tests the validity of topological constraints of membrane helices. There also exists a script to build a dataset from TOPDB fasta files which contain topological annotation.

This release isn't designed for redistribution or reapplication, but rather for data replication. Because the scripts have been developed with this in mind, the documentation may be hard to follow. If you are having trouble running the scripts, or if there is something unclear, do feel free to open an issue.

The file structure includes the scripts that make the datasets at the top level. Once a dataset is made, this can be copied to the `tables_and_figures` folder for analysis.

<!--This needs updating-->
```
├── tables_and_figures
│   ├── Figures.py
├── topdb_to_table(fasta).py
├── *.txt #The database files
└── Uniprot_to_table.py
  ```


## Features
- Mines Uniprot files into tables with easier to handle information about their transmembrane domain and neighboring residue sequences.
	- Detects annotation for transmembrane regions.
	- Detects annotation for orientation.
	- Detects annotation for singlepass and multipass transmembrane proteins.
    - Detects for overlapping residues in the flanks and evenly reassigns the flanking residues to the respective TMH.
- Once a dataset is prepared, it can be copied to `tables_and_figures`. There, scripts exist to analyse statistical variance distributions of the residues. Please refer to the paper for a more detailed explanation for what each figure and table represents about the data.

## Installation
These scripts require Biopython, numpy, and [python 2.7](https://www.python.org/downloads/). To install Biopython and numpy quickly on unix or unix like (OSX) system, run the following commands.
1. Open a terminal.
2. Run `sudo easy_install pip; sudo pip install numpy; sudo pip install Biopython`  and enter your password as appropriate.

If you come across any errors it is probably because python is not installed in the default locations, or a package has already been installed before you did these commands.

## Usage

Whilst the scripts are running, do not move or change files in the directory. This will cause problems in the current run, and the next run.

## v0.2 Patch notes

- Fixed a bug where overlapping flanking residues were being added twice to the dataset in many but not all cases. The flank length restriction now kicks in correctly when flanks overlap with other flanks. For comparison, files were also made that have no flank overlap removal.
- Added comments throughout the scripts.
- Numbering system in dataset files now reflect the respective dataset nomenclature rather than the index value used for extraction.
- Scripts now conform much more to PEP8 style guidelines.
- The `readme` is now more informative.
