### Supplementary materials to:

# Charged residues next to transmembrane regions revisited: "Positive-inside rule" is complemented by the "negative inside depletion/outside enrichment rule"

## Abstract from manuscript

Transmembrane helices (TMHs) are frequently occurring among protein architectures as means for protein attachment to or embedding into biological membranes. Physical constraints such as the membrane's hydrophobicity and electrostatic potential apply uniform requirements to TMHs and their flanking regions; consequently, they are mirrored in their sequence patterns (besides TMHs being a span of generally hydrophobic residues) on top of variations enforced by the specific protein's biological functions. With statistics derived from a large body of protein sequences, we demonstrate that, in addition to the positive charge preference at the cytoplasmic inside ("positive-inside rule"), negatively charged residues preferentially occur or are even enriched at the non-cytoplasmic flank (negative outside rule) or, at least, they are suppressed at the cytoplasmic flank (a "negative-notinside rule"). As negative residues are generally rare within or near TMHs, the statistical significance is sensitive with regard to details of TMH alignment and residue frequency normalization besides to dataset size and, therefore, this trend was obscured in previous work. We observe variations among taxa as well as for organelles along the secretory pathway. The effect is most pronounced for TMHs from single-pass transmembrane (bitopic) proteins compared to those with multiple TMH (polytopic proteins) and especially for the class of simple TMHs that evolved for the sole role as membrane anchors. The charged residue flank bias is only one of the TMH sequence features with a role in the anchorage mechanisms, others apparently being the leucine intra-helix propensity skew towards the cytoplasmic side, tryptophan flanking as well as the cysteine and tyrosine inside preference. These observations will stimulate new prediction methods for TMHs and protein topology from sequence as well as new engineering designs for artificial membrane proteins.

--------------------------------------------------------------------------------

## Supplementary to methods and results section.

### About

Provided here are the files and scripts used to perform the studies outlined in the abstract. These scripts can be used to mine database files from uniprot and TopDB into transmembrane and perform sequence distribution data analysis.

This release isn't designed for redistribution or reapplication but is provided 'as is'. If you are having trouble running the scripts, or if there is something unclear, do feel free to [open an issue on GitHub](https://github.com/jbkr/nni_supplementary_material/issues/new) or contact the authors directly.

### Download

The downloaded zip includes the original files downloaded from the respective databases and the Uniprot non-redundant datasets. [Click here to visit the TopDB site](http://topdb.enzim.hu/?m=download&mid=2). The zip also contains the python scripts used to generate the datasets, tables, figures, as well as the parsed .csv files of each dataset. Within the .csv files is the ID for the respective databasse, the full protein sequence, the TMH sequences, the flank sequences (each file for each dataset has a different cut-off flank length: 5, 10 or 20), the number of TMHs in the given protein, and the orientation of the TMH.

 - [Click here to download](download.zip).

The file structure includes the scripts that make the datasets at the top level. Once a dataset is made, this can be copied to the `tables_and_figures` folder for analysis.

<!-- This needs updating -->

```
.
├── ExpAll90%ID_list.txt
├── readme.md
├── representative_list_topdb_to_uniprot.py
├── tables_and_figures
│   ├── FigureX.py
│   ├── TableX.py
│   ├── List_of_Complex_IDs.txt
│   ├── List_of_Simple_IDs.txt
│   ├── List_of_Twilight_IDs.txt
│   ├── old_datasets
│   │   ├── TopDB_10_flanklength.csv
│   │   ├── TopDB_20_flanklength.csv
│   │   ├── TopDB_5_flanklength.csv
│   │   ├── UniArch_10_flanklength.csv
│   │   ├── UniArch_20_flanklength.csv
│   │   ├── UniArch_5_flanklength.csv
│   │   ├── UniBacilli_10_flanklength.csv
│   │   ├── UniBacilli_20_flanklength.csv
│   │   ├── UniBacilli_5_flanklength.csv
│   │   ├── UniCress_10_flanklength.csv
│   │   ├── UniCress_20_flanklength.csv
│   │   ├── UniCress_5_flanklength.csv
│   │   ├── UniEcoli_10_flanklength.csv
│   │   ├── UniEcoli_20_flanklength.csv
│   │   ├── UniEcoli_5_flanklength.csv
│   │   ├── UniER_10_flanklength.csv
│   │   ├── UniER_20_flanklength.csv
│   │   ├── UniER_5_flanklength.csv
│   │   ├── UniFungi_10_flanklength.csv
│   │   ├── UniFungi_20_flanklength.csv
│   │   ├── UniFungi_5_flanklength.csv
│   │   ├── UniGolgi_10_flanklength.csv
│   │   ├── UniGolgi_20_flanklength.csv
│   │   ├── UniGolgi_5_flanklength.csv
│   │   ├── UniHuman_10_flanklength.csv
│   │   ├── UniHuman_20_flanklength.csv
│   │   ├── UniHuman_5_flanklength.csv
│   │   ├── UniPM_10_flanklength.csv
│   │   ├── UniPM_20_flanklength.csv
│   │   └── UniPM_5_flanklength.csv
│   ├── top_all_10_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── top_all_10_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── top_all_10_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── top_all_10_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── top_all_20_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── top_all_20_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── top_all_20_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── top_all_20_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── top_all_5_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── top_all_5_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── top_all_5_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── top_all_5_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── TopDB_10_flanklength_flankclashFalse.csv
│   ├── TopDB_10_flanklength_flankclashTrue.csv
│   ├── TopDB_20_flanklength_flankclashFalse.csv
│   ├── TopDB_20_flanklength_flankclashTrue.csv
│   ├── TopDB_5_flanklength_flankclashFalse.csv
│   ├── TopDB_5_flanklength_flankclashTrue.csv
│   ├── UniArch_10_flanklength_flankclashFalse.csv
│   ├── UniArch_10_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniArch_10_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniArch_10_flanklength_flankclashTrue.csv
│   ├── UniArch_10_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniArch_10_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniArch_20_flanklength_flankclashFalse.csv
│   ├── UniArch_20_flanklength_flankclashFalse_logged_lengthexclusionIDs.txt
│   ├── UniArch_20_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniArch_20_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniArch_20_flanklength_flankclashTrue.csv
│   ├── UniArch_20_flanklength_flankclashTrue_logged_lengthexclusionIDs.txt
│   ├── UniArch_20_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniArch_20_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniArch_5_flanklength_flankclashFalse.csv
│   ├── UniArch_5_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniArch_5_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniArch_5_flanklength_flankclashTrue.csv
│   ├── UniArch_5_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniArch_5_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniBacilli_10_flanklength_flankclashFalse.csv
│   ├── UniBacilli_10_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniBacilli_10_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniBacilli_10_flanklength_flankclashTrue.csv
│   ├── UniBacilli_10_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniBacilli_10_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniBacilli_20_flanklength_flankclashFalse.csv
│   ├── UniBacilli_20_flanklength_flankclashFalse_logged_lengthexclusionIDs.txt
│   ├── UniBacilli_20_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniBacilli_20_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniBacilli_20_flanklength_flankclashTrue.csv
│   ├── UniBacilli_20_flanklength_flankclashTrue_logged_lengthexclusionIDs.txt
│   ├── UniBacilli_20_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniBacilli_20_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniBacilli_5_flanklength_flankclashFalse.csv
│   ├── UniBacilli_5_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniBacilli_5_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniBacilli_5_flanklength_flankclashTrue.csv
│   ├── UniBacilli_5_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniBacilli_5_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniCress_10_flanklength_flankclashFalse.csv
│   ├── UniCress_10_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniCress_10_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniCress_10_flanklength_flankclashTrue.csv
│   ├── UniCress_10_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniCress_10_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniCress_20_flanklength_flankclashFalse.csv
│   ├── UniCress_20_flanklength_flankclashFalse_logged_lengthexclusionIDs.txt
│   ├── UniCress_20_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniCress_20_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniCress_20_flanklength_flankclashTrue.csv
│   ├── UniCress_20_flanklength_flankclashTrue_logged_lengthexclusionIDs.txt
│   ├── UniCress_20_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniCress_20_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniCress_5_flanklength_flankclashFalse.csv
│   ├── UniCress_5_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniCress_5_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniCress_5_flanklength_flankclashTrue.csv
│   ├── UniCress_5_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniCress_5_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniEcoli_10_flanklength_flankclashFalse.csv
│   ├── UniEcoli_10_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniEcoli_10_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniEcoli_10_flanklength_flankclashTrue.csv
│   ├── UniEcoli_10_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniEcoli_10_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniEcoli_20_flanklength_flankclashFalse.csv
│   ├── UniEcoli_20_flanklength_flankclashFalse_logged_lengthexclusionIDs.txt
│   ├── UniEcoli_20_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniEcoli_20_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniEcoli_20_flanklength_flankclashTrue.csv
│   ├── UniEcoli_20_flanklength_flankclashTrue_logged_lengthexclusionIDs.txt
│   ├── UniEcoli_20_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniEcoli_20_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniEcoli_5_flanklength_flankclashFalse.csv
│   ├── UniEcoli_5_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniEcoli_5_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniEcoli_5_flanklength_flankclashTrue.csv
│   ├── UniEcoli_5_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniEcoli_5_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniER_10_flanklength_flankclashFalse.csv
│   ├── UniER_10_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniER_10_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniER_10_flanklength_flankclashTrue.csv
│   ├── UniER_10_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniER_10_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniER_20_flanklength_flankclashFalse.csv
│   ├── UniER_20_flanklength_flankclashFalse_logged_lengthexclusionIDs.txt
│   ├── UniER_20_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniER_20_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniER_20_flanklength_flankclashTrue.csv
│   ├── UniER_20_flanklength_flankclashTrue_logged_lengthexclusionIDs.txt
│   ├── UniER_20_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniER_20_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniER_5_flanklength_flankclashFalse.csv
│   ├── UniER_5_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniER_5_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniER_5_flanklength_flankclashTrue.csv
│   ├── UniER_5_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniER_5_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniFungi_10_flanklength_flankclashFalse.csv
│   ├── UniFungi_10_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniFungi_10_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniFungi_10_flanklength_flankclashTrue.csv
│   ├── UniFungi_10_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniFungi_10_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniFungi_20_flanklength_flankclashFalse.csv
│   ├── UniFungi_20_flanklength_flankclashFalse_logged_lengthexclusionIDs.txt
│   ├── UniFungi_20_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniFungi_20_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniFungi_20_flanklength_flankclashTrue.csv
│   ├── UniFungi_20_flanklength_flankclashTrue_logged_lengthexclusionIDs.txt
│   ├── UniFungi_20_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniFungi_20_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniFungi_5_flanklength_flankclashFalse.csv
│   ├── UniFungi_5_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniFungi_5_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniFungi_5_flanklength_flankclashTrue.csv
│   ├── UniFungi_5_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniFungi_5_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniGolgi_10_flanklength_flankclashFalse.csv
│   ├── UniGolgi_10_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniGolgi_10_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniGolgi_10_flanklength_flankclashTrue.csv
│   ├── UniGolgi_10_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniGolgi_10_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniGolgi_20_flanklength_flankclashFalse.csv
│   ├── UniGolgi_20_flanklength_flankclashFalse_logged_lengthexclusionIDs.txt
│   ├── UniGolgi_20_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniGolgi_20_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniGolgi_20_flanklength_flankclashTrue.csv
│   ├── UniGolgi_20_flanklength_flankclashTrue_logged_lengthexclusionIDs.txt
│   ├── UniGolgi_20_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniGolgi_20_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniGolgi_5_flanklength_flankclashFalse.csv
│   ├── UniGolgi_5_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniGolgi_5_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniGolgi_5_flanklength_flankclashTrue.csv
│   ├── UniGolgi_5_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniGolgi_5_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniHuman_10_flanklength_flankclashFalse.csv
│   ├── UniHuman_10_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniHuman_10_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniHuman_10_flanklength_flankclashTrue.csv
│   ├── UniHuman_10_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniHuman_10_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniHuman_20_flanklength_flankclashFalse.csv
│   ├── UniHuman_20_flanklength_flankclashFalse_logged_lengthexclusionIDs.txt
│   ├── UniHuman_20_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniHuman_20_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniHuman_20_flanklength_flankclashTrue.csv
│   ├── UniHuman_20_flanklength_flankclashTrue_logged_lengthexclusionIDs.txt
│   ├── UniHuman_20_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniHuman_20_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniHuman_5_flanklength_flankclashFalse.csv
│   ├── UniHuman_5_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniHuman_5_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniHuman_5_flanklength_flankclashTrue.csv
│   ├── UniHuman_5_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniHuman_5_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniPM_10_flanklength_flankclashFalse.csv
│   ├── UniPM_10_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniPM_10_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniPM_10_flanklength_flankclashTrue.csv
│   ├── UniPM_10_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniPM_10_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniPM_20_flanklength_flankclashFalse.csv
│   ├── UniPM_20_flanklength_flankclashFalse_logged_lengthexclusionIDs.txt
│   ├── UniPM_20_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniPM_20_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniPM_20_flanklength_flankclashTrue.csv
│   ├── UniPM_20_flanklength_flankclashTrue_logged_lengthexclusionIDs.txt
│   ├── UniPM_20_flanklength_flankclashTrue_only_full_flanks.csv
│   ├── UniPM_20_flanklength_flankclashTrue_only_half_flanks.csv
│   ├── UniPM_5_flanklength_flankclashFalse.csv
│   ├── UniPM_5_flanklength_flankclashFalse_only_full_flanks.csv
│   ├── UniPM_5_flanklength_flankclashFalse_only_half_flanks.csv
│   ├── UniPM_5_flanklength_flankclashTrue.csv
│   ├── UniPM_5_flanklength_flankclashTrue_only_full_flanks.csv
│   └── UniPM_5_flanklength_flankclashTrue_only_half_flanks.csv
├── top_all.txt #Download this externally
├── topdb_to_table(fasta).py
├── UniArch.txt
├── UniBacilli.txt
├── UniCress.txt
├── UniEcoli.txt
├── UniER.txt
├── UniFungi.txt
├── UniGolgi.txt
├── UniHuman.txt
├── UniPM.txt
└── Uniprot_to_table.py

```

## Features

- Mines Uniprot files into tables that have easier to handle information about their transmembrane domain and neighboring residue sequences in csv format.

  - Detects annotation for transmembrane regions.
  - Detects annotation for orientation.
  - Detects annotation for singlepass and multipass transmembrane proteins.
  - Detects for overlapping residues in the flanks and evenly reassigns the flanking residues to the respective TMH.

- Once a dataset is prepared, it can be copied to `tables_and_figures`. There, scripts exist to analyse statistical variance distributions of the residues. Please refer to the paper for a more detailed explanation for what each figure and table represents about the data.

## Installation and usage

These scripts require Biopython, numpy, and [python 2.7](https://www.python.org/downloads/). To install Biopython and numpy quickly on unix or unix like (OSX) system, run the following commands.

1. Open a terminal.

2. Run `sudo easy_install pip; sudo pip install numpy; sudo pip install Biopython` and enter your password as appropriate. If you come across any errors it is probably because python is not installed in the default locations, or a package has already been installed before you did these commands.

3. Run the database builder scripts using `python topdb_to_table(fasta).py` and `python Uniprot_to_table.py`.

4. Move the files you want to produce tables and figures into the `tables_and_figures` folder.

5. Run the respective scripts. Some may require moving vectors and numbers outputted into a graph producing software like Excel.
