# STPP-UP data analysis

This repository contains workflows for the analysis of STPP-UP experiments (https://pubmed.ncbi.nlm.nih.gov/37742922/). The prerequisites are that label-free quantification was used during data acquisition and that the raw mass spectrometry files were analyzed using DIA-NN (https://github.com/vdemichev/DiaNN)

## Getting started

Getting started is easy. Download the R script and analysis parameters file to your working directory. To do this:
- Click on the big green button in the top right that says ```<> Code``` and then click on ```Download ZIP```.
- Place the files in your working directory
- Fill in the analysis parameter file accordingly.
- Open the R script with RStudio, select all and [run it](https://www.youtube.com/watch?v=w6QGe-pXgdI).

**Note:** - Installing the DEP package for the first time may take a while. Should be smooth sailing after that

## Parameters explained

Explanations are also present in the Excel file, when hovering over the parameter boxes.

```filtering_type``` - Whether a protein should have a value in all samples ("complete") or all samples in at least one condition/group ("condition").

```common_name_value_columns``` - This should be a common element present in the names of the .raw files. It will be used to select the value columns.

```comparison``` - The type of comparison that will be tested. This can be all possible pairwise comparisons ("all"), a manual selection ("manual") or limited to the comparisons versus the control ("control").

```control``` - In case you choose "control" as your comparison, add the name of the control sample, e.g. DMSO, IgG, untreated, KO.

```contrasts``` - In case you choose "manual" as your comparison, add the comparisons you would like to make, separated by a semicolon. E.g. C1_vs_DMSO;DMSO_vs_IgG.

```unique_peptides``` - The minimum amount of unique peptides a protein should be quantified with.

```volcano``` - Specifiy which proteins to highlight in the volcano plot. Can be either a list of supplied protein names ("protein list") or the top N most differential proteins below a specified p-value ("TopN").

```protein_list``` - list of protein names to be highlighted, separated by a semicolon (;).

```p_value``` - p-value cutoff for significance. **REQUIRED**

```log2_FC``` - log2 FC value cutoff for significance. **REQUIRED**

```TopN``` - TopN number of proteins to highlight.

## Required files

Make sure your .raw files and output are correctly named, as outlined below. Otherwise the script may not run succesfully.

```.raw mass spec files``` - Please ensure that the abundance column names (and thus the .raw files) contain the name of the mass spectrometer used (e.g. Astral, Exploris) and the following information, in order and separated by an underscore:
- heating conditions (e.g. Carrier)
- drug/compound treatment (e.g. C1/DMSO)
- Replicate (e.g. 1)
             
Example: ```20250205_Astral_AUR_DZ114_Carrier_C1_1.raw```

```DiaNN``` - Output files should end on ```*pg_matrix.tsv``` and ```*pr_matrix.tsv```.

### Credit

This script uses a lot of the DEP package, so credit to them for developing the package [(Zhang X, Smits A, van Tilburg G, Ovaa H, Huber W, Vermeulen M (2018). “Proteome-wide identification of ubiquitin interactions using UbIA-MS.” Nature Protocols, 13, 530–550)](https://www.nature.com/articles/nprot.2017.147).
