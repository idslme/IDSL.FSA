# IDSL.FSA <img src='FSA_educational_files/Figures/IDSL.FSA-logo.PNG' width="250px" align="right" />

<!-- badges: start -->
[![Developed-by](https://img.shields.io/badge/Developed_by-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
[![CRAN status](https://www.r-pkg.org/badges/version/IDSL.FSA)](https://cran.r-project.org/package=IDSL.FSA)
![](http://cranlogs.r-pkg.org/badges/IDSL.FSA?color=orange)
![](http://cranlogs.r-pkg.org/badges/grand-total/IDSL.FSA?color=brightgreen)
[![Dependencies](https://tinyverse.netlify.com/badge/IDSL.FSA)](https://cran.r-project.org/package=IDSL.FSA)


[![DOI](https://zenodo.org/badge/140601694.svg)](https://zenodo.org/record/7530397#.Y8Byuv7MK70)
<!-- badges: end -->

**Fragmentation Spectra Analysis (FSA)** by the [**Integrated Data Science Laboratory for Metabolomics and Exposomics (IDSL.ME)**](https://www.idsl.me) is an R package designed to annotate ***.msp*** (mass spectra format) and ***.mgf*** (Mascot generic format) files using a combination of mass spectral entropy similarity, dot product (cosine) similarity, and normalized Euclidean mass error (NEME) criteria followed by intelligent pre-filtering steps for rapid searches. IDSL.FSA also provides a number of modules to convert and manipulate ***.msp*** and ***.mgf*** files.

## Table of Contents

- [Features of IDSL.FSA](https://github.com/idslme/IDSL.FSA#features-of-idslfsa)
- [Installation](https://github.com/idslme/IDSL.FSA#installation)
- [Workflow](https://github.com/idslme/IDSL.FSA#workflow)
- [Quick Batch Example](https://github.com/idslme/IDSL.FSA#quick-batch-example)
- [Wiki](https://github.com/idslme/IDSL.FSA#wiki)
- [Citation](https://github.com/idslme/IDSL.FSA#citation)

## Features of IDSL.FSA

1) Parameter selection through a user-friendly and well-described [parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.FSA/main/FSA_parameters.xlsx)
2) Analyzing population size untargeted studies (n > 500)
3) Consistency with ***.msp*** and ***.mgf*** file formats
4) Measuring matching variables including mass spectral entropy similarity, dot product (cosine) similarity, and normalized Euclidean mass error (NEME)
5) Spectra annotation regardless of presence of precursor values
6) Spectra annotation using optional criteria including retention time and adduct type
7) Spectra annotation using accurate or nominal mass values
8) Generating batch spectra match figures
9) Parallel processing in Windows and Linux environments

## Installation

	install.packages("IDSL.FSA")
	
**Note:** In some instances, the [readxl](https://cran.r-project.org/package=readxl) R package is also required to read the [FSA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.FSA/main/FSA_parameters.xlsx).

	install.packages("readxl")

## Workflow

IDSL.FSA requires [**Fragmentation Spectra DataBase (FSDB)**](https://github.com/idslme/IDSL.FSA/wiki/Fragmentation-Spectra-DataBase-(FSDB)) objects (*.Rdata*) to rapidly annotate mass spectrometry data (***.msp*** and ***.mgf***). IDSL.FSA is also able to generate **FSDB** objects from several reference *.msp* files even with inconsistent settings. When the reference *.msp* files or **FSDB** objects are available, to annotate *.msp* and *.mgf* files, download the [FSA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.FSA/main/FSA_parameters.xlsx) and select the parameters accordingly and then use this spreadsheet as the input for the IDSL.FSA workflow:

	library(IDSL.FSA)
	FSA_workflow("Address of the FSA parameter spreadsheet")

## Quick Batch Example

Follow these steps for a quick msp annotation of an ***.msp*** file

1. Download [Kynurenine_Kynurenic_acid.msp](https://github.com/idslme/IDSL.FSA/blob/main/IDSL.FSA/inst/extdata/Kynurenine_Kynurenic_acid.msp). You should also ensure this file ends up with an ***.msp*** extension after downloading. You may also need to rename this file to have an ***.msp*** extension.

2. Download the positive mode [FSDB](https://zenodo.org/record/7530397#.Y8yAdkHMK71)

3. IDSL.FSA requires parameters for three tabs for a full scale run. For this study, use default parameter values presented in the [FSA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.FSA/main/FSA_parameters.xlsx). Then, provide information for 
	
	3.1. **FSA0002** for the *Address of the MS/MS library (FSDB) generated using `FSDB` tab*
	
	3.2. **FSA0003** for *Location of the INPUT sample of .msp and/or .mgf files*
	
	3.3. **FSA0004** for *Location to save OUTPUT results*
		
	3.4. You may also adjust the number of processing threads using **SPEC0001** in the `SpectraSimilarity` tab according to your computational power

4. Run this command in the R/Rstudio console or terminal.

```
library(IDSL.FSA)
FSA_workflow("Address of the FSA parameter spreadsheet")
```

5. You see the results in the address you provided for **FSA0004** including:

	5.1. individual *SpectraAnnotationTable* for each ***.msp*** file in the *annotated_spectra_tables* directory in the *.Rdata* and *.csv* formats
	
	5.2. if you had selected an integer greater than 0 for **SPEC0019**, match spectra figures are also available in the *plot* folder for each MSP block.

## [**Wiki**](https://github.com/idslme/IDSL.FSA/wiki)

1. [**Fragmentation Spectra DataBase (FSDB)**](https://github.com/idslme/IDSL.FSA/wiki/Fragmentation-Spectra-DataBase-(FSDB))
2. [**MSP files management**](https://github.com/idslme/IDSL.FSA/wiki/MSP-Files-Management)
3. [**MSP Standardizations**](https://github.com/idslme/IDSL.FSA/wiki/MSP-Standardizations)
4. [**Spectral Similarity Metrics**](https://github.com/idslme/IDSL.FSA/wiki/Spectral-Similarity-Metrics)

## Citation

[1] Fakouri Baygi, S., Kumar, Y. Barupal, D.K. [IDSL.CSA: Composite Spectra Analysis for Chemical Annotation of Untargeted Metabolomics Datasets](https://doi.org/10.1021/acs.analchem.3c00376). *Analytical Chemistry*, **2023**, *95(25)*, 9480–9487.
