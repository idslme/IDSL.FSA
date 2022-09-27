# IDSL.FSA <img src='FSA_educational_files/Figures/IDSL.FSA-logo.PNG' width="250px" align="right" />

<!-- badges: start -->
[![Maintainer](https://img.shields.io/badge/maintainer-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
<!-- badges: end -->

[**Fragmentation Spectra Analysis (FSA)**](https://www.fsa.idsl.me/) by the [**Integrated Data Science Laboratory for Metabolomics and Exposomics (IDSL.ME)**](https://www.idsl.me/) is an R package designed to annotate ***.msp*** (mass spectra format) and ***.mgf*** (Mascot generic format) files using a combination of mass spectral entropy similarity, dot product (cosine) similarity, and normalized Euclidean mass error (NEME) criteria. IDSL.FSA also provides a number of modules to convert and manipulate *.msp* and *.mgf* files.

	install.packages("IDSL.FSA") # IDSL.FSA package is set to release on CRAN by the end of October

## Workflow
IDSL.FSA requires [**Fragmentation Spectra DataBase (FSDB)**](https://github.com/idslme/IDSL.FSA/wiki/Fragmentation-Spectra-DataBase-(FSDB)) objects (*.Rdata*) to rapidly annotate mass spectrometry data (***.msp*** and ***.mgf***). IDSL.FSA is also able to generate **FSDB** objects from several reference *.msp* files with inconsistent settings. When the reference *.msp* files or **FSDB** objects are available, to annotate *.msp* and *.mgf* files, download the [FSA parameter spreadsheet](https://raw.githubusercontent.com/idslme/IDSL.FSA/main/FSA_parameters.xlsx) and select the parameters accordingly and then use this spreadsheet as the input for the IDSL.FSA workflow:

	library(IDSL.FSA)
	FSA_workflow("Address of the FSA parameter spreadsheet")


Visit [**wiki**](https://github.com/idslme/IDSL.FSA/wiki) for the detailed documentations and tutorials on specific functions for [**Fragmentation Spectra DataBase (FSDB)**](https://github.com/idslme/IDSL.FSA/wiki/Fragmentation-Spectra-DataBase-(FSDB)), [**MSP files management**](https://github.com/idslme/IDSL.FSA/wiki/MSP-files-management), and [**Spectral Similarity Metrics in IDSL.FSA**](https://github.com/idslme/IDSL.FSA/wiki/Spectral-Similarity-Metrics-in-IDSL.FSA).

Visit https://fsa.idsl.me/ for the detailed documentation and tutorial.