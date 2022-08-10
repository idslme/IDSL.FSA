# IDSL.FSA <img src='FSA_educational_files/Figures/IDSL.FSA-logo.png' width="250px" align="right" />

<!-- badges: start -->
[![Maintainer](https://img.shields.io/badge/maintainer-Sadjad_Fakouri_Baygi-blue)](https://github.com/sajfb)
<!-- badges: end -->

The Fragmentation Spectra Analysis (IDSL.FSA) package was designed to annotate standard .msp (mass spectra format) and .mgf (Mascot generic format) files using mass spectral entropy similarity, dot product (cosine) similarity, and normalized Euclidean mass error (NEME). IDSL.FSA also provides a number of modules to convert and manipulate MSP and MGF files.

# FSDB (Fragmentation Spectra DataBase)
FSDB objects are R readable objects from one or several msp files. FSDB objects facilitate msp annotation and msp data parsing in the R environment. In the *FSDB (Fragmentation Spectra DataBase)* folder, we uploaded FSDB files for GNPS and MoNA public libraries.

## `msp2FSdb`
The `msp2FSdb` module of the IDSL.FSA is used to generate FSDB objects. This module was designed to be consistent with various msp structures particularly from NIST, GNPS, MoNA, IDSL.CSA libraries. The `msp2FSdb` module generally can work for any msp file as long as *Num Peaks* lines are available in the msp file.

## FSDB objects:
FSDB objects are R lists consisting of seven primary objects including:

**logFSdb** : parameters used to create the FSDB object

**PrecursorMZ** : A vector of precursor m/z values

**Retention Time** : A vector of retention time values

**Num Peaks** : A vector of num peaks values indicating number of ions for each fragment spectra

**Spectral Entropy** : A vector of spectral entropy values

**FragmentList** : A list of fragment ions

**MSPLibraryParameters** : A dataframe of tabulated headers and their values for each msp block

##
Visit https://fsa.idsl.me/ for the detailed documentation and tutorial.
