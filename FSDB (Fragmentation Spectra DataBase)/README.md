## FSDB (Fragmentation Spectra DataBase)
FSDB objects are R readable objects from one or several msp files. FSDB objects facilitate msp annotation and msp data parsing in the R environment. In this page we are presenting FSDB objects for GNPS and MoNA public libraries.

# msp2FSdb
`msp2FSdb` module of the IDSL.FSA is used to generate FSDB objects. This module was designed to be consistent with various msp structures particularly from NIST, GNPS, MoNA, IDSL.CSA libraries. The msp2FSdb module generally can work for any msp file as long as `Num Peaks` lines are available in the msp file.

# FSDB objects:
FSDB objects are R lists consisting of seven primary objects

logFSdb: parameters used to create the FSDB object

PrecursorMZ: A vector of Precursor m/z values

Retention Time: A vector of retention time values

Num Peaks: A vector of num peaks values indicating number of ions for each fragment spectra

Spectral Entropy: A vector of spectral entropy values

FragmentList: A list of fragment ions

MSPLibraryParameters: A dataframe of tabulated headers and their values for each msp block

