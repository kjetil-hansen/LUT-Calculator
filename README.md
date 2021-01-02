## Peptide collision energy calculator for IMS-HDX-MS

*last updated by KJ Hansen on 02/01/2021*

The 'HDX-MS LUT calculator' is a python-based program to be used with hydrogen deuterium exchange (HDX) mass spectrometry (MS) experiments. The program handles large peptide .csv files containing a list of peptides and associated MS data. The .csv files (produced by PLGS) are combined and linear regression performed between the theoretical collision energy (CE) and experimental mobility. Ultimately, the optimal CE for each peptide is calculated and a look up table (LUT) produced. This LUT can be utilised on Waters MS instruments with ion mobility for increased HDX peptide detection.

Reference: https://doi.org/10.1021/jasms.9b00133

---

#### LUT calculator requirements:

- Python 3.0+
- ProteinLynx Global Server (Waters PLGS)

#### Flag Usage:

-f 	folder/folder/folder	'filepath of PLGS ion accounting files (optional)'
-s	protein_name		'name of protein (optional)'
-l  'do not filter peptides (optional)'
