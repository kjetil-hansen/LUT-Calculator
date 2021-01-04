## Peptide collision energy calculator for IMS-HDX-MS

*last updated by KJ Hansen on 04/01/2021*

The 'HDX-MS LUT calculator' is a python-based program to be used with hydrogen deuterium exchange (HDX) mass spectrometry (MS) experiments. The program handles large peptide .csv files containing a list of peptides and associated MS data. The .csv files (produced by PLGS) are combined and linear regression performed between the theoretical collision energy (CE) and experimental mobility. Ultimately, the optimal CE for each peptide is calculated and a look up table (LUT) produced. This LUT can be utilised on Waters MS instruments with ion mobility for increased HDX peptide detection.

Reference: https://doi.org/10.1021/jasms.9b00133

---
### Usage

#### LUT calculator requirements:

- Python 3.0+
- ProteinLynx Global Server (Waters PLGS)

#### Flag Usage:

-f 	path/folder	'filepath of PLGS ion accounting files (optional)'<br>
-s	protein_name		'name of protein (optional)'<br>
-l  'do not filter peptides (optional)'<br>

PLGS output files should be placed within the same parent folder of the python script. Alternatively, an -f flag can be used.

---

### Output

After processing PLGS ion accounting files, the LUT calculator will output 4 results files:

- LUT text file (.csv): to be input into the Tune page of a Waters instrument.
- Table of processed peptides (.csv) with calculated CE values.
- Image (.png) of LUT CE gradient.
- Image (.png) of linear regression for processed peptides.

#### Example Outputs
An example of the linear regression of processed peptides:

![Linear regression of peptides from LUTCalc](https://github.com/kjetil-hansen/LUT-Calculator/blob/master/Example_Output/LUT_2Charge_regression.png)

An example of a LUT CE gradient:

![Example LUT CE gradient](https://github.com/kjetil-hansen/LUT-Calculator/blob/master/Example_Output/LUT_2Charge_gradient.png)

Example input and output files are provided.
