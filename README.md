# xq-processing
Code to convert xQuest output to xiVIEW input


## Requisite packages
- pandas
- numpy
- matplotlib
- xml
- pyteomics


## Validate_XL_v2.py
This code is a modified version of ValidateXL available in this repository: https://github.com/ThalassinosLab/ValidateXL
The edits performed allow to:
  - Generate csv files for looplinks (still need to be optimised);
  - Generate csv files for monolinks (still need to be optimised);
  - Add the spectrum number for both light and heavy labelled version of a crosslink;
  - Add the retention time for both light and heavy labelled version of a crosslink;
  - Add the m/z of the light and heavy labelled version of a crosslink.
Refer to the original repository for the information.


## mgf_converter.py
This code allows to generate a csv file compatible with xiVIEW input constraints as well as a MGF file from several mzxml files used in the xQuest search.
3 parameters has to be added as arguments:
  - XXXX.csv
  This parameter is the Validate_XL_v2.py output file to be processed.
  - mzxml_folder XXXX/YYYY/
  This parameter is the folder where all the mzxml used for xQuest 
  
  
## Acknowledgement
We thank Charly Menneteau for his help in the debugging/optimisation of the mgf_converter.py
