### To run preprocessing to analysis 
```
source deactivate
source activate stapled_peptide
cd data
python 1_get_smiles.py
python 2_get_RdkitFeatures.py
python 3_analysis.py
source deactivate
```
### Report is in ./reports

### To view analysis
```
cd data
source deactivate
source activate stapled_peptide
cd data
jupyter notebook 3_analysis.ipynb
source deactivate
```
