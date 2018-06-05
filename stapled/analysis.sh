echo 'Ensure you are in anaconda environment'
source deactivate
source activate stapled_peptide
cd data
python 1_get_smiles.py
python 2_get_RdkitFeatures.py
python 3_analysis.py


