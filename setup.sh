conda env create -f environment.yml
conda create -n stapled_peptide python=3.5  --file spec-file.txt
python -m ipykernel install --user --name stapled_peptide --display-name "stapled_peptide"

