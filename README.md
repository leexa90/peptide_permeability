My stapled_peptide_permeability analysis 

install anaconda from 
https://www.anaconda.com/download/#linux

Create conda environment 
$ conda create -n stapled_peptide python=3.5  --file spec-file.txt

Activate Conda environment
$ source activate stapled_peptide

Install other packages
$ pip install -r requirements.txt

Add environment to jupyter notebook
$ python -m ipykernel install --user --name stapled_peptide --display-name "stapled_peptide"

Deactivate environment
$ source deactivate stapled_peptide
