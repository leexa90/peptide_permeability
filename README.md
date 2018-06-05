## My peptide_permeability analysis. 
The code below reproduces my analysis. The instructions are for ubuntu  14/16 os

### Install git and clone repo
```
sudo apt install git-all
git clone git@github.com:leexa90/peptide_permeability.git
```

### Install anaconda from 
```
https://www.anaconda.com/download/#linux
```

### Create conda environment 
```
  conda create -n stapled_peptide python=3.5  --file spec-file.txt
```

### Activate Conda environment
```
  source activate stapled_peptide
```

### Install other essential packages
```
  pip install -r requirements.txt
```

### Add environment to jupyter notebook
```
  python -m ipykernel install --user --name stapled_peptide --display-name "stapled_peptide"
```
### Perform analysis (./data)
```
bash analysis.sh
```
### Deactivate environment
```
  source deactivate stapled_peptide
```
