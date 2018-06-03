import pandas as pd
train = pd.read_csv('./interim/stapled_200.csv')
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import  numpy as np
import scipy
from scipy import stats


for i in range(1,15):
    print ("train['PEOE_VSA%s'] = list(map( lambda smile : Descriptors.PEOE_VSA%s(smile),train['SMILES'].values))" %(i,i))
for i in range(1,11):
    print ("train['SMR_VSA%s'] = list(map( lambda smile : Descriptors.SMR_VSA%s(smile),train['SMILES'].values))" %(i,i))
for i in range(1,12):
    print ("train['SlogP_VSA%s'] = list(map( lambda smile : Descriptors.SlogP_VSA%s(smile),train['SMILES'].values))" %(i,i))
for i in range(1,11):
    print ("train['VSA_EState%s'] = list(map( lambda smile : Descriptors.VSA_EState%s(smile),train['SMILES'].values))" %(i,i))

def make_info(x,num=None):
    print ('making atom: ' , x)
    dictt = { 1.0 : 0 , 1.5 : 1 , 2.0:2}
    mol0 = Chem.MolFromSmiles(x)
    mol = Chem.AddHs(mol0)
    AllChem.Compute2DCoords(mol)
    atoms = list(map(lambda x:x.GetAtomicNum(),list(mol.GetAtoms())))
    AllChem.ComputeGasteigerCharges(mol)
    charges = [float(mol.GetAtomWithIdx(x).GetProp('_GasteigerCharge')) for x in range(len(atoms))]
    if num == None:
        return  [mol,None,atoms,mol.GetPropsAsDict(),charges]
    else:
        return [mol,temp.astype(np.float32),atoms,mol.GetPropsAsDict(),charges][num]

train['all']= train['SMILES'].apply(lambda x : make_info(x))
train['atoms']= train['all'].apply(lambda x : x[2])
if True:
    train['charges'] = train['all'].apply(lambda x : x[-1])
    train['charge_skew'] = train['charges'].apply(stats.skew)
    train['charge_kurtosis'] = train['charges'].apply(stats.kurtosis)
    train['charge_std']= train['charges'].apply(np.std)
    train['SMILES_temp'] = train['SMILES']
    train['SMILES'] = train['all'].apply(lambda x : x[0])
    for i in set(np.concatenate(train['atoms'].values,0)):
        train['atom_%s'%i] = train['atoms'].apply(lambda x : np.sum((np.array(x) ==i)))
        train['atom_sum_%s'%i] = train['atoms'].apply(lambda x : np.mean((np.array(x) ==i)))
       
    train['TPSA'] = list(map( lambda smile : Descriptors.TPSA(smile),
                              train['SMILES'].values))
    train['TPSA_norm'] = train['TPSA']/train['atoms'].apply(len) 
    train['LabuteASA'] = list(map( lambda smile : Descriptors.LabuteASA(smile),
                              train['SMILES'].values))
    train['LabuteASA_norm'] = train['LabuteASA']/train['atoms'].apply(len) 
    train['PEOE_VSA1'] = list(map( lambda smile : Descriptors.PEOE_VSA1(smile),train['SMILES'].values))
    train['PEOE_VSA2'] = list(map( lambda smile : Descriptors.PEOE_VSA2(smile),train['SMILES'].values))
    train['PEOE_VSA3'] = list(map( lambda smile : Descriptors.PEOE_VSA3(smile),train['SMILES'].values))
    train['PEOE_VSA4'] = list(map( lambda smile : Descriptors.PEOE_VSA4(smile),train['SMILES'].values))
    train['PEOE_VSA5'] = list(map( lambda smile : Descriptors.PEOE_VSA5(smile),train['SMILES'].values))
    train['PEOE_VSA6'] = list(map( lambda smile : Descriptors.PEOE_VSA6(smile),train['SMILES'].values))
    train['PEOE_VSA7'] = list(map( lambda smile : Descriptors.PEOE_VSA7(smile),train['SMILES'].values))
    train['PEOE_VSA8'] = list(map( lambda smile : Descriptors.PEOE_VSA8(smile),train['SMILES'].values))
    train['PEOE_VSA9'] = list(map( lambda smile : Descriptors.PEOE_VSA9(smile),train['SMILES'].values))
    train['PEOE_VSA10'] = list(map( lambda smile : Descriptors.PEOE_VSA10(smile),train['SMILES'].values))
    train['PEOE_VSA11'] = list(map( lambda smile : Descriptors.PEOE_VSA11(smile),train['SMILES'].values))
    train['PEOE_VSA12'] = list(map( lambda smile : Descriptors.PEOE_VSA12(smile),train['SMILES'].values))
    train['PEOE_VSA13'] = list(map( lambda smile : Descriptors.PEOE_VSA13(smile),train['SMILES'].values))
    train['PEOE_VSA14'] = list(map( lambda smile : Descriptors.PEOE_VSA14(smile),train['SMILES'].values))
    train['SMR_VSA1'] = list(map( lambda smile : Descriptors.SMR_VSA1(smile),train['SMILES'].values))
    train['SMR_VSA2'] = list(map( lambda smile : Descriptors.SMR_VSA2(smile),train['SMILES'].values))
    train['SMR_VSA3'] = list(map( lambda smile : Descriptors.SMR_VSA3(smile),train['SMILES'].values))
    train['SMR_VSA4'] = list(map( lambda smile : Descriptors.SMR_VSA4(smile),train['SMILES'].values))
    train['SMR_VSA5'] = list(map( lambda smile : Descriptors.SMR_VSA5(smile),train['SMILES'].values))
    train['SMR_VSA6'] = list(map( lambda smile : Descriptors.SMR_VSA6(smile),train['SMILES'].values))
    train['SMR_VSA7'] = list(map( lambda smile : Descriptors.SMR_VSA7(smile),train['SMILES'].values))
    train['SMR_VSA8'] = list(map( lambda smile : Descriptors.SMR_VSA8(smile),train['SMILES'].values))
    train['SMR_VSA9'] = list(map( lambda smile : Descriptors.SMR_VSA9(smile),train['SMILES'].values))
    train['SMR_VSA10'] = list(map( lambda smile : Descriptors.SMR_VSA10(smile),train['SMILES'].values))
    train['SlogP_VSA1'] = list(map( lambda smile : Descriptors.SlogP_VSA1(smile),train['SMILES'].values))
    train['SlogP_VSA2'] = list(map( lambda smile : Descriptors.SlogP_VSA2(smile),train['SMILES'].values))
    train['SlogP_VSA3'] = list(map( lambda smile : Descriptors.SlogP_VSA3(smile),train['SMILES'].values))
    train['SlogP_VSA4'] = list(map( lambda smile : Descriptors.SlogP_VSA4(smile),train['SMILES'].values))
    train['SlogP_VSA5'] = list(map( lambda smile : Descriptors.SlogP_VSA5(smile),train['SMILES'].values))
    train['SlogP_VSA6'] = list(map( lambda smile : Descriptors.SlogP_VSA6(smile),train['SMILES'].values))
    train['SlogP_VSA7'] = list(map( lambda smile : Descriptors.SlogP_VSA7(smile),train['SMILES'].values))
    train['SlogP_VSA8'] = list(map( lambda smile : Descriptors.SlogP_VSA8(smile),train['SMILES'].values))
    train['SlogP_VSA9'] = list(map( lambda smile : Descriptors.SlogP_VSA9(smile),train['SMILES'].values))
    train['SlogP_VSA10'] = list(map( lambda smile : Descriptors.SlogP_VSA10(smile),train['SMILES'].values))
    train['SlogP_VSA11'] = list(map( lambda smile : Descriptors.SlogP_VSA11(smile),train['SMILES'].values))
    train['VSA_EState1'] = list(map( lambda smile : Descriptors.VSA_EState1(smile),train['SMILES'].values))
    train['VSA_EState2'] = list(map( lambda smile : Descriptors.VSA_EState2(smile),train['SMILES'].values))
    train['VSA_EState3'] = list(map( lambda smile : Descriptors.VSA_EState3(smile),train['SMILES'].values))
    train['VSA_EState4'] = list(map( lambda smile : Descriptors.VSA_EState4(smile),train['SMILES'].values))
    train['VSA_EState5'] = list(map( lambda smile : Descriptors.VSA_EState5(smile),train['SMILES'].values))
    train['VSA_EState6'] = list(map( lambda smile : Descriptors.VSA_EState6(smile),train['SMILES'].values))
    train['VSA_EState7'] = list(map( lambda smile : Descriptors.VSA_EState7(smile),train['SMILES'].values))
    train['VSA_EState8'] = list(map( lambda smile : Descriptors.VSA_EState8(smile),train['SMILES'].values))
    train['VSA_EState9'] = list(map( lambda smile : Descriptors.VSA_EState9(smile),train['SMILES'].values))
    train['VSA_EState10'] = list(map( lambda smile : Descriptors.VSA_EState10(smile),train['SMILES'].values))
    train['HBA'] = list(map( lambda smile : Chem.rdMolDescriptors.CalcNumHBA(smile),train['SMILES'].values))
    train['HBD'] = list(map( lambda smile : Chem.rdMolDescriptors.CalcNumHBD(smile),train['SMILES'].values))
    train['Rotatable_num'] = list(map( lambda smile : Descriptors.NumRotatableBonds(smile),train['SMILES'].values))
    train['HBA_norm'] = train['HBA']/train['atoms'].apply(len) 
    train['HBD_norm'] = train['HBD']/train['atoms'].apply(len)
    train['SMILES'] = train['SMILES_temp']
    train['atom_len'] = train['atoms'].apply(len)
    train['atom_size'] = train['atoms'].apply(sum)
    del train['atoms']
    #del train['bond']
    del train['SMILES_temp']
    del train['charges']
    train['permability'] = train['Mean_Flr_intensity']
    train.to_csv('./processed/stapled_peptide_permability_features.csv',index=0)
