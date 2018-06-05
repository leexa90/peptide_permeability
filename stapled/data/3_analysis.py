
# coding: utf-8

# # here is the analysiss of the processfile to find important residues

# In[1]:


import pandas as pd
pd.options.display.max_rows = 20
import numpy as np
import matplotlib.pyplot as plt
import collections
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from scipy import histogram, digitize, stats, mean, std
import sys
import statsmodels.api as sm
from sklearn import linear_model,preprocessing ,decomposition
from sklearn.decomposition import PCA
#from sklearn.preprocessing import StandardScaler #used custom one
from sklearn.linear_model import Lasso
from sklearn.metrics import r2_score


# In[2]:


test = pd.read_csv('./processed/stapled_peptide_permability_features.csv')
test = test[test['Peptide_type'] != 1].reset_index(drop=True)
test = test[test['sequence'] != '-AC-AHL-R8-LCLEKL-S5-GLV-(K-PEG1--FITC-)'].reset_index(drop=True)
test['Lg_Flr'] = np.log10(test['Mean_Flr_intensity'])
del test['permability']


# In[3]:


# function process sequence into list
def process_string(str,len=len):
    result = []
    special_resi = False
    for i in str:
        if i =='(':
            return len(result)
        elif special_resi and i !='-': #inside special residue
            temp += i
        elif i =='-':
            if special_resi : #closing of special residue
                result += [temp,]
                special_resi = False
                temp  = ''
            else: #starting of special residue
                special_resi = True 
                temp = ''
        else: #normal residue
            result += [i,]
            #if i in ['Q','N']:
                #result += [i,]
    return len(result)


# In[4]:


# more features 
if True:
    test['len']=test['sequence'].apply(process_string)
    test['res_list']=test['sequence'].apply(lambda x : process_string(x,list))
    test['list']=test['sequence'].apply(lambda x : np.array(process_string(x,list)))
    test['res_list_QN']=test['res_list'].apply(lambda x : collections.Counter(x)['Q']+collections.Counter(x)['N'])
    test['res_list']=test['res_list'].apply(len)
    dictt = collections.Counter(np.concatenate(test['list'].values))
    test['Aro_Ccycle_num'] = test['SMILES'].apply(lambda x :        Descriptors.NumAromaticCarbocycles(Chem.MolFromSmiles(x)))
    test['Aro_Hcycle_num'] = test['SMILES'].apply(lambda x :        Descriptors.NumAromaticHeterocycles(Chem.MolFromSmiles(x)))
    test['Aro_Ring_num'] = test['SMILES'].apply(lambda x :        Descriptors.NumAromaticRings(Chem.MolFromSmiles(x)))
    test['Ali_Ccycle_num'] = test['SMILES'].apply(lambda x :        Descriptors.NumAliphaticCarbocycles(Chem.MolFromSmiles(x)))
    test['Ali_Hcycle_num'] = test['SMILES'].apply(lambda x :        Descriptors.NumAliphaticHeterocycles(Chem.MolFromSmiles(x)))
    test['Ali_Ring_num'] = test['SMILES'].apply(lambda x :        Descriptors.NumAliphaticRings(Chem.MolFromSmiles(x)))
    for i in [x for x in test.keys() if ('Ali' in x or 'Aro' in x)]:
        test[i[:-4]+'_norm'] = test[i]/test['res_list']
    test['TPSA-DIV-LabuteASA'] = test['TPSA']/test['LabuteASA'] #ratio of polar SA


def get_surface_area(smile):
    #print (smile[-25:])
    mol0 = Chem.MolFromSmiles(smile)
    mol = Chem.AddHs(mol0)
    AllChem.Compute2DCoords(mol)
    adj = (Chem.GetDistanceMatrix(mol)==1)*1
    adj2 = (Chem.GetDistanceMatrix(mol)==2)*1
    molMMFF = AllChem.MMFFGetMoleculeProperties(mol)
    # Chem.MolSurf._LabuteHelper(mol) indiv contribution of surface area
    atoms = list(
                            map(lambda x: molMMFF.GetMMFFAtomType(x),
                                range(len(mol.GetAtoms()))
                                 )
                            )
    AllChem.ComputeGasteigerCharges(mol)
    charges = np.array([float(mol.GetAtomWithIdx(x).GetProp('_GasteigerCharge')) for x in range(len(atoms))])
    surf= np.array(Chem.MolSurf._LabuteHelper(mol))
    return (charges,surf[1:],atoms)
test['charge_surf'] = test['SMILES'].apply( get_surface_area)
test['charge_num'] = test['charge_surf'].apply(lambda x : np.sum(x[0]))
test['charge_mean'] = test['charge_surf'].apply(lambda x : np.mean(x[0]))

#combine similar amino acids together
dictt_vol = {'FITC': 163.32732360530238, 'PEG2': 75.972, 'pff': 84.691, 'Y': 68.658, 'I': 47.901,
             'E': 50.492, 'A': 28.806, 'P': 40.74, 'L': 47.901, 'AC': 16.891, 'B8': 71.665, 'B88': 146.863,
             'S5': 53.418, 'BAla': 28.806, 'PEG': 58.128, 'S': 33.6, 'D': 44.127, 'EEA': 93.299, 'F': 63.863,
             'R8': 72.513, 'T': 39.965, 'R5': 50.236, 'NL': 47.901, 'N': 44.673, 'Q': 51.038, 'R88': 115.722,
             'S8': 72.513, 'C': 39.961, 'B5': 71.665, 'G': 22.441, 'K': 53.241, 'R': 63.476, 'PEG1': 58.128,
             'W': 79.754, 'M': 51.659, 'V': 35.171, 'PEG5': 135.867, 'H': 56.292}
for i in [('W','F','Y','pff'),('N','Q'),('L','I','V','NL'),('D','E'),
          ('S8','R8','B8'),('R5','S5'),('PEG2','PEG5','PEG1'),
          ('R','K'),('H','R','K'),('S','T'),('S8','R8'),('S8','B8'),('R8','B8'),
          ('W','F','pff'),('F','pff'),('F','Y')]:
    dictt[i] = 999

#feature for each aa groups
for res in dictt.keys():
        test['%s_num' %str(res)] = test['list'].apply(lambda x : sum(np.in1d(x,np.array(res))))
        test['%s_norm' %str(res)] = 1.0*test['%s_num'%str(res)]/test['res_list']
        test['%s_normSA' %str(res)] = test['list'].apply(lambda y: np.array(list(map(lambda x : dictt_vol[x],y)))).apply(np.sum) #total area
        test['%s_normSA' %str(res)] = test['list'].apply(lambda y : sum(np.array(list(map(lambda x : dictt_vol[x],y)))[np.in1d(y,np.array(res))]))                                      /test['list'].apply(lambda y: np.array(list(map(lambda x : dictt_vol[x],y)))).apply(np.sum)

print ('frequencies of each amino acid：',dictt)


# In[5]:


# Import pairwise2 module
from Bio import pairwise2

# Import format_alignment method
from Bio.pairwise2 import format_alignment
dictt_counter = {}
counter =0
# convert the amino acids into single letter for sequence anaylsis,as some amino acids contain  > 1 letters
string='abcdefghijklmnopqrstuvwxyz1234567890!@#$%^&*()'
for i in sorted(dictt,key=lambda x : dictt[x])[::-1]:
  if type(i) != tuple:  
    dictt_counter[i] = string[counter]
    counter += 1
def func(x):
    result = ''
    for i in x:
        result += dictt_counter[i]
    return result
test['string'] = test['list'].apply(func)


# In[6]:


corr_mat = np.zeros((len(test),len(test)))
# Define two sequences to be aligned
for i in range(len(test)):
    for j in range(i,len(test)):
        alignments = pairwise2.align.globalms(test.iloc[i]['string'],
                                             test.iloc[j]['string'],
                                              2, -1, -0.5, -0.1)
        alignments = sorted(alignments,key= lambda x : x[-3])[-1]
        score = 2*np.sum(np.array(list(map(lambda x : x,alignments[0])))==np.array(list(map(lambda x : x,alignments[1]))))
        corr_mat[i,j] = score/(len(test.iloc[i]['string'])+len(test.iloc[j]['string']))
        corr_mat[j,i] = score/(len(test.iloc[i]['string'])+len(test.iloc[j]['string']))


# In[7]:


if True:
    import scipy.spatial.distance as ssd
    heat_map = abs(corr_mat-1)
    import scipy.cluster.hierarchy as sch
    from scipy.cluster.hierarchy import fcluster
    Y = sch.linkage(ssd.squareform(heat_map),method='centroid')
    Z = sch.dendrogram(Y,orientation='right')
    index = Z['leaves']
    D = heat_map[index,:]
    D = D[:,index]
    sorted_keys = [test.iloc[x]['sequence'] for x in index]
    plt.close()
    plt.xticks(range(len(index)),sorted_keys,rotation='vertical', fontsize=2)
    plt.yticks(range(len(index)),sorted_keys, fontsize=2)
    plt.imshow(D);plt.savefig('../reports/figures/cluster_peptide.png',dpi=300, bbox_inches='tight');#plt.show()
    plt.close()
    Z = sch.dendrogram(Y,orientation='right')
    plt.yticks(np.array(range(1,1+len(index)))*10-5,sorted_keys, fontsize=2)
    plt.savefig('../reports/figures/cluster_peptide.png',dpi=300)
plt.show()


# In[8]:


#plot each residue vs log flourscnecne
import scipy
for i in list(test.keys()):
    if 'num' in i and '(' in i :
        plt.plot(test[i],test['Lg_Flr'],'go',markersize=.5)
        plt.xlabel(i)
        plt.ylabel('log flr')
        pval = scipy.stats.pearsonr(test[i],test['Lg_Flr'])
        if pval[1] > 0.05:
            continue
        plt.title('corr coefficient %s \npval %s'%(pval[0],pval[1]))
        plt.savefig('../reports/figures/'+i+'.png',dpi=300)
        plt.show()
        plt.close()


# In[9]:


class StandardScaler(object):
    #custom weighted scale function
    def fit(self,X,weights):
        if type(weights)==type(None):
            mean = np.mean(X, 0).values
            std = np.std(X, 0).values
        else:
            weights= np.array([weights,]*X.shape[1]).T
            mean = np.sum(np.array(weights)*np.array(X),0)/np.sum(np.array(weights),0)
            std = np.sum(np.array(weights)*((np.array(X)-mean)**2),0)/np.sum(np.array(weights),0)
            std = std**.5
        self.mean = mean
        self.std = std
        return self
    def transform(self,X):
        return (X.values-self.mean)/self.std


# In[10]:


def plot_LR_coefficient(results,cutoff,test):
    reg_coef = pd.DataFrame(results,columns=                             ['name','simple LR model R2','corr coef','Se','CI','present in ids']).                            sort_values('simple LR model R2').reset_index(drop=True)
    features = []
    if True:
        plt.close()
        counter = 1
        plots = []
        for i in [x for x in sorted(results,key =  lambda x : x[-4])[-250:]]:
                if i[2] < -0.01 and 'num' in i[0] and '.(' not in i[0]:
                        features += [i[0],]
                        plots += [i,]
                        plt.errorbar([i[-4],],[counter,]*1,xerr=i[-3],fmt='r')
                        plt.plot([i[-4],],[counter,]*1,'ro')
                        if i[-4]-1.19*i[-3] > 0 or  i[-4]+1.19*i[-3] < 0:
                            txt = i[0]+'***'
                        elif i[-4]-i[-3] > 0 or i[-4]+i[-3] < 0:
                            txt = i[0]+'**'
                        else:
                            txt = i[0]+'*'
                        plt.text(i[-4],counter+.25,txt,horizontalalignment='center',fontsize=5)
                        counter += 1
        try:
            i = [x for x in sorted(results,key =  lambda x : x[-4]) if 'res_list' in x[0]][0]
            plots += [i,]
            plt.errorbar([i[-4],],[counter,]*1,xerr=i[-3],fmt='g')
            plt.plot([i[-4],],[counter,]*1,'go')
            if i[-4]-1.19*i[-3] > 0 or  i[-4]+1.19*i[-3] < 0:
                txt = i[0]+'***'
            elif i[-4]-i[-3] > 0 or i[-4]+i[-3] < 0:
                txt = i[0]+'**'
            else:
                txt = i[0]+'*'
            plt.text(i[-4],counter+.25,'Number of residues**',horizontalalignment='center',fontsize=5)
            counter += 1
        except : None
        for i in [x for x in sorted(results,key =  lambda x : x[-4])[-250:]]:
                if i[2] > 0.01 and 'num' in i[0] and '.(' not in i[0]:
                        features += [i[0],]
                        plots += [i,]
                        plt.errorbar([i[-4],],[counter,]*1,xerr=i[-3],fmt='b')
                        plt.plot([i[-4],],[counter,]*1,'bo')
                        if i[-4]-1.19*i[-3] > 0 or  i[-4]+1.19*i[-3] < 0:
                            txt = i[0]+'***'
                        elif i[-4]-i[-3] > 0 or i[-4]+i[-3] < 0:
                            txt = i[0]+'**'
                        else:
                            txt = i[0]+'*'
                        plt.text(i[-4],counter+.25,txt,horizontalalignment='center',fontsize=5)
                        counter += 1
        plt.title('Number Features, cutoff: '+str(cutoff))
        plt.xlabel('Coefficient and 95% CI of feature')
        plt.ylabel('Number Features')
        plt.yticks([],[])
        plt.ylim([0,counter])
        plt.savefig('../reports/figures/B_coef_cutoff:_'+str(cutoff)+'.png',dpi=300),plt.show()
        plt.close()
        counter = 1
        plots = []
        for i in [x for x in sorted(results,key =  lambda x : x[-4])[-250:]]:
                if i[2] < -0.01 and 'num' not in i[0] and 'list' not in i[0] and 'normSA' not in i[0] and 'VSA' not in i[0]:
                        features += [i[0],]
                        plots += [i,]
                        plt.errorbar([i[-4],],[counter,]*1,xerr=i[-3],fmt='r')
                        plt.plot([i[-4],],[counter,]*1,'ro')
                        if i[-4]-1.19*i[-3] > 0 or  i[-4]+1.19*i[-3] < 0:
                            txt = i[0]+'***'
                        elif i[-4]-i[-3] > 0 or i[-4]+i[-3] < 0:
                            txt = i[0]+'**'
                        else:
                            txt = i[0]+'*'
                        plt.text(i[-4],counter+.25,txt,horizontalalignment='center',fontsize=5)
                        counter += 1
        for i in [x for x in sorted(results,key =  lambda x : x[-4])[-250:]]:
                if i[2] > 0.01 and 'num' not in i[0] and 'list' not in i[0] and 'normSA' not in i[0] and 'VSA' not in i[0]:
                        features += [i[0],]
                        plots += [i,]
                        plt.errorbar([i[-4],],[counter,]*1,xerr=i[-3],fmt='b')
                        plt.plot([i[-4],],[counter,]*1,'bo')
                        if i[-4]-1.19*i[-3] > 0 or  i[-4]+1.19*i[-3] < 0:
                            txt = i[0]+'***'
                        elif i[-4]-i[-3] > 0 or i[-4]+i[-3] < 0:
                            txt = i[0]+'**'
                        else:
                            txt = i[0]+'*'
                        plt.text(i[-4],counter+.25,txt,horizontalalignment='center',fontsize=5)
                        counter += 1
        plt.ylim([0,counter])
        plt.title('Normalized (by seq length) Features, cutoff: '+str(cutoff))
        plt.xlabel('Coefficient and 95% CI of feature')
        plt.ylabel('Percentage Features')
        plt.yticks([],[])
        plt.savefig('../reports/figures/B_norm_cutoff:_'+str(cutoff)+'.png',dpi=300),plt.show()    


    dataA = []
    all_pred = []
    for seed in range(1,3):
        if seed ==1 :
            a,b = 0,1
        else:
            a,b = 1,0
        replace = False
        test1=test.sort_values(['weight0','Mean_Flr_intensity']).iloc[a::2].reset_index(drop=True)
        test2=test.sort_values(['weight0','Mean_Flr_intensity']).iloc[b::2].reset_index(drop=True)
        X,Y1,Y2 = [],[],[]
        for alpha in np.linspace(-5,-2,5)[-2:-1]:# [0.03,0.1,0.3,1.0,3.0]:
            alpha  = -7 #made small, smillar to OWLS
            alpha =10**alpha
            from scipy.optimize import minimize
            temptrain  = StandardScaler().fit(test1[features],weights=test1['weight']).transform(test1[features])
            X = np.concatenate([temptrain[:,0:1]*0+1,temptrain],1)
            B = np.zeros((X.shape[1],))
            W= test1[['weight',]].values
            y = test1[['Lg_Flr',]].values
            def weighted_lasso(B,X=X,y=y,W=W,alpha=alpha):
                B = np.reshape(B,(B.shape[0],1))
                score = W*(y-np.matmul(X,B))**2
                score = np.sum(score) /np.sum(W)
                score = score + np.sum(np.abs(B*alpha))
                return score
            B = minimize(weighted_lasso,B)['x']
            temptrain  = StandardScaler().fit(test1[features],weights=test1['weight']).transform(test1[features])
            X = np.concatenate([temptrain[:,0:1]*0+1,temptrain],1)
            preds = np.matmul(X,np.reshape(B,(1+len(features),1)))[:,0]
            temptrain  = StandardScaler().fit(test1[features],weights=test1['weight']).transform(test2[features])
            X = np.concatenate([temptrain[:,0:1]*0+1,temptrain],1)
            preds2 = np.matmul(X,np.reshape(B,(1+len(features),1)))[:,0]
            a,b=np.mean((preds -np.log10(test1['Mean_Flr_intensity']))**2)**.5,np.mean((preds2 -np.log10(test2['Mean_Flr_intensity']))**2)**.5
            #print (alpha,a,b)
            X += [alpha,]
            Y1 += [a,]
            Y2 += [b,]    
            r1 = 1-np.mean((preds -test1['Lg_Flr'])**2)/np.mean((np.mean(test1['Lg_Flr']) -test1['Lg_Flr'])**2)
            r2 = r2 = str(np.corrcoef(preds2 ,np.log10(test2['Mean_Flr_intensity']))[0,1])[:4]#str(np.mean((preds2 -np.log10(test2['Mean_Flr_intensity']))**2)**.5)[:4]#str(1-np.mean((10**preds2 -test2['Mean_Flr_intensity'])**2)/np.mean((np.mean(test2['Mean_Flr_intensity']) -test2['Mean_Flr_intensity'])**2))[:4]
            dataA += [r2,]
            #print (r2)
            #print (alpha,a,b,r1,r2)
            if True:
                plt.plot([2.4,3.8],[2.4,3.8],'g',label='y = x');
                plt.plot(preds,np.log10(test1['Mean_Flr_intensity']),'ro',label='train predictions');
                plt.ylim([2.4,3.8]);plt.xlim([2.4,3.8]);
                plt.title('Combined_model\corr coref = %s'%r2);
                plt.xlabel('predictions of log flourescence');plt.ylabel('ground truth of log flourescence');
                plt.plot(preds2,np.log10(test2['Mean_Flr_intensity']),'bo',label='test predictions');#plt.legend()
                plt.savefig('../reports/figures/Combined_model_cutoff:_'+str(cutoff)+'.png');
                plt.legend()
                plt.show()


# In[11]:


for cutoff in [0.01,0.1,0.15,0.2,0.3]:
    #weights for each cluster == 1     
    cluster = fcluster(Y, cutoff, criterion='distance')
    test = test.set_value(test.index,'weight0',cluster)
    test['weight'] = 1.0/test['weight0'].map(collections.Counter(fcluster(Y, cutoff, criterion='distance')))#1/(1+rho*(len(y)/(len(pd.unique(cluster))-1)))

    dictt = collections.Counter(np.concatenate(test['list'].values))

    for i in test.keys():
        try:
            if np.std(test[i])==0:
                del test[i]
        except : None
    dictt_name = {}
    if True:
        results = []
        for i in list(test.keys())[15:]:
            if 'weight' not in i and 'Lg_Flr' not in i:
                try:
                    X,Y1,Y2 = [],[],[]
                    clf = linear_model.LinearRegression()
                    temptrain  = StandardScaler().fit(test[[i,]],weights=test['weight']).transform(test[[i,]])
                    ids_non_zero = test[test[i]!=0].index
                    if 'num' in i or 'res_list' in i or i in ['atom_1', 'atom_6', 'atom_7', 'atom_8', 'atom_9', 'atom_16', 'atom_len', 'atom_size']: #do not scale count features. 
                        temptrain = test[[i,]].values #not scaled
                    temptrain2 = np.concatenate([temptrain*0+1,temptrain],1)
                    preds = clf.fit(temptrain,np.log10(test['Mean_Flr_intensity']),sample_weight=test['weight']).predict(temptrain)
                    # equailvalent ~ inv(X.T*W*X)*X.T*W*Y, analytic form of linear regression
                    W=np.diag(test['weight'])
                    X=temptrain2
                    y=np.log10(test['Mean_Flr_intensity'])
                    B=np.matmul(np.linalg.inv(np.matmul(np.matmul(X.T,W),X)) ,np.matmul(np.matmul(X.T,W),y))
                    preds = np.matmul(B,X.T)
                    se = np.sum(
                        (test['weight']*(preds-np.log10(test['Mean_Flr_intensity']))**2)
                        )*np.linalg.inv(
                             np.matmul(np.matmul(temptrain2.T,np.diag(test['weight'])),temptrain2)
                             )[-1,-1]
                    dof = sum(test['weight'])**2/sum(test['weight']**2)
                    se = 2.0*(se/dof)**.5

                    if  (clf.coef_[0]-2*se/2 >0 or clf.coef_[0]+2*se/2  <0) and r2_score(np.log10(test['Mean_Flr_intensity']),preds,sample_weight=test['weight']) < 0.7:
                        dictt_name[i] = (i,r2_score(np.log10(test['Mean_Flr_intensity']),preds,sample_weight=test['weight']),
                                 clf.coef_[0],se,[clf.coef_[0]-se,clf.coef_[0]+se],ids_non_zero)
                        results += [(i,r2_score(np.log10(test['Mean_Flr_intensity']),preds,sample_weight=test['weight']),
                                     clf.coef_[0],se,[clf.coef_[0]-se,clf.coef_[0]+se],ids_non_zero)]

                except  : print (i)
        print ([x[:3] for x in sorted(results,key =  lambda x : x[1])[-20:]])
        plot_LR_coefficient(results,cutoff,test)
        


# In[12]:



for cutoff in (0.01,0.1,0.15,0.2):
    cluster = fcluster(Y, cutoff, criterion='distance')
    test = test.set_value(test.index,'weight0',cluster)
    test['weight'] = 1.0/test['weight0'].map(collections.Counter(fcluster(Y, cutoff, criterion='distance')))#1/(1+rho*(len(y)/(len(pd.unique(cluster))-1)))

    def makePlot(i) :
            plt.close()
            temp = test[['Lg_Flr',i,'weight0']].groupby('weight0')[i].apply(np.mean).reset_index()
            temp = pd.merge(temp, test[['Lg_Flr',i,'weight0']].groupby('weight0')['Lg_Flr'].apply(np.mean).reset_index(),
                            on = ['weight0'])
            temptrain  = StandardScaler().fit(test[[i,]],test['weight']).transform(test[[i,]])
            if 'num' in i or 'res_list' in i or i in ['atom_1', 'atom_6', 'atom_7', 'atom_8', 'atom_9', 'atom_16', 'atom_len', 'atom_size']: #do not scale count features. 
                temptrain = test[[i,]].values
                print  (i)
            temptrain2 = np.concatenate([temptrain*0+1,temptrain],1)
            W=np.diag(test['weight'])
            X=temptrain2
            y=np.log10(test['Mean_Flr_intensity'])
            B=np.matmul(np.linalg.inv(np.matmul(np.matmul(X.T,W),X)) ,np.matmul(np.matmul(X.T,W),y))
            preds = np.matmul(B,X.T)
            se = np.sum(
                (test['weight']*(preds-np.log10(test['Mean_Flr_intensity']))**2)
                )*np.linalg.inv(
                     np.matmul(np.matmul(temptrain2.T,np.diag(test['weight'])),temptrain2)
                     )
            dof = sum(test['weight'])**2/sum(test['weight']**2)
            se = 2.0*(se/dof)
            
            err2=np.sum(
                (test['weight']*(preds-np.log10(test['Mean_Flr_intensity']))**2)
                )/(sum(test['weight'])-2)
            err2 = err2**.5
            err=np.diag(np.matmul(np.matmul(temptrain2,2*se),temptrain2.T))**.5
            B = np.round(B,3)
            plt.plot(test[i],preds,'r',label='y = %sX + %s'%(B[1],B[0]));
            plt.plot(temp[i],temp['Lg_Flr'],'bo',label='data points from clusters');
            plt.ylabel('Log flourescence');
            #plt.errorbar(test[i],preds,color='green',yerr=err2*2,label='regresion sigma');
            plt.errorbar(test[i],preds,color='red',yerr=err*2);
            axes = plt.axes()
            xlim = axes.get_xlim()
            # example of how to zoomout by a factor of 0.1
            factor = 0.1
            plt.title(i+' cutoff: '+str(cutoff),fontsize=24)
            new_xlim = (xlim[0] + xlim[1])/2 + np.array((-0.5, 0.5)) * (xlim[1] - xlim[0]) * (1 + factor) 
            axes.set_xlim(new_xlim)
            plt.xlabel(i);plt.legend();plt.savefig('../reports/figures/%s.png'%(str(i)+str(cutoff)),bbox_inches='tight' );plt.show()

    makePlot('R_num')
    makePlot("('W', 'F', 'pff')_num")
    makePlot("('L', 'I', 'V', 'NL')_num")
    makePlot("len")


# In[13]:


sum(test['weight'])**2/sum(test['weight']**2)

