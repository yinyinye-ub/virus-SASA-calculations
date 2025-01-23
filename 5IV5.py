# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 14:39:41 2024

@author: yinyinye
This one is to calculate SASA from mmCIF file
"""

# %% final version
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.MMCIFParser import MMCIFParser
import pandas as pd
import numpy as np

#parse the MMCIF structure
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("5IV5", "5IV5_clean.cif") #assmebly.cif

model = structure[0]


chain_ids = []
for chain in model:
    chain_ids.append(chain.get_id())
print(chain_ids)

sr = ShrakeRupley()
sr = ShrakeRupley(n_points=960)
sr.compute(model, level="R")

my_list = []
for chain in model:
    for res in chain:
        my_list.append((chain.get_id(),res.get_id()[1],res.get_resname(),round(res.sasa,2)))

df = pd.DataFrame(my_list, columns = ['UniProtid','AAid','AA','SASA'])
#df.to_csv('5IV5_raw.csv', index=False)
# match auth chain - protien uniprot ID
arr = np.array(df['UniProtid'])

#update the following with matched auth chain # and uniprotID
arr = np.where(np.isin(arr, ['A','B','BH','BI','EA','EB','GD','GE','X','Y','u','v']), 'P19060', arr)
arr = np.where(np.isin(arr, ['BJ','C','EC','GF','Z','w']), 'P19061', arr)
arr = np.where(np.isin(arr, ['CA','CB','D','E','ED','EE','GG','GH','a','b','x','y']), 'P19062', arr)
arr = np.where(np.isin(arr, ['AA','AB','CC','CD','CE','EF','EG','EH','F','G','GI','GJ','H','HA','c','d','e','z']), 'P10927', arr)
arr = np.where(np.isin(arr, ['AC','AD','AE','CF','CG','CH','EI','EJ','FA','HB','HC','HD','I','J','K','f','g','h']), 'P10928', arr)
arr = np.where(np.isin(arr, ['AF','AG','AH','CI','CJ','DA','FB','FC','FD','HE','HF','HG','L','M','N','i','j','k']), 'P10929', arr)
arr = np.where(np.isin(arr, ['AI','AJ','BA','DB','DC','DD','FE','FF','FG','HH','HI','HJ','O','P','Q','l','m','n']), 'P10930', arr)
arr = np.where(np.isin(arr, ['BB','BC','DE','DF','FH','FI','IA','IB','R','S','o','p']), 'P13333', arr)
arr = np.where(np.isin(arr, ['BD','DG','FJ','IC','T','q']), 'P09425', arr)
arr = np.where(np.isin(arr, ['BE','DH','GA','ID','U','r']), 'P13339', arr)
arr = np.where(np.isin(arr, ['BF','DI','GB','IE','V','s']), 'P16011', arr)
arr = np.where(np.isin(arr, ['BG','DJ','GC','IF','W','t']), 'P13341', arr)
arr = np.where(np.isin(arr, ['YA','YB','YC']), 'P16009', arr)
arr = np.where(np.isin(arr, ['YD','YE','YF']), 'P17172', arr)
arr = np.where(np.isin(arr, ['ZA']), 'P39234', arr)


df['UniProtid'] = arr

g_data = df.groupby(['UniProtid','AAid','AA'])['SASA'].sum().reset_index()
# save file-match cif name and csv name
g_data.to_csv('5IV5_1000.csv', index=False)

