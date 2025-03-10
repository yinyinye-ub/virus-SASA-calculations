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
structure = parser.get_structure("5MUU", "5MUU_clean.cif") #assmebly.cif

model = structure[0]

'''
chain_ids = []
for chain in model:
    chain_ids.append(chain.get_id())
print(chain_ids)
'''

sr = ShrakeRupley()
sr = ShrakeRupley(n_points=960)
sr.compute(model, level="R")

my_list = []
for chain in model:
    for res in chain:
        my_list.append((chain.get_id()[0],res.get_id()[1],res.get_resname(),round(res.sasa,2)))

df = pd.DataFrame(my_list, columns = ['UniProtid','AAid','AA','SASA'])
#df.to_csv('5IV5_raw.csv', index=False)
# match auth chain - protien uniprot ID
arr = np.array(df['UniProtid'])

#update the following with matched auth chain # and uniprotID
arr = np.where(np.isin(arr, ['A','B']), 'P11126', arr)
arr = np.where(np.isin(arr, ['C']), 'P11125', arr)
arr = np.where(np.isin(arr, ['D','E','F','G','H','I','J','K','L','M']), 'P07579', arr)

df['UniProtid'] = arr

g_data = df.groupby(['UniProtid','AAid','AA'])['SASA'].sum().reset_index()
# save file-match cif name and csv name
g_data.to_csv('7VS5_1000.csv', index=False)

