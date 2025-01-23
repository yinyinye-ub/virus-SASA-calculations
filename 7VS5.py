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
structure = parser.get_structure("7VS5", "7VS5_clean.cif") #assmebly.cif

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
        my_list.append((chain.get_id()[0:2],res.get_id()[1],res.get_resname(),round(res.sasa,2)))

df = pd.DataFrame(my_list, columns = ['UniProtid','AAid','AA','SASA'])
#df.to_csv('5IV5_raw.csv', index=False)
# match auth chain - protien uniprot ID
arr = np.array(df['UniProtid'])

#update the following with matched auth chain # and uniprotID
arr = np.where(np.isin(arr, ['aa','ab','ac','ad','ae','af','ag','ah','ai','aj','ak','al','am','an','ao','ap','aq','ar','as','at','au','av','aw','ax','ay','az','ba','bb','bc','bd','be','bf','bg','bh','bi','bj','bk','bl','bm','bn','bo','bp','bq','br','bs','bt','bu','bv','bw','bx','by','bz','ca','cb','cc','cd','ce','cf','cg','ch','ci','cj','ck','cl','cm','cn','co','cp','cq','cr','cs','ct','cu','cv','cw','cx','cy','cz','da','db','dc','dd','de','df','dg','dh','di','dj','dk','dl','dm','dn','do','dp','dq','dr','ds','dt','du','dv','dw','dx','dy','dz','ea','eb','ec','ed','ee','ef','eg','eh','ei','ej','ek','el','em','en','eo','ep','eq','er','es','et','eu','ev','ew','ex','ey','ez','fa','fb','fc','fd','fe','ff','fg','fh','fi','fj','fk','fl','fm','fn','fo','fp','fq','fr','fs','ft','fu','fv','fw','fx','fy','fz','ga','gb','gc','gd','ge','gf','gg','gh','gi','gj','gk','gl','gm','gn','go','gp','gq','gr','gs','gt','gu','gv','gw','gx','gy','gz','ha','hb','hc','hd']), 'P04535', arr)
arr = np.where(np.isin(arr, ['he','hf','hg','hh','hi','hj','hk','hl','hm','hn','ho']), 'P19896', arr)
arr = np.where(np.isin(arr, ['hp','hq','hr','hs','ht','hu','hv','hw','hx','hy','hz','ia','ib','ic','id','ie','if','ig','ih','ii','ij','ik','il','im','in','io','ip','iq','ir','is','it','iu','iv','iw','ix','iy','iz','ja','jb','jc','jd','je','jf','jg','jh','ji','jj','jk','jl','jm','jn','jo','jp','jq','jr','js','jt','ju','jv','jw','jx','jy','jz','ka','kb','kc','kd','ke','kf','kg','kh','ki','kj','kk','kl','km','kn','ko','kp','kq','kr','ks','kt','ku','kv','kw','kx','ky','kz','la','lb','lc','ld','le','lf','lg','lh','li','lj','lk','ll','lm','ln','lo','lp','lq','lr','ls','lt','lu','lv','lw','lx','ly','lz','ma','mb','mc','md','me','mf','mg','mh','mi','mj','mk','ml','mm','mn','mo','mp','mq','ms','mt','mu','mv','mw','mx','my','mz','na','nb','nc','nd','ne','nf','ng','nh','ni','nj','nk','nl','nm','nn','no','np','nq','ns','nt','nu','nv','nw','nx','ny','nz','oa','ob','oc','od','oe','of','og']), 'P03715', arr)


df['UniProtid'] = arr

g_data = df.groupby(['UniProtid','AAid','AA'])['SASA'].sum().reset_index()
# save file-match cif name and csv name
g_data.to_csv('7VS5_1000.csv', index=False)

