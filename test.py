import sys, os
import pycbh
from rdkit import Chem

print_cbh=True
#print_cbh=False
if print_cbh:
  from pycbh import cbh_print
else:
  cbh_print = lambda c : print(['FAIL','PASS','PASS (batch)'][len(c[0:2])])

fn_smi='sample_inputs/example.smi'
fn_xyz='sample_inputs/example.xyz'
fn_mol='sample_inputs/example.mol'

#cbh_print(pycbh.smi2cbh('NC(Cc1ccc(O)cc1)C(=O)O', 2))
#sys.exit()
'''
coordinates file(s) (.xyz) --> xyz2cbh --> cbh-n
'''
cbh=pycbh.xyz2cbh(fn_xyz, 2)
cbh_print(cbh)
input('\npress enter to continue...\n')

'''
SMILES string(s) --> smi2cbh --> cbh-n
'''
smi_str = 'O=Cc1ccc(O)c(OC)c1'
cbh=pycbh.smi2cbh(smi_str, [0,1,2,3,4,5,6])
cbh_print(cbh)

input('\npress enter to continue...\n')
'''
RDKit .mol object(s) --> mol2cbh --> cbh-n
'''
m1 = Chem.MolFromMolFile(fn_mol)
m2 = Chem.MolFromSmiles(pycbh.smifromsmifile(fn_smi))
m3 = Chem.RemoveHs(pycbh.molfromxyz('sample_inputs/example_zwit.xyz')[-1])


#for m in [m1,m2,m3]:
#  print(Chem.MolToSmiles(m))
#sys.exit()

cbh=pycbh.mol2cbh([m1,m2,m3], [0,1])
cbh_print(cbh)
input('\npress enter to continue...\n')
'''
General conversion from files:
  filename(s) (.smi | .xyz | .mol) --> files2cbh --> cbh-n
'''
mixed_fns = [ fn_smi, fn_xyz, fn_mol ]

cbh=pycbh.files2cbh(mixed_fns, 0)
cbh_print(cbh)

 
