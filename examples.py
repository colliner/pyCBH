import sys, os
import pycbh
from rdkit import Chem
import numpy as np


fn_smi='sample_inputs/example.smi'
fn_xyz='sample_inputs/example.xyz'
fn_mol='sample_inputs/example.mol'


def test_xyzfile2cbh():
  print('coordinates file(s) (.xyz) --> xyz2cbh --> cbh-n')
  cbh=pycbh.xyz2cbh(fn_xyz, 2)
  pycbh.cbh_print(cbh)

def test_xyzfile2cbh_batch():
  print('specifying multiple files and multiple rungs')
  cbh=pycbh.xyz2cbh([fn_xyz,fn_xyz,fn_xyz], [0,1,2])
  pycbh.cbh_print(cbh)

def test_smistr2cbh():
  print('SMILES string(s) --> smi2cbh --> cbh-n')
  smi_str = 'O=Cc1ccc(O)c(OC)c1'
  cbh=pycbh.smi2cbh(smi_str, 0)
  pycbh.cbh_print(cbh)

def test_smistr2cbh_batch():
  print('specifying multiple rungs of CBH')
  smi_str = 'O=Cc1ccc(O)c(OC)c1'
  cbh=pycbh.smi2cbh([smi_str], [0,1,2,3,4,5,6])
  pycbh.cbh_print(cbh)

def test_molfile2cbh():
  print('RDKit .mol object file(s) --> mol2cbh --> cbh-n')
  cbh=pycbh.mol2cbh(fn_mol, 0)
  pycbh.cbh_print(cbh)

def test_mixedinp2cbh():
  print('General conversion from files:')
  print('  filename(s) (.smi | .xyz | .mol) --> files2cbh --> cbh-n')
  mixed_fns = [ fn_smi, fn_xyz, fn_mol ]
  rungs = [0,1,2]
  cbh=pycbh.files2cbh(mixed_fns, rungs)
  pycbh.cbh_print(cbh)

def test_cbhlookup():
  print('cbh reaction then energy lookup')
  smi_str = 'O=Cc1ccc(O)c(OC)c1'
  cbh=pycbh.smi2cbh(smi_str, [0,1,2,3])
  pycbh.cbh_print(cbh)

  methods=[['g4(0k)','zpe'],'pbe-d3',['g4(0k)','pbe-d3','zpe']]

  key_fn="fragment_lookup/keys.txt"
  energy_fn="fragment_lookup/lookup_energy.txt"

  cbh_e = pycbh.cbh_store2energy(cbh, key_fn=key_fn, energy_fn=energy_fn, levels_of_theory=methods)

  print('\nMethods : {}'.format(' '.join([x if type(x)!=list else '-'.join(x) for x in methods])))

  for key in sorted(cbh_e.keys()):
    if key != 'methods':
      print('{} {}'.format(key,' '.join([str(x) for x in cbh_e[key]])))

print('-'*50)
test_xyzfile2cbh()
input('-'*50)
test_xyzfile2cbh_batch()
input('-'*50)
test_smistr2cbh()
input('-'*50)
test_smistr2cbh_batch()
input('-'*50)
test_molfile2cbh()
input('-'*50)
test_mixedinp2cbh()
input('-'*50)
test_cbhlookup()
print('-'*50)
