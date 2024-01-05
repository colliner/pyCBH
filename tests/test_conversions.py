import sys, os
import pycbh
from rdkit import Chem
import pytest
import numpy as np

fn_smi = 'sample_inputs/example.smi'
fn_xyz = 'sample_inputs/example.xyz'
fn_mol = 'sample_inputs/example.mol'


def test_xyzfile2cbh():
    '''
  coordinates file(s) (.xyz) --> xyz2cbh --> cbh-n
  '''
    cbh = pycbh.xyz2cbh(fn_xyz, 2)
    assert np.array(cbh).shape == (1, 2)


def test_xyzfile2cbh_batch():
    cbh = pycbh.xyz2cbh([fn_xyz, fn_xyz, fn_xyz], [0, 1, 2])
    assert np.array(cbh).shape == (9, 2)


def test_smistr2cbh():
    '''
  SMILES string(s) --> smi2cbh --> cbh-n
  '''
    smi_str = 'O=Cc1ccc(O)c(OC)c1'
    cbh = pycbh.smi2cbh(smi_str, 0)
    assert np.array(cbh).shape == (1, 2)


def test_smistr2cbh_batch():
    smi_str = 'O=Cc1ccc(O)c(OC)c1'
    cbh = pycbh.smi2cbh([smi_str], [0, 1, 2, 3, 4, 5, 6])
    assert np.array(cbh).shape == (7, 2)


def test_molfile2cbh():
    '''
  RDKit .mol object file(s) --> mol2cbh --> cbh-n
  '''
    cbh = pycbh.mol2cbh(fn_mol, 0)
    assert np.array(cbh).shape == (1, 2)


def test_molobj2cbh_batch():
    '''
  RDKit .mol object(s) --> mol2cbh --> cbh-n
  '''
    m1 = Chem.MolFromMolFile(fn_mol)
    m2 = Chem.MolFromSmiles(pycbh.smifromsmifile(fn_smi))
    m3 = Chem.RemoveHs(pycbh.molfromxyz('sample_inputs/example_zwit.xyz')[-1])
    cbh = pycbh.mol2cbh([m1, m2, m3], [0, 1])
    assert np.array(cbh).shape == (6, 2)


def test_mixedinp2cbh():
    '''
  General conversion from files:
    filename(s) (.smi | .xyz | .mol) --> files2cbh --> cbh-n
  '''
    mixed_fns = [fn_smi, fn_xyz, fn_mol]
    rungs = [0, 1, 2]
    cbh = pycbh.files2cbh(mixed_fns, rungs)
    assert np.array(cbh).shape == (9, 2)
