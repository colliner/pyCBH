#!/usr/bin/python

from pycbh.utils import *

import pickle
from rdkit import Chem
import os, re, glob, sys

import numpy as np
import math
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem.Descriptors import MolWt

from xyz2mol import xyz2mol

import sys
import numpy as np

np.set_printoptions(threshold=sys.maxsize)
from collections import Counter

from pycbh.iep import iep


def load_frags(key_fn):
    with open(key_fn, "a"):
        pass
    frag_keys = dict()
    for y in [
            x.split(" ")
            for x in list(filter(None,
                                 open(key_fn, "r").read().split("\n")))
    ]:
        if y[0] in frag_keys:
            frag_keys[y[0]].append(y[1::])
        else:
            frag_keys[y[0]] = [y[1::]]
    return frag_keys


def cbh_store2vec(cbh_store):
    if type(cbh_store) != list:
        cbh_store = [cbh_store]
    if type(cbh_store[0]) != list:
        cbh_store = [cbh_store]
    vec = list()
    for cbh_s in cbh_store:
        #print(cbh_s)
        label = cbh_s[0] + ' ' + cbh_s[1].split('\n')[0]
        v = [x for x in cbh_s[1].split(' ') if len(x) > 0][2:]
        left, right = list(), list()
        '''
    ['+', '3', '[H][H]', '-->', 'CF', '+', 'CBr', '+', 'CCl']
    '''
        v.append('+')
        if len(v) > 1 and v[0] == '+':
            v = v[1:]
        if len(v) > 1 and v[-1] == '+':
            v = v[:-1]
        v = [[y
              for y in x.split(' + ')
              if len(y) > 0]
             for x in ' '.join(v).split('-->')]
        if len(v) > 1:

            for x in v[0]:
                #print('left',x)
                x = [y for y in x.split(' ') if len(y) > 0]
                if len(x) > 1:
                    for i in range(int(x[0])):
                        left.append(x[-1])
                else:
                    left.append(x[-1])
        for x in v[-1]:
            #print('right',x)
            x = [y for y in x.split(' ') if len(y) > 0]
            if len(x) > 1:
                for i in range(int(x[0])):
                    right.append(x[-1])
            else:
                right.append(x[-1])
        vec.append([label, left, right])
    return vec


def cbh_store2fragvec(cbh_store, keys=[None]):
    if type(cbh_store) != list:
        cbh_store = [cbh_store]
    if type(cbh_store[0]) != list:
        cbh_store = [cbh_store]
    vec_ls = list()
    cbh_store = cbh_store2vec(cbh_store)
    for cbh_s in cbh_store:
        vec = list()
        left, right = cbh_s[0], cbh_s[-1]
        #print('reac {}\nprod {}'.format(left,right))
        sum_l, sum_r = 0, 0
        for k in keys:
            if k in left:
                vec.append(-1 * left.count(k))
                sum_l += left.count(k)
            elif k in right:
                vec.append(right.count(k))
                sum_r += right.count(k)
            else:
                vec.append(0)
        #print('sum_r {} {}\nsum_l {} {}'.format(sum_r, len(right),sum_l,len(left)))
        if sum_r == len(right) and sum_l == len(left):
            vec_ls.append(vec)
        else:
            vec_ls.append([])
    return vec_ls


def smi2cbh(smi, rung, key_fn='fragment_lookup/keys.txt', graph_only=False):
    frag_keys = load_frags(key_fn)
    if type(smi) != list:
        smi = [smi]
    cbh_store = list()
    for s in smi:
        print('\nsmi2cbh({}):'.format(str(s)), end=' ')
        graph, cg_graph = smi2graph(s)
        mol_idx, mol = molfromsmi(s)
        mol, graph, cbh_store = molgraph2cbh(str(s), mol, graph, rung,
                                             frag_keys, cbh_store)
    #pprint_graph(graph)
    return cbh_store


def mol2cbh(mol, rung, key_fn='fragment_lookup/keys.txt', graph_only=False):
    frag_keys = load_frags(key_fn)
    if type(mol) != list:
        mol = [mol]
    cbh_store = list()
    for m in mol:
        print('\nmol2cbh({}):'.format(str(m)), end=' ')
        graph, cg_graph = molfn2graph(m)
        try:
            m = Chem.MolFromMolFile(m)
        except:
            pass
        m, graph, cbh_store = molgraph2cbh(str(m), m, graph, rung, frag_keys,
                                           cbh_store)
    return cbh_store


def xyz2cbh(fn, rung, key_fn='fragment_lookup/keys.txt', graph_only=False):
    frag_keys = load_frags(key_fn)
    if type(fn) != list:
        fn = [fn]
    cbh_store = list()
    for f in fn:
        print('\nxyz2cbh({}):'.format(str(f)), end=' ')
        graph, cg_graph = fn2graph(f)
        #pprint_graph(graph)
        mol_idx, mol = molfromxyz(f)
        mol, graph, cbh_store = molgraph2cbh(str(f), mol, graph, rung,
                                             frag_keys, cbh_store)
    return cbh_store


def files2cbh(fns,
              rung_ls,
              save_graph=False,
              key_fn='fragment_lookup/keys.txt',
              graph_only=False,
              fully_connected=False,
              coarse_grain=False):
    if type(fns) != list:
        fns = [fns]
    frag_keys = load_frags(key_fn)
    graphs_ls = list()
    cbh_store = list()
    for fn in fns:
        print('\nfiles2cbh({}):'.format(fn), end=' ')
        try:
            if '.smi' in fn:  #smiles:
                graph, cg_graph = smi2graph(fn)  #smifromsmifile(fn))
            elif '.mol' in fn:
                graph, cg_graph = molfn2graph(fn)
            else:
                graph, cg_graph = fn2graph(fn, fully_connected=fully_connected)
            #print('G:',graph)
            #sys.exit()
            #graphs_ls.append(graph)
            if '.smi' in fn:
                mol_idx, mol = molfromsmi(fn)
            elif '.mol' in fn:
                mol_idx = 0
                mol = Chem.MolFromMolFile(fn)
            else:
                mol_idx, mol = molfromxyz(fn)
            mol, graph, cbh_store = molgraph2cbh(fn, mol, graph, rung_ls,
                                                 frag_keys, cbh_store)
            graphs_ls.append(graph)

            if save_graph:
                #pprint_graph(graph)
                if graph_only:
                    cbh_store[-1] = [cbh_store[-1][0], graph]
                else:
                    cbh_store[-1].append(graph)
            #pprint_graph(graph)
        except:
            print(' FAILED pycbh/pycbh.py :172 ', end='')
            pass
    #print('\nLoaded {} graphs'.format(len(graphs_ls)))
    #for c in cbh_store:
    #  print('\n'+'\n'.join(c))

    with open(key_fn, 'w') as new_key_fn:
        i = 0
        for key in sorted(frag_keys):
            for frag in frag_keys[key]:
                new_key_fn.write('{} {}\n'.format(
                    key, " ".join([str(x) for x in frag])))
                i += 1
    #print('\n{} updated to {} fragments'.format(key_fn,i))
    return cbh_store


def molgraph2cbh(fn, mol, graph, rung_ls, frag_keys, cbh_store):
    if type(rung_ls) != list:
        try:
            if '+' in str(rung_ls):
                rung_ls = list(
                    range(1 +
                          int(str(rung_ls).replace('+', '').replace(' ', ''))))
            else:
                rung_ls = [rung_ls]
        except:
            rung_ls = [rung_ls]
    #Chem.Kekulize(mol, clearAromaticFlags=True)
    #mol = Chem.RemoveHs(mol)
    #print('{}'.format(Chem.MolToSmiles(Chem.RemoveHs(mol))))
    #print(rung_ls)
    for rung in rung_ls:
        #print('rung : CBH-{}'.format(rung))
        cbh_p_, cbh_r = calc_cbh(rung, graph)
        #f_dict = iep(cbh_p_)
        cbh_p__ = list()
        for frag in cbh_p_:
            f = list()
            for i in frag:
                f.extend(graph['nodes'][i][2])
            cbh_p__.append(f)
        f_dict = iep(cbh_p__)
        #print(f_dict)
        cbh_p, cbh_r = list(), list()
        for key in f_dict:
            if key % 2:
                cbh_p.extend(f_dict[key])
            else:
                cbh_r.extend(f_dict[key])
        p_test = Counter([str(x) for x in cbh_p])
        r_test = Counter([str(x) for x in cbh_r])
        #print('Products {} -> {}\nReactants {} -> {}'.format(len(cbh_p),len(p_test.keys()),len(cbh_r),len(r_test.keys())))
        #print('primary = {}'.format(cbh_p))
        #print('overlap = {}'.format(cbh_r))
        cbh_dict = dict()
        r_atoms = mol2formula(Chem.AddHs(mol), incl_H=True)
        p_atoms = dict()
        #print('  pdt start')
        for frag in p_test:
            f = str2list(frag)
            #f = list()
            #for i in str2list(frag):
            #  f.extend(graph['nodes'][i][2])
            #print('sending f ({}) to fraginc2smi'.format(f))
            try:
                smi, atoms, frag_fn, frag_keys = fraginc2smi(
                    f, mol, frag_keys, frag_type='primary', kekulize=True)
            except:
                #sys.exit()
                #pprint_graph(graph)
                print('    fraginc2smi(kek=False)')
                sys.exit()
                smi, atoms, frag_fn, frag_keys = fraginc2smi(
                    f, mol, frag_keys, frag_type='primary', kekulize=False)
            #sys.exit()
            atoms.update((x, y * p_test[frag]) for x, y in atoms.items())
            p_atoms = combine_dicts(p_atoms, atoms)
            '''
            atom_or_bond == 0 for atom centric
            atom_or_bond == 1 for bond centric
      '''
            if f in cbh_p_:
                atom_or_bond = rung % 2
                which_idx = [f[0]]
                if atom_or_bond:
                    try:
                        which_idx.append(f[1])
                    except:
                        pass
                #print('  info:',rung,which_idx,smi,frag_fn)
                graph = add_attr2graph(graph, rung, which_idx, smi, frag_fn)
            if '.' in smi:
                smi_ls = smi.split('.')
            else:
                smi_ls = [smi]
            for smi in smi_ls:
                if smi in cbh_dict:
                    cbh_dict[smi] += p_test[frag]
                else:
                    cbh_dict[smi] = p_test[frag]
        #print('  pdt done')
        #print('  rnt start')
        for frag in r_test:
            f = str2list(frag)
            #f = list()
            #for i in str2list(frag):
            #  f.extend(graph['nodes'][i][2])
            try:
                smi, atoms, frag_fn, frag_keys = fraginc2smi(
                    f, mol, frag_keys, frag_type='overlap', kekulize=True)
            except:
                #print('kekulize')
                smi, atoms, frag_fn, frag_keys = fraginc2smi(
                    f, mol, frag_keys, frag_type='overlap', kekulize=False)
            atoms.update((x, y * r_test[frag]) for x, y in atoms.items())
            r_atoms = combine_dicts(r_atoms, atoms)
            if '.' in smi:
                smi_ls = smi.split('.')
            else:
                smi_ls = [smi]
            for smi in smi_ls:
                if smi in cbh_dict:
                    cbh_dict[smi] -= r_test[frag]
                else:
                    cbh_dict[smi] = -1 * r_test[frag]
        #print('  rnt done')
        if '[HH]' in cbh_dict:
            cbh_dict['[H][H]'] = cbh_dict.pop('[HH]')
        if rung == 0:
            '''
      bond_deg = 0.0
      for bond in graph['edges']:
        if bond[0] != 1 and bond[1] != 1:
          bond_deg+=bond[2]
      cbh_dict['[H][H]'] = -1*int(0.5*bond_deg)
      '''
            if 'H' not in p_atoms:
                p_atoms['H'] = 0
            if 'H' not in r_atoms:
                r_atoms['H'] = 0
            if not p_atoms['H'] == r_atoms['H']:
                net_H = abs(p_atoms['H'] - r_atoms['H'])
                cbh_dict['[H][H]'] = -1 * int(0.5 * net_H)
                r_atoms['H'] += net_H
        #print('print_cbh')
        cbh = print_cbh(mol, cbh_dict, rung)

        if not sorted(p_atoms.items()) == sorted(r_atoms.items()):
            print("  WARNING: Atom counts dont match")
            net_atoms = dict()
            for key in p_atoms:
                if key in r_atoms:
                    net_atoms[key] = p_atoms[key] - r_atoms[key]
            #print('  Product:  {}\n  Reactant: {}'.format(sorted(p_atoms.items()),sorted(r_atoms.items())))
            print('  Atom counts: {}'.format(net_atoms))
            cbh_store.append([
                fn, cbh, "  WARNING: Atom counts dont match",
                "  Atom counts: {}".format(net_atoms)
            ])
            #sys.exit(cbh_store)
        else:
            cbh_store.append([fn, cbh])  #,"  (Atom counts are balanced)"])
            #print("\n  (Atom counts are balanced)")
        #print('  Atom counts: {}, {}'.format(p_atoms,r_atoms))
        #if len(graph['globals']) < 2:
        #  graph['globals'].append(Chem.MolToSmiles(Chem.RemoveHs(mol),kekuleSmiles=True,canonical=True))
    return mol, graph, cbh_store
