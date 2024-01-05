#!/usr/bin/python

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
import subprocess
import pprint

from rdkit import RDLogger

from pycbh.cg_utils import *
#from pycbh.iep import iep


def vprint(*argu, end='\n'):
    #global verbose
    if verbose == 1:
        for arg in argu:
            print(arg, end=end)
        print()
    return


def str2list(st):
    return [
        int(x) for x in st.replace(' ', '').replace('[', '').replace(
            ']', '').split(',')
    ]


def molfromxyz(fn):
    atoms, charge, xyz_coordinates = xyz2mol.read_xyz_file(fn)
    try:
        mol_idx = int(fn.replace('xyz/qm7b_', '').replace('.xyz', '')) - 1
    except:
        mol_idx = int(0)
    mol = xyz2mol.xyz2mol(atoms, xyz_coordinates, charge,
                          embed_chiral=True)  #False)
    return mol_idx, mol


def atom2label(label):
    atoms = [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
        'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
        'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
        'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',
        'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
        'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At',
        'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
        'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
        'Ds', 'Rg', 'Cn', 'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo'
    ]
    if label in atoms:
        return atoms.index(label) + 1
    else:
        return label


def xyzfromsmi(smi):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    #print(smi, mol.GetNumAtoms())
    try:
        mol_ = Chem.RemoveHs(mol)
        with open('fragment_lookup/tmp.mol', "w") as FILE:
            FILE.write(Chem.MolToMolBlock(mol_))
        xyz_coordinates = list()
        #print(Chem.MolToMolBlock(mol))
        bashCommand = 'obabel -imol fragment_lookup/tmp.mol -oxyz --gen3d -xb'
        process = subprocess.Popen(bashCommand.split(),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        output, error = process.communicate()
        output = [x.split(' ') for x in output.decode("utf-8").split("\n")[2::]]
        for i, x_ in enumerate(output):
            #vprint(x_,len(x_))
            if len(x_) > 3:
                xyz_coordinates.append(
                    [float(x) for x in x_[1::] if len(x) > 0])
    except:
        xyz_coordinates = list()
        #print(Chem.MolToMolBlock(mol))
        with open('tmp', 'w') as FILE:
            FILE.write(smi)
        bashCommand = 'obabel -ismi tmp -oxyz --gen3d -xb '
        process = subprocess.Popen(bashCommand.split(),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        output, error = process.communicate()
        #print(output)
        output = [x.split(' ') for x in output.decode("utf-8").split("\n")[2::]]
        for i, x_ in enumerate(output):
            #print(x_,len(x_))
            if len(x_) > 3:
                xyz_coordinates.append(
                    [float(x) for x in x_[1::] if len(x) > 0])
    #print(coords)
    atoms = [atom2label(atom.GetSymbol()) for atom in mol.GetAtoms()]
    #print(atoms)
    return atoms, xyz_coordinates


def smifromsmifile(fn):
    return str(list(filter(None, open(fn, "r").read().split("\n")))[0])


def molfromsmi(fn):
    try:
        smi = smifromsmifile(fn)
    except:
        smi = fn
    mol = Chem.MolFromSmiles(smi)  #fromsmifile(fn))
    return 0, mol


def xyzfromfn(fn, heavy_only=False):
    atoms, charge, xyz_coordinates = xyz2mol.read_xyz_file(fn)
    heavy_atoms, heavy_xyz = list(), list()
    for idx, atom in enumerate(atoms):
        #print(idx,atom,xyz_coordinates[idx])
        if atom != 1 or not heavy_only:
            heavy_atoms.append(atom)
            heavy_xyz.append([round(float(x), 8) for x in xyz_coordinates[idx]])
    return heavy_atoms, heavy_xyz


def atom2onehot(atomtypes, atom):
    return [int(atom == x) for x in atomtypes]


def gaussian_expansion(d,
                       n_centers=20,
                       sigma=0.5,
                       min_d=0.0,
                       max_d=4.0,
                       centers=[None]):
    if None in centers:
        centers = np.linspace(min_d, max_d, n_centers)
        #print('{} centers between {} and {} :\n  {}'.format(n_centers,min_d,max_d,centers))
    return [round(math.exp(-(d - x)**2 / sigma**2), 8) for x in centers
           ], centers


def transform_edges(edge,
                    verbose=False,
                    use_GauExp=True,
                    use_bondinverse=False,
                    use_z=True,
                    use_onehotatom=True):
    '''
  input: edge (list) : [ bond_len, atom_num1, atom_num2 ]
  '''
    atomtypes = [6, 7, 8, 9, 16, 17]
    #print(atomtypes)
    new_edge = list()
    if use_z:
        new_edge.append(edge[1])
        new_edge.append(edge[2])
    new_edge.append(edge[0])
    if use_GauExp:
        new_edge.extend(gaussian_expansion(edge[0])[0])
    elif use_bondinverse:
        new_edge.append(1 / edge[0])
    if use_onehotatom:
        new_edge.extend(atom2onehot(atomtypes, edge[1]))
        new_edge.extend(atom2onehot(atomtypes, edge[2]))
    #if verbose:
    #print(atomtypes)
    #print('\nold edge: {}'.format(edge))
    #print('\nnew edge: {}'.format(new_edge))
    return [float(x) for x in new_edge]


def get_dist(x1, x2):
    return round(
        ((x1[0] - x2[0])**2 + (x1[1] - x2[1])**2 + (x1[2] - x2[2])**2)**0.5, 8)


def pprint_graph(graph):
    print("graph = {")
    print("  'nodes':", end=' ')
    for idx, a in enumerate(graph['nodes']):
        if idx == 0:
            print(a)
        else:
            print('           {}'.format(a))
    print("  'edges':", end=' ')
    for idx, e in enumerate(graph['edges']):
        if idx == 0:
            print(e)
        else:
            print('           {}'.format(e))
    print("  'senders':   {}".format(graph['senders']))
    print("  'receivers': {}".format(graph['receivers']))
    print("  'globals': {}".format(graph['globals']))
    print("}")
    print('graph[nodes] : {}'.format(np.asarray(graph['nodes']).shape))
    print('graph[edges] : {}'.format(np.asarray(graph['edges']).shape))
    print('graph[senders] : {}'.format(np.asarray(graph['senders']).shape))
    print('graph[receivers] : {}'.format(np.asarray(graph['receivers']).shape))
    print('graph[globals] : {}'.format(np.asarray(graph['globals']).shape))
    return


def molneighbors_old(idx, mol):
    incl = list()
    for bond in mol.GetBonds():
        idx1 = bond.GetBeginAtomIdx()
        idx2 = bond.GetEndAtomIdx()
        if idx in [idx1, idx2]:
            if idx1 not in incl:
                incl.append(idx1)
            if idx2 not in incl:
                incl.append(idx2)
    return incl


def coarse_grain_old(idx, atoms, mol, cg_opts):
    if len(cg_opts) > 0:
        print('Coarse graining')
        incl = list()
    else:
        incl = [[idx]]
    if 'rings' in cg_opts:
        print('  rings')
        if mol.GetAtomWithIdx(idx).IsInRing():
            ssr = Chem.GetSymmSSSR(mol)
            for r in ssr:
                if idx in r:
                    incl.append([idx] + [x for x in r])
    if len(incl) < 1:
        incl = [[idx]]
    #incl = [x for i, x in enumerate(incl) if incl[i::].count(x) < 2]
    if 'halogen' in cg_opts or 'halogens' in cg_opts:
        print('  halogens')
        halogens = [9, 17, 35, 53, 85, 117]
        #if int(atoms[idx]) in halogens:
        for idx_, inc in enumerate(incl):
            for i in inc:
                #incl=molneighbors(idx,mol)
                if int(atoms[i]) in halogens:
                    for x in molneighbors(i, mol):
                        if x not in inc:
                            inc.append(x)
                            #incl[jdx].append(x)
                    #incl.extend(molneighbors(i,mol))
                else:
                    incl_ = molneighbors(i, mol)
                    #incl_=molneighbors(idx,mol)
                    for jdx in incl_:
                        if int(atoms[jdx]) in halogens:
                            if jdx not in incl:
                                #incl.append([idx]+jdx)
                                inc.append(jdx)
            incl[idx_] = inc
    if len(cg_opts) > 0:
        print('  incl:{}'.format(incl))
    return incl


def float_len(x):
    return float(len(x))


def smi2graph(smi, simple_graph=False):
    try:
        smi = smifromsmifile(smi)
    except:
        pass
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    graph, cg_graph = mol2graph(smi, 0, mol, simple_graph=simple_graph)
    return graph, cg_graph


def molfn2graph(mol, simple_graph=False):
    try:
        mol = Chem.MolFromMolFile(mol)
    except:
        pass
    mol = Chem.AddHs(mol)
    graph, cg_graph = mol2graph('', 0, mol, simple_graph=simple_graph)
    return graph, cg_graph


def fn2graph(fn, simple_graph=False, fully_connected=False):
    mol_idx, mol = molfromxyz(fn)
    #print(mol)
    graph, cg_graph = mol2graph(fn,
                                mol_idx,
                                mol,
                                simple_graph=simple_graph,
                                fully_connected=fully_connected)
    return graph, cg_graph


def mol2graph(fn,
              mol_idx,
              mol,
              simple_graph=False,
              cg_opts=[],
              fully_connected=False):
    smi = Chem.MolToSmiles(mol, kekuleSmiles=True, canonical=True)
    mol_ = Chem.RemoveHs(mol)
    smi = Chem.MolToSmiles(mol_)
    if '.xyz' in fn:
        atoms, charge, xyz_coordinates = xyz2mol.read_xyz_file(fn)
    else:
        with open('fragment_lookup/tmp.mol', "w") as FILE:
            FILE.write(Chem.MolToMolBlock(mol_))
        xyz_coordinates = list()
        #print(Chem.MolToMolBlock(mol))
        bashCommand = 'obabel -imol fragment_lookup/tmp.mol -oxyz --gen3d -xb'
        process = subprocess.Popen(bashCommand.split(),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        output, error = process.communicate()
        output = [x.split(' ') for x in output.decode("utf-8").split("\n")[2::]]
        for i, x_ in enumerate(output):
            #vprint(x_,len(x_))
            if len(x_) > 3:
                xyz_coordinates.append(
                    [float(x) for x in x_[1::] if len(x) > 0])
        #vprint(xyz_coordinates)
        atoms = [atom2label(atom.GetSymbol()) for atom in mol.GetAtoms()]
    #print(len(atoms),atoms)
    #print(len(xyz_coordinates),xyz_coordinates)
    #sys.exit()
    '''
  try:
    atoms, xyz_coordinates = xyzfromfn(fn)
  except:
    try:
      atoms, xyz_coordinates = xyzfromsmi(fn)
    except:
      try:
        atoms, xyz_coordinates = xyzfromsmi(smi)
      except:
        sys.exit('failed')
  '''
    #atomtypes = [6, 7, 8, 9, 16, 17]
    graph = {
        'nodes': list(),
        'edges': list(),
        'senders': list(),
        'receivers': list(),
        'globals': list()
    }
    cg_graph = {
        'nodes': list(),
        'edges': list(),
        'senders': list(),
        'receivers': list(),
        'globals': list()
    }
    #for atom in mol.GetAtoms():
    #  print(atom.GetExplicitValence())
    #print('here')
    include_H = True
    cg_incld = list()
    #mol = Chem.MolFromMolFile('fragment_lookup/tmp.mol')
    use_coords = False
    oxy_ls = [idx for idx, x in enumerate(atoms) if x == 8]
    #print('oxy list: {}'.format(oxy_ls))
    for idx, atom in enumerate(atoms):
        if idx in cg_incld + oxy_ls:
            #print('{} already incl'.format(idx))
            pass
        elif atom != 1.0 or include_H:
            cg_opts = []
            cg_incl_ls = coarse_grain(idx, atoms, mol, cg_opts=cg_opts)
            #print('idx: {}'.format(idx))
            atom_obj = mol.GetAtoms()[idx]
            #print('  ',idx, float(atom), atom_obj.GetExplicitValence(), xyz_coordinates[idx])
            for cg_incl in cg_incl_ls:
                if use_coords:
                    try:
                        graph['nodes'].append([
                            idx,
                            float(atom), cg_incl,
                            atom_obj.GetExplicitValence(),
                            [
                                round(x, 4)
                                for x in xyz_coordinates[idx].tolist()
                            ]
                        ])
                    except:
                        #print('{} failed'.format(idx))
                        graph['nodes'].append([
                            idx,
                            float(atom), cg_incl,
                            atom_obj.GetExplicitValence(),
                            [round(x, 4) for x in xyz_coordinates[idx]]
                        ])
                else:
                    graph['nodes'].append([idx, float(atom), cg_incl])
                if idx not in cg_incld:
                    cg_graph['nodes'].append([idx, float(atom), cg_incl])
                    #print('added {}'.format([idx, float(atom), cg_incl]))
            if idx not in cg_incld:
                for cg_incl in cg_incl_ls:
                    cg_incld.extend(cg_incl)
            #print(cg_incld)
    atom = 8.0
    for idx in oxy_ls:
        if idx not in cg_incld:
            cg_incl_ls = coarse_grain(idx, atoms, mol, cg_opts=cg_opts)
            atom_obj = mol.GetAtoms()[idx]
            for cg_incl in cg_incl_ls:
                if use_coords:
                    try:
                        graph['nodes'].append([
                            idx,
                            float(atom), cg_incl,
                            atom_obj.GetExplicitValence(),
                            [
                                round(x, 4)
                                for x in xyz_coordinates[idx].tolist()
                            ]
                        ])
                    except:
                        graph['nodes'].append([
                            idx,
                            float(atom), cg_incl,
                            atom_obj.GetExplicitValence(),
                            [round(x, 4) for x in xyz_coordinates[idx]]
                        ])
                else:
                    graph['nodes'].append([idx, float(atom), cg_incl])
                if idx not in cg_incld:
                    cg_graph['nodes'].append([idx, float(atom), cg_incl])
            if idx not in cg_incld:
                for cg_incl in cg_incl_ls:
                    cg_incld.extend(cg_incl)
    graph['nodes'] = sorted(graph['nodes'])
    #pprint_graph(graph)
    #pprint.pprint(invar)
    #print('before kek')
    Chem.Kekulize(mol, clearAromaticFlags=True)
    #print('here')
    for bond in mol.GetBonds():
        idx1 = bond.GetBeginAtomIdx()
        idx2 = bond.GetEndAtomIdx()
        inc = True
        if atoms[idx1] == 1.0 or atoms[idx2] == 1.0:
            if include_H:
                inc = True
            else:
                inc = False
        if inc:  #atoms[idx1] != 1.0 and atoms[idx2] != 1.0:
            if use_coords:
                dist = get_dist(xyz_coordinates[idx1], xyz_coordinates[idx2])
                coulomb_inter = round(
                    float(atoms[idx1]) * float(atoms[idx2]) / float(dist), 8)
                edge1 = [
                    atoms[idx1], atoms[idx2],
                    bond.GetBondTypeAsDouble(), dist, coulomb_inter
                ]
                edge2 = [
                    atoms[idx2], atoms[idx1],
                    bond.GetBondTypeAsDouble(), dist, coulomb_inter
                ]
            else:
                edge1 = [atoms[idx1], atoms[idx2], bond.GetBondTypeAsDouble()]
                edge2 = [atoms[idx2], atoms[idx1], bond.GetBondTypeAsDouble()]
            graph['edges'].append(edge1)
            graph['senders'].append(idx1)
            graph['receivers'].append(idx2)
            graph['edges'].append(edge2)
            graph['senders'].append(idx2)
            graph['receivers'].append(idx1)
            for idx, node in enumerate(cg_graph['nodes']):
                if idx1 in node[2]:
                    if idx2 not in node[2]:
                        for jdx, jnode in enumerate(cg_graph['nodes']):
                            create_edge = True
                            if idx2 in jnode[2]:
                                for kdx, sender in enumerate(
                                        cg_graph['senders']):
                                    if idx == sender:
                                        if jdx == cg_graph['receivers'][kdx]:
                                            create_edge = False
                                            break
                                if create_edge:
                                    cg_graph['edges'].append(edge1)
                                    cg_graph['senders'].append(idx)
                                    cg_graph['receivers'].append(jdx)
                                    cg_graph['edges'].append(edge2)
                                    cg_graph['senders'].append(jdx)
                                    cg_graph['receivers'].append(idx)
    if fully_connected:
        for idx1 in range(len(cg_graph['nodes']) - 1):
            for idx2 in range(idx1 + 1, len(cg_graph['nodes'])):
                if not is_connected(idx1, idx2, cg_graph):
                    dist = get_dist(xyz_coordinates[idx1],
                                    xyz_coordinates[idx2])
                    coulomb_inter = round(
                        float(atoms[idx1]) * float(atoms[idx2]) / float(dist),
                        8)
                    edge1 = [atoms[idx1], atoms[idx2], 0.0, dist, coulomb_inter]
                    edge2 = [atoms[idx2], atoms[idx1], 0.0, dist, coulomb_inter]
                    cg_graph['edges'].append(edge1)
                    cg_graph['senders'].append(idx1)
                    cg_graph['receivers'].append(idx2)
                    cg_graph['edges'].append(edge2)
                    cg_graph['senders'].append(idx2)
                    cg_graph['receivers'].append(idx1)

    graph['globals'] = [fn, smi]
    cg_graph['globals'] = [fn, smi]
    #pprint_graph(graph)
    #pprint_graph(cg_graph)
    #sys.exit()
    return cg_graph, graph  #, cg_graph


def is_connected(idx1, idx2, graph):
    for i, idx in enumerate(graph['senders']):
        if idx == idx1:
            if graph['receivers'][i] == idx2:
                return True
    return False


def check_dup(f, f_ls):
    #return f_ls.count(f)
    f = sorted(f)
    f_ls = [sorted(x) for x in f_ls]
    return f_ls.count(f)


def overlap(fragments, rung, graph):
    '''
  atom_or_bond == 0 for atom centric
  atom_or_bond == 1 for bond centric
  '''
    atom_or_bond = rung % 2
    overlaps = list()
    for idx1, frag1 in enumerate(fragments):
        for idx2, frag2 in enumerate(fragments):
            calc_overlap = False
            if atom_or_bond:
                if len([value for value in frag1[0:2] if value in frag2[0:2]
                       ]) == 1:
                    calc_overlap = True
            else:
                #print('is_connected(frag1[0],frag2[0])=is_connected({},{})'.format(frag1,frag2))
                if is_connected(frag1[0], frag2[0], graph):
                    calc_overlap = True
            if calc_overlap and idx1 < idx2:
                new_frag = [value for value in frag1 if value in frag2]
                #print('overlap between {} and {} = {}'.format(frag1,frag2,new_frag))
                if len(new_frag) > 0:
                    overlaps.append(new_frag)
    if rung == 0:
        bond_deg = 0.0
        for bond in graph['edges']:
            if bond[0] != 1 and bond[1] != 1:
                bond_deg += bond[2]

    #print('overlaps = {}'.format(overlaps))
    return overlaps


def collect_neighbors(idx, rung, graph):
    '''
  atom_or_bond == 0 for atom centric
  atom_or_bond == 1 for bond centric
  '''
    atom_or_bond = rung % 2
    incl_ls = list()
    if atom_or_bond:
        l = 1
        found = False
        for i, sender in enumerate(graph['senders']):
            if sender == idx:
                idx2 = graph['receivers'][i]
                connected_node = graph['nodes'][idx2]
                if connected_node[1] != 1.0:
                    if idx2 > idx and graph['edges'][i][2] > 0.0:
                        #incl=[idx,idx2]+graph['nodes'][idx][2]+graph['nodes'][idx2][2]
                        #incl_ls.append(incl)
                        incl_ls.append([idx, idx2])
                        found = True
        if not found:
            incl_ls.append([idx])
    else:
        l = 0
        #incl=[idx]+graph['nodes'][idx][2]
        #incl_ls.append(incl)
        incl_ls.append([idx])
    while l < rung:
        for j, incl in enumerate(incl_ls):
            new_incl = list()
            for i, sender in enumerate(graph['senders']):
                if sender in incl:
                    idx2 = graph['receivers'][i]
                    connected_node = graph['nodes'][idx2]
                    if connected_node[1] != 1.0:
                        if idx2 not in incl and graph['edges'][i][2] > 0.0:
                            new_incl.append(idx2)
                            #new_incl.extend(graph['nodes'][idx2][2])
            incl_ls[j].extend(new_incl)
        l += 2
    return incl_ls


def fix_cbh(f, o):
    done = False
    while not done:
        done = True
        which = list()
        for idx, x in enumerate(o):
            if check_dup(x, o) > 2:
                which.append(idx)
                done = False
                break
        for i in which:
            del o[i]
        which = list()
        for idx, x in enumerate(f):
            if check_dup(x, f) > 2:
                which.append(idx)
                done = False
                break
        for i in which:
            del f[i]
    return f, o


def calc_cbh(rung, graph):
    fragments = list()
    for idx, node in enumerate(graph['nodes']):
        if node[1] != 1.0:
            fragments.extend(collect_neighbors(idx, rung, graph))
    #print('primary fragments:',fragments)
    overlaps = overlap(fragments, rung, graph)
    if rung > 0:
        return fix_cbh(fragments, overlaps)
    else:
        bond_deg = 0.0
        for bond in graph['edges']:
            if bond[0] != 1 and bond[1] != 1:
                bond_deg += bond[2]
        return fragments, overlaps


def mol2formula(mol, incl_H=False):
    atom_types = dict()
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol != 'H' or incl_H:
            if symbol in atom_types:
                atom_types[symbol] += 1
            else:
                atom_types[symbol] = 1
    return atom_types


def moldict2hill(atom_types):
    #print('atom types: {}'.format(atom_types))
    hill = list()
    if 'C' in atom_types:
        hill.append('C')
        num = atom_types['C']
        if num > 1:
            hill.append('{}'.format(atom_types['C']))
    if 'H' in atom_types:
        hill.append('H')
        num = atom_types['H']
        if num > 1:
            hill.append('{}'.format(atom_types['H']))
    for key in sorted(atom_types):
        if key not in ['C', 'H']:
            hill.append(key)
            num = atom_types[key]
            if num > 1:
                hill.append('{}'.format(atom_types[key]))
    return ''.join(hill)


def combine_dicts(d1, d2):
    for key in d2:
        if key in d1:
            d1[key] += d2[key]
        else:
            d1[key] = d2[key]
    return d1


def chunks(data, size=5000):
    for i in range(0, len(data), size):
        #yield {k:data[k] for k in islice(it, size)}
        yield data[i:i + size]


def formula2ha_old(formula):
    ha_count = 0
    spl = list(formula)
    for i in range(len(spl)):
        try:
            c = int(spl[i]) - 1
            if spl[i - 1] != 'H':
                ha_count += c
        except:
            if spl[i] != 'H':
                ha_count += 1
            continue
    return ha_count


def formula2ha(f):
    ha_count = 0
    spl = list(f)
    k, v = '', '0'
    num_mode = True
    for x in spl:
        if x.isalpha():
            if num_mode:
                if k != 'H':
                    ha_count += int(v)
                k, v = '', ''
                num_mode = False
            elif x.isupper():
                ha_count += 1
                k, v = '', ''
            k += x
        else:
            try:
                int(x)
                v += str(x)
                num_mode = True
            except:
                pass
    if not num_mode:
        ha_count += 1
    elif k != 'H':
        ha_count += int(v)
    return ha_count


def get_frag_fn(formula, smi, smi2, frag_keys):
    '''
  frag_keys structure:
  dict() = {formula : [idx, smi, smi2, fn]}
  example:
      'C2H6' : [[0, 'CC',  '[H][C]([H])([H])[C]([H])([H])[H]H', 'f2_C2H6_000']]
  '''
    if formula in frag_keys:
        for f in frag_keys[formula]:
            if smi == f[1]:
                #with open('coeff.txt','a') as FILE:
                #  FILE.write(' '.join([str(x) for x in f])+'\n')
                return f[-1], False, frag_keys
        idx = len(frag_keys[formula])
        fn = 'f{}_{}_{:03d}'.format(formula2ha(formula), formula, idx)
        frag_keys[formula].append([idx, smi, smi2, fn])
        #with open('coeff.txt','a') as FILE:
        #  FILE.write(' '.join([str(x) for x in f])+'\n')
    else:
        fn = 'f{}_{}_{:03d}'.format(formula2ha(formula), formula, 0)
        frag_keys[formula] = [[0, smi, smi2, fn]]
        #with open('coeff.txt','a') as FILE:
        #  FILE.write(' '.join([str(x) for x in [0, smi, smi2, fn]])+'\n')
    return fn, True, frag_keys


def fraginc2smi(f, mol, frag_keys, frag_type=None, kekulize=False):
    RDLogger.DisableLog('rdApp.*')
    smi = Chem.MolToSmiles(mol)
    #print('{:02d}'.format(f[0]), end=' ')
    if kekulize:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    mw = Chem.RWMol(mol)
    numatoms = mw.GetNumAtoms()
    total_deg = [atom.GetTotalValence() for atom in mw.GetAtoms()]
    for i in range(numatoms):
        idx = numatoms - 1 - i
        if idx not in f:
            mw.RemoveAtom(idx)
    numatoms = mw.GetNumAtoms()
    #if len(Chem.GetSymmSSSR(mw)) < 1:
    mw = Chem.RWMol(Chem.AddHs(mw))
    #print(total_deg)
    #print('a : {}'.format([atom.GetAtomicNum() for atom in mol.GetAtoms()]))
    #print('f : {}'.format(f))
    for idx, val in enumerate(total_deg):
        if idx in f:
            idx2 = sorted(list(set(f))).index(idx)
            atom = mw.GetAtomWithIdx(idx2)
            if atom.GetAtomicNum() != 1:
                if atom.GetTotalValence() != val:
                    #print('{}({})'.format(idx,atom.GetAtomicNum()))
                    #print('VALENCE DOES NOT MATCH {} -> {}'.format(atom.GetTotalValence(),val))
                    #print('numatoms : {}'.format(mw.GetNumAtoms()))
                    for _ in range(val - atom.GetTotalValence()):
                        idx_h = mw.AddAtom(Chem.Atom(1))
                        #print('added H {}'.format(idx_h))
                        mw.AddBond(idx2, idx_h, Chem.BondType.SINGLE)
                    #print('numatoms : {}'.format(mw.GetNumAtoms()))
                    #sys.exit()
    idx_rings = list()
    for r in Chem.GetSymmSSSR(mw):
        for x in r:
            if x not in idx_rings:
                idx_rings.append(x)
    #print(idx_rings)
    if True:
        for idx, atom in enumerate(mw.GetAtoms()):
            if idx not in idx_rings:
                atom.SetIsAromatic(False)
    if len(Chem.GetSymmSSSR(mw)) < 1:
        try:
            Chem.Kekulize(mw, clearAromaticFlags=True)
            smi = Chem.MolToSmiles(mw, kekuleSmiles=True, canonical=True)
        except:
            print('Cannot kekulize mw')
            smi = Chem.MolToSmiles(mw)

    #smi = Chem.MolToSmiles(mw,kekuleSmiles=True,canonical=True)
    else:
        smi = Chem.MolToSmiles(mw)
    mol = Chem.MolFromSmiles(smi)
    if mol == None:
        if 'n' in smi:
            smi = smi.replace('n', '[nH]')
        elif ':O:' in smi:
            smi = smi.replace(':O:', '[O]')
        mol = Chem.MolFromSmiles(smi)
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
        smi = Chem.MolToSmiles(mol, kekuleSmiles=True)
    except:
        pass
    #smi = Chem.MolToSmiles(mol,kekuleSmiles=True)
    #print(smi, mol)
    #mw = Chem.AddHs(mw)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    #print(smi, mol.GetNumAtoms())
    formula = moldict2hill(mol2formula(mol, incl_H=True))
    smi2 = Chem.MolToSmiles(mol, allHsExplicit=True, allBondsExplicit=False)
    mol.SetProp("_Name", "{}  {}  {}".format(formula, smi, smi2))
    if '.' in smi:
        smi_ls = smi.split('.')
        for s in smi_ls:
            mol_s = Chem.MolFromSmiles(s)
            mol_s = Chem.AddHs(mol_s)
            AllChem.EmbedMolecule(mol_s)
            s2 = Chem.MolToSmiles(mol_s,
                                  allHsExplicit=True,
                                  allBondsExplicit=False)
            formula = moldict2hill(mol2formula(mol_s, incl_H=True))
            mol.SetProp("_Name", "{}  {}  {}".format(formula, s, s2))
            frag_fn, make_mol, frag_keys = get_frag_fn(formula, s, s2,
                                                       frag_keys)
            if make_mol:
                with open('fragment_lookup/' + frag_fn + '.mol', "w") as fn:
                    fn.write(Chem.MolToMolBlock(mol_s))
                #print('Written to fragment_lookup/{}.mol'.format(frag_fn))
    else:
        frag_fn, make_mol, frag_keys = get_frag_fn(formula, smi, smi2,
                                                   frag_keys)
        #print(Chem.MolToMolBlock(mol))
        if make_mol:
            with open('fragment_lookup/' + frag_fn + '.mol', "w") as fn:
                fn.write(Chem.MolToMolBlock(mol))
            #print('Written to fragment_lookup/{}.mol'.format(frag_fn))
    return smi, mol2formula(mol, incl_H=True), frag_fn, frag_keys


def cbh_print(cbh):
    print()
    for c in cbh:
        if type(c[-1]) == dict:
            c = c[1:-1]
        else:
            c = c[1:]
        print('\n' + '\n'.join([' ' + x for x in c]))
    pass


def print_cbh(mol, cbh_dict, rung):
    react_ls, prod_ls = list(), list()
    #print('\nCBH-{}:'.format(rung))
    mol = Chem.RWMol(mol)
    mol = Chem.RemoveHs(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    #print('  '+Chem.MolToSmiles(mol,kekuleSmiles=True,canonical=True),end=' ')
    for key in cbh_dict:
        if cbh_dict[key] < 0:
            if cbh_dict[key] == -1:
                react_ls.append(' {} '.format(key))
            else:
                react_ls.append(' {} {} '.format(-1 * cbh_dict[key], key))
        elif cbh_dict[key] > 0:
            if cbh_dict[key] == 1:
                prod_ls.append(' {} '.format(key))
            else:
                prod_ls.append(' {} {} '.format(cbh_dict[key], key))
    cbh_str = 'CBH-{}:\n'.format(rung) + '  ' + Chem.MolToSmiles(
        mol, kekuleSmiles=True, canonical=True) + ' '
    if len(react_ls) > 0:
        #print('+',end='')
        cbh_str = cbh_str + '+'
    #print('+'.join(react_ls),end=' --> ')
    #print('+'.join(prod_ls))
    #print(cbh_dict)
    return cbh_str + '+'.join(react_ls) + ' --> ' + '+'.join(prod_ls)


def add_attr2graph(graph, rung, which_idx, smi, frag_fn):
    if len(which_idx) == 1:
        if rung % 2 == 0:
            graph['nodes'][which_idx[0]].append(
                ['CBH-{}'.format(rung), smi, frag_fn])
    else:
        for idx, idx1 in enumerate(graph['senders']):
            idx2 = graph['receivers'][idx]
            if idx1 == which_idx[0] and idx2 == which_idx[1]:
                graph['edges'][idx].append(
                    ['CBH-{}'.format(rung), smi, frag_fn])
            elif idx2 == which_idx[0] and idx1 == which_idx[1]:
                graph['edges'][idx].append(
                    ['CBH-{}'.format(rung), smi, frag_fn])
    return graph


if __name__ == '__main__':
    if len(sys.argv[1::]) < 2:
        sys.exit('Usage python3 fp_gen.py CBH_rung FILENAME(s)')

    fns = sys.argv[2::]
    try:
        rung_ls = [int(sys.argv[1])]
    except:
        if sys.argv[1] == 'all':
            rung_ls = [0, 1, 2, 3, 4]
        else:
            sys.exit(
                'Usage: python3 fp_gen.py CBH_rung FILENAME(s)\n  CBH rung must be an integer'
            )
    key_fn = 'fragment_lookup/key.txt'
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
    verbose = False
    graph = True
    save_graph = False
    #rung=2
    #smiles=True
    if graph:
        graphs_ls = list()
        cbh_store = list()
        for fn in fns:
            print('{} :'.format(fn), end=' ')
            try:
                if '.smi' in fn:  #smiles:
                    graph = smi2graph(smifromsmifile(fn))
                else:
                    graph = fn2graph(fn)
                #print('G:',graph)
                #graphs_ls.append(graph)
                if '.smi' in fn:
                    mol_idx, mol = molfromsmi(fn)
                else:
                    mol_idx, mol = molfromxyz(fn)
                Chem.Kekulize(mol, clearAromaticFlags=True)
                print('{}'.format(
                    Chem.MolToSmiles(Chem.RemoveHs(mol),
                                     kekuleSmiles=True,
                                     canonical=True)))
                for rung in rung_ls:
                    cbh_p_, cbh_r = calc_cbh(rung, graph)
                    f_dict = iep(cbh_p_)
                    print(f_dict)
                    cbh_p, cbh_r = list(), list()
                    for key in f_dict:
                        if key % 2:
                            cbh_p.extend(f_dict[key])
                        else:
                            cbh_r.extend(f_dict[key])
                    print('primary = {}'.format(cbh_p))
                    print('overlap = {}'.format(cbh_r))
                    '''
          if '.smi' in fn:
            mol_idx, mol = molfromsmi(fn)
          else:
            mol_idx, mol = molfromxyz(fn)
          Chem.Kekulize(mol, clearAromaticFlags=True)
          '''
                    cbh_dict = dict()
                    r_atoms = mol2formula(Chem.AddHs(mol), incl_H=True)
                    p_atoms = dict()
                    for f in cbh_p:
                        #print('sending f ({}) to fraginc2smi'.format(f))
                        smi, atoms, frag_fn, frag_keys = fraginc2smi(
                            f, mol, frag_keys, frag_type='primary')
                        p_atoms = combine_dicts(p_atoms, atoms)
                        '''
            atom_or_bond == 0 for atom centric
            atom_or_bond == 1 for bond centric
            '''
                        if f in cbh_p_:
                            atom_or_bond = rung % 2
                            which_idx = [f[0]]
                            if atom_or_bond:
                                which_idx.append(f[1])
                            print('  info:', rung, which_idx, smi, frag_fn)
                            graph = add_attr2graph(graph, rung, which_idx, smi,
                                                   frag_fn)
                        if '.' in smi:
                            smi_ls = smi.split('.')
                        else:
                            smi_ls = [smi]
                        for smi in smi_ls:
                            if smi in cbh_dict:
                                cbh_dict[smi] += 1
                            else:
                                cbh_dict[smi] = 1
                    for f in cbh_r:
                        smi, atoms, frag_fn, frag_keys = fraginc2smi(
                            f, mol, frag_keys, frag_type='overlap')
                        r_atoms = combine_dicts(r_atoms, atoms)
                        if [sorted(x) for x in cbh_r].count(sorted(f)) > 1:
                            print('\n{} matches something else'.format(f))
                            '''
              if smi in cbh_dict:
                cbh_dict[smi]+=0.5
              else:
                cbh_dict[smi]=0.5
            
              atoms.update((x, y*-0.5) for x, y in atoms.items())
              print(atoms)
              r_atoms = combine_dicts(r_atoms, atoms)
              '''
                        if '.' in smi:
                            smi_ls = smi.split('.')
                        else:
                            smi_ls = [smi]
                        for smi in smi_ls:
                            if smi in cbh_dict:
                                cbh_dict[smi] -= 1
                            else:
                                cbh_dict[smi] = -1
                    if rung == 0:
                        '''
            bond_deg = 0.0
            for bond in graph['edges']:
              if bond[0] != 1 and bond[1] != 1:
                bond_deg+=bond[2]
            cbh_dict['[H][H]'] = -1*int(0.5*bond_deg)
            '''
                        if not p_atoms['H'] == r_atoms['H']:
                            net_H = abs(p_atoms['H'] - r_atoms['H'])
                            cbh_dict['[H][H]'] = -1 * int(0.5 * net_H)
                            r_atoms['H'] += net_H
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
                        cbh_store.append([fn, cbh
                                         ])  #,"  (Atom counts are balanced)"])
                        print("\n  (Atom counts are balanced)")
                    print('  Atom counts: {}'.format(net_atoms))
                graphs_ls.append(graph)
                pprint_graph(graph)
            except:
                print(' FAILED pycbh/utils.py : 976')
                pass
        print('\nLoaded {} graphs'.format(len(graphs_ls)))
        for c in cbh_store:
            print('\n' + '\n'.join(c))

        with open(key_fn, 'w') as new_key_fn:
            i = 0
            for key in sorted(frag_keys):
                for frag in frag_keys[key]:
                    new_key_fn.write('{} {}\n'.format(
                        key, " ".join([str(x) for x in frag])))
                    i += 1
        print('\n{} updated to {} fragments'.format(key_fn, i))

        if save_graph:
            graph_chunks = list()
            for item in chunks(graphs_ls, size=len(graphs_ls) // 25):
                graph_chunks.append(item)
            print('  {} chunks made'.format(len(graph_chunks)))
            which_fp = picklefn.replace('.pickle', '')
            if not os.path.exists("qm7b_graphs/" + which_fp):
                os.makedirs("qm7b_graphs/" + which_fp)
                print('Made dir: {}'.format("qm7b_graphs/" + which_fp))

            fn_base = "qm7b_graphs/" + which_fp + "/qm7b_graphs_" + which_fp
            print('\nWriting to files {}_###.npy'.format(fn_base))

            for idx, item in enumerate(graph_chunks):
                fn = fn_base + '_{:03d}'.format(idx) + '.npy'
                np.save(fn, item)
                print('  written to {}'.format(fn))

    else:
        for fn in fns:
            fp = fn2fp(fn, V)
            print('{} : {}'.format(fn, ' '.join([str(x) for x in fp.tolist()])))
