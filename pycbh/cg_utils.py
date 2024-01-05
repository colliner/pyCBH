import os, re, glob, sys
from rdkit import Chem

def molneighbors(idx, mol):
    """
    Returns a list of atom indices that are neighbors to the specified atom index in the molecule.

    Parameters:
    - idx (int): Index of the target atom.
    - mol (rdkit.Chem.Mol): RDKit molecule object.

    Returns:
    - list of int: List of atom indices that are neighbors to the specified atom index.
    """

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


def check_if_aromatic(mol, r):
    """
    Checks if all atoms in the specified ring are aromatic, and if all bonds within the ring are aromatic.

    Parameters:
    - mol (rdkit.Chem.Mol): RDKit molecule object.
    - r (list of int): List of atom indices representing a ring.

    Returns:
    - bool: True if all atoms and bonds in the ring are aromatic, False otherwise.
    """

    atom_cond = True
    for idx in r:
        if not mol.GetAtomWithIdx(idx).GetIsAromatic():
            atom_cond = False
            break
    if atom_cond:
        return True
    else:
        for bond in mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            if idx1 in r and idx2 in r:
                if bond.GetBondType() != 'AROMATIC':
                    return False
    return True


def check_if_oxy(i, mol, atoms, atype):
    """
    Checks if the specified atom is bonded to oxygen atoms based on its neighbors.

    Parameters:
    - i (int): Index of the target atom.
    - mol (rdkit.Chem.Mol): RDKit molecule object.
    - atoms (list): List of atom types.
    - atype (int): Type of the target atom.

    Returns:
    - bool: True if the conditions for oxygen bonding are met, False otherwise.
    """

    o_bonds = [0, 0]
    for x in molneighbors(i, mol):
        if int(atoms[x]) == 8:
            o_bonds[0] += 1
            if mol.GetBondBetweenAtoms(i, x).GetBondTypeAsDouble() > 1.0:
                o_bonds[1] += 1
    if atype == 6:
        return o_bonds[0] > 1 and o_bonds[1] > 0
    elif atype == 7:
        return o_bonds[0] > 0
    else:
        return o_bonds[0] > 1 or o_bonds[1] > 0


def coarse_grain(idx, atoms, mol, cg_opts, cut_bonds_btwn_groups=True):
    """
    Coarse grain module with various options for grouping atoms based on specific criteria.

    Parameters:
    - idx (int): Index of the target atom.
    - atoms (list): List of atom types.
    - mol (rdkit.Chem.Mol): RDKit molecule object.
    - cg_opts (list of str): List of coarse-graining options.
    - cut_bonds_btwn_groups (bool): Flag to cut bonds between groups if True.

    Returns:
    - list of list of int: List of lists, where each inner list represents a group of atom indices based on the coarse-graining options.
    """

    if len(cg_opts) > 0:
        incl = list()
    else:
        return [[idx]]
    ring_sizes, oxygen_groups = list(), list()
    for opt in cg_opts:
        if 'rings' in opt:
            size = opt.replace('rings', '')
            if len(size) > 0:
                try:
                    ring_sizes.append(int(size))
                except:
                    print('invalid ring size: {}'.format(size))
        elif opt == 'nitro':
            oxygen_groups.append(7)
        elif opt == 'sulfo':
            oxygen_groups.append(16)
        elif opt == 'phospho':
            oxygen_groups.append(15)
        elif opt == 'carbo':
            oxygen_groups.append(6)
    if 'aromatic' in cg_opts:
        arom_only = True
    else:
        arom_only = False
    if 'rings' in ''.join(cg_opts) or arom_only:
        if mol.GetAtomWithIdx(idx).IsInRing():
            ssr = Chem.GetSymmSSSR(mol)
            for r in ssr:
                if idx in r:
                    if len(ring_sizes) > 0:
                        if len(r) in ring_sizes:
                            incl.append([idx] + [x for x in r])
                            continue
                    elif arom_only:
                        if check_if_aromatic(mol, r):
                            incl.append([idx] + [x for x in r])
                            continue
                    else:
                        incl.append([idx] + [x for x in r])
    if len(incl) < 1:
        incl = [[idx]]
    if 'halogen' in cg_opts or 'halogens' in cg_opts:
        halogens = [9, 17, 35, 53, 85, 117]
        for idx_, inc in enumerate(incl):
            new_inc = inc.copy()
            for i in inc:
                if int(atoms[i]) in halogens:
                    for x in molneighbors(i, mol):
                        if x not in new_inc:
                            new_inc.append(x)
                else:
                    incl_ = molneighbors(i, mol)
                    for jdx in incl_:
                        if int(atoms[jdx]) in halogens:
                            if jdx not in new_inc:
                                new_inc.append(jdx)
                    incl[idx_] = new_inc
    if len(oxygen_groups) > 0:
        for idx_, inc in enumerate(incl):
            new_inc = inc.copy()
            for i in inc:
                if int(atoms[i]) in oxygen_groups:
                    if check_if_oxy(i, mol, atoms, int(atoms[i])):
                        for x in molneighbors(i, mol):
                            if cut_bonds_btwn_groups:
                                if int(atoms[x]) == 8:
                                    new_inc.append(x)
                            elif x not in new_inc:
                                new_inc.append(x)
                elif not cut_bonds_btwn_groups:
                    incl_ = molneighbors(i, mol)
                    for jdx in incl_:
                        if int(atoms[jdx]) in oxygen_groups:
                            if check_if_oxy(jdx, mol, atoms, int(atoms[i])):
                                for x in [jdx] + molneighbors(jdx, mol):
                                    if cut_bonds_btwn_groups:
                                        if int(atoms[x]) == 8:
                                            new_inc.append(x)
                                            if jdx not in new_inc:
                                                new_inc.append(jdx)
                                    elif x not in new_inc:
                                        new_inc.append(x)
                                        if jdx not in new_inc:
                                            new_inc.append(jdx)
            incl[idx_] = new_inc
    return incl


