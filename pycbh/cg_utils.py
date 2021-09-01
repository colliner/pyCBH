import os,re,glob,sys
from rdkit import Chem


def molneighbors(idx,mol):
  incl=list()
  for bond in mol.GetBonds():
    idx1=bond.GetBeginAtomIdx()
    idx2=bond.GetEndAtomIdx()
    if idx in [idx1, idx2]:
      if idx1 not in incl:
        incl.append(idx1)
      if idx2 not in incl:
        incl.append(idx2)
  return incl

def check_if_aromatic(mol,r):
  atom_cond=True
  for idx in r:
    if not mol.GetAtomWithIdx(idx).GetIsAromatic():
      atom_cond=False
      break
    #print('{} : {}'.format(idx,mol.GetAtomWithIdx(idx).GetIsAromatic()))
  if atom_cond:
    #print('all atoms are aromatic')
    return True
  else:
    for bond in mol.GetBonds():
      idx1=bond.GetBeginAtomIdx()
      idx2=bond.GetEndAtomIdx()
      #print('bond between {} and {} is {}'.format(idx1, idx2, bond.GetBondType()))
      if idx1 in r and idx2 in r:
        if bond.GetBondType() != 'AROMATIC':
          return False
  #print('all bonds are aromatic')
  return True

def check_if_oxy(i,mol,atoms,atype):
  o_bonds=[0,0]
  for x in molneighbors(i,mol):
    if int(atoms[x]) == 8:
      o_bonds[0]+=1
      if mol.GetBondBetweenAtoms(i,x).GetBondTypeAsDouble() > 1.0:
        o_bonds[1]+=1
  #print('  oxygen bonds : {}'.format(o_bonds))
  #print(atype, o_bonds)
  if atype == 6:
    return o_bonds[0] > 1 and o_bonds[1] > 0
  elif atype == 7:
    return o_bonds[0] > 0
  else:
    return o_bonds[0] > 1 or o_bonds[1] > 0

def coarse_grain(idx,atoms,mol,cg_opts,cut_bonds_btwn_groups=True):
  '''
  coarse grain module
  options:
    rings : preserves all rings
    rings# : preserves only #s specified
    aromatic : preserves only aromatic rings
    halogen : preserves any bond to [F, Cl, Br, I, At, Ts]
    nitro : preserves nitro groups
    sulfo : preserves sulfoxide and sulfone groups
    phospho : preserves phosphate groups
    carbo : preserves carboxyl groups
  '''
  if len(cg_opts) > 0:
    #print('Coarse graining')
    #print(cg_opts)
    incl=list()
  else:
    #print('not coarse graining')
    return([[idx]])
  ring_sizes, oxygen_groups = list(), list()
  for opt in cg_opts:
    if 'rings' in opt:
      size = opt.replace('rings','')
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
    '''
    print('  rings',end=' ')
    if arom_only:
      print('(aromatic only)')
    else:
      print()
    if len(ring_sizes) > 0:
      print('    ring sizes to include: {}'.format(ring_sizes))
    '''
    if mol.GetAtomWithIdx(idx).IsInRing():
      ssr=Chem.GetSymmSSSR(mol)
      for r in ssr:
        if idx in r:
          if len(ring_sizes) > 0:
            if len(r) in ring_sizes:
              incl.append([idx]+[x for x in r])
              continue
          elif arom_only:
            if check_if_aromatic(mol,r): 
              incl.append([idx]+[x for x in r])
              continue
          else:
            incl.append([idx]+[x for x in r])
  if len(incl) < 1:
    incl=[[idx]]
  #incl = [x for i, x in enumerate(incl) if incl[i::].count(x) < 2]
  if 'halogen' in cg_opts or 'halogens' in cg_opts:
    #print('  halogens')
    halogens=[9, 17, 35, 53, 85, 117]
    #if int(atoms[idx]) in halogens:
    for idx_, inc in enumerate(incl):
      new_inc=inc.copy()
      for i in inc:
        #incl=molneighbors(idx,mol)
        if int(atoms[i]) in halogens:
          for x in molneighbors(i,mol):
            if x not in new_inc:
              new_inc.append(x)
              #incl[jdx].append(x)
          #incl.extend(molneighbors(i,mol))
        else:
          incl_=molneighbors(i,mol)
          #incl_=molneighbors(idx,mol)
          for jdx in incl_:
            if int(atoms[jdx]) in halogens:
              if jdx not in new_inc:
                #incl.append([idx]+jdx)
                new_inc.append(jdx)
      incl[idx_]=new_inc
  if len(oxygen_groups) > 0:
    '''
    nitro / sulfo / phospho
    '''
    '''
    for opt in ['nitro', 'sulfo', 'phospho']:
      if opt in cg_opts:
        print('  {}'.format(opt))
    '''
    for idx_, inc in enumerate(incl):
      new_inc=inc.copy()
      for i in inc:
        #incl=molneighbors(idx,mol)
        if int(atoms[i]) in oxygen_groups:
          if check_if_oxy(i,mol,atoms,int(atoms[i])):
            for x in molneighbors(i,mol):
              if cut_bonds_btwn_groups:
                if int(atoms[x]) == 8:
                  new_inc.append(x)
              elif x not in new_inc:
                new_inc.append(x)
        elif not cut_bonds_btwn_groups:
          incl_=molneighbors(i,mol)
          #incl_=molneighbors(idx,mol)
          for jdx in incl_:
            if int(atoms[jdx]) in oxygen_groups:
              if check_if_oxy(jdx,mol,atoms,int(atoms[i])):
                for x in [jdx]+molneighbors(jdx, mol):
                  if cut_bonds_btwn_groups:
                    if int(atoms[x]) == 8:
                      new_inc.append(x)
                      if jdx not in new_inc:
                        new_inc.append(jdx)
                  elif x not in new_inc:
                    #incl.append([idx]+jdx)
                    new_inc.append(x)
                    if jdx not in new_inc:
                      new_inc.append(jdx)
        #print(atoms[i],new_inc,[atoms[x] for x in molneighbors(i,mol)])
      incl[idx_]=new_inc
  #if len(cg_opts) > 0:
    #print('  incl:{}'.format(incl))
  return incl
