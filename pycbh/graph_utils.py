import sys, os
import numpy as np
from pycbh.gaussian_expansion import gaussian_expansion
from pycbh.smi2attr import smi2attr
#from pycbh.smi2attr import smi2mol2vec


def graphs_from_npy(fns):
    """
    Reads a list of graph representations from NPY files.

    Parameters:
    - fns (list or str): List of NPY file paths or a single NPY file path.

    Returns:
    - tuple: A tuple containing two lists:
        - list: List of graph representations.
        - list: List of loaded file paths.
    """
    if type(fns) != list:
        fns = [fns]
    graphs_ls, loaded_ls = list(), list()

    np_load_old = np.load

    np.load = lambda *a, **k: np_load_old(*a, allow_pickle=True, **k)

    for fn in fns:
        try:
            g = np.load(fn).tolist()
            if type(g) == list:
                graphs_ls.extend(g)
            else:
                graphs_ls.append(g)
            #pprint_graph(g)
            loaded_ls.append(fn)
        except:
            pass

    #print('Loaded {}/{} graphs'.format(len(graphs_ls),len(fns)))

    np.load = np_load_old

    return graphs_ls, loaded_ls


def get_atomtypes(graphs, bad_types=[None]):
    """
    Extracts unique atom types from the edges of a list of graph representations.

    Parameters:
    - graphs (list): List of graph representations.
    - bad_types (list): List of atom types to be ignored.

    Returns:
    - tuple: A tuple containing two elements:
        - list: Sorted list of unique atom types.
        - list: Modified list of graph representations.
    """
    atomtypes = list()
    for idx1, graph in enumerate(graphs):
        rmd_idx = list()
        for idx2, edge in enumerate(graph['edges']):
            if edge[0] in bad_types or edge[1] in bad_types:
                rmd_idx.append(idx2)
                #pprint.pprint(graph)
            else:
                if edge[0] not in atomtypes:
                    atomtypes.append(edge[0])
                if edge[1] not in atomtypes:
                    atomtypes.append(edge[1])
        if len(rmd_idx) > 0:
            graph_new = dict()
            nodes, edges, senders, receivers = list(), list(), list(), list()
            for idx2, edge in enumerate(graph['edges']):
                if idx2 not in rmd_idx:
                    edges.append(edge)
                    senders.append(graph['senders'][idx2])
                    receivers.append(graph['receivers'][idx2])
            good_nodes = sorted(set(senders + receivers))
            rmd_nodes = list()
            for idx2, node in enumerate(graph['nodes']):
                if idx2 in good_nodes:
                    nodes.append(node)
                else:
                    rmd_nodes.append(idx2)
            rmd_nodes.sort()
            for idx2, sender in enumerate(senders):
                for x in rmd_nodes:
                    if sender > x:
                        senders[idx2] -= 1
            for idx2, receiver in enumerate(receivers):
                for x in rmd_nodes:
                    if receiver > x:
                        receivers[idx2] -= 1
            graph_new = graph
            graph_new['nodes'] = nodes
            graph_new['edges'] = edges
            graph_new['senders'] = senders
            graph_new['receivers'] = receivers
            #print('\n\ngraph_new:')
            #pprint.pprint(graph_new)
            #print('rmd node(s): {}'.format(' '.join([str(x) for x in rmd_nodes])))
            graphs[idx1] = graph_new

    return sorted(atomtypes), graphs


def atom2onehot(atomtypes, atom)
    """
    Converts an atom type into a one-hot encoded representation.

    Parameters:
    - atomtypes (list): List of unique atom types.
    - atom: Atom type to be encoded.

    Returns:
    - list: One-hot encoded representation of the atom type.
    """
    return [int(atom == x) for x in atomtypes]


def transform_edges(graphs,
                    atom_types=None,
                    fp_dict=dict(),
                    bad_types=[None],
                    use_GauExp=False,
                    use_bondinverse=True,
                    use_z=True,
                    use_onehotatom=True,
                    use_bond_type=True,
                    use_onehotbond=True,
                    use_overlapfp=False,
                    use_coulomb=True,
                    use_bond_length=True,
                    which_rung=0,
                    which_fp='default'):
    """
    Transforms edge information in a list of graph representations based on specified options.

    Parameters:
    - graphs (list): List of graph representations.
    - atom_types (list or None): List of atom types to consider.
    - fp_dict (dict): Dictionary to store fingerprint information.
    - bad_types (list): List of atom types to be ignored.
    - use_GauExp (bool): Flag to use Gaussian expansion.
    - use_bondinverse (bool): Flag to use inverse of bond length.
    - use_z (bool): Flag to include atomic numbers in edge information.
    - use_onehotatom (bool): Flag to use one-hot encoded atom types.
    - use_bond_type (bool): Flag to include bond types in edge information.
    - use_onehotbond (bool): Flag to use one-hot encoded bond types.
    - use_overlapfp (bool): Flag to use overlapping fingerprints.
    - use_coulomb (bool): Flag to include Coulombic interaction in edge information.
    - use_bond_length (bool): Flag to include bond length in edge information.
    - which_rung (int): Rung index for fingerprint calculation.
    - which_fp (str or list): Type(s) of fingerprints to use.

    Returns:
    - tuple: A tuple containing two elements:
        - list: Modified list of graph representations.
        - dict: Updated fingerprint dictionary.
    """
    atomtypes, graphs = get_atomtypes(graphs, bad_types=bad_types)
    #atomtypes = [6, 7, 8, 9, 16, 17]
    print(atomtypes, end='')
    if atom_types is not None:
        atomtypes = [x for x in atom_types if x not in bad_types]
        print(' -> {}'.format(atomtypes), end='')
    print()
    '''
  edge:
  [6, 6, 1.0, 1.52618438, 23.58823775, ['CBH-1', 'CC', 'f2_C2H6_000']] 
  '''
    for idx1, graph in enumerate(graphs):
        for idx2, edge in enumerate(graph['edges']):
            bond_type = edge[2]
            bond_length = edge[3]
            coulomb = edge[4]
            new_edge = list()
            if use_z:
                new_edge = [edge[0], edge[1]]
            if use_bond_type:
                new_edge.append(bond_type)
            if use_onehotatom:
                new_edge.extend(atom2onehot(atomtypes, edge[0]))
                new_edge.extend(atom2onehot(atomtypes, edge[1]))
            if use_onehotbond:
                new_edge.extend(
                    atom2onehot([0.0, 1.0, 1.5, 2.0, 3.0, 4.0], bond_type))
            if use_bond_length:
                new_edge.append(bond_length)
            if use_coulomb:
                new_edge.append(coulomb)
            if use_GauExp:
                new_edge.extend(gaussian_expansion(bond_length)[0])
            elif use_bondinverse:
                new_edge.append(1 / bond_length)
            if use_overlapfp:
                smi, fn = None, None
                for x in node:
                    if type(x) is list:
                        if x[0] == 'CBH-{}'.format(which_rung - 1):
                            smi, fn = x[1], x[2]
                fp = None
                if smi is not None or fn is not None:
                    if which_fp in fp_dict:
                        if smi in fp_dict[which_fp]:
                            fp = fp_dict[which_fp][smi]
                        elif fn in fp_dict[which_fp]:
                            fp = fp_dict[which_fp][fn]
                    if fp is None:
                        fp, fp_dict = get_fp(smi=smi,
                                             fn=fn,
                                             which_fp=which_fp,
                                             fp_dict=fp_dict)
                #print('get_fp({})'.format(edge[3]))
                #new_edge.extend(get_fp(fp_dict, edge[3], which_fp='mol2vec'))

            #print('\n\nold edge: {}'.format(edge))
            #print('\nnew edge: {}'.format(new_edge))
            graph['edges'][idx2] = new_edge
        graphs[idx1] = graph
    return graphs, fp_dict


def get_fp(smi=None, fn=None, which_fp='default', fp_dict=dict()):
    """
    Obtains a fingerprint for a given SMILES string or filename and updates the fingerprint dictionary.

    Parameters:
    - smi (str or None): SMILES string.
    - fn (str or None): Filename.
    - which_fp (str): Type of fingerprint to use.
    - fp_dict (dict): Dictionary to store fingerprint information.

    Returns:
    - tuple: A tuple containing two elements:
        - list: Fingerprint.
        - dict: Updated fingerprint dictionary.
    """
    if which_fp not in fp_dict:
        fp_dict[which_fp] = dict()
    if smi is None and fn is None:
        fp = [None]
    elif smi is not None:
        if which_fp == 'default':
            fp = smi2attr(smi)
            '''
    elif which_fp=='mol2vec':
      if 'model' in fp_dict['mol2vec']:
        fp, fp_dict['mol2vec']['model'] = smi2mol2vec(smi, model=fp_dict['mol2vec']['model'])
      else:
        fp, fp_dict['mol2vec']['model'] = smi2mol2vec(smi)
      '''
        elif which_fp in fp_dict:
            fp = [None]
        fp_dict[which_fp][smi] = fp
    elif fn is not None:
        if which_fp in fp_dict:
            fp = [None]
        fp_dict[which_fp][fn] = fp
    return fp, fp_dict


def transform_nodes(graphs,
                    atom_types=None,
                    bad_types=[None],
                    fp_dict=dict(),
                    which_rung=0,
                    use_fp=True,
                    which_fp='default'):
    """
    Transforms node information in a list of graph representations based on specified options.

    Parameters:
    - graphs (list): List of graph representations.
    - atom_types (list or None): List of atom types to consider.
    - bad_types (list): List of atom types to be ignored.
    - fp_dict (dict): Dictionary to store fingerprint information.
    - which_rung (int): Rung index for fingerprint calculation.
    - use_fp (bool): Flag to use fingerprints.
    - which_fp (str or list): Type(s) of fingerprints to use.

    Returns:
    - tuple: A tuple containing two elements:
        - list: Modified list of graph representations.
        - dict: Updated fingerprint dictionary.
    """
    atomtypes, graphs = get_atomtypes(graphs, bad_types=bad_types)
    print(atomtypes, end='')
    if atom_types is not None:
        atomtypes = [x for x in atom_types if x not in bad_types]
        print(' -> {}'.format(atomtypes), end='')
    print()
    '''
  node:
  [0, 6.0, [0], ['CBH-0', 'C', 'f1_CH4_000'], ['CBH-2', 'CC', 'f2_C2H6_000'], ['CBH-4', 'CCC', 'f3_C3H8_000']] 
  '''
    for idx1, graph in enumerate(graphs):
        for idx2, node in enumerate(graph['nodes']):
            new_node = list()
            new_node.append(node[1])
            smi, fn = None, None
            for x in node:
                if type(x) is list:
                    if x[0] == 'CBH-{}'.format(which_rung):
                        smi, fn = x[1], x[2]
            if use_fp:
                if type(which_fp) != list:
                    which_fp = [which_fp]
                for fp_type in which_fp:
                    fp = None
                    if fp_type in fp_dict:
                        if smi in fp_dict[fp_type]:
                            fp = fp_dict[fp_type][smi]
                        elif fn in fp_dict[fp_type]:
                            fp = fp_dict[fp_type][fn]
                    if fp is None:
                        fp, fp_dict = get_fp(smi=smi,
                                             fn=fn,
                                             which_fp=fp_type,
                                             fp_dict=fp_dict)
                    new_node.extend(fp)

            graph['nodes'][idx2] = new_node
        graphs[idx1] = graph
    return graphs, fp_dict


def transform_globals(graphs,
                      atom_types=None,
                      fp_dict=dict(),
                      use_fp=True,
                      which_fp='default'):
    """
    Transforms global information in a list of graph representations based on specified options.

    Parameters:
    - graphs (list): List of graph representations.
    - atom_types (list or None): List of atom types to consider.
    - fp_dict (dict): Dictionary to store fingerprint information.
    - use_fp (bool): Flag to use fingerprints.
    - which_fp (str or list): Type(s) of fingerprints to use.

    Returns:
    - tuple: A tuple containing two elements:
        - list: Modified list of graph representations.
        - dict: Updated fingerprint dictionary.
    """
    for idx1, graph in enumerate(graphs):
        new_global = graph['globals'].copy()
        smi, fn = graph['globals'][1], graph['globals'][0]
        if use_fp:
            if type(which_fp) != list:
                which_fp = [which_fp]
            for fp_type in which_fp:
                fp = None
                if fp_type in fp_dict:
                    if smi in fp_dict[fp_type]:
                        fp = fp_dict[fp_type][smi]
                    elif fn in fp_dict[fp_type]:
                        fp = fp_dict[fp_type][fn]
                if fp is None:
                    fp, fp_dict = get_fp(smi=smi,
                                         fn=fn,
                                         which_fp=fp_type,
                                         fp_dict=fp_dict)
                new_global.extend(fp)
        graphs[idx1]['globals'] = new_global
    return graphs, fp_dict
