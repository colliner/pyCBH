import numpy as np
import sys
import json
from pycbh import pycbh
import pprint


def remove_values_from_list(the_list, val):
    return [value for value in the_list if value != val]


def smifromfn(smi_file, how_many=None):
    """
    Extract SMILES strings from a file.

    Parameters:
    - smi_file (str): File containing SMILES strings.
    - how_many (int): Number of SMILES strings to extract.

    Returns:
    - list: List of SMILES strings.
    """
    smiles_arr = [
        x.split(' ')
        for x in list(filter(None,
                             open(smi_file, "r").read().split("\n")))
    ]
    if how_many == None:
        how_many = len(smiles_arr)
    return smiles_arr[-1 * how_many:]


def smi2cbhvec(smi_ls, cbh_rung):
    """
    Convert SMILES strings to CBH vectors.

    Parameters:
    - smi_ls (list): List of SMILES strings.
    - cbh_rung: CBH rung.

    Returns:
    - tuple: Tuple containing lists of good SMILES and CBH vectors.
    """
    if type(smi_ls[0]) != list:
        smi_ls = [smi_ls]
    good_smi, cbh_ls = list(), list()
    for smi in smi_ls:
        try:
            _, left, right = pycbh.cbh_store2vec(pycbh.smi2cbh(
                smi[1], cbh_rung))[0]
            good_smi.append(smi)
            cbh_ls.append([left, right])
        except:
            print('Failed : {}'.format(smi))
    return good_smi, cbh_ls


def xyz2cbhvec(fn_ls, cbh_rung):
    """
    Convert XYZ files to CBH vectors.

    Parameters:
    - fn_ls (list): List of file names.
    - cbh_rung: CBH rung.

    Returns:
    - tuple: Tuple containing lists of good file names and CBH vectors.
    """
    if type(fn_ls) != list:
        fn_ls = [fn_ls]
    good_fn, cbh_ls = list(), list()
    for fn in fn_ls:
        try:
            _, left, right = pycbh.cbh_store2vec(pycbh.xyz2cbh(fn, cbh_rung))[0]
            good_fn.append(fn)
            cbh_ls.append([left, right])
        except:
            print('Failed : {}'.format(fn))
    return good_fn, cbh_ls


def cbhvec2multihot(cbh_ls):
    """
    Convert CBH vectors to multihot representation.

    Parameters:
    - cbh_ls (list): List of CBH vectors.

    Returns:
    - tuple: Tuple containing the multihot representation and labels.
    """
    left, right = list(), list()
    for l, r in cbh_ls:
        left.extend(x for x in l if x not in left)
        right.extend(x for x in r if x not in right)
    left.sort()
    right.sort()
    cbh_mh = list()
    for cbh in cbh_ls:
        mh = [cbh[1].count(x) for x in right
             ] + [-1 * cbh[0].count(x) for x in left]
        cbh_mh.append(mh)
    return cbh_mh, right + left


def smi2mh(smiles_arr, rung):
    """
    Convert SMILES strings to multihot representation.

    Parameters:
    - smiles_arr (list): List of SMILES strings.
    - rung: CBH rung.

    Returns:
    - tuple: Tuple containing the multihot representation, good SMILES, and labels.
    """
    good_smi, cbh_ls = smi2cbhvec(smiles_arr, rung)
    print('multihot vector')
    cbh_mh, labels = cbhvec2multihot(cbh_ls)
    return cbh_mh, good_smi, labels


def xyzfn2mh(fn_ls, rung):
    """
    Convert XYZ file names to multihot representation.

    Parameters:
    - fn_ls (list): List of file names.
    - rung: CBH rung.

    Returns:
    - tuple: Tuple containing the multihot representation, good file names, and labels.
    """
    good_fn, cbh_ls = xyz2cbhvec(fn_ls, rung)
    print('multihot vector')
    cbh_mh, labels = cbhvec2multihot(cbh_ls)
    return cbh_mh, good_fn, labels


def cbh_store2fndict(cbh_store, key_fn=None):
    """
    Convert CBH store to a dictionary of file names.

    Parameters:
    - cbh_store (list): List of CBH vectors.
    - key_fn (str): Key file name.

    Returns:
    - dict: Dictionary of file names.
    """
    if type(cbh_store) != list:
        cbh_store = [cbh_store]
    if type(cbh_store[0]) != list:
        cbh_store = [cbh_store]
    cbh_store = pycbh.cbh_store2vec(cbh_store)
    frag_dict = dict()
    for frag_vecs in pycbh.load_frags(key_fn).values():
        for frag_vec in frag_vecs:
            frag_dict[frag_vec[1]] = frag_vec[-1]
    '''
  for key, val in frag_dict.items():
    print('{} : {}'.format(key,val))
  sys.exit()
  '''
    fn_dicts = dict()
    for cbh_s in cbh_store:
        fn_dict = dict()
        label, left, right = cbh_s[0], cbh_s[1], cbh_s[2]
        for smi in left:
            if smi not in frag_dict:
                print('error: {} not in frag_dict'.format(smi))
                fn_dict['fragment_error_flag'] = -1
            else:
                fn = frag_dict[smi]
                if fn not in fn_dict:
                    fn_dict[fn] = -1
                else:
                    fn_dict[fn] -= 1
        for smi in right:
            if smi not in frag_dict:
                print('error: {} not in frag_dict'.format(smi))
                fn_dict['fragment_error_flag'] = -1
            else:
                fn = frag_dict[smi]
                if fn not in fn_dict:
                    fn_dict[fn] = 1
                else:
                    fn_dict[fn] += 1
        if label not in fn_dicts:
            fn_dicts[label] = fn_dict
        else:
            print('ERROR: {} already in fn_dicts'.format(label))
    return fn_dicts


if __name__ == '__main__':
    if len(sys.argv[1::]) < 2:
        sys.exit('USAGE: python3 cbhvec_utils.py fn_type FILENAME(s)')
    fn_type = sys.argv[1]
    fn_ls = sys.argv[2::]
    save_vecs = True
    fn_dir = 'data/'
    fn_base = 'MLCBH'
    rung = 2
    if fn_type == 'smi':
        smi_fn = sys.argv[2]
        smiles_arr = smifromfn(smi_file, 100)
        cbh_mh, good_x, labels = smiarr2mh(smiles_arr, rung)
    elif fn_type == 'xyz':
        fn_ls = sys.argv[2::]
        cbh_mh, good_x, labels = xyzfn2mh(fn_ls, rung)
    else:
        sys.exit('File type ({}) not supported'.format(fn_type))
    if save_vecs:
        x_inputs = open(
            fn_dir + fn_base.replace('CBH', 'CBH{}'.format(rung)) +
            '_inputs.txt', "w")
        y_output = open(
            fn_dir + fn_base.replace('CBH', 'CBH{}'.format(rung)) +
            '_outputs.txt', "w")
        for idx, v in enumerate(cbh_mh):
            x_inputs.write(" ".join([str(x) for x in v]) + "\n")
            y_output.write(" ".join([str(x) for x in good_x[idx]]) + "\n")
        x_inputs.close()
        y_output.close()

        print('Length of good smiles : {}'.format(len(good_x)))

        with open(
                fn_dir + fn_base.replace('CBH', 'CBH{}'.format(rung)) +
                '_labels.txt', "w") as fn:
            fn.write(" ".join([str(x) for x in labels]) + "\n")

        print('Length of input vec : {}'.format(len(labels)))
