import numpy as np
import sys
import json
from pycbh import vec_utils
from pycbh import pycbh
import pprint


def cbh_vec2uniq(cbh_vec):
    """
    Converts a dictionary of CBH vectors into a dictionary of unique fragments.

    Parameters:
    - cbh_vec (dict): Dictionary containing CBH vectors.

    Returns:
    - dict: Dictionary with unique fragments grouped by CBH rung.
    """
    uniq_dict = {'all': list()}
    for key, value in cbh_vec.items():
        cbh_rung = key.split(' ')[-1].replace(':', '')
        if cbh_rung not in uniq_dict:
            uniq_dict[cbh_rung] = list()
        for x in value:
            if x not in uniq_dict[cbh_rung]:
                uniq_dict[cbh_rung].append(x)
            if x not in uniq_dict['all']:
                uniq_dict['all'].append(x)
    return uniq_dict


def cbh_store2energy(cbh_store,
                     levels_of_theory=['B3LYP', 'HF', 'MP2', 'CCSD(T)'],
                     key_fn="fragment_lookup/keys.txt",
                     energy_fn="fragment_lookup/lookup_energy.txt"):
    """
    Converts a CBH store to a dictionary containing energy information.

    Parameters:
    - cbh_store (dict): CBH store containing fragment coefficients.
    - levels_of_theory (list): List of methods or levels of theory.
    - key_fn (str): Filepath for the keys.
    - energy_fn (str): Filepath for the energy lookup table.

    Returns:
    - dict: Dictionary containing energy information for CBH vectors.
    """
    energy_fn = "fragment_lookup/lookup_energy.txt"
    cbh_e = {'methods': levels_of_theory}
    cbh_vec = vec_utils.cbh_store2fndict(cbh_store, key_fn=key_fn)
    cbh_lookup = pycbh.load_frags(energy_fn)
    #print(cbh_lookup)
    #print(cbh_vec)
    '''
  for v, frags in cbh_vec.items():
    print('\n{}  {}'.format(v,' '.join(list(frags.keys())))) 
  #sys.exit()
  unique_frags = cbh_vec2uniq(cbh_vec)
  for rung, frags in unique_frags.items():
    if rung == 'all':
      outfile='PDS_allfrags.txt'
      with open(outfile, 'w') as OUTF:
        OUTF.write('\n'.join(sorted(frags)))
        print('written to {}'.format(outfile))
    print('\nUnique {} fragments:'.format(rung))
    print('  {}'.format(sorted(frags)))

  sys.exit()
  '''
    print()
    missing_f, missing_m = list(), list()
    for key, val in cbh_vec.items():
        #print('\n{} \n{}'.format(key,val))
        energy = list()
        for method in levels_of_theory:
            if type(method) != list:
                method = [method]
            e = 0.
            bad_e = False
            for idx1, m in enumerate(method):
                if idx1 == 0:
                    sign = 1
                else:
                    sign = -1

                if m not in cbh_lookup['methods'][0]:
                    #print('  {} not in cbh_lookup'.format(m))
                    if m not in missing_m:
                        missing_m.append(m)
                    bad_e = True
                else:
                    #print('{}:'.format(method))
                    idx = cbh_lookup['methods'][0].index(m)
                    if len(val.keys()) <= 1:
                        #print('Full molecule')
                        e = 'Full'
                    else:
                        for mol, coeff in val.items():
                            if mol not in cbh_lookup:
                                #print('  {} not in cbh_lookup'.format(mol))
                                if mol not in missing_f:
                                    missing_f.append(mol)
                                bad_e = True
                                #break
                            else:
                                e += sign * float(coeff) * float(
                                    cbh_lookup[mol][0][idx])
                                #print(mol, m, cbh_lookup[mol][0][idx])
            if bad_e:
                e = None
            else:
                e = round(e, 9)
            energy.append(e)
        if key not in cbh_e:
            cbh_e[key] = energy
        else:
            print('error: {} already in cbh_e'.format(key))
    #print('cbh_e:\n{}'.format(cbh_e))
    len_m = len(missing_m)
    len_f = len(missing_f)
    if len_m + len_f > 0:
        print('Missing entries written to missing_lookup.txt')
        with open('missing_lookup.txt', 'w') as FILE:
            if len_m > 0:
                FILE.write('Missing Method(s):\n  ' +
                           '\n  '.join(sorted(missing_m)) + '\n')
                print('  Missing {} Methods in {}'.format(len_m, energy_fn))
            if len_f > 0:
                FILE.write('Missing Fragment(s):\n  ' +
                           '\n  '.join(sorted(missing_f)) + '\n')
                print('  Missing {} Fragments in {}'.format(len_f, energy_fn))
    return cbh_e
