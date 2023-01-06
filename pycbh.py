import sys, os
import pycbh
from rdkit import Chem
import numpy as np

print_cbh=True
##print_cbh=False
if print_cbh:
  from pycbh import cbh_print
else:
  cbh_print = lambda c : print(['FAIL','PASS','PASS (batch)'][len(c[0:2])])

if __name__=='__main__':
  if len(sys.argv[1::]) < 2:
    sys.exit('Usage python3 pycbh.py CBH_rung FILENAME(s)')
  save_graph = False
  coarse_grain = False
  fns=list()
  for arg in sys.argv[2:]:
    if arg == '-g':
      save_graph = True
    elif arg == '-cg':
      coarse_grain = True
    else:
      fns.append(arg)
  try:
    rung_ls=[int(sys.argv[1])]
  except:
    if sys.argv[1]=='all':
      rung_ls=list(range(5))
    elif '+' not in sys.argv[1]:
      sys.exit('Usage: python3 pycbh.py CBH_rung FILENAME(s)\n  CBH rung must be an integer')
    else:
      rung_ls=sys.argv[1]
  '''
  General conversion from files:
    filename(s) (.smi | .xyz | .mol) --> files2cbh --> cbh-n
  '''
  cbh=pycbh.files2cbh(fns, rung_ls, save_graph=save_graph,fully_connected=False,coarse_grain=coarse_grain)
  cbh_print(cbh)
  methods=[['g4(0k)','zpe'],'pbe-d3',['g4(0k)','pbe-d3','zpe']]
  key_fn="fragment_lookup/keys.txt"
  energy_fn="fragment_lookup/lookup_energy.txt"
  cbh_e = pycbh.cbh_store2energy(cbh, key_fn=key_fn, energy_fn=energy_fn, levels_of_theory=methods)
  print('\nMethods : {}'.format(' '.join([x if type(x)!=list else '-'.join(x) for x in methods])))
  for key in sorted(cbh_e.keys()):
    if key != 'methods':
      print('{} {}'.format(key,' '.join([str(x) for x in cbh_e[key]])))
  

