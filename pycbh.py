import sys, os
import pycbh
from rdkit import Chem

print_cbh=True
#print_cbh=False
if print_cbh:
  from pycbh import cbh_print
else:
  cbh_print = lambda c : print(['FAIL','PASS','PASS (batch)'][len(c[0:2])])

if __name__=='__main__':
  if len(sys.argv[1::]) < 2:
    sys.exit('Usage python3 pycbh.py CBH_rung FILENAME(s)')
  fns=sys.argv[2::]
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

  cbh=pycbh.files2cbh(fns, rung_ls)
  cbh_print(cbh)

 
