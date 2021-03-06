import numpy as np
import sys
import json
from pycbh import pycbh
import pprint

def remove_values_from_list(the_list, val):
   return [value for value in the_list if value != val]

def smifromfn(smi_file, how_many=None):
  smiles_arr=[x.split(' ') for x in list(filter(None,open(smi_file, "r").read().split("\n")))]
  if how_many == None:
    how_many = len(smiles_arr)
  return smiles_arr[-1*how_many:]

def smi2cbhvec(smi_ls, cbh_rung):
  '''
  should pass
  [ label, smi]
  '''
  if type(smi_ls[0]) != list:
    smi_ls = [smi_ls]
  good_smi, cbh_ls = list(), list()
  for smi in smi_ls:
    try:
      left, right = pycbh.cbh_store2vec(pycbh.smi2cbh(smi[1], cbh_rung))[0]
      good_smi.append(smi)
      cbh_ls.append([left,right])
    except:
      print('Failed : {}'.format(smi))
  return good_smi, cbh_ls

def xyz2cbhvec(fn_ls, cbh_rung):
  if type(fn_ls) != list:
    fn_ls = [fn_ls]
  good_fn, cbh_ls = list(), list()
  for fn in fn_ls:
    try:
      left, right = pycbh.cbh_store2vec(pycbh.xyz2cbh(fn, cbh_rung))[0]
      good_fn.append(fn)
      cbh_ls.append([left,right])
    except:
      print('Failed : {}'.format(fn))
  return good_fn, cbh_ls

def cbhvec2multihot(cbh_ls):
  left, right = list(), list()
  for l, r in cbh_ls:
    left.extend(x for x in l if x not in left)
    right.extend(x for x in r if x not in right)
  left.sort()
  right.sort()
  cbh_mh = list()
  for cbh in cbh_ls:
    mh = [cbh[1].count(x) for x in right] + [-1*cbh[0].count(x) for x in left]
    cbh_mh.append(mh)
    #print(cbh)
    #print(mh)
  return cbh_mh, right+left

def smi2mh(smiles_arr, rung):
  good_smi, cbh_ls = smi2cbhvec(smiles_arr, rung)
  print('multihot vector')
  #print('good_smi : {}'.format(good_smi))
  #print('cbh_ls   : {}'.format(cbh_ls))
  cbh_mh, labels = cbhvec2multihot(cbh_ls)
  return cbh_mh, good_smi, labels

def xyzfn2mh(fn_ls, rung):
  good_fn, cbh_ls = xyz2cbhvec(fn_ls, rung)
  print('multihot vector')
  #print('good_fn : {}'.format(good_fn))
  #print('cbh_ls   : {}'.format(cbh_ls))
  cbh_mh, labels = cbhvec2multihot(cbh_ls)
  return cbh_mh, good_fn, labels

if __name__=='__main__':
  if len(sys.argv[1::]) < 2:
    sys.exit('USAGE: python3 vec_utils.py fn_type FILENAME(s)')
  fn_type = sys.argv[1]
  fn_ls = sys.argv[2::]
  save_vecs = True
  fn_dir  = 'data/'
  fn_base = 'MLCBH'
  rung = 2
  if fn_type == 'smi':
    smi_fn = sys.argv[2]
    smiles_arr=smifromfn(smi_file, 100)
    cbh_mh, good_x, labels = smiarr2mh(smiles_arr, rung)
  elif fn_type == 'xyz':
    fn_ls = sys.argv[2::]
    cbh_mh, good_x, labels = xyzfn2mh(fn_ls, rung)
  else:
    sys.exit('File type ({}) not supported'.format(fn_type))
  if save_vecs:
    x_inputs=open(fn_dir+fn_base.replace('CBH','CBH{}'.format(rung))+'_inputs.txt',"w")
    y_output=open(fn_dir+fn_base.replace('CBH','CBH{}'.format(rung))+'_outputs.txt', "w")
    for idx, v in enumerate(cbh_mh):
      x_inputs.write(" ".join([str(x) for x in v])+"\n")
      y_output.write(" ".join([str(x) for x in good_x[idx]])+"\n")
    x_inputs.close()
    y_output.close()

    print('Length of good smiles : {}'.format(len(good_x)))

    with open(fn_dir+fn_base.replace('CBH','CBH{}'.format(rung))+'_labels.txt',"w") as fn:
      fn.write(" ".join([str(x) for x in labels])+"\n")

    print('Length of input vec : {}'.format(len(labels)))


