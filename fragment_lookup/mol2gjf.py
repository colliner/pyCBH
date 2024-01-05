import sys, os
import subprocess
import re


def mol2gjf(fns):
    if type(fns) != list:
        fns = [fns]
    goodf, badf = list(), list()
    for fn in fns:
        try:
            outf = fn.replace('mol', 'gjf')
            chkf = fn.replace('mol', 'chk')
            inlines = open(fn, 'r').read().split('\n')
            charge = 0
            spin_m = 1
            for line in inlines:
                if re.findall('^M  CHG ', line):
                    charge = line.split(' ')[-1]
                elif re.findall('^M  RAD ', line):
                    spin_m = line.split(' ')[-1]
            bashCommand = 'obabel -imol {} -oxyz --gen3d -xb'.format(fn)
            process = subprocess.Popen(bashCommand.split(),
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            output, error = process.communicate()
            output = output.decode("utf-8").split("\n")
            with open('gjf_files/' + outf, 'w') as FILE:
                FILE.write('%chk={}\n%mem=7GB\n%nproc=2\n'.format(chkf))
                FILE.write('# g4\n\n'
                          )  #opt freq=(scale=0.9854) b3lyp/6-31G(2df,p)\n\n')
                FILE.write('pyCBH fragment file ({}) : {}\n\n'.format(
                    fn, output[1]))
                FILE.write('{} {}\n'.format(charge, spin_m))
                FILE.write('\n'.join(output[2::]) + '\n')
            print('written to gjf_files/{}'.format(outf))
            goodf.append(outf)
        except:
            print('ERROR: {} failed to write'.format(outf))
            badf.append(outf)
            pass
    return [goodf, badf]


if __name__ == '__main__':
    if len(sys.argv[1::]) < 1:
        sys.exit('USAGE: python3 mol2gjf.py FILENAME(s)')
    fns = sys.argv[1::]
    _ = mol2gjf(fns)
