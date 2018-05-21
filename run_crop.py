import copy
from pprint import pprint
import pickle
import os
import sys

import numpy as np

import ROOT
ROOT.gStyle.SetOptFit(111)
ROOT.gStyle.SetOptTitle(1)

root_obj = []

def get_rundata():
  files = ['/'.join(['/media/hdd1/rshara01', 'traindata', name]) for name in os.listdir('/media/hdd1/rshara01/traindata')]
  files = sorted(files)
  i = len(files)

  #print [str(files[a].split('/')[5].split('_')[1].split('.')[0]) for a in xrange(10)]
  run_nums = [str(f.split('/')[5].split('_')[1]) for f in files]
  indexes = [files.index(f) for f in files]
  
  return {'files': files, 'run_nums': run_nums, 'indexes': indexes}


def main(argv=None):
  if argv is None:
    argv=sys.argv[1:]

  ROOT.gROOT.SetBatch(1)
  
  rundata = get_rundata()

  files = rundata['files']
  idx   = rundata['indexes']
  #outfiles = ['_'.join(['/media/hdd1/rshara01/larflow/train_cropped/larcv_test', '%03d'%(i)]) for i in idx ]
  outfiles = ['_'.join(['/media/hdd1/rshara01/larflow/cropper/files/larcv_test', '%03d'%(i)]) for i in idx ]
  outfiles = ['.'.join([out, 'root']) for out in outfiles]

  for a in [0]:
    args = "%s %s"%(files[a],outfiles[a])
    print args
    os.system("./bin/larflow_imgcrop %s"%(args))
    
if __name__ == '__main__':
  main()


