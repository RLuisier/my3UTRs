#!/usr/bin/python

import argparse;
import numpy as np;
from ctypes import *;
import os;
import random;
#import statsmodels.api as sm;
from array import array;
import linecache as lc;
from numpy import diff;
import time
from scipy import weave


def import_file(full_path_to_module):
    try:
        import os
        module_dir, module_file = os.path.split(full_path_to_module)
        module_name, module_ext = os.path.splitext(module_file)
        save_cwd = os.getcwd()
        os.chdir(module_dir)
        module_obj = __import__(module_name)
        module_obj.__file__ = full_path_to_module
        globals()[module_name] = module_obj
        os.chdir(save_cwd)
    except:
        raise ImportError

import_file('./scripts/DetectRegions')

t1 = time.time()

parser=argparse.ArgumentParser(description='Find regions along the chromosome which have reads above a certain threshold for each position at certain percetage P');
parser.add_argument('--Cov','-i', help='Input Coverage file (use - from standard input)');
parser.add_argument('--winS','-s',type=int,help='WindowSize');
parser.add_argument('--Tr','-t',type=int,help='Number of reads minimum per nucleotide');
parser.add_argument('--Comments', '-c',type=int, help='If 1 then will print the comments');
parser.add_argument('--out_range','-o', help='Write as an output.');

args=parser.parse_args();


#Input variables from input:
lim       = float(args.Tr)/args.winS;
winS      = int(args.winS);
output    = args.out_range;
toRead    = args.Cov;
comment   = int(args.Comments);
Lim       = int(args.Tr);


out = DetectRegions.DetectRegions(lim,Lim,winS,toRead,comment)

#Write output:
np.savetxt(output, out, delimiter=' ', header=('Start End'), comments='')
t2 = time.time()
deltaT=t2-t1;
print 'Processing time:'+str(deltaT);
