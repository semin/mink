#!/usr/local/stackless-python/bin/python

import os
import os.path
import sys
import subprocess
import stackless
from optparse import OptionParser

def run_minkoe(index, options):
    while index < len(files):
        pdb = os.path.join(idir, files[index])
        mnk = os.path.join(odir, os.path.splitext(os.path.basename(pdb))[0] + '.mnk')
        arg = options.mink + ' ' + pdb + ' > ' + mnk
        print 'Processing ' + arg
        pop = subprocess.Popen(arg, cwd=home, shell=True)
        while pop.poll() == None:
            stackless.schedule()
        print 'Finished', pdb, " (%d/%d)" % (index+1, len(files))
        index += options.nthreads

argparser = OptionParser()
argparser.add_option('-x', dest="mink", type='string',
        help="Mink executable path", default='mink')
argparser.add_option('-l', dest="pdblist", type='string',
        help="Filename of list of PDBs to process (pdblist.txt)", default=None)
argparser.add_option('-i', dest="idir", type='string',
        help="Directory where PDB files are found", default='./')
argparser.add_option('-o', dest="odir", type='string',
        help="Directory where MINK files are created", default='./mink')
argparser.add_option('-n', dest="nthreads", type='int',
        help="Number of threads to run simultaneously", default=5)
argparser.add_option('-e', dest="expfactor", type='float',
        help="Number of threads to run simultaneously", default=1.0)
argparser.add_option('-r', dest="maxdil", type='int',
        help="Number of threads to run simultaneously", default=100)

(options, args) = argparser.parse_args()

if not options.pdblist:
    files = os.listdir(options.idir)
else:
    mapfile = open(options.pdblist)
    files   = [ line.split()[-1]+'.pdb' for line in mapfile if line.split()[-1].isdigit() ]
    mapfile.close()
    
home = os.getcwd()
idir = os.path.join(home, options.idir)
odir = os.path.join(home, options.odir)

if os.path.isdir(odir):
    pass
else:
    os.mkdir(odir)
    print odir + " directory created."

count       = 0
fullcount   = 0

for i in range(options.nthreads):
    stackless.tasklet(run_minkoe)(i, options)
    
stackless.run()
