#!/usr/local/stackless-python/bin/python

import os
import os.path
import sys
import subprocess
import stackless
from optparse import OptionParser

def run_minkoe(index, options):
    while index<len(files):
        filename = options.dir + files[index]
        args = home+options.minkoe + ' -f ' + filename + ' -e ' + str(options.expfactor) + ' -r ' + str(options.maxdil) + ' -b'
        print 'Running', ''.join(args)
        mink = subprocess.Popen(args,cwd=home,shell=True)
        while mink.poll()==None:
            stackless.schedule()
        subprocess.call('pwd', shell=True)
        mnk = os.path.splitext(filename)[0] + '.mnk ' 
        subprocess.call('mv '+mnk+results, shell=True)
        print 'Finished', filename
        index += options.nthreads

argparser = OptionParser()
argparser.add_option('-x', dest="minkoe", type='string', help="MinkOE executable path", default='minkoe')
argparser.add_option('-l', dest="pdblist", type='string', help="Filename of list of PDBs to process (master.map)", default=None)
argparser.add_option('-d', dest="dir", type='string', help="Directory where PDB files are found", default='./')
argparser.add_option('-n', dest="nthreads", type='int', help="Number of threads to run simultaneously", default=5)
argparser.add_option('-e', dest="expfactor", type='float', help="Number of threads to run simultaneously", default=1.0)
argparser.add_option('-r', dest="maxdil", type='int', help="Number of threads to run simultaneously", default=100)

(options, args) = argparser.parse_args()

if not options.pdblist:
    files = os.listdir(options.dir)
else:
    mapfile = open(options.pdblist)
    files = [ line.split()[-1]+'.pdb' for line in mapfile if line.split()[-1].isdigit() ]
    mapfile.close()
    
home = os.getcwd() + '/'

if os.path.isdir("results"):
    pass
else:
    os.mkdir("results")
    print home+"/results/ directory created."

results = home+"results/"
count = 0
fullcount = 0

for i in range(options.nthreads):
    stackless.tasklet(run_minkoe)(i,options)
    
stackless.run()
 
    
    
