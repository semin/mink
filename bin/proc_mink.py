#!/usr/bin/env python
import os
import sys

import numpy as np

l = os.listdir(os.getcwd())

files = []
for idf in l:
    name, ext = os.path.splitext(idf)
    if ext == ".mnk":
        files.append(idf)

out_file = open("vectors.dat","w")

i = 0
for file in files:
    name, ext = os.path.splitext(file)
    in_file = open(file, "r")
    lines = in_file.readlines()
    in_file.close()
    
    p = []
    a = []
    b = []
    c = []
    d = []
    for line in lines:
        tokens = line.split()
        try:
            p.append(int(tokens[0]))
            a.append(int(tokens[1])/1000)
            b.append(int(tokens[2])/1000)
            c.append(int(tokens[3]))
            d.append(int(tokens[4]))
        except ValueError:
            #print "Values: ", tokens[0], tokens[1], tokens[2], tokens[3], tokens[5]
            print line

    p = np.array(p)
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)
    d = np.array(d)
    
    try:
        # Area graph
        area_a = np.trapz(a)
        for i in xrange(10, 200):
            new_area = np.trapz(a[:i])
            if new_area >= area_a/2.0:
                r_half_a= i
                break
        
        std_a = np.std(a)
        
        
        # Perimeter graph
        area_p = np.trapz(b)
        for i in xrange(10, 200):
            new_area = np.trapz(b[:i])
            if new_area >= area_p/2.0:
                r_half_p= i
                break
        
        std_p = np.std(b)
        
        
        # Mean breadth graph
        n = c.shape[0]
        mean = np.mean(c)
        std_mb  = np.std(c)
        variance = np.var(c)
        kurtosis = np.sum((c - mean)**4)/(n*variance**2) - 3.0
        skewness = np.sum((c - mean)**3)/(n*variance**(3/2))
        
        
        # Euler number graph
        prev = d[0]
        Si = 0
        Di = 0
        Pi = 0
        for i in xrange(1, 200):
            curr = d[i]
            diff =  curr - prev
            if diff > 0:
                Si += 1
            elif diff < 0:
                Di += 1
            else:
                Pi += 1
            prev = curr
            
        Si = float(Si)
        Di = float(Di)
        Pi = float(Pi)
        IS = Si/(Si+Di+Pi)
        
        pol   = np.poly1d(np.polyfit(p, d, 15))
        dpol  = pol(p)
        ddiff = d - dpol
        
        area_e = np.trapz(dpol)
        std_e = np.std(ddiff)
    
        print >> out_file, name, area_a, r_half_a, std_a, area_p, r_half_p, std_p, mean, std_mb, kurtosis, skewness, area_e, std_e, IS
    except:
        print >> out_file, name, " error"
    i += 1
    if i%20 == 0:
        out_file.flush()

out_file.close()
