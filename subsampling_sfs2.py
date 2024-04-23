#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 15:18:32 2022

@author: cxs6193
"""

"""
Need an output directory
python3 subsampling_sfs.py sfs.txt H name_out

"""

import sys
from sys import argv
import math
import numpy as np
import matplotlib.pyplot as plt
import copy as cp

#argv=["","output/test.sfs.txt","10","Subsampling_sfs.test"]

sfs_file = argv[1]
H = argv[2]
name_out = argv[3]

#Reading the sfs file and putting the values into a matrix
Sn=np.empty((0, 2))
with open(sfs_file, "rt") as sfs:
    for line in sfs:
        fields = line.strip().split()
        Sn = np.vstack([Sn,[int(k) for k in fields]])
        
def binom(k,n):
    if k==0 or k==n:
        return 1
    elif k>n:
        return 0
    else:
        return int((math.factorial(n)/(math.factorial(k)*math.factorial(n-k))))

#S if the originial SFS of sample size n, H is the new sample size (H<n)
#S is a matrix first column is the count of derived allele and the second colunm is the number of alleles exhibiting the corresponding nb of derived alleles in the pop
def subsamp(Sn,H):
    H=int(H)
    Sh = np.zeros(((H+1),2))
    #print(Sh)
    n=len(Sn)-1
    print(n)
    for j in range(0,H+1):
        print('j=',j)
        for i in range(j,n+1):
            print('i=',i)
            sn=Sn[(i),1]
            x=binom(j,i)*binom((H-j),(n-i))/binom(H,n)
            #print(x)
            #print(sn)
            Sh[(j),:]=[j,float(Sh[(j),1]+sn*x)]
    return(Sh)


Sh = subsamp(Sn,H)

##OUTPUT

#txt file
output = open("output/"+name_out+".subsampled.txt", "a")
for i in range(len(Sh)):
    output.write(str(Sh[i][0])+"\t"+str(Sh[i][1])+"\n")
output.close()

#plot 
#plt.figure()
#plt.bar(Sh[:,0], Sh[:,1])
#plt.title(name_out)
#plt.xlabel('Derived allele count')
#plt.ylabel('Number of sites')
#plt.savefig("output/"+name_out+".png")
