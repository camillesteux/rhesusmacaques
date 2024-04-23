#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 13:30:49 2022

@author: cxs6193
"""

import sys
from sys import argv
import os
from itertools import chain
import matplotlib.pyplot as plt


vcf_file = argv[1]
name_out = argv[2]

#vcf_file = "test.python.vcf"
#name_out = "test"

pos_list=[]
all_freq=[]
with open(vcf_file, "rt") as vcf:
    for line in vcf:
        fields = line.strip().split()
        if fields[0] == "#CHROM":
            nind = len(fields[9:])    
        elif "chr" in fields[0][0:3]:
            pos = "\t".join(str(e) for e in fields[0:2])
            pos_list.append(pos)
            gen = fields[9:]
            count = sum([int(i) for i in list(chain.from_iterable([g.split(":")[0].split("|") for g in gen])) if i != "."])
            all_freq.append(count)
            #print(count)


##Explicit allele count

output = open("output/"+name_out+".allelecount.txt", "a")
for i in range(len(pos_list)):
    output.write(str(pos_list[i])+"\t"+str(all_freq[i])+"\n")
output.close()



##SFS

sfs = {x:all_freq.count(x) for x in all_freq}

#Remove monomorphic sites in the pop
#if 0 in [x[0] for x in sorted(sfs.items())]:
    #sfs.pop(0)

#if 2*nind in [x[0] for x in sorted(sfs.items())]:
    #sfs.pop(2*nind)    

#Need to add the missing alleles counts
for i in range(0,(2*nind+1)):
    if i not in [x[0] for x in sorted(sfs.items())]:
        sfs[i]=0

freq_list = [x[0] for x in sorted(sfs.items())]
count_list = [x[1] for x in sorted(sfs.items())]

##OUTPUT

#txt file
output = open("output/"+name_out+".sfs.txt", "a")
for i in range(len(freq_list)):
    output.write(str(freq_list[i])+"\t"+str(count_list[i])+"\n")
output.close()

#plot
#plt.figure()
#plt.bar(freq_list, count_list)
#plt.title(name_out)
#plt.xlabel('Derived allele count')
#plt.ylabel('Number of sites')
#plt.savefig("output/"+name_out+".png")