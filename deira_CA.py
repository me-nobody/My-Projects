#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 15:34:12 2019

@author: teetul
"""

import glob
import re

def read_dir():
    lst=[]
    for file in glob.glob(r"/home/teetul/Documents/Deira_G4/*"):
        lst.append(file)
    return lst      



def get_seqs():
    seqlist=[]
    filenames=read_dir()
    for file in filenames:
        fh=open(file,"r")    
        data=fh.read()
        seqlist.append(data)
    print (len(seqlist))
    
get_seqs()    

def pattern_compile():
    lst=[]
    patterns=["AAA[ATGC]{1,7}AAA[ATGC]{1,7}AAA[ATGC]{1,7}AAA","GGG[ATGC]{1,7}GGG[ATGC]{1,7}GGG[ATGC]{1,7}GGG"\
              "CCC[ATGC]{1,7}CCC[ATGC]{1,7}CCC[ATGC]{1,7}CCC","TTT[ATGC]{1,7}TTT[ATGC]{1,7}TTT[ATGC]{1,7}TTT"]
    for pattern in patterns:
         regex=re.compile(pattern)
         lst.append(regex)
    return patterns,lst
    

def pattern_search():
    dna=get_seqs()
    pat_lst,pat_relist=pattern_compile()
    count=0
    for item in pat_relist:        
        print(count," ",pat_lst[count]," length is ",len(re.findall(item,dna)))
        count=count+1

#pattern_search()
