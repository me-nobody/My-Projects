#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 15:34:12 2019

@author: teetul
"""

import re
def get_seqs():
    with open(r"/home/teetul/Documents/Deira_G4/deiraDNA.fna","r") as fh:
        data=fh.read()
    return data

def pattern_compile():
    lst=[]
    patterns=["AAA[ATGC]{1,7}AAA[ATGC]{1,7}AAA[ATGC]{1,7}AAA","GGG[ATGC]{1,7}GGG[ATGC]{1,7}GGG[ATGC]{1,7}GGG"]
    for pattern in patterns:
         regex=re.compile(pattern)
         
         lst.append(regex)
    return lst, patterns

def pattern_search():
    dna=get_seqs()
    pat_lst=pattern_compile()
    for pattern in pat_lst:
        print(pattern," length is ",len(re.findall(pattern,dna)))

pattern_search()
