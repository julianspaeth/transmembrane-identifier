#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 10:06:54 2018

@author: Juli
"""

import pickle
import numpy as np
import collections
import random
from sklearn.svm import SVC


#TODO try different kernels, find optimal parameters with cv


def getAAlst(helix):
    c = collections.Counter(helix)
    AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'H', 'I', 'G', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    lst = []
    for aa in AA:
        lst.append(c[aa])
    return lst

if __name__ == "__main__":
    #generate data
    #this should be cleaned for U's
    trans_hel = pickle.load(open("pdbtmhel", "rb"))
    mol_hel = pickle.load(open("pdbhel", "rb"))
    
    #makes sure that i can interpret the accuracy
    trans_hel = trans_hel[:5000]       
    
    data = []
    
    for helix in trans_hel:
        temp = getAAlst(helix)
        temp.append(1)
        data.append(temp)
        
    for helix in mol_hel:
        temp = getAAlst(helix)
        temp.append(-1)
        data.append(temp)
        
    random.shuffle(data)       
    npdata = np.array(data)
    
    split = 2*len(npdata)//3
    
    train = npdata[:split]
    test = npdata[split:]
    
    X_train = train[:,0:20]
    Y_train = np.ravel(train[:,20:21])
    X_test = test[:,0:20]
    Y_test = np.ravel(test[:,20:21])
    
    
    #svm
    sv = SVC()
    
    sv.fit(X_train, Y_train)
    print(sv.score(X_test, Y_test))


