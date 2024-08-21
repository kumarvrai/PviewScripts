#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

fNameIn  = 'NekLineData_2D.csv'
fNameOut = 'NekAvgData_2D.csv'
fid = open(fNameIn,'r'); 
header_str = fid.readline();
header_str += 'nx=8 \n'
header_str += 'ny=1000'
fid.close()
word = "#time"
count = 0
with open(fNameIn, 'r') as f: 
    for line in f: 
        words = line.split("=") 
        for i in words: 
            if(i==word): 
                count=count+1
print("Time snapshots: ", count)
nt = count
#inputs
data = np.loadtxt(fNameIn, delimiter=',', skiprows=1);
nxy,nvars = np.shape(data);
print("Data shape: ", nxy, nvars, nt)
data = np.reshape(data,(nt,int(nxy/nt),nvars));
data = np.mean(data,axis=0);
#data = np.permute(data,(1,2,0));
print("Reduced Data shape: ", np.shape(data))
np.savetxt(fNameOut,data,delimiter = ',',header = header_str)
