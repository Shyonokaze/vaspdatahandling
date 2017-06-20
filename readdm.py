# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 21:37:17 2017

@author: pyh
"""

import numpy as np

fid = open('dm','rt')

nm = [[3,3],[9,10]] #原子序列+原子质量
layer=[[1,3,6],[2,4,5]] #分层

N=3*sum(nm[0])
mass=[None]*N
FCX = [[0 for i in range(len(layer))] for i in range(len(layer)) ]
FCY = [[0 for i in range(len(layer))] for i in range(len(layer)) ]
FCZ = [[0 for i in range(len(layer))] for i in range(len(layer)) ]

count = 0
for i in range(len(nm[0])):
    for j in range(3*nm[0][i]):
        mass[j+count]=[nm[1][i]]*N
    count += 3*nm[0][i]

mass = np.sqrt(np.multiply(np.mat(mass),np.transpose(np.mat(mass))))

fc = np.mat(np.zeros((N,N)))

for i in range(N):
    fc[i,:]=np.mat(fid.readline())

D2 = np.multiply(fc,mass)


for i in range(len(layer)):
    for j in range(len(layer)):
        add = 0
        for k in layer[i]:
            for l in layer[j]:
                FCX[i][j] = FCX[i][j] + D2[k*3-3,l*3-3]
                FCY[i][j] = FCY[i][j] + D2[k*3-2,l*3-2]
                FCZ[i][j] = FCZ[i][j] + D2[k*3-1,l*3-1]
