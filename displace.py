# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 09:45:43 2017

@author: pyh
"""
def readvasprun():
     import numpy as np
     import re
     fid = open('vasprun.xml','rt')
     lc = np.mat(np.zeros((3,3)))
     fir = True
     while fir == True:
         line = fid.readline()
         if 'NSW' in line:
             line = re.findall('(?<=> ).*(?=<)',line)
             nstep = int(line[0])
         if '<atoms>' in line:
             line = re.findall('(?<=> ).*(?=<)',line)
             natom = int(line[0])
             fir=False        
         if 'basis' in line:
             for direct in range(3):
                 line = fid.readline()
                 line = line.replace('<v>','')
                 line = line.replace('</v>','')
                 lc[direct,:] = np.mat(line)                      
     posa=[]
     pos=np.mat(np.zeros((natom,3)))
         
     step=0
     atom=0        
     line = fid.readline()
     while line:
         if 'name="positions"' in line:
             for atom in range(natom):
                 line = fid.readline()
                 line = line.replace('<v>','')
                 line = line.replace('</v>','')
                 pos[atom,:] = np.mat(line)
             posa.append(pos.copy())
             step+=1
         line = fid.readline()
     fid.close()
     return nstep,lc,posa


def nearby(lc,pos1,pos2): 
    import numpy as np
    new_pos=pos2
    dp1=pos1*lc
    dp2=pos2*lc
    dis=(dp1-dp2)*np.transpose(dp1-dp2)
    old_dis=0
    n_pos=new_pos.copy()
    while np.all(old_dis != dis):
        old_dis=dis
        for i in range(0,3):
            for j in range(0,3):
                for k in range(0,3):
                    dp2=(new_pos+np.mat([i-1,j-1,k-1]))*lc
                    n_dis=(dp1-dp2)*np.transpose(dp1-dp2)
                    if np.all(n_dis < dis):
                        n_pos=new_pos+np.mat([i-1,j-1,k-1])
                        dis=n_dis
        new_pos=n_pos.copy()
    return new_pos
          

def displacement(posa,lc,arg1,step):
    import numpy as np
    import math
    arg1 -=1
    posa[0][arg1,:]=nearby(lc,posa[step][arg1,:],posa[0][arg1,:])
    dp1=posa[step][arg1,:]*lc
    dp2=posa[0][arg1,:]*lc
    dis=(dp1-dp2)*np.transpose(dp1-dp2)
    dis=math.sqrt(dis)
    return dis   

NSW,lc,posa = readvasprun()
dis=[0]*(NSW-1)
atom=8

fid=open('record','wt')

for i in range(1,len(posa)):
    dis[i-1]=displacement(posa,lc,atom,i)
#    print(dis[i-1])
    print("%d\t%.16f"%(i,dis[i-1]),file=fid)
