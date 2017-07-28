# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 10:13:13 2017

@author: pyh
"""

def nearby(lp_old,lp_new,pos):
    dpos=pos*lp_old
    pos_new=dpos*lp_new.I
    dis=dpos*np.transpose(dpos)
    dis_old=0
    pos_n=pos_new
    while dis_old !=dis:
        dis_old=dis
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    dpos_new=(pos_new+np.mat([i-1,j-1,k-1]))*lp_new
                    dis_new=dpos_new*np.transpose(dpos_new)
                    if dis > dis_new:
                        dis = dis_new
                        pos_n=pos_new+np.mat([i-1,j-1,k-1])
        pos_new=pos_n       
    return pos_new

def mul(arr1,arr2):
    return [arr1[1]*arr2[2]-arr1[2]*arr2[1],
            arr1[2]*arr2[0]-arr1[0]*arr2[2],
            arr1[0]*arr2[1]-arr1[1]*arr2[0]]
    

def callp(lp_real_old,lp_real_new):
    lp_old=np.mat(np.zeros((3,3)))
    lp_new=np.mat(np.zeros((3,3)))
    v_old=float(np.mat(lp_real_old[0])*np.transpose(np.mat(mul(lp_real_old[1],lp_real_old[2]))))
    v_new=float(np.mat(lp_real_new[0])*np.transpose(np.mat(mul(lp_real_new[1],lp_real_new[2]))))

    lp_old[0]=np.mat(mul(lp_real_old[1],lp_real_old[2]))/v_old
    lp_old[1]=np.mat(mul(lp_real_old[2],lp_real_old[0]))/v_old
    lp_old[2]=np.mat(mul(lp_real_old[0],lp_real_old[1]))/v_old

    lp_new[0]=np.mat(mul(lp_real_new[1],lp_real_new[2]))/v_new
    lp_new[1]=np.mat(mul(lp_real_new[2],lp_real_new[0]))/v_new
    lp_new[2]=np.mat(mul(lp_real_new[0],lp_real_new[1]))/v_new
    
    return lp_old,lp_new
    
    
if __name__=='__main__':
    
    import numpy as np
    "请输入实空间折叠前、折叠后的晶格常数矩阵、折叠前的坐标"
#==============================================================================

    lp_real_old=[[1,0,0],
                 [-0.5,np.sqrt(3)/2,0],
                 [0,0,30]]
    lp_real_new=[[1,0,0],
                 [0,np.sqrt(3),0],
                 [0,0,30]]
    
    pos=np.mat([1/2,0,0])
#==============================================================================

    print(nearby(callp(lp_real_old,lp_real_new)[0],callp(lp_real_old,lp_real_new)[1],pos))
