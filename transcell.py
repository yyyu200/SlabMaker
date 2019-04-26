# -*- coding: utf-8 -*-
""" @author: yyyu200@163.com """

import numpy as np


class CELL(object):
    def __init__(self, fnam):
        fi=open(fnam)
        ll=fi.readlines() 
        self.system=ll[0]
        self.alat=float(ll[1])
        self.cell=np.zeros([3,3])
        for i in range(3):
            for j in range(3):
                self.cell[i,j]=float(ll[2+i].split()[j])
        self.ntyp=len(ll[5].split())
        self.typ_name=ll[5].split()
        self.typ_num=np.zeros([self.ntyp],dtype=np.int32)
        for i in range(self.ntyp):
            self.typ_num[i]=int(ll[6].split()[i])
        self.coordsystem=ll[7]
        self.nat=int(sum(self.typ_num))
        self.atpos=np.zeros([self.nat,3])
        for i in range(self.nat):
            for j in range(3):
                self.atpos[i,j]=float(ll[8+i].split()[j])
        
        fi.close() 

    def mbc2prim(self,newCELL):
        P=np.mat([[0.5,-0.5,0],[0.5,0.5,0],[0,0,1.0]],dtype=np.float64)
        Q=np.mat([[1,1,0],[-1,1,0],[0,0,1]],dtype=np.float64)
        newCELL.cell=(np.mat(self.cell).T*P).T
        newCELL.nat=self.nat
        for i in range(newCELL.nat):
            newCELL.atpos[i]=np.array(Q*(np.mat(self.atpos[i]).T)).flatten()
        #TODO: 去掉重复的原子 

    def printcell(self,fnam):
        fo=open(fnam,"w")
        fo.write("system: "+' '.join(self.typ_name)+"\n")
        fo.write("1.0\n")
        for i in range(3):
            for j in range(3):
                fo.write(" %f" % self.cell[i,j])
            fo.write("\n")
        for i in range(self.ntyp):
            fo.write(self.typ_name[i]+" ")
        fo.write("\n")
        for i in range(self.ntyp):
            fo.write(str(self.typ_num[i])+" ")
        fo.write("\n")
        fo.write("Direct\n")
        for i in range(self.nat):
            for j in range(3):
                #dig=np.modf(self.atpos[i][j])[0]
                #if dig < 0:
                #    dig=dig+1.0
                fo.write(" %f" % (np.modf(self.atpos[i][j])[0]))
            fo.write("\n")
    
    def findfour():
        pass    

    def findprim():
        findfour()
        pass

    def vasp2pw():
        pass

c1=CELL("mC.vasp")
import copy
prim=copy.deepcopy(c1)
c1.mbc2prim(prim)

prim.printcell("prim.vasp")
