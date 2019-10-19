#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" @author: yyyu200@163.com """

import numpy as np
import copy

def dist(a,b):
    return np.linalg.norm(a-b)


class CELL(object):
    close_thr=1.0e-4
    def __init__(self, fnam, fmt='POSCAR'):
        '''
        init from POSCAR
        '''
        if fmt=='POSCAR':
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
            self.coordsystem=ll[7] # Direct only, no 'Selective Dynamics' line expected
            self.nat=int(sum(self.typ_num))
            self.atpos=np.zeros([self.nat,3],dtype=np.float64)
            for i in range(self.nat):
                for j in range(3):
                    self.atpos[i,j]=float(ll[8+i].split()[j])
            
            fi.close() 
        elif fmt=='QE':
            pass
        else:
            raise NotImplementedError

    def __str__(self):
        return "cell\n"+self.cell.__str__() + "\natom positions:\n"+ self.atpos.__str__()

    def mbc2prim(self,newCELL):
        P=np.mat([[0.5,-0.5,0],[0.5,0.5,0],[0,0,1.0]],dtype=np.float64)
        Q=np.mat([[1,1,0],[-1,1,0],[0,0,1]],dtype=np.float64)
        newCELL.cell=(np.mat(self.cell).T*P).T
        newCELL.nat=self.nat
        for i in range(newCELL.nat):
            newCELL.atpos[i]=np.array(Q*(np.mat(self.atpos[i]).T)).flatten()
        #TODO: 去掉重复的原子
        newCELL.tidy_up()

    @staticmethod
    def hex2rh(hexcell):
        rhcell=copy.deepcopy(hexcell)
        P=np.mat([[2/3.0,-1/3.0,-1/3.0],[1/3.0,1/3.0,-2/3.0],[1/3.0,1/3.0,1/3.0]],dtype=np.float64)
        Q=np.mat([[1,0,1],[-1,1,1],[0,-1,1]],dtype=np.float64)
      
        rhcell.cell=(np.mat(hexcell.cell).T*P).T
        rhcell.nat=hexcell.nat
        for i in range(rhcell.nat):
            rhcell.atpos[i]=np.array(Q*(np.mat(hexcell.atpos[i]).T)).flatten()
        rhcell.tidy_up()
        #TODO: 去掉重复的原子
        rhcell.merge_duplicate()
 
        return rhcell

    def merge_duplicate(self):
        pass

    def tidy_up(self): # tanslate to [0-1), sort by z axis
        for i in range(self.nat):
            for j in range(3):
                tmp_f, tmp_i=np.modf(self.atpos[i][j]) # the fractional part , the interger part
                if tmp_f < 0.0:
                    self.atpos[i][j]=tmp_f+1.0
                else:
                    self.atpos[i][j]=tmp_f
        self.atpos=self.atpos[self.atpos[:,2].argsort(kind='mergesort')]

    def print_poscar(self,fnam):
        '''
        print to POSCAR
        '''
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
                fo.write(" %.12f" % (np.modf(self.atpos[i][j])[0]))
            fo.write("\n")
    
    def findfour():
        pass    

    def findprim():
        findfour()
        pass

    def vasp2pw():
        pass


if __name__ == '__main__':
    #c1=CELL("mC.vasp")
    #prim=copy.deepcopy(c1)
    #c1.mbc2prim(prim)
    #prim.print_poscar("prim.vasp")

    hexcell=CELL("Al2O3.vasp")
    rhcell=CELL.hex2rh(hexcell)
    rhcell.print_poscar("rh.vasp")
    print(rhcell)
