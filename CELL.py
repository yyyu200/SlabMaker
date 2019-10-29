#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" @author: yyyu200@163.com """

import numpy as np
import copy

def ext_euclid(a, b):
     # find root of ax+by=gcd(a,b) =q
     if b == 0:
         return 1, 0, a
     else:
         x, y, q = ext_euclid(b, a % b) # q = gcd(a, b) = gcd(b, a%b)
         x, y = y, (x - (a // b) * y)
         return x, y, q

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
            self.attyp=np.zeros([self.nat],dtype=np.int32)
            k=0
            for i in range(self.ntyp):
                for j in range(self.typ_num[i]):
                    self.attyp[k]=i
                    k+=1

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

    @staticmethod
    def dist(a,b):
        return np.linalg.norm(a-b)

    def unique(self):
        n=self.nat
        i_kind=np.zeros([n],dtype=np.int32)
        for i in range(n):
            i_kind[i]=i

        for i in range(n):
            for j in range(i+1,n):
                if CELL.dist(self.atpos[i],self.atpos[j])<self.close_thr:
                    assert(self.attyp[i]==self.attyp[j])
                    i_kind[j]=i_kind[i]
        uniq=[]
        typ_of_uniq=[]
        for i in np.unique(i_kind):
            uniq.append(self.atpos[i])
            typ_of_uniq.append(self.attyp[i])
        
        #print(i_kind)
        self.nat=len(uniq)
        self.atpos=np.array(uniq)
        self.attyp=np.array(typ_of_uniq)

        
        for i in range(self.ntyp):
            self.typ_num[i]=0
            for j in range(self.nat):
                if self.attyp[j]==i:
                    self.typ_num[i]+=1

    def tidy_up(self): # tanslate to [0-1), sort by z axis
        for i in range(self.nat):
            for j in range(3):
                tmp_f, tmp_i=np.modf(self.atpos[i][j]) # the fractional part , the interger part
                if tmp_f < 0.0:
                    self.atpos[i][j]=tmp_f+1.0
                else:
                    self.atpos[i][j]=tmp_f

        self.at_sort()

    def at_sort(self): # sort by element type
        tmp=self.attyp[:].argsort(kind='mergesort')
        self.attyp=self.attyp[tmp]
        self.atpos=self.atpos[tmp]
        #tmp=self.atpos[:,2].argsort(kind='mergesort')


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
                fo.write(" %.12f" % ( self.atpos[i][j]))
            fo.write("\n")
    
    def findfour():
        pass    

    def findprim():
        findfour()
        pass

    def get_volume(self):
        self.volume=np.dot(np.cross(self.cell[0], self.cell[1]),self.cell[2])
        assert self.volume>0.0
        return self.volume

    def get_rec(self):
        '''
            get reciprocal lattice, crystallographer's definition, without factor of 2 \pi
        '''
        self.rec=np.zeros([3,3],dtype=np.float64)
        self.volume=self.get_volume()
        for i in range(3):
            self.rec[i]=np.cross(self.cell[(i+1)%3], self.cell[(i+2)%3])/self.volume;
         
        return self.rec

    def cart2direct(self, a):
        b=np.matmul(self.get_rec(), a)
        return b

    def append(self, atpos, attyp):
        atpos=atpos.reshape(1,3)
        attyp=np.array([attyp])
        self.atpos=np.concatenate((self.atpos, atpos), axis=0)
        self.attyp=np.concatenate((self.attyp, attyp), axis=0)
        self.nat+=1

    @staticmethod
    def cell2supercell(cell, P):
        '''
        use a cell to fill in new cell by translate of base vectors
        '''
        
        supercell=copy.deepcopy(cell)
        supercell.cell=(np.mat(cell.cell).T*P).T
        supercell.cell=np.array(supercell.cell) # mat to array

        assert np.linalg.det(supercell.cell)>=0

        Q=np.linalg.inv(P)
        for i in range(cell.nat):
            supercell.atpos[i]=np.array(Q*(np.mat(cell.atpos[i]).T)).flatten()

        trans=np.zeros([3,3], dtype=np.float64)
        for i in range(3):
            cell_i_frac=cell.cart2direct(cell.cell[i]) # one hot
            trans[i]=np.array(Q*(np.mat(cell_i_frac).T)).flatten()
        
        L,M,N=-1,-1,-1
        i=1
        while L<0 or M<0 or N<0:
            if ((trans[0]*i)[:]>=1.0).any() and L<0:
                L=i+1
            if ((trans[1]*i)[:]>=1.0).any() and M<0:
                M=i+1
            if ((trans[2]*i)[:]>=1.0).any() and N<0:
                N=i+1
            if i>1000:
                break
            i+=1
        #print(L,M,N)
         
        for n in range(supercell.nat):
            for i in range(L):
                for j in range(M):
                    for k in range(N):
                        supercell.append(supercell.atpos[n]+i*trans[0]+j*trans[1]+k*trans[2], supercell.attyp[n])

        supercell.tidy_up()
        supercell.unique()

        return supercell

    @staticmethod
    def unit2prim(unitcell, ibrav):
        primcell=copy.deepcopy(unitcell)

        if ibrav in [2,10]: # face center
            P=np.mat([[0.0, 0.5, 0.5],[0.5, 0.0, 0.5],[0.5, 0.5, 0.0]],dtype=np.float64)
            Q=np.mat([[-1,1,1],[1,-1,1],[1,1,-1]],dtype=np.float64)
        elif ibrav in [3,7,11]: # body center
            P=np.mat([[-0.5, 0.5, 0.5],[0.5, -0.5, 0.5],[0.5, 0.5, -0.5]],dtype=np.float64)
            Q=np.mat([[0,1,1],[1,0,1],[1,1,0]],dtype=np.float64)
        elif ibrav==5: # rhombohedral
            P=np.mat([[2/3.0,-1/3.0,-1/3.0],[1/3.0,1/3.0,-2/3.0],[1/3.0,1/3.0,1/3.0]],dtype=np.float64)
            Q=np.mat([[1,0,1],[-1,1,1],[0,-1,1]],dtype=np.float64)
        elif ibrav in [9,13]: # base center
            P=np.mat([[0.5,-0.5,0],[0.5,0.5,0],[0,0,1.0]],dtype=np.float64)
            Q=np.mat([[1,1,0],[-1,1,0],[0,0,1]],dtype=np.float64)
        else:
            print("unit cell is primitive for ibrav: ", ibrav)
            P=np.mat([[1,0,0],[0,1,0],[0,0,1]], dtype=np.float64)
            Q=np.mat([[1,0,0],[0,1,0],[0,0,1]], dtype=np.float64)
        
        primcell.cell=(np.mat(unitcell.cell).T*P).T
        primcell.cell=np.array(primcell.cell)
        assert np.linalg.det(primcell.cell)>=0
        primcell.nat=unitcell.nat
        for i in range(primcell.nat):
            primcell.atpos[i]=np.array(Q*(np.mat(unitcell.atpos[i]).T)).flatten()
        primcell.tidy_up()
        primcell.unique()
        # assert 变换之后填满了原胞，不需要继续填充

        return primcell

    def makeslab(self, miller_index, length=-1.0, layer=-1, origin_shift=0.0, vacuum=15.0):
        P=np.mat(np.eye(3, dtype=np.float64))
        h,k,l=miller_index
        e=np.gcd(h, np.gcd(k,l))
        if e>1:
            h=h//e
            k=k//e
            l=l//e
            print("find miller index ( ",h,k,l," ) instead.")

        if h==0 and k==0 and l==0:
            print("miller_index cannot be 0 0 0!")
            raise Exception
        elif h==0 and k==0:
            P[2,2]=layer
        elif k==0 and l==0:
            P[0,0]=layer
        elif h==0 and l==0:
            P[1,1]=layer
        else:
            p,q,_=ext_euclid(k,l)
            P[0]=[p*k+q*l,-p*h,-q*h]
            P[1]=[0, l, -k]
            P[2]=np.cross(P[0], P[1])
        
        print(P)    
        slab=CELL.cell2supercell(self,P) #TODO  
        
        return slab

if __name__ == '__main__':
    c1=CELL("Al2O3.vasp")
    prim=CELL.unit2prim(c1,5)
    prim.print_poscar("rh.vasp")

    P=np.mat([[2,0,0],[0,2,0],[0,0,1]], dtype=np.float64)
    c2=CELL.cell2supercell(c1,P)
    c2.print_poscar("./test/c2.vasp")

    slab=c1.makeslab([1,1,1], layer=3)
    slab.print_poscar("./test/slab.vasp")

