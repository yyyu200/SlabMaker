#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" @author: yyyu200@163.com """

import numpy as np
import copy
import re

def ext_euclid(a, b):
     # find root of ax+by=gcd(a,b) =q
     if b == 0:
         return 1, 0, a
     else:
         x, y, q = ext_euclid(b, a % b) # q = gcd(a, b) = gcd(b, a%b)
         x, y = y, (x - (a // b) * y)
         return x, y, q

def ext_gcd(a, b):
    if b == 0:
        return 1, 0
    elif a % b == 0:
        return 0, 1
    else:
        x, y = ext_gcd(b, a % b)
        return y, x - y * (a // b)

def parse_lines(key, lines):
    res=None
    for l in lines:
        if not res:
            res=parse_str(key,l)

    return res

def parse_str(key, line): 
    # multiple key in a line must seperate by ',', not ';' or ' '
    findkey=re.search(key, line)
    if findkey:
        r1=line.split('!')[0].split(',')
        for s in r1:
            r2=re.search(key, s)
            if r2:
                return s.split('=')[1]
    else:
        return None

class CELL(object):
    close_thr=1.0e-4
    def __init__(self, fnam, fmt='POSCAR'):
        '''
        init from POSCAR, test only, use with caution
        '''
        if fmt=='POSCAR':
            fi=open(fnam)
            ll=fi.readlines() 
            self.system=ll[0]
            self.alat=float(ll[1])
            self.cell=np.zeros([3,3])
            for i in range(3):
                for j in range(3):
                    self.cell[i,j]=float(ll[2+i].split()[j])*self.alat
            self.ntyp=len(ll[5].split())
            self.typ_name=ll[5].split()
            self.typ_num=np.zeros([self.ntyp],dtype=np.int32)
            for i in range(self.ntyp):
                self.typ_num[i]=int(ll[6].split()[i])
            self.coordsystem=ll[7] # Direct only, no 'Selective Dynamics' line expected
            assert self.coordsystem[0]=='D' or self.coordsystem[0]=='d'

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
        elif fmt=='QE': # experimental
            fi=open(fnam)
            ll=fi.readlines()
            ibrav=parse_lines('ibrav',ll)
            self.nat=parse_lines('nat', ll)
            self.ntyp=parse_lines("ntyp",ll)

            # assert not None
            ibrav, self.nat, self.ntyp=int(ibrav), int(self.nat), int(self.ntyp)

            self.cell=np.zeros([3,3], dtype=np.float64)
            if ibrav==0:
                i=0
                for l in ll:
                    cpara=re.search("CELL_PARAMETERS",l)
                    i+=1
                    if cpara:
                        is_alat=re.search('alat',l)
                        is_angstrom=re.search('angstrom',l)
                        is_bohr=re.search('bohr',l)
                        for j in range(3):
                            self.cell[j,0]=np.float64(ll[j+i].split()[0])
                            self.cell[j,1]=np.float64(ll[j+i].split()[1])
                            self.cell[j,2]=np.float64(ll[j+i].split()[2])
                        break
                    else:
                        continue
                assert cpara

                if is_alat:
                    cd1=parse_lines("celldm\(1\)",ll)
                    A=parse_lines("A",ll) #TODO: only test
                    if cd1:
                        cd1=np.float64(cd1)
                        self.cell=self.cell*cd1*0.52917720859
                    elif A:
                        A=np.float64(A)
                        self.cell=self.cell*A
                    else:
                        raise Exception
                elif is_bohr:
                    self.cell=self.cell*0.52917720859
                elif is_angstrom:
                    self.cell=self.cell
                        
            elif ibrav==1: # sc
                cd1=parse_lines("celldm\(1\)",ll)
                cd1=np.float64(cd1)
                self.cell=cd1*np.eye(3, dtype=np.float64)*0.52917720859
            elif ibrav==2: # fcc
                cd1=parse_lines("celldm\(1\)",ll)
                cd1=np.float64(cd1)
                self.cell=cd1*0.5*np.array([[-1,0,1],[0,1,1],[-1,1,0]])*0.52917720859
            elif ibrav==3: # bcc
                cd1=parse_lines("celldm\(1\)",ll)
                cd1=np.float64(cd1)
                self.cell=cd1*0.5*np.array([[1,1,1],[-1,1,1],[-1,-1,1]])*0.52917720859
            elif ibrav==4: # hex
                cd1=parse_lines("celldm\(1\)",ll)
                cd1=np.float64(cd1)
                cd3=parse_lines("celldm\(3\)",ll)
                cd3=np.float64(cd3)
                self.cell=cd1*np.array([[1,0,0],[-1.0/2,np.sqrt(3)/2,0],[0,0,cd3]])*0.52917720859
            elif ibrav==6: # tetra P
                cd1=parse_lines("celldm\(1\)",ll)
                cd1=np.float64(cd1)
                cd3=parse_lines("celldm\(3\)",ll)
                cd3=np.float64(cd3)
                self.cell=cd1*np.array([[1,0,0],[0,1,0],[0,0,cd3]])*0.52917720859
            elif ibrav==8: # orth P
                cd1=parse_lines("celldm\(1\)",ll)
                cd1=np.float64(cd1)
                cd2=parse_lines("celldm\(2\)",ll)
                cd2=np.float64(cd2)
                cd3=parse_lines("celldm\(3\)",ll)
                cd3=np.float64(cd3)
                self.cell=cd1*np.array([[1,0,0],[0,cd2,0],[0,0,cd3]])*0.52917720859
            else:
                # other ibrav
                raise NotImplementedError
        else:
            # other code? maybe Abinit, castep
            raise NotImplementedError

    def __str__(self):
        return "cell\n"+self.cell.__str__() + "\natom positions:\n"+ self.atpos.__str__()

    @staticmethod
    def dist(a,b): # TODO: a and b are fractional coordinates
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
                # consider when tmp_i is small negative, e.g. -1e-17 , -1e-17+1=1
                if tmp_f < 0.0- self.close_thr:
                    self.atpos[i][j]=tmp_f+1.0
                elif np.fabs(tmp_f)<=self.close_thr or np.fabs(tmp_f-1.0)<=self.close_thr:
                    self.atpos[i][j]=0.0
                else:
                    self.atpos[i][j]=tmp_f
                assert self.atpos[i][j]<1.0 and self.atpos[i][j]>=0.0- self.close_thr

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

        fo.close()
    
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
    def is_inside(A, B, La,Ma,Na, Lb,Mb,Nb):
        ''' 
            cell A is inside cell B repeat by La~Lb, Ma~Mb, Na~Nb, A and B have the same origin
            the 8 corners of A in inside repeated B
        '''
        if Lb<=La or Mb<=Ma or Nb<=Na:
            return False

        P=np.mat(np.eye(3,dtype=np.float64))
        P[0,0]=Lb-La
        P[1,1]=Mb-Ma
        P[2,2]=Nb-Na
        Brepeat=copy.deepcopy(B)
        Brepeat.cell=np.array((np.mat(B.cell).T*P).T)
        Acell=np.zeros([3,3],dtype=np.float64)
        for i in range(3):
            Acell[i]=Brepeat.cart2direct(A.cell[i])
            #print("$",Acell[i])
    
        for i in range(8):
            corner=np.array([0.0 if np.fabs(La)<CELL.close_thr else -1.0/La, 
                             0.0 if np.fabs(Ma)<CELL.close_thr else -1.0/Ma, 
                             0.0 if np.fabs(Na)<CELL.close_thr else -1.0/Na])
            for j in range(3):
                corner+=(i>>j)%2*Acell[j]
            if (corner>1.0+CELL.close_thr).any() or ( corner<0.0-CELL.close_thr).any():
                return False
      
        return True 

    @staticmethod
    def cell2supercell(cell, P):
        '''
        use a cell to fill in new cell by translate of base vectors
        '''
        
        supercell=copy.deepcopy(cell)
        supercell.cell=((np.mat(cell.cell).T)*P).T
        supercell.cell=np.array(supercell.cell) # mat to array

        #print("\ncell.cell=\n", np.mat(cell.cell),"\nP=\n",P,"\nsupercell.cell=\n", supercell.cell)

        assert np.linalg.det(supercell.cell)>0

        Q=np.linalg.inv(P)
        for i in range(cell.nat):
            supercell.atpos[i]=np.array(Q*(np.mat(cell.atpos[i]).T)).flatten()
        
        # trans: fractional coordinates of cell vectors in supercell
        # atoms translate by trans is allowed in supercell
        trans=np.zeros([3,3], dtype=np.float64)
        for i in range(3):
            cell_i_frac=cell.cart2direct(cell.cell[i]) # one hot
            trans[i]=np.array(Q*(np.mat(cell_i_frac).T)).flatten()
       
        La,Ma,Na=0,0,0
        Lb,Mb,Nb=1,1,1
        while not CELL.is_inside(supercell,cell,La,Ma,Na,Lb,Mb,Nb):
            La-=1;Ma-=1;Na-=1
            Lb+=1;Mb+=1;Nb+=1

        while CELL.is_inside(supercell,cell,La+1,Ma,Na,Lb,Mb,Nb):
            La+=1
        while CELL.is_inside(supercell,cell,La,Ma+1,Na,Lb,Mb,Nb):
            Ma+=1
        while CELL.is_inside(supercell,cell,La,Ma,Na+1,Lb,Mb,Nb):
            Na+=1
        while CELL.is_inside(supercell,cell,La,Ma,Na,Lb-1,Mb,Nb):
            Lb-=1
        while CELL.is_inside(supercell,cell,La,Ma,Na,Lb,Mb-1,Nb):
            Mb-=1
        while CELL.is_inside(supercell,cell,La,Ma,Na,Lb,Mb,Nb-1):
            Nb-=1

        #print(La,Ma,Na, Lb,Mb,Nb)
         
        for n in range(supercell.nat):
            for i in range(La,Lb+1):
                for j in range(Ma,Mb+1):
                    for k in range(Na,Nb+1):
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

    def set_cell(self,cell):
        self.cell=np.mat(cell)

    def cell_redefine(self):
        '''Re-define unit cell to have the x-axis parallel with a surface vector and z perpendicular to the surface
         
        lattice: Atoms object
            The cell is rotated.
        '''
        from numpy.linalg import norm
        a1, a2, a3 = self.cell
        self.set_cell([a1, a2, np.cross(a1, a2) * np.dot(a3, np.cross(a1, a2)) /norm(np.cross(a1, a2))**2])
    
        a1, a2, a3 = np.array(self.cell)
        self.set_cell([[norm(a1), 0, 0], 
                          [np.dot(a1, a2) / norm(a1), 
                          np.sqrt(norm(a2)**2 - (np.dot(a1, a2) / norm(a1))**2), 0],
                          [0, 0, norm(a3)]] )

    def makeslab(self, miller_index, length=-1.0, layer=-1, method="point-group", origin_shift=0.0, vacuum=15.0):
        '''
        self is unit-cell
        '''
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
        elif k==0 and l==0: # slab not along z
            P[0,0]=layer
        elif h==0 and l==0:
            P[1,1]=layer
        else:
            p,q, _ =ext_euclid(k,l)
            c1,c2,c3=self.cell
    
            k1=np.dot(p*(k*c1-h*c2)+q*(l*c1-h*c3), l*c2-k*c3)
            k2=np.dot(l*(k*c1-h*c2)-k*(l*c1-h*c3), l*c2-k*c3)
            if np.fabs(k2) > self.close_thr:
                i=-int(round(k1/k2))
                p,q=p+i*l,q-i*k
            
            P[0]=(p * k + q * l, -p * h, -q * h)
            P[1]=np.array((0, l, -k)) // abs(np.gcd(l, k))
            P[2]=(h,k,l)
            P[2]*=layer
            P=P.T
      
        print("P = ",P) 
        slab=CELL.cell2supercell(self,P)
        print("slab cell \n", slab.cell)
        slab.cell_redefine()
        zmax=np.max(slab.atpos[:,2])
        zmin=np.min(slab.atpos[:,2])

        oldC=np.linalg.norm(slab.cell[2])
        newC=oldC*(zmax-zmin)+vacuum
        print(zmin,zmax, oldC, newC)
        slab.cell[2,2]=newC
        print("slab cell \n", slab.cell)
        for i in range(slab.nat):
            slab.atpos[i][2]=(vacuum/2.0+(slab.atpos[i][2]-zmin)*oldC)/newC

        print("slab cell \n", slab.cell)

        return slab

if __name__ == '__main__':
    c1=CELL("Al2O3.vasp")
    prim=CELL.unit2prim(c1,5)
    prim.print_poscar("rh.vasp")

    P=np.mat([[1,0,1],[-1,1,1],[0,-1,1]],dtype=np.float64)
    c2=CELL.cell2supercell(prim,P)
    c2.print_poscar("./test/c2.vasp")
    print(c2.get_volume(), prim.get_volume())

    P=np.mat([[2/3.0,-1/3.0,-1/3.0],[1/3.0,1/3.0,-2/3.0],[1/3.0,1/3.0,1/3.0]],dtype=np.float64)    
    prim=CELL.cell2supercell(c1,P)
    prim.print_poscar("./test/rh2.vasp")

    slab=c1.makeslab([1,1,0], layer=2)
    slab.print_poscar("./test/slab.vasp")

