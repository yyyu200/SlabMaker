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

def parse_lines(key, lines):
    res=None
    for l in lines:
        if not res:
            res=parse_str(key,l)

    return res

def dist2(a,b=[0,0,0]): # TODO: a and b are fractional coordinates
        return (a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2

def parse_str(key, line): 
    # multiple key in a line must seperated by ',', not ';' nor ' ' as in QE native code
    findkey=re.search(key, line)
    if findkey:
        r1=line.split('!')[0].split(',')
        for s in r1:
            r2=re.search(key, s)
            if r2:
                return s.split('=')[1]
    else:
        return None

def find_common_min(vecs, vecs_frac):
    #vecs: cartesian coords

    nat=len(vecs)

    # nv: number of inplane vec for each atom
    nv=np.zeros([nat],dtype=np.int64)
    nv[0]=vecs[0].shape[0]

    is_common0=np.ones([nv[0]],dtype=np.int64)
    for i in range(nat):
        #print("vecs_frac",i,vecs_frac[i])
        nv[i]=vecs[i].shape[0]
        for k in range(nv[0]):
            if is_common0[k]==0:
                continue
            for j in range(nv[i]):
                if dist2(vecs[i][j], vecs[0][k])<CELL.close_thr:
                    break
                elif j==nv[i]-1:
                    is_common0[k]=0
                else:
                    pass
                   
    com_vec=[]
    com_vec_frac=[]
    for i in range(nv[0]):
        if is_common0[i]==1:
            com_vec.append(vecs[0][i])
            com_vec_frac.append(vecs_frac[0][i])

    N=len(com_vec)
    #print(N, " com_vec_frac: ",com_vec_frac)
   
    # generate areas of vector pairs
    vol=np.zeros([N,N],dtype=np.float64)
    tmp=np.zeros([3],dtype=np.float64)
    for i in range(N):
        for j in range(N):
            tmp=np.cross(com_vec[i],com_vec[j])
            vol[i,j]=np.sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2])
    
    # choose from candidates, most small: minarea
    minarea=99999999.0
    for i in range(N):
        for j in range(N):
            if vol[i,j]<minarea-CELL.close_thr and vol[i,j]>CELL.close_thr:
                minarea=vol[i,j]
   
    # which pairs are minarea
    cij=[]
    for i in range(N):
        for j in range(N):
            if abs(vol[i,j]-minarea)<CELL.close_thr:
                cij.append([i,j])

    #  choose from most small, most closest to 90 degree
    most_close_to_rect=89.9
    mi=-1
    mj=-1
    for ij in cij:
            tmpang=fan(com_vec[ij[0]],com_vec[ij[1]])
            if abs(tmpang-90)< most_close_to_rect-0.01: # close thr for angle
                most_close_to_rect=abs(tmpang-90)
                mi=ij[0]
                mj=ij[1]

    #TODO:  exist some randomness in choice mi,mj, can be safely disregarded
    #print("inplane vectors: ", "\nu:\n",com_vec_frac[mi], "\nv:\n", com_vec_frac[mj]) 

    P=np.mat(np.eye(3,dtype=np.float64))
    
    P[0,0]=com_vec_frac[mi][0]
    P[1,0]=com_vec_frac[mi][1]
    P[2,0]=com_vec_frac[mi][2]
    P[0,1]=com_vec_frac[mj][0]
    P[1,1]=com_vec_frac[mj][1]
    P[2,1]=com_vec_frac[mj][2]

    # keep the right-hand system
    if np.cross(com_vec[mi],com_vec[mj])[2]<-CELL.close_thr:
        P[2,2]=-1
    else:
        P[2,2]=1

    print("P2 = ",P)
    return P

def fan(v1,v2):
    c=np.dot(v1,v2)/np.sqrt(v1.dot(v1))/np.sqrt(v2.dot(v2))
    return np.arccos(c)*180/np.pi

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
            self.alat=float(ll[1]) # Angstrom as unit
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
                        self.alat=cd1*0.52917720859
                    elif A:
                        A=np.float64(A)
                        self.cell=self.cell*A
                        self.alat=A
                    else:
                        raise Exception
                elif is_bohr:
                    self.cell=self.cell*0.52917720859
                    self.alat=np.sqrt(dist2(self.cell[0]))
                elif is_angstrom:
                    self.cell=self.cell
                    self.alat=np.sqrt(dist2(self.cell[0]))
                        
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

            if ibrav!=0:
                self.alat=cd1*0.52917720859
         
            # parse atomic type names
            i=0
            self.typ_name=[]
            self.typ_num=np.zeros([self.ntyp], dtype=np.int32)
            for l in ll:
                is_line_atspecies=re.search("ATOMIC_SPECIES",l)
                i+=1
                if is_line_atspecies:
                    for j in range(self.ntyp):
                        tmp=ll[j+i].split()[0]
                        self.typ_name.append(tmp)

            # parse atomic positions
            i=0
            self.atpos=np.zeros([self.nat,3], dtype=np.float64)
            self.attyp=np.zeros([self.nat], dtype=np.int32)
            for l in ll:
                is_line_atpos=re.search("ATOMIC_POSITIONS",l)
                i+=1
                if is_line_atpos: # do not support Capital, such as 'ALAT'
                    is_atpos_alat=re.search('alat',l)
                    is_atpos_angstrom=re.search('angstrom',l)
                    is_atpos_bohr=re.search('bohr',l)
                    is_atpos_crystal_sg=re.search('crystal_sg',l)
                    if not is_atpos_crystal_sg:
                        is_atpos_crystal=re.search('crystal',l)

                    assert is_atpos_alat or is_atpos_angstrom or is_atpos_bohr or is_atpos_crystal_sg or is_atpos_crystal

                    for j in range(self.nat):
                        tmp_element_name=ll[j+i].split()[0]
                        tmp_pos0=np.float64(ll[j+i].split()[1])
                        tmp_pos1=np.float64(ll[j+i].split()[2])
                        tmp_pos2=np.float64(ll[j+i].split()[3])
                    
                        self.attyp[j]=self.element_name2attyp(tmp_element_name)
                        self.typ_num[self.attyp[j]]+=1
                        if is_atpos_crystal:
                            self.atpos[j,0]=tmp_pos0
                            self.atpos[j,1]=tmp_pos1
                            self.atpos[j,2]=tmp_pos2
                        elif is_atpos_crystal_sg:
                            raise NotImplementedError
                        else:
                            tmp_c=np.array([tmp_pos0, tmp_pos1, tmp_pos2], dtype=np.float64)
                            if is_atpos_alat:
                                self.atpos[j]=self.cart2direct(tmp_c*self.alat)
                            elif is_atpos_angstrom:
                                self.atpos[j]=self.cart2direct(tmp_c)
                            elif is_atpos_bohr:
                                self.atpos[j]=self.cart2direct(tmp_c*0.52917720859)

                    break
                else:
                    continue

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
                dd=self.atpos[i]-self.atpos[j]
                if abs(dd[0])>self.close_thr:
                    continue
                if abs(dd[1])>self.close_thr:
                    continue
                if abs(dd[2])>self.close_thr:
                    continue
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

    def at_sort(self): # sort by element type
        tmp=self.attyp[:].argsort(kind='mergesort')
        self.attyp=self.attyp[tmp]
        self.atpos=self.atpos[tmp]

    def tidy_up(self): # tanslate to [0-1), sort by element
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

    def print_poscar(self,fnam):
        '''
        print to POSCAR
        '''
        fo=open(fnam,"w")
        fo.write("system: "+' '.join(self.typ_name)+"\n")
        fo.write("1.0\n")
        for i in range(3):
            for j in range(3):
                fo.write(" %15.10f" % self.cell[i,j])
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

    def element_name2attyp(self, element_name):
        for i in range(self.ntyp):
            if element_name == self.typ_name[i]:
                return i

    def direct2cart(self, a):
        b=np.matmul(self.cell.T, a)
        return b

    def cart2direct(self, a):
        b=np.matmul(self.get_rec(), a)
        return b

    def append(self, atpos, attyp, update_typ_num=False):# TODO: number of each typ (typ_num) not added
        atpos=atpos.reshape(1,3)
        attyp=np.array([attyp])
        self.atpos=np.concatenate((self.atpos, atpos), axis=0)
        self.attyp=np.concatenate((self.attyp, attyp), axis=0)
        self.nat+=1

    def unique_append(self, atpos, attyp):
        for i in range(self.nat):
            d=atpos-self.atpos[i]
            if abs(d[0])+abs(d[1])+abs(d[2])<3*self.close_thr:
                return
        attyp=np.array([attyp])
        atpos=atpos.reshape(1,3)
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
        
        # fraction coords of A basis vectors in Brepeat basis
        Acell=np.zeros([3,3],dtype=np.float64)
        for i in range(3):
            Acell[i]=Brepeat.cart2direct(A.cell[i])
    
        corner_o=np.array([0.0 if np.fabs(La)<CELL.close_thr else -La/(Lb-La), 
                           0.0 if np.fabs(Ma)<CELL.close_thr else -Ma/(Mb-Ma), 
                           0.0 if np.fabs(Na)<CELL.close_thr else -Na/(Nb-Na)])
        
        for i in range(8): # fraction coords of eight corners of A
            corner=corner_o.copy()
            for j in range(3):
                corner+=(i>>j)%2*Acell[j]
            if (corner>1.0+CELL.close_thr).any() or (corner<0.0-CELL.close_thr).any():
                return False
      
        #print(La,Ma,Na,Lb,Mb,Nb,"T",A.cell, Brepeat.cell, Acell)
        return True 

    @staticmethod
    def cell2supercell(cell, P, checkunique=True):
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
        Lb,Mb,Nb=0,0,0
        while not CELL.is_inside(supercell,cell,La,Ma,Na,Lb,Mb,Nb):
            La-=1;Ma-=1;Na-=1
            Lb+=1;Mb+=1;Nb+=1
        assert CELL.is_inside(supercell,cell,La,Ma,Na,Lb,Mb,Nb)

        while CELL.is_inside(supercell,cell,La+1,Ma,Na,Lb,Mb,Nb):
            La+=1
        while CELL.is_inside(supercell,cell,La,Ma,Na,Lb-1,Mb,Nb):
            Lb-=1
        while CELL.is_inside(supercell,cell,La,Ma+1,Na,Lb,Mb,Nb):
            Ma+=1
        while CELL.is_inside(supercell,cell,La,Ma,Na,Lb,Mb-1,Nb):
            Mb-=1
        while CELL.is_inside(supercell,cell,La,Ma,Na+1,Lb,Mb,Nb):
            Na+=1
        while CELL.is_inside(supercell,cell,La,Ma,Na,Lb,Mb,Nb-1):
            Nb-=1
        
        N=supercell.nat
        for n in range(N):
            for i in range(La,Lb):
                for j in range(Ma,Mb):
                    for k in range(Na,Nb):
                        #print(n,i,j,k)
                        supercell.append(supercell.atpos[n]+i*trans[0]+j*trans[1]+k*trans[2], supercell.attyp[n])
        #print("final:",La,Ma,Na, Lb,Mb,Nb, supercell.nat)

        supercell.tidy_up()
        if checkunique:
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
        primcell.cell=np.array(primcell.cell)
        assert np.linalg.det(primcell.cell)>=0
        primcell.nat=unitcell.nat
        for i in range(primcell.nat):
            primcell.atpos[i]=np.array(Q*(np.mat(unitcell.atpos[i]).T)).flatten()
        primcell.tidy_up()
        primcell.unique()
        # assert transform fully covered the primcell, no need for filling
        return primcell

    def set_cell(self,cell):
        self.cell=np.array(cell)

    def cell_redefine(self):
        '''Re-define slab cell to have the x-axis parallel with a surface vector 
           and z perpendicular to the surface. Only keep the inplane periodicity.
           First c is projected along z, fractional coords are changed in tis step.
           Then the cell is rotated.
        '''
        from numpy.linalg import norm
        a1, a2, a3 = self.cell

        carts=np.zeros([self.nat,3],dtype=np.float64)
        for i in range(self.nat):
            carts[i]=self.direct2cart(self.atpos[i])
        
        self.set_cell([a1, a2, np.cross(a1, a2) * np.dot(a3, np.cross(a1, a2)) /norm(np.cross(a1, a2))**2])
        for i in range(self.nat):
            self.atpos[i]=self.cart2direct(carts[i])
        self.tidy_up()

        a1, a2, a3 = np.array(self.cell)
        self.set_cell([[norm(a1), 0, 0], 
            [np.dot(a1, a2) / norm(a1), np.sqrt(norm(a2)**2 - (np.dot(a1, a2) / norm(a1))**2), 0],
            [0, 0, norm(a3)]] )

    def is_inlattice(self,i,by_s_tau):
        tmp=self.atpos[i]+by_s_tau
        for k in range(3):
            tf,ti=np.modf(tmp[k])
            if tf < 0.0- self.close_thr:
                tf=tf+1.0
            elif np.fabs(tf)<=self.close_thr or np.fabs(tf-1.0)<=self.close_thr:
                tf=0.0
            tmp[k]=tf

        for k in range(self.nat):
            if dist2(self.atpos[k],tmp)<=self.close_thr and self.attyp[k]==self.attyp[i]:
                return True

        return False

    def find_inplane(self,i,search_range=6):
        # norm vector, after redefine, n_ is [0,0,c]
        n_=np.cross(self.cell[0],self.cell[1])
       
        # find d in plane-equation:hx+ky+lz=d
        d=[]
        for k in range(self.nat):
            d.append(np.dot(self.atpos[k].reshape(3,),n_.reshape(3,)))
        
        # ikind, index of atoms inplane with atom i
        ikind=[]
        for k in range(self.nat):
            if k==i:
                continue
            if abs(d[k]-d[i])<self.close_thr and self.attyp[i]==self.attyp[k]:
                tau=self.atpos[k]-self.atpos[i]
                #print(i,k,tau)
                is_inplane_and_trans=True
                for s in range(1,search_range):
                    by_s_tau=s*tau
                    if not self.is_inlattice(i,by_s_tau):
                        is_inplane_and_trans=False
                        #print(i,by_s_tau,s)
                        break
                if is_inplane_and_trans:
                    ikind.append(k)

        iin=len(ikind) # number of inplane atoms with the i-th atom
        ir=1 # neighbor -ir, -ir+1,..., ir-1, ir
        inplane_range=2*ir+1 
        Nin=inplane_range**2
        re=np.zeros([iin*Nin,3],dtype=np.float64)

        lat_neighbor=np.zeros([Nin,3],dtype=np.float64)
        for il in range(Nin):
            rx=il//inplane_range-ir
            ry=il%inplane_range-ir
            lat_neighbor[il]=np.array([rx,ry,0],dtype=np.float64)
        
        for k in range(iin):
            tmp=self.atpos[ikind[k]]-self.atpos[i]
            for il in range(Nin):
                re[k*Nin+il]=tmp+lat_neighbor[il]

        re=np.concatenate((re,lat_neighbor),axis=0)
        
        re2=np.zeros([re.shape[0],re.shape[1]],dtype=np.float64)
        for k in range(re.shape[0]):
            re2[k]=self.direct2cart(re[k])
        
        #print("re",re.shape) 
        return re2, re

    @staticmethod
    def reduce_slab(slab,axis='c'):
        # find 2d minimum cell for the slab, assume c is surface norm
        vecs=[]
        vecs_frac=[]
        for i in range(slab.nat):
            tmp_cart, tmp_frac=slab.find_inplane(i)
            vecs.append(tmp_cart)
            vecs_frac.append(tmp_frac)

        P=find_common_min(vecs, vecs_frac)

        reduced=CELL.cell2supercell(slab,P)

        return reduced

    def makeslab(self, miller_index, length=-1.0, layer=-1, method="bf", origin_shift=0.0, vacuum=15.0):
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
            raise AssertionError
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
      
        if np.linalg.det(P)<0:
            P[0,2]*=-1
            P[1,2]*=-1
            P[2,2]*=-1
        print("P1 = ",P) 
        slab=CELL.cell2supercell(self,P)
        slab.cell_redefine()
        zmax=np.max(slab.atpos[:,2])
        zmin=np.min(slab.atpos[:,2])

        oldC=np.linalg.norm(slab.cell[2])
        newC=oldC*(zmax-zmin)+vacuum
        slab.cell[2,2]=newC
        for i in range(slab.nat):
            slab.atpos[i,2]=(vacuum/2.0+(slab.atpos[i,2]-zmin)*oldC)/newC

        reduced_slab=CELL.reduce_slab(slab)
        print("reduced slab cell \n", reduced_slab.cell)
        ang_B=fan(reduced_slab.cell[0],reduced_slab.cell[1])
        edg_a=np.linalg.norm(reduced_slab.cell[0])
        edg_c=np.linalg.norm(reduced_slab.cell[1])
        print("reduced slab No. of atoms: ", reduced_slab.nat)
        print("slab and vacuum length: ", newC-vacuum, vacuum, "Ang.")
        print("inplane edge and angle: ", edg_a, edg_c, ang_B," degree.")
        print("!reduced slab cell area: ", np.sin(ang_B/180*np.pi)*edg_a*edg_c, " Ang^2.")

        return reduced_slab

if __name__ == '__main__':
    import os

    if not os.path.exists("tmp"):
        os.makedirs("tmp")

    unit=CELL("Al2O3.vasp")

    slab=unit.makeslab([1,1,0], layer=1)
    slab.print_poscar("./tmp/slab.vasp")

