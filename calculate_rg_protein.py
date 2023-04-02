from Bio import PDB
import numpy as np
import codecs

def calcum(A):#transfer the str in float mass
    if(A=='C'):
        m=12
    elif(A=='O'):
        m=16
    elif(A=='N'):
        m=14
    elif(A=='H'):
        m=1
    elif(A=='FE'):
        m=56
    elif(A=='P'):
        m=31
    elif(A=='S'):
        m=32
    elif(A=='CU'):
        m=64
    elif(A=='CA'):
        m=40
    else:
        m=0
    return m

def calcuRc(x,m,M,dim):#calculate the radio of the center of mass
    X=0
    for i in range(dim):
        X+=m[i]*x[i]/M
    return X

def calcuRg(x,y,z,m,M,dim):#calculate the Rg
    X=calcuRc(x,m,M,dim)
    Y=calcuRc(y,m,M,dim)
    Z=calcuRc(z,m,M,dim)
    Rg2=0
    for i in range(dim):
        Rg2+=m[i]*((x[i]-X)**2+(y[i]-Y)**2+(z[i]-Z)**2)
    return np.sqrt(Rg2/M)

parser=PDB.PDBParser()#make a subjection
io=PDB.PDBIO()#make an io

structure=parser.get_structure('1tub','1tub.pdb')#get structure, put in the pdb you like.

class Select():#select some special atom, maybe sometime feel useful
    def accept_model(self, model):
        if model.get_id()==0:#just select the model id =0, or name=1
            return True
        else:
            return False
    def accept_chain(self, chain):
        return True
    def accept_residue(self, residue):
        return True
    def accept_atom(self, atom):
        return True

io.set_structure(structure)        
io.save('1tubatom.pdb',Select())#save in a new file

A=codecs.open('1tubatom.pdb',mode='r',encoding='utf-8')# read the contents of new file
line=A.readline()
list1=[]
list2=[]

while line:
    a=line.split()
    if a[11:12]!=[]:
        b=a[6:7]
        c=a[7:8]
        d=a[8:9]
        e=a[11:12]
    elif a[6:7]!=[]:#this situation is just for some str link to each other, so you may modulate some profiles, passing somethings just like the "TER" in the midium of the pdb file.
        if len(a[9:10][0])<5:#for the 2:3 3:4 link together
            b=a[5:6]
            c=a[6:7]
            d=a[7:8]
            e=a[10:11]
        if len(a[9:10][0])>5:#for the 9:10,10:11 link together
            b=a[6:7]
            c=a[7:8]
            d=a[8:9]
            e=a[10:11]
    list1.append(b)
    list1.append(c)
    list1.append(d)
    list2.append(e)
    line=A.readline()
    
A.close()
"""print(list2[0][0]=='N')"""

dim=int(len(list1)/3-2)
'''print(dim)'''
atom_M=np.array([0.0]*dim)#initialize the arry of mass and coordinates
atom_x=np.array([0.0]*dim)
atom_y=np.array([0.0]*dim)
atom_z=np.array([0.0]*dim)

for i in range(dim):
    atom_x[i]=float(list1[i*3+0][0])
    atom_y[i]=float(list1[i*3+1][0])
    atom_z[i]=float(list1[i*3+2][0])
'''for i in range(dim//100+1):
    for j in range(100*i,np.minimum(100*(i+1),dim)):
        mm=calcum(list1[j*4+3][0])
        atom_M[j]=float(mm)
print(atom_M[230])'''
total_M=0
for j in range(dim):
    atom_M[j]=float(calcum(list2[j][0]))
    total_M+=atom_M[j]

rg=calcuRg(atom_x,atom_y,atom_y,atom_M,total_M,dim)#calculate the Rg

print(rg)