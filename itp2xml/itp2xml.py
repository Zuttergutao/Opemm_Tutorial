import re
from datetime import datetime
from xml.etree.ElementTree import Element
from xml.etree.ElementTree import SubElement
from xml.etree.ElementTree import ElementTree
import math
import sys

# 将生成的itp文件转换为xml
with open(sys.argv[1]) as f:
    lines=f.readlines()

# elemnt为传进来的Elment类，参数indent用于缩进，newline用于换行   
def prettyXml(element, indent, newline, level = 0): 
    # 判断element是否有子元素
    if element:
        # 如果element的text没有内容      
        if element.text == None or element.text.isspace():     
            element.text = newline + indent * (level + 1)      
        else:    
            element.text = newline + indent * (level + 1) + element.text.strip() + newline + indent * (level + 1)    
    # 此处两行如果把注释去掉，Element的text也会另起一行 
    #else:     
        #element.text = newline + indent * (level + 1) + element.text.strip() + newline + indent * level    
    temp = list(element) # 将elemnt转成list    
    for subelement in temp:    
        # 如果不是list的最后一个元素，说明下一个行是同级别元素的起始，缩进应一致
        if temp.index(subelement) < (len(temp) - 1):     
            subelement.tail = newline + indent * (level + 1)    
        else:  # 如果是list的最后一个元素， 说明下一行是母元素的结束，缩进应该少一个    
            subelement.tail = newline + indent * level   
        # 对子元素进行递归操作 
        prettyXml(subelement, indent, newline, level = level + 1)     
        

ForceField = Element('ForceField')
info = SubElement(ForceField, 'Info')
DataGenerated=SubElement(info,"DataGenerated")
DataGenerated.text=str(datetime.today())
Source =SubElement(info,"Source")
Source.text=lines[0].strip(";")
Reference=SubElement(info,"Reference")
Reference.text="itp2xml (version 1.1) was developed by Casea; 2023.03.03"

namelists={}
nameindex=[]
count=0
for index,line in enumerate(lines):
    if line.strip().startswith("["):
        name=re.findall(re.compile(r"^\[([\sa-zA-Z]+)\]"),line)
        try:
            if line.split(";")[1].strip() == "impropers":
                name[0]="imdiherdrals"
            elif line.split(";")[1].strip() == "propers":
                name[0]="diherdrals"
        except:
            pass
        namelists[name[0].strip()]=count
        nameindex.append(index)
        count=count+1
namelists["empty"]=count
nameindex.append(-1)

# atoms
atoms=[]
for line in lines[nameindex[namelists["atoms"]]+1:nameindex[namelists["atoms"]+1]]:
    if line.startswith(";") or line.isspace():
        pass
    else:
        atoms.append(line.split())

# bonds
bonds=[]
for line in lines[nameindex[namelists["bonds"]]+1:nameindex[namelists["bonds"]+1]]:
    if line.startswith(";") or line.isspace():
        pass
    else:
        bonds.append(line.split())


# angles
angles=[]
for line in lines[nameindex[namelists["angles"]]+1:nameindex[namelists["angles"]+1]]:
    if line.startswith(";") or line.isspace():
        pass
    else:
        angles.append(line.split())

# diherdrals
diherdrals=[]
for line in lines[nameindex[namelists["diherdrals"]]+1:nameindex[namelists["diherdrals"]+1]]:
    if line.startswith(";") or line.isspace():
        pass
    else:
        diherdrals.append(line.split())

# pair
pairs=[]
for line in lines[nameindex[namelists["pairs"]]+1:nameindex[namelists["pairs"]+1]]:
    if line.startswith(";") or line.isspace():
        pass
    else:
        pairs.append(line.split())

# imdihedral
imdiherdrals=[]
for line in lines[nameindex[namelists["imdiherdrals"]]+1:nameindex[namelists["imdiherdrals"]+1]]:
    if line.startswith(";") or line.isspace():
        pass
    else:
        imdiherdrals.append(line.split())

# Atomtypes
atomtypes=[]
for line in lines[nameindex[namelists["atomtypes"]]+1:nameindex[namelists["atomtypes"]+1]]:
    if line.startswith(";") or line.isspace():
        pass
    else:
        atomtypes.append(line.split())

AtomTypes=SubElement(ForceField,"AtomTypes")

atomtypeslist={}
for i in range(len(atomtypes)):
    atomtypeslist[atomtypes[i][0]]=atomtypes[i][2]


atlist={}
for i in range(len(atoms)):
    atlist["atom_"+str(i)]= SubElement(AtomTypes,"Type")
    atlist["atom_"+str(i)].set("class",atoms[i][1])
    atlist["atom_"+str(i)].set("element",atoms[i][1][0].upper())
    atlist["atom_"+str(i)].set("mass",atoms[i][7])
    atlist["atom_"+str(i)].set("name","LIG-"+atoms[i][4])

epsilonlist={}
for i in range(len(atomtypes)):
    epsilonlist[atomtypes[i][0]]=atomtypes[i][5]

sigmalist={}
for i in range(len(atomtypes)):
    sigmalist[atomtypes[i][0]]=atomtypes[i][6]

# Residue
Residues=SubElement(ForceField,"Residues")
Residue=SubElement(Residues,"Residue")
Residue.set("name",atoms[0][3])
atomslist={}
for i in range(len(atoms)):
    atomslist["atom_"+str(i)]= SubElement(Residue,"Atom")
    atomslist["atom_"+str(i)].set("charge",atoms[i][6])
    atomslist["atom_"+str(i)].set("name",atoms[i][4])
    atomslist["atom_"+str(i)].set("type","LIG-"+atoms[i][4])

for i in range(len(bonds)):
    atomslist["bond_"+str(i)]= SubElement(Residue,"Bond")
    atomslist["bond_"+str(i)].set("atomName1",atoms[int(bonds[i][0])-1][4])
    atomslist["bond_"+str(i)].set("atomName2",atoms[int(bonds[i][1])-1][4])

# HarmonicBondForce
HarmonicBondForce=SubElement(ForceField,"HarmonicBondForce")
Bondlist={}
for i in range(len(bonds)):
    Bondlist["bond_"+str(i)]=SubElement(HarmonicBondForce,"Bond")
    Bondlist["bond_"+str(i)].set("k",str(eval(bonds[i][4])))
    Bondlist["bond_"+str(i)].set("length",str(eval(bonds[i][3])))
    Bondlist["bond_"+str(i)].set("type1","LIG-"+atoms[int(bonds[i][0])-1][4])
    Bondlist["bond_"+str(i)].set("type2","LIG-"+atoms[int(bonds[i][1])-1][4])

# HarmonicAngleForce
HarmonicAngleForce=SubElement(ForceField,"HarmonicAngleForce")
Anglelist={}
for i in range(len(angles)):
    Anglelist["angle_"+str(i)]=SubElement(HarmonicAngleForce,"Angle")
    Anglelist["angle_"+str(i)].set("angle",str(math.radians(float(angles[i][4]))))
    Anglelist["angle_"+str(i)].set("k",str(eval(angles[i][5])))
    Anglelist["angle_"+str(i)].set("type1","LIG-"+atoms[int(angles[i][0])-1][4])
    Anglelist["angle_"+str(i)].set("type2","LIG-"+atoms[int(angles[i][1])-1][4])
    Anglelist["angle_"+str(i)].set("type3","LIG-"+atoms[int(angles[i][2])-1][4])

# PeriodicTorsionForce
PeriodicTorsionForce=SubElement(ForceField,"PeriodicTorsionForce")
PeriodicTorsionForce.set("ordering","amber")

dihehrallist={}
for i in range(len(diherdrals)):
    dihehrallist["proper_"+str(i)]=SubElement(PeriodicTorsionForce,"Proper")
    dihehrallist["proper_"+str(i)].set("k1",diherdrals[i][6])
    dihehrallist["proper_"+str(i)].set("periodicity1",diherdrals[i][7])
    dihehrallist["proper_"+str(i)].set("phase1",str(math.radians(float(diherdrals[i][5]))))
    dihehrallist["proper_"+str(i)].set("type1","LIG-"+atoms[int(diherdrals[i][0])-1][4])
    dihehrallist["proper_"+str(i)].set("type2","LIG-"+atoms[int(diherdrals[i][1])-1][4])
    dihehrallist["proper_"+str(i)].set("type3","LIG-"+atoms[int(diherdrals[i][2])-1][4])
    dihehrallist["proper_"+str(i)].set("type4","LIG-"+atoms[int(diherdrals[i][3])-1][4])

for i in range(len(imdiherdrals)):
    dihehrallist["improper_"+str(i)]=SubElement(PeriodicTorsionForce,"Improper")
    dihehrallist["improper_"+str(i)].set("k1",imdiherdrals[i][6])
    dihehrallist["improper_"+str(i)].set("periodicity1",imdiherdrals[i][7])
    dihehrallist["improper_"+str(i)].set("phase1",str(math.radians(float(imdiherdrals[i][5]))))
    dihehrallist["improper_"+str(i)].set("type1","LIG-"+atoms[int(imdiherdrals[i][0])-1][4])
    dihehrallist["improper_"+str(i)].set("type2","LIG-"+atoms[int(imdiherdrals[i][1])-1][4])
    dihehrallist["improper_"+str(i)].set("type3","LIG-"+atoms[int(imdiherdrals[i][2])-1][4])
    dihehrallist["improper_"+str(i)].set("type4","LIG-"+atoms[int(imdiherdrals[i][3])-1][4])

NonbondedForce=SubElement(ForceField,"NonbondedForce")
NonbondedForce.set("coulomb14scale","0.8333333333333334")
NonbondedForce.set("lj14scale","0.5")
# UseAttributeFromResidue=SubElement(NonbondedForce,"UseAttributeFromResidue")
# UseAttributeFromResidue.set("name","charge")
atomtypeslist={}
for i in range(len(atoms)):
    atomtypeslist["Atom_"+str(i)]= SubElement(NonbondedForce,"Atom")
    atomtypeslist["Atom_"+str(i)].set("epsilon",str(eval(epsilonlist[atoms[i][1]])))
    atomtypeslist["Atom_"+str(i)].set("sigma",str(eval(sigmalist[atoms[i][1]])))
    atomtypeslist["Atom_"+str(i)].set("charge",atoms[i][6])
    atomtypeslist["Atom_"+str(i)].set("type","LIG-"+atoms[i][4])

tree = ElementTree(ForceField)
root = tree.getroot()                  #得到根元素，Element类    
prettyXml(root, '\t', '\n')            #执行美化方法    
# write out xml data
tree.write(sys.argv[2], encoding = 'utf-8')
print("itp file has been converted into xml file")