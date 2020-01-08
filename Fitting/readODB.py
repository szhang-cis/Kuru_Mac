# ABAQUS/PYTHON POST PROCESSING SCRIPT
# Run using abaqus pyton / abaqus viewer -noGUI / abaqus cae -noGUI

# Set of some variables
work_dir='/home/local/EMSE2000/joan.laubrie/Documents/Fitting'

# Import of modules
from os import chdir
from odbAccess import *

# Working directory
chdir(work_dir)
# choose node to analisys
nodei=779
nodeo=775
#call ODB file
odb=openOdb('tube.odb')
fileCo=open('coordo.dat','w')
fileCi=open('coordi.dat','w')
fileT=open('thickn.dat','w')
fileS=open('stress.dat','w')
# Initial geometrie
Coordo=odb.rootAssembly.instances['PART-1-1'].nodes[nodeo].coordinates
Coordi=odb.rootAssembly.instances['PART-1-1'].nodes[nodei].coordinates
# Definition of arrays
Displo=[]
Displi=[]
stress=[]
inc=len(odb.steps['Step-1'].frames)
# Reading from ODB file
# First, taking first frame from every step in time
for k in range(inc):
    Displo.append(odb.steps['Step-1'].frames[k].fieldOutputs['U'].values[nodeo].data)
    Displi.append(odb.steps['Step-1'].frames[k].fieldOutputs['U'].values[nodei].data)
    stress.append(odb.steps['Step-1'].frames[k].fieldOutputs['S'].values[nodei].data)
# Output of some important variables
for k in range(inc):
    coopdisoX = Coordo[0] + Displo[k][0]
    coopdisoY = Coordo[1] + Displo[k][1]
    coopdisoZ = Coordo[2] + Displo[k][2]
    fileCo.write('%3d %f %f %f\n'%(k,coopdisoX,coopdisoY,coopdisoZ))
    coopdisiX = Coordi[0] + Displi[k][0]
    coopdisiY = Coordi[1] + Displi[k][1]
    coopdisiZ = Coordi[2] + Displi[k][2]
    fileCi.write('%3d %f %f %f\n'%(k,coopdisiX,coopdisiY,coopdisiZ))
    thickness = ((coopdisoX - coopdisiX)**2 + (coopdisoY - coopdisiY)**2)**0.5
    fileT.write('%3d %f\n'%(k,thickness))
    fileS.write('%3d %f %f %f %f %f %f\n'%(k,stress[k][0],stress[k][1],stress[k][2],stress[k][3],stress[k][4],stress[k][5]))

odb.close()
fileS.close()
fileCo.close()
fileCi.close()
fileT.close()
