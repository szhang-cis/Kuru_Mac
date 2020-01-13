# System libraries and interaction with user
import sys, os
# Build a path for python to Florence
Path = os.path.dirname(os.path.realpath(__file__))
#Path = Path + '/2nd_GR'

Files = os.listdir(Path)

for filename in Files:
    os.rename(filename,filename.replace('_quantity_0_increment_3',''))
