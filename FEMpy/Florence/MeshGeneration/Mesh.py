#from __future__ import print_function, division
import os, sys, warnings, platform
#from time import time
import numpy as np
#if "PyPy" not in platform.python_implementation():
#    from scipy.io import loadmat, savemat
#from Florence.Tensor import makezero, itemfreq, unique2d, in2d
#from Florence.Utils import insensitive
#from .vtk_writer import write_vtu
#try:
#    import meshpy.triangle as triangle
#    has_meshpy = True
#except ImportError:
#    has_meshpy = False
#from .HigherOrderMeshing import *
from .NodeArrangement import *
#from .GeometricPath import *
#from warnings import warn
#from copy import deepcopy



"""
Mesh class providing most of the pre-processing functionalities of the Core module

Roman Poya - 13/06/2015
"""



class Mesh(object):

    """Mesh class provides the following functionalities:
    1. Generating higher order meshes based on a linear mesh, for tris, tets, quads and hexes
    2. Generating linear tri and tet meshes based on meshpy back-end
    3. Generating linear tri meshes based on distmesh back-end
    4. Finding bounary edges and faces for tris and tets, in case they are not provided by the mesh generator
    5. Reading Salome meshes in binary (.dat/.txt/etc) format
    6. Reading gmsh files .msh
    7. Checking for node numbering order of elements and fixing it if desired
    8. Writing meshes to unstructured vtk file format (.vtu) in xml and binary formats,
        including high order elements
    """

    def __init__(self, element_type=None):
        super(Mesh, self).__init__()
        # self.faces and self.edges ARE BOUNDARY FACES
        # AND BOUNDARY EDGES, RESPECTIVELY

        self.degree = None
        self.ndim = None
        self.edim = None
        self.nelem = None
        self.nnode = None

        self.elements = None
        self.points = None
        self.corners = None
        self.edges = None
        self.faces = None
        self.element_type = element_type

        self.face_to_element = None
        self.edge_to_element = None
        self.boundary_edge_to_element = None
        self.boundary_face_to_element = None
        self.all_faces = None
        self.all_edges = None
        self.interior_faces = None
        self.interior_edges = None
        # TYPE OF BOUNDARY FACES/EDGES
        self.boundary_element_type = None

        # FOR GEOMETRICAL CURVES/SURFACES
        self.edge_to_curve = None
        self.face_to_surface = None


        self.spatial_dimension = None
        self.reader_type = None
        self.reader_type_format = None
        self.reader_type_version = None
        self.writer_type = None

        self.filename = None
        # self.has_meshpy = has_meshpy

    def GetFaces(self):
        assert self.element_type is not None
        if self.element_type == "tet":
            self.GetFacesTet()
        elif self.element_type == "hex":
            self.GetFacesHex()
        elif self.element_type=="tri" or self.element_type=="quad":
            raise ValueError("2D mesh does not have faces")
        else:
            raise ValueError('Type of element not understood')
        return self.all_faces

    def GetFacesHex(self):
        """Find all faces (surfaces) in the hexahedral mesh (boundary & interior).
            Sets all_faces property and returns it

        returns:

            arr:            numpy ndarray of all faces

        """

        # DETERMINE DEGREE
        p = self.InferPolynomialDegree()

        # DO NOT COMPUTE IF ALREADY COMPUTED
        if isinstance(self.all_faces,np.ndarray):
            if self.all_faces.shape[0] > 1:
                # IF LINEAR VERSION IS COMPUTED, DO COMPUTE HIGHER VERSION
                if self.all_faces.shape[1] == 4 and p > 1:
                    pass
                else:
                    return self.all_faces

        node_arranger = NodeArrangementHex(p-1)[0]
        fsize = int((p+1)**3)

        # GET ALL FACES FROM THE ELEMENT CONNECTIVITY
        faces = np.concatenate((np.concatenate((
                np.concatenate((np.concatenate((np.concatenate((self.elements[:,node_arranger[0,:]],
                self.elements[:,node_arranger[1,:]]),axis=0),self.elements[:,node_arranger[2,:]]),axis=0),
                self.elements[:,node_arranger[3,:]]),axis=0),self.elements[:,node_arranger[4,:]]),axis=0),
                self.elements[:,node_arranger[5,:]]),axis=0).astype(np.int64)

        # REMOVE DUPLICATES
        self.all_faces, idx = unique2d(faces,consider_sort=True,order=False,return_index=True)

        face_to_element = np.zeros((self.all_faces.shape[0],2),np.int64)
        face_to_element[:,0] =  idx % self.elements.shape[0]
        face_to_element[:,1] =  idx // self.elements.shape[0]

        self.face_to_element = face_to_element

        return self.all_faces

    def Read(self, filename=None, element_type="tri", reader_type=None, reader_type_format=None,
        reader_type_version=None, order=0, read_surface_info=False, **kwargs):
        """Convenience mesh reader method to dispatch call to subsequent apporpriate methods"""

        if not isinstance(filename,str):
            raise ValueError("filename must be a string")
            return
        if reader_type is not None:
            if not isinstance(filename,str):
                raise ValueError("filename must be a string")
                return

        if reader_type is None:
            if filename.split('.')[-1] == "msh":
                reader_type = "gmsh"
            elif filename.split('.')[-1] == "obj":
                reader_type = "obj"
            elif filename.split('.')[-1] == "unv":
                reader_type = "unv"
            elif filename.split('.')[-1] == "fro":
                reader_type = "fro"
            elif filename.split('.')[-1] == "dat":
                for key in kwargs.keys():
                    inkey = insensitive(key)
                    if "connectivity" in inkey and "delimiter" not in inkey:
                        reader_type = "read_separate"
                        break
            if reader_type is None:
                raise ValueError("Mesh file format was not undertood. Please specify it using reader_type keyword")


        self.filename = filename
        self.reader_type = reader_type
        self.reader_type_format = reader_type_format
        self.reader_type_version = reader_type_version

        if self.reader_type is 'salome':
            self.ReadSalome(filename, element_type=element_type, read_surface_info=read_surface_info)
        elif reader_type is 'GID':
            self.ReadGIDMesh(filename, element_type, order)
        elif self.reader_type is 'gmsh':
            self.ReadGmsh(filename, element_type=element_type, read_surface_info=read_surface_info)
        elif self.reader_type is 'obj':
            self.ReadOBJ(filename, element_type=element_type, read_surface_info=read_surface_info)
        elif self.reader_type is 'fenics':
            self.ReadFenics(filename, element_type)
        elif self.reader_type is 'vtu':
            self.ReadVTK(filename)
        elif self.reader_type is 'unv':
            self.ReadUNV(filename, element_type)
        elif self.reader_type is 'fro':
            self.ReadFRO(filename, element_type)
        elif self.reader_type is 'read_separate':
            # READ MESH FROM SEPARATE FILES FOR CONNECTIVITY AND COORDINATES
            from Florence.Utils import insensitive
            # return insensitive(kwargs.keys())
            for key in kwargs.keys():
                inkey = insensitive(key)
                if "connectivity" in inkey and "delimiter" not in inkey:
                    connectivity_file = kwargs.get(key)
                if "coordinate" in insensitive(key) and "delimiter" not in inkey:
                    coordinates_file = kwargs.get(key)

            self.ReadSeparate(connectivity_file,coordinates_file,element_type,
                delimiter_connectivity=',',delimiter_coordinates=',')
        elif self.reader_type is 'ReadHDF5':
            self.ReadHDF5(filename)

        self.nnode = self.points.shape[0]
        # MAKE SURE MESH DATA IS CONTIGUOUS
        self.points   = np.ascontiguousarray(self.points)
        self.elements = np.ascontiguousarray(self.elements)
        return

    def ReadGmsh(self, filename, element_type, read_surface_info=False):
        """Read gmsh (.msh) file"""

        try:
            fid = open(filename, "r")
        except IOError:
            print("File '%s' not found." % (filename))
            sys.exit()

        msh_version = None
        # CHECK MSH FILE VERSION
        if "MeshFormat" in fid.readline():
            msh_version = int(np.floor(float(fid.readline().split(" ")[0])))
            if 4 != msh_version and 2 != msh_version:
                raise IOError("Only ASCII version 2 and 4 (>=4.1) .msh file formats are supported")
        if 4 != msh_version and 2 != msh_version:
            raise IOError("Only ASCII version 2 and 4 (>=4.1) .msh file formats are supported")
        fid.close()

        if self.elements is not None and self.points is not None:
            self.__reset__()

        self.filename = filename

        bel = -1
        if element_type == "line":
            el = 1
        elif element_type == "tri":
            el = 2
            bel = 2
        elif element_type == "quad":
            el = 3
            bel = 3
        elif element_type == "tet":
            el = 4
            bel = 2
        elif element_type == "hex":
            el = 5
            bel = 3
        else:
            raise ValueError("Element type not understood")


        # NEW FAST READER
        var = 0 # for old gmsh versions - needs checks
        node_blocks, elem_blocks, face_blocks = None, None, None
        rem_nnode, rem_nelem, rem_faces = int(1e09), int(1e09), int(1e09)
        face_counter = 0
        for line_counter, line in enumerate(open(filename)):
            item = line.rstrip()
            plist = item.split()
            if plist[0] == "Dimension":
                self.ndim = plist[1]
            elif plist[0] == "Vertices":
                rem_nnode = line_counter+1
                continue
            elif plist[0] == "$Nodes":
                rem_nnode = line_counter+1
                continue
            elif plist[0] == "Triangles":
                rem_faces = line_counter+1
                continue
            elif plist[0] == "Tetrahedra":
                rem_nelem = line_counter+1
                continue
            elif plist[0] == "$Elements":
                rem_nelem = line_counter+1
                var = 1
                continue

            if msh_version == 2:
                if rem_nnode == line_counter:
                    self.nnode = int(plist[0])
                if rem_faces == line_counter:
                    face_counter = int(plist[0])
                if rem_nelem == line_counter:
                    self.nelem = int(plist[0])
                    break
            else:
                if rem_nnode == line_counter:
                    node_blocks, self.nnode = int(plist[0]), int(plist[1])
                if rem_faces == line_counter:
                    face_blocks, face_counter = int(plist[0]), int(plist[1])
                if rem_nelem == line_counter:
                    elem_blocks, self.nelem = int(plist[0]), int(plist[1])
                    break

        points, elements, faces, face_to_surface = [],[], [], []
        if msh_version == 2:
            # RE-READ
            ns = self.InferNumberOfNodesPerElement(p=1,element_type=element_type)
            for line_counter, line in enumerate(open(filename)):
                item = line.rstrip()
                plist = item.split()
                if var == 0:
                    if line_counter > rem_nnode and line_counter < self.nnode+rem_nnode+1:
                        points.append([float(i) for i in plist[:3]])
                    if line_counter > rem_nelem and line_counter < self.nelem+rem_nelem+1:
                        elements.append([int(i) for i in plist[:4]])
                elif var == 1:
                    if line_counter > rem_nnode and line_counter < self.nnode+rem_nnode+1:
                        points.append([float(i) for i in plist[1:]])
                    if line_counter > rem_nelem and line_counter < self.nelem+rem_nelem+1:
                        if int(plist[1]) == el:
                            elements.append([int(i) for i in plist[-ns:]])

                        # READ SURFACE INFO - CERTAINLY ONLY IF SURFACE ELEMENT TYPE IS QUADS/TRIS
                        if read_surface_info:
                            if int(plist[1]) == bel:
                                faces.append([int(i) for i in plist[5:]])
                                face_to_surface.append(int(plist[4]))


        elif msh_version == 4:
            # RE-READ
            fid = open(filename)
            content = fid.readlines()

            # READ NODES
            nodes_content = content[rem_nnode+1:2*self.nnode+node_blocks+rem_nnode+1]
            incrementer, line_number = 0, 0
            # LOOP OVER BLOCKS
            for i in range(node_blocks):
                incrementer = int(nodes_content[line_number].rstrip().split()[3])
                # LOOP OVER NODES OF EACH BLOCK
                for j in range(line_number+1, line_number+2*incrementer+1):
                    plist = nodes_content[j].rstrip().split()
                    if len(plist) == 1:
                        continue
                    points.append([float(plist[k]) for k in range(0,len(plist))])
                line_number += 2*incrementer + 1

            # READ ELEMENTS
            elems_content = content[rem_nelem+1:self.nelem+elem_blocks+rem_nelem+1]
            incrementer, line_number = 0, 0
            # LOOP OVER BLOCKS
            for i in range(elem_blocks):
                incrementer = int(elems_content[line_number].rstrip().split()[3])
                if el == int(elems_content[line_number].rstrip().split()[2]):
                    # LOOP OVER ELEMENTS OF EACH BLOCK
                    for j in range(line_number+1, line_number+incrementer+1):
                        plist = elems_content[j].rstrip().split()
                        elements.append([int(plist[k]) for k in range(1,len(plist))])
                line_number += incrementer + 1

            if read_surface_info:
                # READ FACES
                incrementer, line_number = 0, 0
                # LOOP OVER BLOCKS
                for i in range(elem_blocks):
                    incrementer = int(elems_content[line_number].rstrip().split()[3])
                    surface_tag = int(elems_content[line_number].rstrip().split()[1])
                    if bel == int(elems_content[line_number].rstrip().split()[2]):
                        # LOOP OVER FACES OF EACH BLOCK
                        for j in range(line_number+1, line_number+incrementer+1):
                            plist = elems_content[j].rstrip().split()
                            faces.append([int(plist[k]) for k in range(1,len(plist))])
                            face_to_surface.append(surface_tag)
                    line_number += incrementer + 1


        self.points = np.array(points,copy=True)
        self.elements = np.array(elements,copy=True) - 1
        # CORRECT
        self.nelem = self.elements.shape[0]
        self.nnode = self.points.shape[0]
        if self.nelem == 0:
            raise ValueError("msh file does not contain {} elements".format(element_type))

        if read_surface_info:
            self.faces = np.array(faces,copy=True) - 1
            self.face_to_surface = np.array(face_to_surface, dtype=np.int64, copy=True).flatten()
            self.face_to_surface -= 1
            # CHECK IF FILLED
            if isinstance(self.face_to_surface,list):
                if not self.face_to_surface:
                    self.face_to_surface = None
            elif isinstance(self.face_to_surface,np.ndarray):
                if self.face_to_surface.shape[0]==0:
                    self.face_to_surface = None

        if self.points.shape[1] == 3:
            if np.allclose(self.points[:,2],0.):
                self.points = np.ascontiguousarray(self.points[:,:2])

        self.element_type = element_type
        if self.element_type == "tri" or self.element_type == "quad":
            self.GetEdges()
            self.GetBoundaryEdges()
        elif self.element_type == "tet" or self.element_type == "hex":
            self.GetFaces()
            self.GetBoundaryFaces()
            self.GetBoundaryEdges()

        return

    def InferPolynomialDegree(self):
        """Infer the degree of interpolation (p) based on the shape of
            self.elements

            returns:        [int] polynomial degree
            """

        assert self.element_type is not None
        assert self.elements is not None

        if self.degree is not None:
            if isinstance(self.degree,np.ndarray):
                self.degree = np.asscalar(self.degree)
            i = self.degree
            if self.element_type == "tet" and (i+1)*(i+2)*(i+3)/6==self.elements.shape[1]:
                return self.degree
            if self.element_type == "tri" and (i+1)*(i+2)/2==self.elements.shape[1]:
                return self.degree


        p = 0
        if self.element_type == "tet":
            for i in range(100):
                if (i+1)*(i+2)*(i+3)/6==self.elements.shape[1]:
                    p = i
                    break

        elif self.element_type == "tri":
            for i in range(100):
                if (i+1)*(i+2)/2==self.elements.shape[1]:
                    p = i
                    break

        elif self.element_type == "hex":
            for i in range(100):
                if int((i+1)**3)==self.elements.shape[1]:
                    p = i
                    break

        elif self.element_type == "quad":
            for i in range(100):
                if int((i+1)**2)==self.elements.shape[1]:
                    p = i
                    break

        elif self.element_type == "line":
            for i in range(100):
                if int(i+1)==self.elements.shape[1]:
                    p = i
                    break

        elif self.element_type == "pent":
            if 5==self.elements.shape[1]:
                p = 1
            else:
                raise NotImplementedError("High order pentagonal elements are not supported yet")

        self.degree = p
        return p

    def InferNumberOfNodesPerElement(self, p=None, element_type=None):
        """Infers number of nodes per element. If p and element_type are
            not None then returns the number of nodes required for the given
            element type with the given polynomial degree"""

        if p is not None and element_type is not None:
            if element_type=="line":
                return int(p+1)
            elif element_type=="tri":
                return int((p+1)*(p+2)/2)
            elif element_type=="quad":
                return int((p+1)**2)
            elif element_type=="tet":
                return int((p+1)*(p+2)*(p+3)/6)
            elif element_type=="hex":
                return int((p+1)**3)
            else:
                raise ValueError("Did not understand element type")

        assert self.elements.shape[0] is not None
        return self.elements.shape[1]

