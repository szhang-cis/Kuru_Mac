from __future__ import print_function, division
import os, sys, warnings, platform
from time import time
import numpy as np
#if "PyPy" not in platform.python_implementation():
#    from scipy.io import loadmat, savemat
from Kuru.Tensor import unique2d, itemfreq, in2d, makezero
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
from warnings import warn
from copy import deepcopy



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

    @property
    def Bounds(self):
        """Returns bounds of a mesh i.e. the minimum and maximum coordinate values
            in every direction
        """
        assert self.points is not None

        if self.points.shape[1] == 3:
            bounds = np.array([[np.min(self.points[:,0]),
                        np.min(self.points[:,1]),
                        np.min(self.points[:,2])],
                        [np.max(self.points[:,0]),
                        np.max(self.points[:,1]),
                        np.max(self.points[:,2])]])
            makezero(bounds)
            return bounds
        elif self.points.shape[1] == 2:
            bounds = np.array([[np.min(self.points[:,0]),
                        np.min(self.points[:,1])],
                        [np.max(self.points[:,0]),
                        np.max(self.points[:,1])]])
            makezero(bounds)
            return bounds
        elif self.points.shape[1] == 1:
            bounds = np.array([[np.min(self.points[:,0])],
                        [np.max(self.points[:,0])]])
            makezero(bounds)
            return bounds
        else:
            raise ValueError("Invalid dimension for mesh coordinates")

    def GetBoundaryEdges(self):
        assert self.element_type is not None
        if self.element_type == "tri":
            self.GetBoundaryEdgesTri()
        elif self.element_type == "quad":
            self.GetBoundaryEdgesQuad()
        elif self.element_type == "pent":
            self.GetBoundaryEdgesPent()
        elif self.element_type == "tet":
            self.GetBoundaryEdgesTet()
        elif self.element_type == "hex":
            self.GetBoundaryEdgesHex()
        else:
            raise ValueError('Type of element not understood')
        return self.edges

    def GetBoundaryEdgesHex(self):
        """Find boundary edges (lines) of hexahedral mesh.
        """

        p = self.InferPolynomialDegree()
        # DO NOT COMPUTE IF ALREADY COMPUTED
        if isinstance(self.edges,np.ndarray):
            if self.edges.shape[0] > 1:
                # IF LINEAR VERSION IS COMPUTED, DO COMPUTE HIGHER VERSION
                if self.edges.shape[1] == 2 and p > 1:
                    pass
                else:
                    return


        # FIRST GET BOUNDARY FACES
        if not isinstance(self.faces,np.ndarray):
            self.GetBoundaryFacesHex()

        # BUILD A 2D MESH
        tmesh = Mesh()
        tmesh.element_type = "quad"
        tmesh.elements = self.faces
        tmesh.nelem = tmesh.elements.shape[0]
        del tmesh.faces
        del tmesh.points

        # ALL THE EDGES CORRESPONDING TO THESE BOUNDARY FACES ARE BOUNDARY EDGES
        self.edges =  tmesh.GetEdgesQuad()

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

    def GetBoundaryFaces(self):
        assert self.element_type is not None
        if self.element_type == "tet":
            self.GetBoundaryFacesTet()
        elif self.element_type == "hex":
            self.GetBoundaryFacesHex()
        elif self.element_type=="tri" or self.element_type=="quad":
            raise ValueError("2D mesh does not have faces")
        else:
            raise ValueError('Type of element not understood')
        return self.faces

    def GetBoundaryFacesHex(self):
        """Find boundary faces (surfaces) of a hexahedral mesh"""

        p = self.InferPolynomialDegree()

        # DO NOT COMPUTE IF ALREADY COMPUTED
        if isinstance(self.faces,np.ndarray):
            if self.faces.shape[0] > 1:
                # IF LINEAR VERSION IS COMPUTED, DO COMPUTE HIGHER VERSION
                if self.faces.shape[1] == 4 and p > 1:
                    pass
                else:
                    return

        node_arranger = NodeArrangementHex(p-1)[0]

        # CONCATENATE ALL THE FACES MADE FROM ELEMENTS
        all_faces = np.concatenate((np.concatenate((
                np.concatenate((np.concatenate((np.concatenate((self.elements[:,node_arranger[0,:]],
                self.elements[:,node_arranger[1,:]]),axis=0),self.elements[:,node_arranger[2,:]]),axis=0),
                self.elements[:,node_arranger[3,:]]),axis=0),self.elements[:,node_arranger[4,:]]),axis=0),
                self.elements[:,node_arranger[5,:]]),axis=0).astype(np.int64)
        # GET UNIQUE ROWS
        uniques, idx, inv = unique2d(all_faces,consider_sort=True,order=False,return_index=True,return_inverse=True)

        # ROWS THAT APPEAR ONLY ONCE CORRESPOND TO BOUNDARY FACES
        freqs_inv = itemfreq(inv)
        faces_ext_flags = freqs_inv[freqs_inv[:,1]==1,0]
        # NOT ARRANGED
        self.faces = uniques[faces_ext_flags,:]

        # DETERMINE WHICH FACE OF THE ELEMENT THEY ARE
        boundary_face_to_element = np.zeros((faces_ext_flags.shape[0],2),dtype=np.int64)

        # FURTHER RE-ARRANGEMENT / ARANGE THE NODES BASED ON THE ORDER THEY APPEAR
        # IN ELEMENT CONNECTIVITY
        # THIS STEP IS NOT NECESSARY INDEED - ITS JUST FOR RE-ARANGMENT OF FACES
        all_faces_in_faces = in2d(all_faces,self.faces,consider_sort=True)
        all_faces_in_faces = np.where(all_faces_in_faces==True)[0]

        # boundary_face_to_element = np.zeros((all_faces_in_faces.shape[0],2),dtype=np.int64)
        boundary_face_to_element[:,0] = all_faces_in_faces % self.elements.shape[0]
        boundary_face_to_element[:,1] = all_faces_in_faces // self.elements.shape[0]

        # ARRANGE FOR ANY ORDER OF BASES/ELEMENTS AND ASSIGN DATA MEMBERS
        self.faces = self.elements[boundary_face_to_element[:,0][:,None],node_arranger[boundary_face_to_element[:,1],:]]
        self.faces = self.faces.astype(np.uint64)
        self.boundary_face_to_element = boundary_face_to_element

    def GetEdgesQuad(self):
        """Find the all edges of a quadrilateral mesh.
            Sets all_edges property and returns it
        returns:
            arr:            numpy ndarray of all edges"""

        p = self.InferPolynomialDegree()

        # DO NOT COMPUTE IF ALREADY COMPUTED
        if isinstance(self.all_edges,np.ndarray):
            if self.all_edges.shape[0] > 1:
                # IF LINEAR VERSION IS COMPUTED, DO COMPUTE HIGHER VERSION
                if self.all_edges.shape[1]==2 and p > 1:
                    pass
                else:
                    return self.all_edges

        node_arranger = NodeArrangementQuad(p-1)[0]

        # GET ALL EDGES FROM THE ELEMENT CONNECTIVITY
        edges = np.concatenate((self.elements[:,node_arranger[0,:]],self.elements[:,node_arranger[1,:]],
            self.elements[:,node_arranger[2,:]],self.elements[:,node_arranger[3,:]]),axis=0).astype(np.uint64)

        # REMOVE DUPLICATES
        edges, idx = unique2d(edges,consider_sort=True,order=False,return_index=True)

        edge_to_element = np.zeros((edges.shape[0],2),np.int64)
        edge_to_element[:,0] =  idx % self.elements.shape[0]
        edge_to_element[:,1] =  idx // self.elements.shape[0]

        self.edge_to_element = edge_to_element

        # DO NOT SET all_edges IF THE CALLER FUNCTION IS GetBoundaryEdgesHex
        import inspect
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)[1][3]

        if calframe != "GetBoundaryEdgesHex":
            self.all_edges = edges

        return edges

    def GetBoundaryEdgesQuad(self):
        """Find boundary edges (lines) of a quadrilateral mesh"""

        p = self.InferPolynomialDegree()

        # DO NOT COMPUTE IF ALREADY COMPUTED
        if isinstance(self.edges,np.ndarray):
            if self.edges.shape[0] > 1:
                # IF LINEAR VERSION IS COMPUTED, DO COMPUTE HIGHER VERSION
                if self.edges.shape[1] == 2 and p > 1:
                    pass
                else:
                    return

        node_arranger = NodeArrangementQuad(p-1)[0]

        # GET ALL EDGES FROM THE ELEMENT CONNECTIVITY
        all_edges = np.concatenate((self.elements[:,node_arranger[0,:]],self.elements[:,node_arranger[1,:]],
            self.elements[:,node_arranger[2,:]],self.elements[:,node_arranger[3,:]]),axis=0).astype(np.uint64)

        # GET UNIQUE ROWS
        uniques, idx, inv = unique2d(all_edges,consider_sort=True,order=False,return_index=True,return_inverse=True)

        # ROWS THAT APPEAR ONLY ONCE CORRESPOND TO BOUNDARY EDGES
        freqs_inv = itemfreq(inv)
        edges_ext_flags = freqs_inv[freqs_inv[:,1]==1,0]
        # NOT ARRANGED
        self.edges = uniques[edges_ext_flags,:]

        # DETERMINE WHICH FACE OF THE ELEMENT THEY ARE
        boundary_edge_to_element = np.zeros((edges_ext_flags.shape[0],2),dtype=np.int64)

        # FURTHER RE-ARRANGEMENT / ARANGE THE NODES BASED ON THE ORDER THEY APPEAR
        # IN ELEMENT CONNECTIVITY
        # THIS STEP IS NOT NECESSARY INDEED - ITS JUST FOR RE-ARANGMENT OF EDGES
        all_edges_in_edges = in2d(all_edges,self.edges,consider_sort=True)
        all_edges_in_edges = np.where(all_edges_in_edges==True)[0]

        boundary_edge_to_element[:,0] = all_edges_in_edges % self.elements.shape[0]
        boundary_edge_to_element[:,1] = all_edges_in_edges // self.elements.shape[0]

        # ARRANGE FOR ANY ORDER OF BASES/ELEMENTS AND ASSIGN DATA MEMBERS
        self.edges = self.elements[boundary_edge_to_element[:,0][:,None],node_arranger[boundary_edge_to_element[:,1],:]]
        self.edges = self.edges.astype(np.uint64)
        self.boundary_edge_to_element = boundary_edge_to_element

        return self.edges

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

    def GetHighOrderMesh(self,p=1, silent=True, **kwargs):
        """Given a linear tri, tet, quad or hex mesh compute high order mesh based on it.
        This is a static method linked to the HigherOrderMeshing module"""

        if not isinstance(p,int):
            raise ValueError("p must be an integer")
        else:
            if p < 1:
                raise ValueError("Value of p={} is not acceptable. Provide p>=1.".format(p))

        if self.degree is None:
            self.InferPolynomialDegree()

        C = p-1
        if 'C' in kwargs.keys():
            if kwargs['C'] != p - 1:
                raise ValueError("Did not understand the specified interpolation degree of the mesh")
            del kwargs['C']

        # DO NOT COMPUTE IF ALREADY COMPUTED FOR THE SAME ORDER
        if self.degree == None:
            self.degree = self.InferPolynomialDegree()
        if self.degree == p:
            return

        # SITUATIONS WHEN ANOTHER HIGH ORDER MESH IS REQUIRED, WITH ONE HIGH
        # ORDER MESH ALREADY AVAILABLE
        if self.degree != 1 and self.degree - 1 != C:
            dum = self.GetLinearMesh(remap=True)
            self.__dict__.update(dum.__dict__)

        if not silent:
            print('Generating p = '+str(C+1)+' mesh based on the linear mesh...')
        t_mesh = time()
        # BUILD A NEW MESH BASED ON THE LINEAR MESH
        if self.element_type == 'line':
            nmesh = HighOrderMeshLine(C,self,**kwargs)
        if self.element_type == 'tri':
            if self.edges is None:
                self.GetBoundaryEdgesTri()
            # nmesh = HighOrderMeshTri(C,self,**kwargs)
            nmesh = HighOrderMeshTri_SEMISTABLE(C,self,**kwargs)
        elif self.element_type == 'tet':
            # nmesh = HighOrderMeshTet(C,self,**kwargs)
            nmesh = HighOrderMeshTet_SEMISTABLE(C,self,**kwargs)
        elif self.element_type == 'quad':
            if self.edges is None:
                self.GetBoundaryEdgesTri()
            nmesh = HighOrderMeshQuad(C,self,**kwargs)
        elif self.element_type == 'hex':
            nmesh = HighOrderMeshHex(C,self,**kwargs)

        self.points = nmesh.points
        self.elements = nmesh.elements.astype(np.uint64)
        if isinstance(self.corners,np.ndarray):
            # NOT NECESSARY BUT GENERIC
            self.corners = nmesh.corners.astype(np.uint64)
        if isinstance(self.edges,np.ndarray):
            self.edges = nmesh.edges.astype(np.uint64)
        if isinstance(self.faces,np.ndarray):
            if isinstance(nmesh.faces,np.ndarray):
                self.faces = nmesh.faces.astype(np.uint64)
        self.nelem = nmesh.nelem
        self.nnode = self.points.shape[0]
        self.element_type = nmesh.info
        self.degree = C+1

        self.ChangeType()

        if not silent:
            print('Finished generating the high order mesh. Time taken', time()-t_mesh,'sec')

    def Rectangle(self,lower_left_point=(0,0), upper_right_point=(2,1),
        nx=5, ny=5, element_type="tri"):
        """Creates a quad/tri mesh of a rectangle"""

        if element_type != "tri" and element_type != "quad":
            raise ValueError("Element type should either be tri or quad")

        if self.elements is not None and self.points is not None:
            self.__reset__()

        if (lower_left_point[0] > upper_right_point[0]) or \
            (lower_left_point[1] > upper_right_point[1]):
            raise ValueError("Incorrect coordinate for lower left and upper right vertices")

        nx, ny = int(nx), int(ny)
        if nx <= 0 or ny <= 0:
            raise ValueError("Number of discretisation cannot be zero or negative: nx={} ny={}".format(nx,ny))



        from scipy.spatial import Delaunay

        x=np.linspace(lower_left_point[0],upper_right_point[0],nx+1)
        y=np.linspace(lower_left_point[1],upper_right_point[1],ny+1)

        X,Y = np.meshgrid(x,y)
        coordinates = np.dstack((X.ravel(),Y.ravel()))[0,:,:]

        if element_type == "tri":
            tri_func = Delaunay(coordinates)
            self.element_type = "tri"
            self.elements = tri_func.simplices
            self.nelem = self.elements.shape[0]
            self.points = tri_func.points
            self.nnode = self.points.shape[0]
            self.GetBoundaryEdgesTri()

        elif element_type == "quad":

            self.nelem = int(nx*ny)
            elements = np.zeros((self.nelem,4),dtype=np.int64)

            dum_0 = np.arange((nx+1)*ny)
            dum_1 = np.array([(nx+1)*i+nx for i in range(ny)])
            col0 = np.delete(dum_0,dum_1)
            elements[:,0] = col0
            elements[:,1] = col0 + 1
            elements[:,2] = col0 +  nx + 2
            elements[:,3] = col0 +  nx + 1

            self.nnode = int((nx+1)*(ny+1))
            self.element_type = "quad"
            self.elements = elements
            self.points = coordinates
            self.nnode = self.points.shape[0]
            self.GetBoundaryEdgesQuad()
            self.GetEdgesQuad()

    def GetNodeCommonality(self):
        """Finds the elements sharing a node.
            The return values are linked lists [list of numpy of arrays].
            Each numpy array within the list gives the elements that contain a given node.
            As a result the size of the linked list is nnode

            outputs:
                els:                        [list of numpy arrays] element numbers containing nodes
                pos:                        [list of numpy arrays] elemental positions of the nodes
                res_flat:                   [list of numpy arrays] position of nodes in the
                                            flattened element connectivity.
        """

        self.__do_essential_memebers_exist__()

        elements = self.elements.ravel()
        idx_sort = np.argsort(elements)
        sorted_elements = elements[idx_sort]
        vals, idx_start = np.unique(sorted_elements, return_index=True)

        # Sets of indices
        flat_pos = np.split(idx_sort, idx_start[1:])
        els = np.split(idx_sort // int(self.elements.shape[1]), idx_start[1:])
        pos = np.split(idx_sort %  int(self.elements.shape[1]), idx_start[1:])

        # In case one wants to return only the duplicates i.e. filter keeping only items occurring more than once
        # vals, idx_start, count = np.unique(sorted_elements, return_counts=True, return_index=True)
        # vals = vals[count > 1]
        # res = filter(lambda x: x.size > 1, res)

        return els, pos, flat_pos

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
            #self.ReadSalome(filename, element_type=element_type, read_surface_info=read_surface_info)
            raise ValueError("Reader not implemented yet")
        elif reader_type is 'GID':
            #self.ReadGIDMesh(filename, element_type, order)
            raise ValueError("Reader not implemented yet")
        elif self.reader_type is 'gmsh':
            self.ReadGmsh(filename, element_type=element_type, read_surface_info=read_surface_info)
        elif self.reader_type is 'obj':
            #self.ReadOBJ(filename, element_type=element_type, read_surface_info=read_surface_info)
            raise ValueError("Reader not implemented yet")
        elif self.reader_type is 'fenics':
            #self.ReadFenics(filename, element_type)
            raise ValueError("Reader not implemented yet")
        elif self.reader_type is 'vtu':
            #self.ReadVTK(filename)
            raise ValueError("Reader not implemented yet")
        elif self.reader_type is 'unv':
            #self.ReadUNV(filename, element_type)
            raise ValueError("Reader not implemented yet")
        elif self.reader_type is 'fro':
            #self.ReadFRO(filename, element_type)
            raise ValueError("Reader not implemented yet")
        elif self.reader_type is 'read_separate':
            # READ MESH FROM SEPARATE FILES FOR CONNECTIVITY AND COORDINATES
            raise ValueError("Reader not implemented yet")
            from Florence.Utils import insensitive
            # return insensitive(kwargs.keys())
            #for key in kwargs.keys():
            #    inkey = insensitive(key)
            #    if "connectivity" in inkey and "delimiter" not in inkey:
            #        connectivity_file = kwargs.get(key)
            #    if "coordinate" in insensitive(key) and "delimiter" not in inkey:
            #        coordinates_file = kwargs.get(key)

            #self.ReadSeparate(connectivity_file,coordinates_file,element_type,
            #    delimiter_connectivity=',',delimiter_coordinates=',')
        elif self.reader_type is 'ReadHDF5':
            #self.ReadHDF5(filename)
            raise ValueError("Reader not implemented yet")

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

    def ChangeType(self):
        """Change mesh data type from signed to unsigned"""

        self.__do_essential_memebers_exist__()
        self.points = np.ascontiguousarray(self.points.astype(np.float64))
        if isinstance(self.elements,np.ndarray):
            self.elements = np.ascontiguousarray(self.elements.astype(np.uint64))
        if hasattr(self, 'edges'):
            if isinstance(self.edges,np.ndarray):
                self.edges = np.ascontiguousarray(self.edges.astype(np.uint64))
        if hasattr(self, 'faces'):
            if isinstance(self.faces,np.ndarray):
                self.faces = np.ascontiguousarray(self.faces.astype(np.uint64))

    def IsHighOrder(self):
        is_high_order = False
        if self.InferPolynomialDegree() > 1:
            is_high_order = True
        return is_high_order

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

    def InferNumberOfNodesPerLinearElement(self, element_type=None):
        """Infers number of nodes per element. If element_type are
            not None then returns the number of nodes required for the given
            element type"""

        if element_type is None and self.element_type is None:
            raise ValueError("Did not understand element type")
        if element_type is None:
            element_type = self.element_type

        tmp = self.element_type
        if element_type != self.element_type:
            self.element_type = element_type

        nodeperelem = None
        if element_type=="line":
            nodeperelem = 2
        elif element_type=="tri":
            nodeperelem = 3
        elif element_type=="quad":
            nodeperelem = 4
        elif element_type=="tet":
            nodeperelem = 4
        elif element_type=="hex":
            nodeperelem = 8
        else:
            raise ValueError("Did not understand element type")

        self.element_type = tmp

        return nodeperelem

    def InferSpatialDimension(self):
        """Infer the spatial dimension of the mesh"""

        assert self.points is not None
        # if self.points.shape[1] == 3:
        #     if self.element_type == "tri" or self.element_type == "quad":
        #         print("3D surface mesh of ", self.element_type)

        return self.points.shape[1]

    def InferElementType(self):

        if self.element_type is not None:
            return self.element_type

        assert self.elements is not None
        assert self.points is not None

        ndim = self.InferSpatialDimension()
        nodeperelem = self.InferNumberOfNodesPerElement()
        nn = 20

        if ndim==3:
            if nodeperelem in [int((i+1)*(i+2)*(i+3)/6) for i in range(1,nn)]:
                self.element_type = "tet"
            elif nodeperelem in [int((i+1)**3) for i in range(1,nn)]:
                self.element_type = "hex"
            else:
                if nodeperelem in [int((i+1)*(i+2)/2) for i in range(1,nn)]:
                    self.element_type = "tri"
                elif nodeperelem in [int((i+1)**2) for i in range(1,nn)]:
                    self.element_type = "quad"
                else:
                    raise ValueError("Element type not understood")
        elif ndim==2:
            if nodeperelem in [int((i+1)*(i+2)/2) for i in range(1,nn)]:
                self.element_type = "tri"
            elif nodeperelem in [int((i+1)**2) for i in range(1,nn)]:
                self.element_type = "quad"
            else:
                raise ValueError("Element type not understood")
        elif ndim==1:
            self.element_type = "line"
        else:
            raise ValueError("Element type not understood")

        # IF POINTS ARE CO-PLANAR THEN IT IS NOT TET BUT QUAD
        if ndim == 3 and self.element_type == "tet":
            a = self.points[self.elements[:,0],:]
            b = self.points[self.elements[:,1],:]
            c = self.points[self.elements[:,2],:]
            d = self.points[self.elements[:,3],:]

            det_array = np.dstack((a-d,b-d,c-d))
            # FIND VOLUME OF ALL THE ELEMENTS
            volume = 1./6.*np.linalg.det(det_array)
            if np.allclose(volume,0.0):
                self.element_type = "quad"

        return self.element_type

    def InferBoundaryElementType(self):

        self.InferElementType()
        if self.element_type == "hex":
            self.boundary_element_type = "quad"
        elif self.element_type == "tet":
            self.boundary_element_type = "tri"
        elif self.element_type == "quad" or self.element_type == "tri":
            self.boundary_element_type = "line"
        elif self.element_type == "line":
            self.boundary_element_type = "point"
        else:
            raise ValueError("Could not understand element type")

        return self.boundary_element_type

    def CreateDummyLowerDimensionalMesh(self):
        """Create a dummy lower dimensional mesh that would have some specific mesh attributes at least.
            The objective is that the lower dimensional mesh should have the same element type as the
            boundary faces/edges of the actual mesh and be the same order"""


        sys.stdout = open(os.devnull, "w")
        p = self.InferPolynomialDegree()
        mesh = Mesh()
        if self.element_type == "tet":
            mesh.Rectangle(nx=1,ny=1, element_type="tri")
            mesh.GetHighOrderMesh(p=p)
        elif self.element_type == "hex":
            mesh.Rectangle(nx=1,ny=1, element_type="quad")
            mesh.GetHighOrderMesh(p=p)
        elif self.element_type == "tri" or self.element_type == "quad":
            mesh.Line(n=1, p=p)
        sys.stdout = sys.__stdout__

        return mesh

    def GetLinearMesh(self, solution=None, remap=False):
        """Returns the linear mesh from a high order mesh. If mesh is already linear returns the same mesh.
            Also maps any solution vector/tensor of high order mesh to the linear mesh, if supplied.
            For safety purposes, always makes a copy"""

        self.__do_essential_memebers_exist__()

        ndim = self.InferSpatialDimension()
        if ndim==2:
            if self.element_type == "tri" or self.element_type == "quad":
                assert self.edges is not None
        elif ndim==3:
            if self.element_type == "tet" or self.element_type == "hex":
                assert self.faces is not None


        if self.IsHighOrder is False:
            if solution is not None:
                return deepcopy(self), deepcopy(solution)
            return deepcopy(self)
        else:
            if not remap:
                # WORKS ONLY IF THE FIST COLUMNS CORRESPOND TO
                # LINEAR CONNECTIVITY
                lmesh = Mesh()
                lmesh.element_type = self.element_type
                lmesh.degree = 1
                if self.element_type == "tri":
                    lmesh.elements = np.copy(self.elements[:,:3])
                    lmesh.edges = np.copy(self.edges[:,:2])
                    lmesh.nnode = int(np.max(lmesh.elements)+1)
                    lmesh.points = np.copy(self.points[:lmesh.nnode,:])
                elif self.element_type == "tet":
                    lmesh.elements = np.copy(self.elements[:,:4])
                    lmesh.faces = np.copy(self.faces[:,:3])
                    lmesh.nnode = int(np.max(lmesh.elements)+1)
                    lmesh.points = np.copy(self.points[:lmesh.nnode,:])
                elif self.element_type == "quad":
                    lmesh.elements = np.copy(self.elements[:,:4])
                    lmesh.edges = np.copy(self.edges[:,:2])
                    lmesh.nnode = int(np.max(lmesh.elements)+1)
                    lmesh.points = np.copy(self.points[:lmesh.nnode,:])
                elif self.element_type == "hex":
                    lmesh.elements = np.copy(self.elements[:,:8])
                    lmesh.faces = np.copy(self.faces[:,:4])
                    lmesh.nnode = int(np.max(lmesh.elements)+1)
                    lmesh.points = np.copy(self.points[:lmesh.nnode,:])
                lmesh.nelem = lmesh.elements.shape[0]

                if solution is not None:
                    solution = solution[np.unique(lmesh.elements),...]
                    return lmesh, solution

            else:
                # WORKS FOR ALL CASES BUT REMAPS (NO MAPPING BETWEEN LOW AND HIGH ORDER)
                nodeperelem = self.InferNumberOfNodesPerLinearElement()
                lmesh = Mesh()
                lmesh.element_type = self.element_type
                lmesh.nelem = self.nelem
                unnodes, inv = np.unique(self.elements[:,:nodeperelem], return_inverse=True)
                aranger = np.arange(lmesh.nelem*nodeperelem)
                lmesh.elements = inv[aranger].reshape(lmesh.nelem,nodeperelem)
                lmesh.points = self.points[unnodes,:]
                if lmesh.element_type == "hex" or lmesh.element_type == "tet":
                    lmesh.GetBoundaryFaces()
                    lmesh.GetBoundaryEdges()
                elif lmesh.element_type == "quad" or lmesh.element_type == "tri":
                    lmesh.GetBoundaryEdges()

                if solution is not None:
                    solution = solution[unnodes,...]
                    return lmesh, solution

        return lmesh

    def __do_essential_memebers_exist__(self):
        """Check if essential members exist"""
        assert self.element_type is not None
        assert self.elements is not None
        assert self.points is not None


    def __update__(self,other):
        self.__dict__.update(other.__dict__)


    def __reset__(self):
        """Class resetter. Resets all elements of the class
        """

        for i in self.__dict__.keys():
            self.__dict__[i] = None
