module Mesh
include("NodeArrangement.jl")
#function __reset__(self)
#   """Class resetter. Resets all elements of the class
#   """
#   for i in self.__dict__.keys():
#      self.__dict__[i] = None
#   end
#end

function InferPolynomialDegree(element_type,nodes_per_element)
   """Infer the degree of interpolation (p) based on the shape of
      self.elements

      returns:        [int] polynomial degree
   """
   #@assert element_type!=nothing
   #@assert elements!=nothing

   # check if degree is already asign
   #if degree!=nothing
   #   if isa(degree,Array)
   #      degree = np.asscalar(self.degree)
   #   end
   #   i = degree
   #   if element_type == "tet" && (i+1)*(i+2)*(i+3)/6==size(elements,2)
   #      return degree
   #   end
   #   if element_type == "tri" && (i+1)*(i+2)/2==size(elements,2)
   #      return degree
   #   end
   #end

   p = 0
   if element_type == "tet"
      for i=1:100
         if (i+1)*(i+2)*(i+3)/6==nodes_per_element
            p = i
            break
         end
      end
   elseif element_type == "tri"
      for i=1:100
         if (i+1)*(i+2)/2==nodes_per_element
            p = i
            break
         end
      end
   elseif element_type == "hex"
      for i=1:100
         if (i+1)^3==nodes_per_element
            p = i
            break
         end
      end
   elseif element_type == "quad"
      for i=1:100
         if (i+1)^2==nodes_per_element
            p = i
            break
         end
      end
   elseif element_type == "line"
      for i=1:100
         if i+1==nodes_per_element
            p = i
            break
         end
      end
   elseif element_type == "pent"
      if 5==nodes_per_element
         p = 1
      else
         ErrorException("High order pentagonal elements are not supported yet")
      end
   end
   #degree = p self.
   return p
end

function GetFacesHex(elements)
   """Find all faces (surfaces) in the hexahedral mesh (boundary & interior).
      Sets all_faces property and returns it

      returns:
      arr:         array of all faces
   """
   # DETERMINE DEGREE
   nodes_per_element = size(elements,2)
   p = InferPolynomialDegree("Hex",nodes_per_element)

   # DO NOT COMPUTE IF ALREADY COMPUTED
   #if isa(all_faces,Array)
   #if size(all_faces,1) > 1
      # IF LINEAR VERSION IS COMPUTED, DO COMPUTE HIGHER VERSION
   #   if size(all_faces,2) == 4 && p > 1
   #      pass
   #   else
   #      return all_faces
   #   end
   #end
   #end
   node_arranger = NodeArrangementHex(p)[1]
   fsize = ((p+1)^3)::Int

   # GET ALL FACES FROM THE ELEMENT CONNECTIVITY
   faces = np.concatenate((np.concatenate((
         np.concatenate((np.concatenate((np.concatenate((self.elements[:,node_arranger[0,:]],
         self.elements[:,node_arranger[1,:]]),axis=0),self.elements[:,node_arranger[2,:]]),axis=0),
         self.elements[:,node_arranger[3,:]]),axis=0),self.elements[:,node_arranger[4,:]]),axis=0),
         self.elements[:,node_arranger[5,:]]),axis=0).astype(np.int64)

   # REMOVE DUPLICATES
   all_faces, idx = unique2d(faces,consider_sort=True,order=false,return_index=true)

   face_to_element = zeros(Int,size(all_faces,1),2)
   face_to_element[:,1] =  idx % size(elements,1)
   face_to_element[:,2] =  idx // size(elements,1)

   #self.face_to_element = face_to_element

   return all_faces
end

function GetFaces(element_type)
   @assert element_type!=nothing
   #if element_type == "tet"
   #   self.GetFacesTet()
   if element_type == "hex"
      all_faces = GetFacesHex()
   elseif element_type=="tri" || element_type=="quad"
      ArgumentError("2D mesh does not have faces")
   else
      ArgumentError("Type of element not understood")
   end
   return all_faces
end

function InferNumberOfNodesPerElement(;p=nothing, element_type=nothing)
   """Infers number of nodes per element. If p and element_type are
      not None then returns the number of nodes required for the given
      element type with the given polynomial degree"""

   if p!=nothing && element_type!=nothing
      if element_type=="line"
         return (p+1)::Int
      elseif element_type=="tri"
         return ((p+1)*(p+2)/2)::Int
      elseif element_type=="quad"
         return ((p+1)^2)::Int
      elseif element_type=="tet"
         return ((p+1)*(p+2)*(p+3)/6)::Int
      elseif element_type=="hex"
         return ((p+1)^3)::Int
      else
         ArgumentError("Did not understand element type")
      end
   end
   @assert size(elements,1)!=nothing "Element list is empty" #self.
   return size(elements,2) #self.
end

function ReadGmsh(filename, element_type; read_surface_info=false)
   """Read gmsh (.msh) file"""
   try
      fid = open(filename, "r")
   catch #IOError
      println("File $filename not found.")
      ErrorException("Error at opening the file")
   end

   #if self.elements is not None and self.points is not None:
   #   self.__reset__()
   #end

   bel = -1
   if element_type == "line"
      el = 1
   elseif element_type == "tri"
      el = 2
      bel = 2
   elseif element_type == "quad"
      el = 3
      bel = 3
   elseif element_type == "tet"
      el = 4
      bel = 2
   elseif element_type == "hex"
      el = 5
      bel = 3
   else
      ArgumentError("Element type not understood")
   end

   # NEW FAST READER
   var = 0  # for old gmsh versions - needs checks
   rem_nnode, rem_nelem, rem_faces = convert(Int,1e09), convert(Int,1e09), convert(Int,1e09)
   face_counter = 0
   ndim = 0
   for (line_counter,line) in enumerate(readlines(open(filename)))
      plist = split(line)
      if plist[1] == "Dimension"
         ndim = parse(Int,plist[2]) #self.
      elseif plist[1] == "Vertices"
         rem_nnode = line_counter+1
         continue
      elseif plist[1] == "\$Nodes"
         rem_nnode = line_counter+1
         continue
      elseif plist[1] == "Triangles"
         rem_faces = line_counter+1
         continue
      elseif plist[1] == "Tetrahedra"
         rem_nelem = line_counter+1
         continue
      elseif plist[1] == "\$Elements"
         rem_nelem = line_counter+1
         var = 1
         continue
      end
      if rem_nnode == line_counter
         nnode = parse(Int,plist[1]) #self.
      end
      if rem_nnode+1 == line_counter && ndim == 0
         ndim = size(plist,1)-1
      end
      if rem_faces == line_counter
         face_counter = parse(Int,plist[1])
      end
      if rem_nelem == line_counter
         nelem = parse(Int,plist[1]) #self.
         break
      end
   end

   ne = InferNumberOfNodesPerElement(p=1,element_type=element_type)
   # Re-read
   points, elements, faces, face_to_surface = Float64[],Int[], Int[], Int[]
   elem_counter = 0 #to correct the number of elements
   face_counter = 0
   for (line_counter,line) in enumerate(readlines(open(filename)))
      global nf
      plist = split(line)
      if var==0
         if line_counter>rem_nnode && line_counter<nnode+rem_nnode+1
            append!(points,[parse(Float64,i) for i in plist[1:3]])
         end
         if line_counter>rem_nelem && line_counter<nelem+rem_nelem+1
            append!(elements,[parse(Int,i) for i in plist[1:4]])
         end
      elseif var==1
         if line_counter>rem_nnode && line_counter<nnode+rem_nnode+1
            append!(points,[parse(Float64,i) for i in plist[2:end]])
         end
         if line_counter>rem_nelem && line_counter<nelem+rem_nelem+1
            if parse(Int,plist[2]) == el
               elem_counter += 1
               append!(elements,[parse(Int,i) for i in plist[end-ne+1:end]])
            end
            # WRITE SURFACE INFO - CERTAINLY ONLY IF ELEMENT TYPE IS QUADS/TRIS
            if read_surface_info
               if parse(Int,plist[2]) == bel
                  nf = size(plist,1)-5
                  face_counter += 1
                  append!(faces,[parse(Int,i) for i in plist[6:end]])
                  append!(face_to_surface,parse(Int,plist[5]))
               end
            end
         end
      end
   end

   points = reshape(points,ndim,nnode) #self.
   points = transpose(points)
   points = copy(points)
   elements = reshape(elements,ne,elem_counter) #self.
   elements = transpose(elements)
   elements = copy(elements)
   # CORRECT
   nelem = size(elements,1) #self.
   nnode = size(points,1) #self.
   if nelem == 0 #self.
      ArgumentError("msh file does not contain $element_type elements")
   end
   if read_surface_info
      faces = reshape(faces,nf,face_counter) #self.
      faces = transpose(faces)
      faces = copy(faces)
      #self.face_to_surface = face_to_surface
      # CHECK IF FILLED
      if size(face_to_surface,1)==0
         face_to_surface = nothing
      end
   end
   # Reduce dimension if Z coordinates are zero
   if size(points,2) == 3 #self.
      if isapprox(points[:,3],zeros(Float64,size(points,1)))
         points = points[:,1:2]
      end
   end
   # Get boundary data of the mesh
   #if element_type == "tri" || element_type == "quad"
   #   self.GetEdges()
   #   self.GetBoundaryEdges()
   if element_type == "tet" || element_type == "hex"
      all_faces = GetFaces(element_type)
      GetBoundaryFaces()
      GetBoundaryEdges()
   end
   return
end

function Read(;filename=nothing, element_type="tri", reader_type=nothing, reader_type_format=nothing,
   reader_type_version=nothing, order=0, read_surface_info=false)
   """Convenience mesh reader method to dispatch call to subsequent apporpriate methods"""

   if !isa(filename,String)
      ArgumentError("filename must be a string")
      return
   end
   if reader_type != nothing
      if !isa(filename,String)
         ArgumentError("filename must be a string")
         return
      end
   end

   if reader_type==nothing
      if endswith(filename,"msh")
         reader_type = "gmsh"
      elseif endswith(filename,"obj")
         reader_type = "obj"
      elseif endswith(filename,"unv")
         reader_type = "unv"
      elseif endswith(filename,"fro")
         reader_type = "fro"
      elseif endswith(filename,"dat")
         ErrorException("A reader for this mesh is still not implemented")
         #for key in kwargs.keys()
         #   inkey = insensitive(key)
         #   if "connectivity" in inkey and "delimiter" not in inkey
         #      reader_type = "read_separate"
         #      break
         #   end
         #end
      else
         ArgumentError("Mesh file format was not understood. Please specify it using reader_type keyword")
      end
   end

   if reader_type=="salome"
      ErrorException("A reader for this mesh is still not implemented")
      #ReadSalome(filename, element_type=element_type, read_surface_info=read_surface_info)
   elseif reader_type=="GID"
      ErrorException("A reader for this mesh is still not implemented")
      #ReadGIDMesh(filename, element_type, order)
   elseif reader_type=="gmsh"
      ReadGmsh(filename, element_type, read_surface_info=read_surface_info)
   elseif reader_type=="obj"
      ErrorException("A reader for this mesh is still not implemented")
      #ReadOBJ(filename, element_type=element_type, read_surface_info=read_surface_info)
   elseif reader_type=="fenics"
      ErrorException("A reader for this mesh is still not implemented")
      #ReadFenics(filename, element_type)
   elseif reader_type=="vtu"
      ErrorException("A reader for this mesh is still not implemented")
      #ReadVTK(filename)
   elseif reader_type=="unv"
      ErrorException("A reader for this mesh is still not implemented")
      #ReadUNV(filename, element_type)
   elseif reader_type=="fro"
      ErrorException("A reader for this mesh is still not implemented")
      #ReadFRO(filename, element_type)
   elseif reader_type=="read_separate"
      ErrorException("A reader for this mesh is still not implemented")
      # READ MESH FROM SEPARATE FILES FOR CONNECTIVITY AND COORDINATES
      #from Florence.Utils import insensitive
      # return insensitive(kwargs.keys())
      #for key in kwargs.keys()
      #   inkey = insensitive(key)
      #   if "connectivity" in inkey and "delimiter" not in inkey
      #      connectivity_file = kwargs.get(key)
      #   end
      #   if "coordinate" in insensitive(key) and "delimiter" not in inkey
      #      coordinates_file = kwargs.get(key)
      #   end
      #end
      #self.ReadSeparate(connectivity_file,coordinates_file,element_type,
      #   delimiter_connectivity=",",delimiter_coordinates=",")
   elseif reader_type=="ReadHDF5"
      ErrorException("A reader for this mesh is still not implemented")
      #ReadHDF5(filename)
   end

   #self.nnode = self.points.shape[0]
   # MAKE SURE MESH DATA IS CONTIGUOUS
   #self.points   = np.ascontiguousarray(self.points)
   #self.elements = np.ascontiguousarray(self.elements)
   return
end #function

end #module
