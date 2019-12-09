# Call FEMme library
push!(LOAD_PATH, homedir() * "/femme/src/")
import FEMme: Mesh
# Function of Homogenized Constrained Mixture Theory in explicit scheme
function homogenized_CMT_explicit()
   """A hyperelastic implicit static example using ArterialMixture model
      of a cylinder under pression with hexahedral elements
   """
   println("Hello FEMme")
   # Mesh processing
   ProblemPath = @__DIR__
   mesh_file = ProblemPath * "/Quarter_Ring.msh"
   println(mesh_file)
   Mesh.Read(filename=mesh_file, reader_type="gmsh", element_type="hex", read_surface_info=true)
end

homogenized_CMT_explicit()
