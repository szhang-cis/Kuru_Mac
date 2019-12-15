module FEMme

#from .Base import Base
#from .QuadratureRules import QuadratureRule
#from .FunctionSpace import FunctionSpace
#-------------------------------------------------
# Get Mesh modules
include("MeshGeneration/GeometricPath.jl")
include("MeshGeneration/Mesh.jl")
using Mesh
include("MeshGeneration/CustomMesher.jl")
#-------------------------------------------------
#from .MaterialLibrary import *
#from .VariationalPrinciple import *
#from .BoundaryCondition import BoundaryCondition
#from .Solver import *
#from .Utils import PWD, RSWD
#from .PostProcessing import *
#from .FiniteElements import AssembleMass, AssembleForm

greet() = print("Hello World!")

end # module
