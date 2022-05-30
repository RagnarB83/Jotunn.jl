"""
Jotunn: a simple QM program in Julia by Ragnar Bjornsson
Integrals via the GaussianBasis.jl interface to Libcint
"""
module Jotunn

#From standard library
using LinearAlgebra
using Statistics
#Extra packages
using Printf
using DelimitedFiles
using PrettyTables
using Crayons
using Molecules #FermiQC
using GaussianBasis #FermiQC
using PyCall

#Linear Algebra speed-up packages
#using LoopVectorization
#using TensorOperations
#using Tullio


#Source-code files
include("fragment.jl")
include("jSCF.jl")
include("Fock.jl")

include("scf_convergers.jl")
include("basis_sets_and_integrals.jl")
include("mol_properties.jl")
include("pop_ana.jl")
include("printing.jl")
include("tools.jl")
include("num_integration.jl")

#Export the callable functions of Jotunn
export jSCF, create_fragment

end
