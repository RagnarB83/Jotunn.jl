"""
Jotunn: a simple QM program in Julia by Ragnar Bjornsson
Integrals via the GaussianBasis.jl interface to Libcint
"""
module Jotunn

using LinearAlgebra
using Statistics
using Printf
using DelimitedFiles
using PrettyTables
using Crayons

#Linear Algebra speed-up packages
using LoopVectorization
#using TensorOperations
#using Tullio

#From FermiQC
using Molecules
using GaussianBasis

#Source-code files
include("fragment.jl")
include("jHF.jl")
include("Fock.jl")
include("basis_sets_and_integrals.jl")
include("scf_convergers.jl")
include("mol_properties.jl")
include("pop_ana.jl")
include("printing.jl")
include("tools.jl")

#Export the callable functions of Jotunn
export jHF, create_fragment

end
