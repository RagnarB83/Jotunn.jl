using Jotunn
using Test

@testset verbose=true "DFT-UKS-tests" begin
    include("test_UKS_H2.jl")
    include("test_UKS_He.jl")
    include("test_UKS_LDAX-atoms.jl")
end

