using Jotunn
using Test


@testset verbose=true "HF-RHF-tests" begin
    include("test_RHF_H2.jl")
    include("test_RHF_H2O.jl")
    include("test_RHF_He.jl")
end

