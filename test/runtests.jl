using Jotunn
using Test

@testset "Jotunn-tests" begin
    include("test_fragments.jl")
    include("test_RHF.jl")
end