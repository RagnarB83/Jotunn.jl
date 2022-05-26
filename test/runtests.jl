using Jotunn
using Test

@testset verbose=true "Jotunn-tests" begin
    include("test_fragments.jl")
    include("test_RHF_H2.jl")
    include("test_UHF_H2.jl")
    include("test_UHF_NO.jl")
    include("test_RHF_H2O.jl")
end
