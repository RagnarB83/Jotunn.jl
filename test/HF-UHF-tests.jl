using Jotunn
using Test


@testset verbose=true "HF-UHF-tests" begin
    include("test_UHF_H2.jl")
    include("test_UHF_NO.jl")
    include("test_UHF-atoms.jl")
end

