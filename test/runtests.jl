using Jotunn
using Test

@testset verbose=true "Jotunn-tests" begin
    include("test_fragments.jl")
    include("test_scfaids_RHF.jl")
    include("HF-RHF-tests.jl")
    include("HF-UHF-tests.jl")
    include("DFT-RKS-tests.jl")
    include("DFT-UKS-tests.jl")
end
