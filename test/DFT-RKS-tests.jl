using Jotunn
using Test

@testset verbose=true "DFT-RKS tests" begin
    include("test_RKS_He.jl")
    include("test_RKS_H2_sto3g-xc-infrastructure.jl")
    include("test_RKS_H2.jl")
    include("test_RKS_LDAX-sto3g-raregas.jl")
    include("test_RKS_Be.jl")
    include("test_RKS_C_singlet.jl")
    include("test_RKS_Zn.jl")
    include("test_RKS_LDAX-Ne-def2svp.jl")    
end

