using Jotunn
using Test
H2 = create_fragment(coords_string="""
    H 0.0 0.0 0.0
    H 0.0 0.0 0.74
    """, charge=0, mult=1)

@testset "UKS/STO-3G on H2 (sparse4c,default)" begin
    #Simple call
    result= jSCF(H2, "sto-3g"; WFtype="UKS", functional="NoXC", printlevel=0, maxiter=20)
    #UKS/STO-3G with no XC from ORCA
    ref_energy=-0.44200338
    @test isapprox(result["energy"],ref_energy)
end
