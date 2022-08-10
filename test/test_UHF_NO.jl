using Jotunn
using Test
    
NO = create_fragment(coords_string="""
N 0.0 0.0 0.0
O 0.0 0.0 1.15
""", charge=0, mult=2)


@testset "UHF/STO-3g on NO (4c)" begin
    #UHF on a doublet radical
    result= jSCF(NO, "sto-3g"; WFtype="UHF", printlevel=0, maxiter=30, tei_type="4c")
    #UHF/STO-3G NO (r=0.1.15 A) from ORCA 5.0.3
    ref_energy=-127.530134374969
    @test isapprox(result["energy"],ref_energy)
end


@testset "UHF/STO-3g on NO (sparse4c)" begin
    #UHF on a doublet radical
    result= jSCF(NO, "sto-3g"; WFtype="UHF", printlevel=0, maxiter=30, tei_type="sparse4c")
    #UHF/STO-3G NO (r=0.1.15 A) from ORCA 5.0.3
    ref_energy=-127.530134374969
    @test isapprox(result["energy"],ref_energy)
end
