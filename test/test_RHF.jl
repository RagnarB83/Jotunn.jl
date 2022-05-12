H2 = create_fragment(coords_string="""
    H 0.0 0.0 0.0
    H 0.0 0.0 0.74
    """, charge=0, mult=1)

@testset "RHF/STO-3G on H2" begin
    #Simple call
    result= jHF(H2, "sto-3g")
    #RHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end