H2 = create_fragment(coords_string="""
    H 0.0 0.0 0.0
    H 0.0 0.0 0.74
    """, charge=0, mult=1)

@testset "RHF/STO-3G on H2 (no SCF aids)" begin
    #Simple no SCF aids
    result= jSCF(H2, "sto-3g", diis=false, levelshift=false, damping=false)
    #RHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "RHF/STO-3G on H2 (Damping)" begin
    #Simple call
    result= jSCF(H2, "sto-3g", diis=false, levelshift=false, damping=true)
    #RHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end


@testset "RHF/STO-3G on H2 (Levelshift)" begin
    #Simple call
    result= jSCF(H2, "sto-3g", diis=false, levelshift=true, damping=false)
    #RHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "RHF/STO-3G on H2 (DIIS)" begin
    #Simple call
    result= jSCF(H2, "sto-3g", diis=true, levelshift=false, damping=false)
    #RHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end


