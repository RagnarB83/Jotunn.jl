H2 = create_fragment(coords_string="""
    H 0.0 0.0 0.0
    H 0.0 0.0 0.74
    """, charge=0, mult=1)


@testset "UHF/STO-3G on H2 (sparse4c,default)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; HFtype="UHF", printlevel=0)
    #UHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "UHF/STO-3G on H2 (4c)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; HFtype="UHF", printlevel=0, tei_type="4c")
    #UHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "UHF/STO-3G on H2 (no SCF-aids)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; HFtype="UHF", printlevel=0, diis=false, damping=false, levelshift=false)
    #UHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "UHF/STO-3G on H2 (diis)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; HFtype="UHF", printlevel=0, diis=true, damping=false, levelshift=false)
    #UHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "UHF/STO-3G on H2 (damping)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; HFtype="UHF", printlevel=0, diis=false, damping=true, levelshift=false)
    #UHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "UHF/STO-3G on H2 (levelshift)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; HFtype="UHF", printlevel=0, diis=false, damping=false, levelshift=true)
    #UHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "UHF/STO-3G on H2 (diis+damp)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; HFtype="UHF", printlevel=0, diis=true, damping=true, levelshift=false)
    #UHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "UHF/STO-3G on H2 (diis+levelshift)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; HFtype="UHF", printlevel=0, diis=true, damping=false, levelshift=true)
    #UHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "UHF/STO-3G on H2 (damping+levelshift)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; HFtype="UHF", printlevel=0, diis=false, damping=true, levelshift=true)
    #UHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "UHF/STO-3G on H2 (diis+damping+levelshift)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; HFtype="UHF", printlevel=0, diis=true, damping=true, levelshift=true)
    #UHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end
