using Jotunn
using Test
H2 = create_fragment(coords_string="""
    H 0.0 0.0 0.0
    H 0.0 0.0 0.74
    """, charge=0, mult=1)

@testset "RHF/STO-3G on H2 (sparse4c,default)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; printlevel=0)
    #RHF/STO-3G H2 (r=0.74) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "RHF/STO-3G on H2 (4c)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; printlevel=0, tei_type="4c")
    #RHF/STO-3G H2 (r=0.74) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "RHF/STO-3G on H2 (no SCF-aids)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; printlevel=0, diis=false, damping=false, levelshift=false)
    #RHF/STO-3G H2 (r=0.74) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "RHF/STO-3G on H2 (diis)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; printlevel=0, diis=true, damping=false, levelshift=false)
    #RHF/STO-3G H2 (r=0.74) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "RHF/STO-3G on H2 (damping)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; printlevel=0, diis=false, damping=true, levelshift=false)
    #RHF/STO-3G H2 (r=0.74) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "RHF/STO-3G on H2 (levelshift)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; printlevel=0, diis=false, damping=false, levelshift=true)
    #RHF/STO-3G H2 (r=0.74) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "RHF/STO-3G on H2 (diis+damp)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; printlevel=0, diis=true, damping=true, levelshift=false)
    #RHF/STO-3G H2 (r=0.74) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "RHF/STO-3G on H2 (diis+levelshift)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; printlevel=0, diis=true, damping=false, levelshift=true)
    #RHF/STO-3G H2 (r=0.74) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "RHF/STO-3G on H2 (damping+levelshift)" begin
    #Simple call
    result= jHF(H2, "sto-3g"; printlevel=0, diis=false, damping=true, levelshift=true)
    #RHF/STO-3G H2 (r=0.74) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "RHF/STO-3G on H2 (diis+damping+levelshift)" begin
   #Simple call
    result= jHF(H2, "sto-3g"; printlevel=0, diis=true, damping=true, levelshift=true)
    #RHF/STO-3G H2 (r=0.74) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end
