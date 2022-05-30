using Jotunn
using Test

#This set of tests uses 1 system, 1 WFtype and 1 basis set for everything but then different jSCF parameters
#Comparison to ORCA energy

#System
He = create_fragment(coords_string="""
    He 0.0 0.0 0.0
    """, charge=0, mult=1)

#Test global settings
mol=He
WFtype="RHF"
basis1="sto-3g"
#Reference energy: RHF/STO-3G He energy from ORCA
ref_energy=-2.807783957535

#Actual tests with different parameters
@testset "$WFtype/$basis1 on He (sparse4c,default)" begin
    #Simple call
    result= jSCF(mol, basis1; printlevel=0)
    @test isapprox(result["energy"],ref_energy)
end

@testset "$WFtype/$basis1 on He (4c)" begin
    #Simple call
    result= jSCF(mol, basis1; printlevel=0, tei_type="4c")
    @test isapprox(result["energy"],ref_energy)
end

@testset "$WFtype/$basis1 on He (no SCF-aids)" begin
    #Simple call
    result= jSCF(mol, basis1; printlevel=0, diis=false, damping=false, levelshift=false)
    @test isapprox(result["energy"],ref_energy)
end

@testset "$WFtype/$basis1 on He (diis)" begin
    #Simple call
    result= jSCF(mol, basis1; printlevel=0, diis=true, damping=false, levelshift=false)
    @test isapprox(result["energy"],ref_energy)
end

@testset "$WFtype/$basis1 on He (damping)" begin
    #Simple call
    result= jSCF(mol, basis1; printlevel=0, diis=false, damping=true, levelshift=false)
    @test isapprox(result["energy"],ref_energy)
end

@testset "$WFtype/$basis1 on He (levelshift)" begin
    #Simple call
    result= jSCF(mol, basis1; printlevel=0, diis=false, damping=false, levelshift=true)
    @test isapprox(result["energy"],ref_energy)
end

@testset "$WFtype/$basis1 on He (diis+damp)" begin
    #Simple call
    result= jSCF(mol, basis1; printlevel=0, diis=true, damping=true, levelshift=false)
    @test isapprox(result["energy"],ref_energy)
end

@testset "$WFtype/$basis1 on He (diis+levelshift)" begin
    #Simple call
    result= jSCF(mol, basis1; printlevel=0, diis=true, damping=false, levelshift=true)
    @test isapprox(result["energy"],ref_energy)
end

@testset "$WFtype/$basis1 on He (damping+levelshift)" begin
    #Simple call
    result= jSCF(mol, basis1; printlevel=0, diis=false, damping=true, levelshift=true)
    @test isapprox(result["energy"],ref_energy)
end

@testset "$WFtype/$basis1 on He (diis+damping+levelshift)" begin
   #Simple call
    result= jSCF(mol, basis1; printlevel=0, diis=true, damping=true, levelshift=true)
    @test isapprox(result["energy"],ref_energy)
end
