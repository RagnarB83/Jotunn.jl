using Jotunn
using Test

H2 = create_fragment(coords_string="""
    H 0.0 0.0 0.0
    H 0.0 0.0 0.74
    """, charge=0, mult=1)

@testset "RHF/STO-3G on H2 (printlevel 0)" begin
    #Simple call
    result= jHF(H2, "sto-3g", printlevel=0)
    #RHF/STO-3G H2 (r=0.74 A) from ORCA
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

#@testset "RHF/STO-3G on H2 (printlevel 1)" begin
#    #Simple call
#    result= jHF(H2, "sto-3g", printlevel=1)
#    #RHF/STO-3G H2 (r=0.74 A) from ORCA
#    ref_energy=-1.116759307204
#    @test isapprox(result["energy"],ref_energy)
#end

#@testset "RHF/STO-3G on H2 (printlevel 2)" begin
#    #Simple call
#    result= jHF(H2, "sto-3g", printlevel=2)
#    #RHF/STO-3G H2 (r=0.74 A) from ORCA
#    ref_energy=-1.116759307204
#    @test isapprox(result["energy"],ref_energy)
#end

#@testset "RHF/STO-3G on H2 (printlevel 3)" begin
#    #Simple call
#    result= jHF(H2, "sto-3g", printlevel=3)
#    #RHF/STO-3G H2 (r=0.74 A) from ORCA
#    ref_energy=-1.116759307204
#    @test isapprox(result["energy"],ref_energy)
#end