using Jotunn

H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)

@testset "RHF/STO-3G on H2O (default)" begin
    #Simple call
    result= jHF(H2O, "sto-3g"; printlevel=0)
    #RHF/STO-3G H2O from ORCA
    ref_energy=-74.964166493793
    @test isapprox(result["energy"],ref_energy)
end

@testset "RHF/def2-SVP on H2O (default)" begin
    #Simple call
    result= jHF(H2O, "def2-svp"; printlevel=0)
    #RHF/STO-3G H2O from ORCA
    ref_energy=-75.960448149155
    @test isapprox(result["energy"],ref_energy)
end

@testset "RHF/def2-TZVPP on H2O (default)" begin
    #Simple call
    result= jHF(H2O, "def2-tzvpp"; printlevel=0)
    #RHF/STO-3G H2O from ORCA
    ref_energy=-76.061779236568
    @test isapprox(result["energy"],ref_energy)
end

@testset "RHF/def2-QZVPP on H2O (default)" begin
    #Simple call
    result= jHF(H2O, "def2-qzvpp"; printlevel=0)
    #RHF/STO-3G H2O from ORCA
    ref_energy=-76.066026165627
    @test isapprox(result["energy"],ref_energy)
end
