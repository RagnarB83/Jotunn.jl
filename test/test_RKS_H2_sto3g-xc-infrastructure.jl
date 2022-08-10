using Jotunn
using Test
H2 = create_fragment(coords_string="""
    H 0.0 0.0 0.0
    H 0.0 0.0 0.74
    """, charge=0, mult=1)

ref_num_electrons=2
ref_energy_LDAX_sto3g=-1.025008086952

#Threshold for passing energy tast
threshold=1e-5
#Threshold for passing integrated electrons test
elec_threshold=1e-4

@testset "RKS/STO-3G (LDA-X via libxc-keyword) on H2" begin
    #Simple call
    result= jSCF(H2, "sto-3g"; WFtype="RKS", libxc_keyword="lda_x", printlevel=0, maxiter=20)
    #RKS/STO-3G with LDA-X from ORCA
    ref_energy=ref_energy_LDAX_sto3g
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

@testset "RKS/STO-3G (No XC approximation) on H2" begin
    #Simple call
    result= jSCF(H2, "sto-3g"; WFtype="RKS", functional="NoXC", printlevel=0, maxiter=20)
    #RKS/STO-3G with no XC from ORCA
    ref_energy=-0.44200338
    #@test isapprox(result["energy"],ref_energy)
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

@testset "RKS/STO-3G (LibXC LDA-Exchange functional) on H2" begin
    #Simple call
    result= jSCF(H2, "sto-3g"; WFtype="RKS", functional="LDA-X(libxc)", printlevel=0, maxiter=20)
    #RKS/STO-3G with LDA-Exchange only from ORCA
    ref_energy=-1.025008086952
    #@test isapprox(result["energy"],ref_energy)
    #Grid sensitivity means isapprox no longer works
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

@testset "RKS/STO-3G (Manual LDA-Exchange functional) on H2" begin
    #Simple call
    result= jSCF(H2, "sto-3g"; WFtype="RKS", functional="manual", manual_func="LDA", printlevel=0, maxiter=20)
    #RKS/STO-3G with LDA-Exchange only from ORCA
    ref_energy=-1.025008086952
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end
