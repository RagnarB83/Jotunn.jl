using Jotunn
using Test

He = create_fragment(coords_string="""He 0.0 0.0 0.0""", charge=0, mult=1)
Ne = create_fragment(coords_string="""Ne 0.0 0.0 0.0""", charge=0, mult=1)
Ar = create_fragment(coords_string="""Ar 0.0 0.0 0.0""", charge=0, mult=1)
Kr = create_fragment(coords_string="""Kr 0.0 0.0 0.0""", charge=0, mult=1)

refE_He_LDAX_sto3g=-2.657311972433
refE_Ne_LDAX_sto3g=-125.389908806637
refE_Ar_LDAX_sto3g=-518.750424668477
refE_Kr_LDAX_sto3g=-2716.835571339233

#Threshold for passing energy tast
threshold=1e-6
#Threshold for passing integrated electrons test
elec_threshold=1e-4

@testset "RKS/STO-3G (LDA-X via libxc) on He" begin
    #Simple call
    result= jSCF(He, "sto-3g"; WFtype="RKS", functional="LDA-X", printlevel=0, maxiter=20)
    #RKS/STO-3G with LDA-X from ORCA
    ref_energy=refE_He_LDAX_sto3g
    ref_num_electrons=He.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

@testset "RKS/STO-3G (NoXC) on Ne" begin
    #Simple call
    result= jSCF(Ne, "sto-3g"; WFtype="RKS", functional="NoXC", printlevel=0, maxiter=20)
    #RKS/STO-3G with LDA-X from ORCA
    ref_energy=-114.147356142907
    ref_num_electrons=Ne.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

@testset "RKS/STO-3G (LDA-X via libxc) on Ne" begin
    #Simple call
    result= jSCF(Ne, "sto-3g"; WFtype="RKS", functional="LDA-X", printlevel=0, maxiter=20, diis=false, damping=false, levelshift=false)
    #RKS/STO-3G with LDA-X from ORCA
    ref_energy=refE_Ne_LDAX_sto3g
    ref_num_electrons=Ne.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

@testset "RKS/STO-3G (LDA-X via libxc) on Ar" begin
    #Simple call
    result= jSCF(Ar, "sto-3g"; WFtype="RKS", functional="LDA-X", printlevel=0, maxiter=20)
    #RKS/STO-3G with LDA-X from ORCA
    ref_energy=refE_Ar_LDAX_sto3g
    ref_num_electrons=Ar.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

#Needs d-functions

@testset "RKS/STO-3G (LDA-X via libxc) on Kr" begin
    #Simple call
    result= jSCF(Kr, "sto-3g"; WFtype="RKS", functional="LDA-X", printlevel=3, maxiter=50, diis=true, damping=true)
    #RKS/STO-3G with LDA-X from ORCA
    ref_energy=refE_Kr_LDAX_sto3g
    ref_num_electrons=Kr.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

