using Jotunn
using Test
Be = create_fragment(coords_string="""Be 0.0 0.0 0.0""", charge=0, mult=1)

ref_num_electrons=4
ref_energy_NoXC_sto3g=-11.583813202694
ref_energy_LDAX_sto3g=-13.98871887
ref_energy_LDAX_def2svp=-14.204253768039
ref_energy_LDAX_def2tzvpp=-14.221657957506
#ref_energy_LDAX_def2tzvpp=-1.043016819359

#Threshold for passing energy tast
threshold=1e-6
#Threshold for passing integrated electrons test
elec_threshold=1e-4


@testset "RKS/STO-3G (NoXC) on Be" begin
    #Simple call
    result= jSCF(Be, "sto-3g"; WFtype="RKS", functional="NoXC", printlevel=0, maxiter=20)
    #RKS/STO-3G with LDA-Exchange only from ORCA
    ref_energy=ref_energy_NoXC_sto3g
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

@testset "RKS/STO-3G (LDA-X) on Be" begin
    #Simple call
    result= jSCF(Be, "sto-3g"; WFtype="RKS", functional="LDA-X", printlevel=1, maxiter=20)
    #RKS/STO-3G with LDA-Exchange only from ORCA
    ref_energy=ref_energy_LDAX_sto3g
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

@testset "RKS/def2-SVP (LDA-X) on Be" begin
    #Simple call
    result= jSCF(Be, "def2-svp"; WFtype="RKS", functional="LDA-X", printlevel=1, maxiter=20)
    #RKS/STO-3G with LDA-Exchange only from ORCA
    ref_energy=ref_energy_LDAX_def2svp
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

#@testset "RKS/def2-TZVPP (LDA-X) on Be" begin
#    #Simple call
#    result= jSCF(Be, "def2-tzvpp"; WFtype="RKS", functional="LDA-X", printlevel=1, maxiter=20)
#    #RKS/STO-3G with LDA-Exchange only from ORCA
#    ref_energy=ref_energy_LDAX_def2tzvpp
#    @test abs(result["energy"] - ref_energy) <  threshold
#    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
#end
