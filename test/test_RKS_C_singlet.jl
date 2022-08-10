using Jotunn
using Test
C_s = create_fragment(coords_string="""C 0.0 0.0 0.0""", charge=0, mult=1)

ref_num_electrons=C_s.numelectrons
ref_energy_NoXC_sto3g=-31.872547699980
ref_energy_LDAX_sto3g=-36.471990350874
ref_energy_LDAX_def2svp=-36.967274927888
ref_energy_LDAX_def2tzvpp=-37.015239223844
#ref_energy_LDAX_def2tzvpp=-1.043016819359

#Threshold for passing energy tast
threshold=1e-6
#Threshold for passing integrated electrons test
elec_threshold=1e-4


@testset "RKS/STO-3G (NoXC) on C_s" begin
    #Simple call
    result= jSCF(C_s, "sto-3g"; WFtype="RKS", functional="NoXC", printlevel=1, maxiter=20)
    #RKS/STO-3G with LDA-Exchange only from ORCA
    ref_energy=ref_energy_NoXC_sto3g
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end


#NOTE: strange orbital energies here
@testset "RKS/STO-3G (LDA-X) on C_s" begin
    #Simple call
    result= jSCF(C_s, "sto-3g"; WFtype="RKS", functional="LDA-X", printlevel=1, maxiter=20)
    #RKS/STO-3G with LDA-Exchange only from ORCA
    ref_energy=ref_energy_LDAX_sto3g
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

