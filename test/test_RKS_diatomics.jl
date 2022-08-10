using Jotunn
using Test
HF = create_fragment(coords_string="""
    H 0.0 0.0 0.0
    F 0.0 0.0 0.91
    """, charge=0, mult=1)
LiH = create_fragment(coords_string="""
    Li 0.0 0.0 0.0
    H 0.0 0.0 1.60
    """, charge=0, mult=1)

ref_num_electrons=2
ref_energy_LDAX_sto3g_HF=-97.529697247330
ref_energy_LDAX_sto3g_LiH=-7.573496109546

#Threshold for passing energy tast
threshold=1e-5
#Threshold for passing integrated electrons test
elec_threshold=1e-4

#@testset "RKS/STO-3G (LibXC LDA-X) on LiH" begin
#    #Simple call
#    result= jSCF(LiH, "sto-3g"; WFtype="RKS", functional="LDA-X", printlevel=1, maxiter=60)
#    #RKS/def2-SVP with LDA-Exchange only from ORCA
#    ref_energy=ref_energy_LDAX_sto3g_LiH
#    ref_num_electrons=LiH.numelectrons
#    @test abs(result["energy"] - ref_energy) <  threshold
#    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
#end

@testset "RKS/STO-3G (LibXC LDA-X) on HF" begin
    #Simple call
    result= jSCF(HF, "sto-3g"; WFtype="RKS", functional="LDA-X", printlevel=2, maxiter=160,
        diis=true, diis_startiter=4, diis_size=10, damping=true, levelshift=true)
    #RKS/def2-SVP with LDA-Exchange only from ORCA
    ref_energy=ref_energy_LDAX_sto3g_HF
    ref_num_electrons=HF.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end
