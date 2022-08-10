using Jotunn
using Test
Zn = create_fragment(coords_string="""Zn 0.0 0.0 0.0""", charge=0, mult=1)

ref_num_electrons=Zn.numelectrons
ref_energy_NoXC_sto3g=-1688.690895164397
ref_energy_LDAX_sto3g=-1752.687861807158
ref_energy_LDAX_def2svp=-1773.565717853807
ref_energy_LDAX_def2tzvpp=-1773.846407267273
#ref_energy_LDAX_def2tzvpp=-1.043016819359

#Threshold for passing energy tast
threshold=1e-6
#Threshold for passing integrated electrons test
elec_threshold=1e-4

#Note: Zn/STO-3G is practically impossible to converge (also by ORCA). Skipping

#Uses f-basis functions but will still pass
@testset "RKS/def2-SVP (LDA-X) on Zn (angmom=f)" begin
    #Simple call
    result= jSCF(Zn, "def2-svp"; WFtype="RKS", functional="LDA-X", printlevel=1, maxiter=60, diis=true, levelshift=true, damping=true)
    #RKS/STO-3G with LDA-Exchange only from ORCA
    ref_energy=ref_energy_LDAX_def2svp
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end


#Big job:
#@testset "RKS/def2-TZVPP (LDA-X) on Zn" begin
#    #Simple call
#    result= jSCF(Zn, "def2-tzvpp"; WFtype="RKS", functional="LDA-X", printlevel=1, maxiter=60)
#    #RKS/STO-3G with LDA-Exchange only from ORCA
#    ref_energy=ref_energy_LDAX_def2tzvpp
#    @test abs(result["energy"] - ref_energy) <  threshold
#    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
#end
