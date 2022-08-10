using Jotunn
using Test
H2 = create_fragment(coords_string="""
    H 0.0 0.0 0.0
    H 0.0 0.0 0.74
    """, charge=0, mult=1)

ref_num_electrons=2
ref_energy_LDAX_sto3g=-1.025008086952
ref_energy_LDAX_def2svp=-1.037405390454
ref_energy_NOXC_def2svp=-0.525123090925
#ref_energy_LDAX_def2tzvpp=-1.043016819359

#Threshold for passing energy tast
threshold=1e-5
#Threshold for passing integrated electrons test
elec_threshold=1e-4


#@testset "RKS/STO-3G (LibXC LDA-X) on H2" begin
#    #Simple call
#    result= jSCF(H2, "sto-3g"; WFtype="RKS", functional="LDA-X", printlevel=0, maxiter=20)
#    #RKS/STO-3G with LDA-Exchange only from ORCA
#    ref_energy=ref_energy_LDAX_sto3g
#    @test abs(result["energy"] - ref_energy) <  threshold
#    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
#end

#@testset "RKS/def2-SVP (LibXC NoXC) on H2" begin
#    #Simple call
#    result= jSCF(H2, "def2-SVP"; WFtype="RKS", functional="NoXC", printlevel=1, maxiter=20,
#        diis=true, diis_size=10, damping=true,levelshift=false)
#    #RKS/def2-SVP with LDA-Exchange only from ORCA
#    ref_energy=ref_energy_NOXC_def2svp
#    @test abs(result["energy"] - ref_energy) <  threshold
#    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
#end

@testset "RKS/def2-SVP (LibXC LDA-X) on H2" begin
    #Simple call
    result= jSCF(H2, "def2-SVP"; WFtype="RKS", functional="LDA-X", printlevel=1, maxiter=30,
        diis=true, diis_size=10, damping=true,levelshift=false)
    #RKS/def2-SVP with LDA-Exchange only from ORCA
    ref_energy=ref_energy_LDAX_def2svp
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

#@testset "RKS/def2-TZVPP (LibXC LDA-Exchange functional) on H2" begin
#   #Simple call
#    result= jSCF(H2, "def2-TZVPP"; WFtype="RKS", functional="LDA-X", printlevel=0, maxiter=40)
#   #RKS/def2-SVP with LDA-Exchange only from ORCA
#    ref_energy=ref_energy_LDAX_def2tzvpp
#    @test abs(result["energy"] - ref_energy) <  threshold
#    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
#end
