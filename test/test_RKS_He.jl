using Jotunn
using Test

He = create_fragment(coords_string="""
    He 0.0 0.0 0.0
    """, charge=0, mult=1)

ref_num_electrons=He.numelectrons

ref_energy_LDAX_sto3g=-2.657311972433
ref_energy_LDAX_def2svp=-2.714655827872
ref_energy_LDAX_def2tzvpp=-2.721023564381
ref_energy_LDAX_def2qzvpp=-2.723479148402

#Threshold for passing energy tast
threshold=1e-6
#Threshold for passing integrated electrons test
elec_threshold=1e-6

@testset "RKS/sto-3g (LDA-X via libxc) on He" begin
    #Simple call
    result= jSCF(He, "sto-3g"; WFtype="RKS", functional="LDA-X", printlevel=0, maxiter=20)
    #RKS/def2-TZVPP with LDA-X from ORCA
    ref_energy=ref_energy_LDAX_sto3g
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end


@testset "RKS/def2-svp (LDA-X via libxc) on He" begin
    #Simple call
    result= jSCF(He, "def2-svp"; WFtype="RKS", functional="LDA-X", printlevel=0, maxiter=20)
    #RKS/def2-TZVPP with LDA-X from ORCA
    ref_energy=ref_energy_LDAX_def2svp
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

#Below needs d-functions

@testset "RKS/def2-tzvpp (LDA-X via libxc) on He" begin
    #Simple call
    result= jSCF(He, "def2-tzvpp"; WFtype="RKS", functional="LDA-X", printlevel=0, maxiter=20)
    #RKS/def2-TZVPP with LDA-X from ORCA
    ref_energy=ref_energy_LDAX_def2tzvpp
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end


#Needs f-functions below

#@testset "RKS/def2-qzvpp (LDA-X via libxc) on He" begin
#    #Simple call
#    result= jSCF(He, "def2-qzvpp"; WFtype="RKS", functional="LDA-X", printlevel=0, maxiter=20)
#    #RKS/def2-TZVPP with LDA-X from ORCA
#    ref_energy=ref_energy_LDAX_def2qzvpp
#    @test abs(result["energy"] - ref_energy) <  threshold
#    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
#end

