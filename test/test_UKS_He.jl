using Jotunn
using Test

He = create_fragment(coords_string="""He 0.0 0.0 0.0""", charge=0, mult=1)

#Threshold for passing energy tast
threshold=1e-6
#Threshold for passing integrated electrons test
elec_threshold=1e-6

system=He
@testset "UKS/STO-3G (LDA-X manual) on $(system.formula)" begin
    #Note: Regualar LDA
    result= jSCF(system, "sto-3g"; WFtype="UKS", functional="manual", manual_func="LDA",
        printlevel=3, maxiter=20, diis=true, damping=false, levelshift=false)
    #UKS/STO-3G with LDA-X from ORCA
    ref_energy=-2.657311972433
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=He
@testset "UKS/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    #Note: Regular LDA
    result= jSCF(system, "sto-3g"; WFtype="UKS", functional="LDA-X",
        printlevel=3, maxiter=20, diis=true, damping=false, levelshift=false)
    #UKS/STO-3G with LDA-X from ORCA
    ref_energy=-2.657311972433
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end
