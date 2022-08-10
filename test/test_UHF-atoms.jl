using Jotunn
using Test

#UHF calculations for atoms H-Ne 
# Checking energies against ORCA: tests SCF, UHF part

#Defining atoms
H = create_fragment(coords_string="""H 0.0 0.0 0.0""", charge=0, mult=2)
He = create_fragment(coords_string="""He 0.0 0.0 0.0""", charge=0, mult=1)
Li = create_fragment(coords_string="""Li 0.0 0.0 0.0""", charge=0, mult=2)
Be = create_fragment(coords_string="""Be 0.0 0.0 0.0""", charge=0, mult=1)
B = create_fragment(coords_string="""B 0.0 0.0 0.0""", charge=0, mult=2)
C = create_fragment(coords_string="""C 0.0 0.0 0.0""", charge=0, mult=3)
N = create_fragment(coords_string="""N 0.0 0.0 0.0""", charge=0, mult=4)
O = create_fragment(coords_string="""O 0.0 0.0 0.0""", charge=0, mult=3)
F = create_fragment(coords_string="""F 0.0 0.0 0.0""", charge=0, mult=2)
Ne = create_fragment(coords_string="""Ne 0.0 0.0 0.0""", charge=0, mult=1)

#Reference energy dict
refE=Dict("refE_H_UHF_sto3g"=> -0.466581849556, "refE_He_UHF_sto3g"=> -2.807783957535, "refE_Li_UHF_sto3g"=> -7.315525981228, 
    "refE_Be_UHF_sto3g" => -14.351880476179, "refE_B_UHF_sto3g" => -24.148988598815, "refE_C_UHF_sto3g" => -37.198392563675,
    "refE_N_UHF_sto3g" => -53.719010162515, "refE_O_UHF_sto3g" => -73.804150233151, "refE_F_UHF_sto3g" => -97.986504958602, 
    "refE_Ne_UHF_sto3g" => -126.604524996637)

#Threshold for passing energy tast
threshold=1e-6
#Threshold for passing integrated electrons test
elec_threshold=1e-4

system=H
@testset "UHF/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UHF",  printlevel=1, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_UHF_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    #@test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=He
@testset "UHF/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UHF",  printlevel=1, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_UHF_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    #@test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=Li
@testset "UHF/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UHF",  printlevel=1, maxiter=20, diis=false, damping=true, levelshift=true)
    ref_energy=refE["refE_$(system.formula)_UHF_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    #@test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=Be
@testset "UHF/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UHF",  printlevel=0, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_UHF_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    #@test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=B
@testset "UHF/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UHF",  printlevel=0, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_UHF_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    #@test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=C
@testset "UHF/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UHF",  printlevel=1, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_UHF_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    #@test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=N
@testset "UHF/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UHF",  printlevel=0, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_UHF_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    #@test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=O
@testset "UHF/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UHF",  printlevel=0, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_UHF_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    #@test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=F
@testset "UHF/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UHF",  printlevel=0, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_UHF_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    #@test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=Ne
@testset "UHF/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UHF",  printlevel=0, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_UHF_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    #@test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end
