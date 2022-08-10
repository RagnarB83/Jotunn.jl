using Jotunn
using Test

#UKS calculations for atoms H-Ne with Slater LDA Exchange.
#Tests Libxc implementation, density generation, grid etc.


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
refE=Dict("refE_H_LDAX_sto3g"=> -0.411379109553, "refE_He_LDAX_sto3g"=> -2.657311972433, "refE_Li_LDAX_sto3g"=> -7.066411545715, 
    "refE_Be_LDAX_sto3g" => -13.988718869619, "refE_B_LDAX_sto3g" => -23.661579845852, "refE_C_LDAX_sto3g" => -36.585336212020,
    "refE_N_LDAX_sto3g" => -52.972280164934, "refE_O_LDAX_sto3g" => -72.907416612476, "refE_F_LDAX_sto3g" => -96.936443159980, 
    "refE_Ne_LDAX_sto3g" => -125.389908806637)

#Threshold for passing energy tast
threshold=1e-6
#Threshold for passing integrated electrons test
elec_threshold=1e-4

system=H
@testset "UKS/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UKS", functional="LDA-X", printlevel=0, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_LDAX_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=He
@testset "UKS/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UKS", functional="LDA-X", printlevel=3, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_LDAX_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=Li
@testset "UKS/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UKS", functional="LDA-X", printlevel=3, maxiter=30, diis=false, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_LDAX_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end


system=Be
@testset "UKS/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UKS", functional="LDA-X", printlevel=0, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_LDAX_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=B
@testset "UKS/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UKS", functional="LDA-X", printlevel=0, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_LDAX_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=C
@testset "UKS/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UKS", functional="LDA-X", printlevel=1, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_LDAX_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=N
@testset "UKS/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UKS", functional="LDA-X", printlevel=0, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_LDAX_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=O
@testset "UKS/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UKS", functional="LDA-X", printlevel=0, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_LDAX_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=F
@testset "UKS/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UKS", functional="LDA-X", printlevel=0, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_LDAX_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

system=Ne
@testset "UKS/STO-3G (LDA-X via libxc) on $(system.formula)" begin
    result= jSCF(system, "sto-3g"; WFtype="UKS", functional="LDA-X", printlevel=0, maxiter=20, diis=true, damping=false, levelshift=false)
    ref_energy=refE["refE_$(system.formula)_LDAX_sto3g"]
    ref_num_electrons=system.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end
