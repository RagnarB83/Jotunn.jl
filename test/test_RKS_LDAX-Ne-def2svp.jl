using Jotunn
using Test

Ne = create_fragment(coords_string="""Ne 0.0 0.0 0.0""", charge=0, mult=1)

refE_Ne_LDAX_def2svp=-127.306729528465

#Threshold for passing energy tast
threshold=1e-6
#Threshold for passing integrated electrons test
elec_threshold=1e-4


#Needs d-functions

@testset "RKS/def2-SVP (LDA-X via libxc) on Ne" begin
    #Simple call
    result= jSCF(Ne, "def2-svp"; WFtype="RKS", functional="LDA-X", printlevel=3, maxiter=20, diis=true, damping=false, levelshift=false)
    #RKS/def2-SVP with LDA-X from ORCA
    ref_energy=refE_Ne_LDAX_def2svp
    ref_num_electrons=Ne.numelectrons
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end
