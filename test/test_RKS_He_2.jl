using Jotunn
using Test
He2 = create_fragment(coords_string="""
    He 0.0 0.0 0.0
    He 0.0 0.0 2.0
    """, charge=0, mult=1)

ref_num_electrons=4
ref_energy_LDAX_def2svp=-5.428816485742

#Threshold for passing energy tast
threshold=1e-5
#Threshold for passing integrated electrons test
elec_threshold=1e-4


@testset "RKS/def2-SVP (LibXC LDA-X) on He2" begin
    #Simple call
    result= jSCF(He2, "def2-SVP"; WFtype="RKS", functional="LDA-X", printlevel=1, maxiter=30,
        diis=true, diis_size=10, damping=true,levelshift=false)
    #RKS/def2-SVP with LDA-Exchange only from ORCA
    ref_energy=ref_energy_LDAX_def2svp
    @test abs(result["energy"] - ref_energy) <  threshold
    @test abs(result["int_electrons"] - ref_num_electrons) <  elec_threshold
end

