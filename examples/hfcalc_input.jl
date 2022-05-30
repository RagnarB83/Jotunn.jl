#using Profile, ProfileView
using Jotunn

##############
#INPUT
##############
#H2 fragment
#H = create_fragment(coords_string="""
#H 0.0 0.0 0.0
#""", charge=0, mult=2)
#H2 fragment
H2 = create_fragment(coords_string="""
H 0.0 0.0 0.0
H 0.0 0.0 0.74
""", charge=0, mult=1)
#O2 fragment
#O2 = create_fragment(coords_string="""
#O 0.0 0.0 0.0
#O 0.0 0.0 1.2
#""", charge=0, mult=3)
#O2 fragment
#NO = create_fragment(coords_string="""
#N 0.0 0.0 0.0
#O 0.0 0.0 1.2
#""", charge=0, mult=2)
#acetaldehyde = create_fragment(xyzfile="acetaldehyde.xyz", charge=0, mult=1)
#Be fragment
#Be = create_fragment(coords_string="""
#Be 0.0 0.0 0.0
#""", charge=0, mult=1)
#Li- fragment
#Li_min = create_fragment(coords_string="""
#Li 0.0 0.0 0.0
#""", charge=-1, mult=1) 
#He fragment
#He = create_fragment(coords_string="""
#He 0.0 0.0 0.0
#""", charge=0, mult=1)
#H2O fragment
#H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)
O3_triplet = create_fragment(xyzfile="o3.xyz", charge=0, mult=3)
#Acetic acid fragment
@time Acetic = create_fragment(xyzfile="acetic.xyz", charge=0, mult=1)

#Call jSCF on fragment with input basis
basisname="sto-3g"

# @time resultA= jSCF(O3_triplet, basisname; maxiter=200, fock_algorithm="turbo", WFtype="UHF",
#     levelshift=false, levelshift_val=0.25, damping=false,tei_type="4c", printlevel=1)
# @time resultX= jSCF(O3_triplet, basisname; maxiter=200, fock_algorithm="turbo", WFtype="UHF",
#     levelshift=true, levelshift_val=0.5, damping=false,tei_type="4c", printlevel=1)
#     @time resultY= jSCF(O3_triplet, basisname; maxiter=200, fock_algorithm="turbo", WFtype="UHF",
#     levelshift=false, damping=true,tei_type="4c", printlevel=1)
#     println("None : ", resultA["finaliter"])
#     println("Levelshift : ", resultX["finaliter"])
#     println("Damping : ", resultY["finaliter"])
# exit()
#@time energy= jSCF(H2O, basisname; maxiter=120, fock_algorithm="loop")#, levelshift=20.0, lshift_thresh=1e-4
@time resultRHF= jSCF(Acetic, basisname; maxiter=200, fock_algorithm="turbo", WFtype="UHF",
    levelshift=false, damping=false,tei_type="4c", printlevel=1)
@time result_l_F_d_F= jSCF(Acetic, basisname; maxiter=200, fock_algorithm="turbo", WFtype="UHF",
    levelshift=false, damping=false,tei_type="4c", printlevel=1)
@time result_l= jSCF(Acetic, basisname; maxiter=200, fock_algorithm="turbo", WFtype="UHF",
levelshift=true, levelshift_val=0.5, lshift_thresh=0.01, damping=false, tei_type="4c", printlevel=1)
@time result_d= jSCF(Acetic, basisname; maxiter=200, fock_algorithm="turbo", WFtype="UHF",
levelshift=false, damping=true, damping_val=0.4, damping_thresh=0.01,tei_type="4c", printlevel=1)
@time result_lf_true= jSCF(Acetic, basisname; maxiter=200, fock_algorithm="turbo", WFtype="UHF",
levelshift=true, damping=true, levelshift_val=0.5, lshift_thresh=0.01, damping_val=0.4, damping_thresh=0.01,tei_type="4c", printlevel=1)
println("Nothing (RHF): ", resultRHF["finaliter"])
println("Nothing : ", result_l_F_d_F["finaliter"])
println("Levelshift : ", result_l["finaliter"])
println("Damping : ", result_d["finaliter"])
println("Damp+Levelshift : ", result_lf_true["finaliter"])
println("---------------------------")
println("Nothing (RHF) : ", resultRHF["energy"])
println("Nothing : ", result_l_F_d_F["energy"])
println("Levelshift : ", result_l["energy"])
println("Damping : ", result_d["energy"])
println("Damp+Levelshift : ", result_lf_true["energy"])
#@time energy= jSCF(H2O, basisname; debugprint=false, maxiter=120, fock_algorithm="tullio")
#@time energy= jSCF(H2O, basisname; debugprint=false, maxiter=120, fock_algorithm="tensor")



#Profile.print()
#ProfileView.view()
#pprof(;webport=58699)

#energy1 = jSCF(elems, coords, charge, mult, basis, basisfile="sto3g-orca.basis")
#energy1 = jSCF(H2, "read_file", basisfile="sto3g-orca.basis")
#energy1b = jSCF(H2, "read_file", basisfile="sto3g-orca-he-normalized.basis")
#energy2= jSCF(H2, "sto-3g")
#energy3= jSCF(H2, "def2-tzvp")
#energy = jSCF(H2, integrals_from_file=true)
