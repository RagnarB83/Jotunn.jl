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
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)
O3_triplet = create_fragment(xyzfile="o3.xyz", charge=0, mult=3)
#Acetic acid fragment
#@time Acetic = create_fragment(xyzfile="acetic.xyz", charge=0, mult=1)

#Call jHF on fragment with input basis
basisname="STO-3g"

#@time energy= jHF(H2O, basisname; maxiter=120, fock_algorithm="loop")#, levelshift=20.0, lshift_thresh=1e-4
@time energy= jHF(H2O, basisname; maxiter=500, fock_algorithm="turbo", HFtype="UHF", levelshift=2.0, lshift_thresh=1e-4,
    tei_type="4c", print_final_matrices=true, debugprint=false, printlevel=1)
#@time energy= jHF(H2O, basisname; debugprint=false, maxiter=120, fock_algorithm="tullio")
#@time energy= jHF(H2O, basisname; debugprint=false, maxiter=120, fock_algorithm="tensor")



#Profile.print()
#ProfileView.view()
#pprof(;webport=58699)

#energy1 = jHF(elems, coords, charge, mult, basis, basisfile="sto3g-orca.basis")
#energy1 = jHF(H2, "read_file", basisfile="sto3g-orca.basis")
#energy1b = jHF(H2, "read_file", basisfile="sto3g-orca-he-normalized.basis")
#energy2= jHF(H2, "STO-3G")
#energy3= jHF(H2, "def2-TZVP")
#energy = jHF(H2, integrals_from_file=true)
