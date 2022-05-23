using Jotunn
using BenchmarkTools

#Geometry used in this tutorial: https://www.chem.fsu.edu/~deprince/programming_projects/diis/

#Specifically the doubled O-H one : https://www.chem.fsu.edu/~deprince/programming_projects/diis/
H2O = create_fragment(xyzfile="h2o_strainedgeo.xyz", charge=0, mult=1)


#More options and more printing
energy= jHF(H2O, "sto-3g"; maxiter=7, diis=true, diis_size=5, levelshift=false, damping=false,damping_thresh=0.001,
    fock_algorithm="loop", HFtype="RHF", printlevel=2)

#energy= jHF(H2O, "sto-3g"; maxiter=150, diis=false, levelshift=false, damping=true, damping_thresh=0.00001,
#    fock_algorithm="loop", HFtype="RHF", printlevel=1)
