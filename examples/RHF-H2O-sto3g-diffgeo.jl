using Jotunn
using BenchmarkTools

#Geometry used in this tutorial: https://www.chem.fsu.edu/~deprince/programming_projects/diis/
H2O = create_fragment(xyzfile="h2o_diffgeo.xyz", charge=0, mult=1)

#More options and more printing
energy= jSCF(H2O, "sto-3g"; maxiter=200, diis=true, levelshift=false, damping=false, 
    fock_algorithm="loop", WFtype="RHF", printlevel=2)
