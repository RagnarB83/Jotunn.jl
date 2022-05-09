using Jotunn
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)

#More options and more printing
energy= jHF(H2O, "def2-qzvpp"; maxiter=200, fock_algorithm="turbo", HFtype="RHF")
