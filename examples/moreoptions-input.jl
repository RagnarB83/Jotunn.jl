using Jotunn
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)

#More options and more printing
result= jSCF(H2O, "sto-3g"; maxiter=200, fock_algorithm="turbo", printlevel=1, 
    WFtype="RHF", levelshift=true, tei_type="4c", 
    print_final_matrices=false, debugprint=false)
