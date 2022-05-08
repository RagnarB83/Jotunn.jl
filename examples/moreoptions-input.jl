using Jotunn
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)

#More keywords
energy= jHF(H2O, "STO-3G"; maxiter=200, fock_algorithm="turbo", 
    HFtype="RHF", levelshift=2.0, lshift_thresh=1e-4, tei_type="4c", 
    print_final_matrices=true, debugprint=true)
