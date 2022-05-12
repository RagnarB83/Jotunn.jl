using Jotunn

H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)

#All features
result = jHF(H2O, "STO-3G"; HFtype="RHF", guess="hcore", basisfile="none", maxiter=120, 
    print_final_matrices=false, rmsDP_threshold=5e-9, maxDP_threshold=1e-7, tei_type="4c",
    energythreshold=1e-8, debugprint=false, fock_algorithm="turbo", 
    levelshift=false, levelshift_val=0.10, lshift_thresh=0.01,
    damping=true, damping_val=0.4, damping_thresh=0.01,
    printlevel=1)


println("Result: $result")
