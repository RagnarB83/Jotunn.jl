using Jotunn
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)
#All features
energy = jHF(H2O, "sto-3g"; HFtype="UHF", guess="hcore", 
    basisfile="none", maxiter=200, 
    print_final_matrices=true, debugprint=true, 
    rmsDP_threshold=5e-9, maxDP_threshold=1e-7, energythreshold=1e-8, 
    tei_type="4c", fock_algorithm="turbo", 
    levelshift=1.0, lshift_thresh=0.001)
