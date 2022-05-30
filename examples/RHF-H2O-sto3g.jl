using Jotunn
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)

#More options and more printing
res_4c = jSCF(H2O, "sto-3g"; maxiter=200, fock_algorithm="loop", WFtype="RHF", tei_type="4c")

res_sp = jSCF(H2O, "sto-3g"; maxiter=200, fock_algorithm="loop", WFtype="RHF", tei_type="sparse4c")


println("Energy (sp): $(res_sp["energy"])")
println("Energy (4c): $(res_4c["energy"])")
