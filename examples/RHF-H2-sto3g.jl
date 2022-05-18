using Jotunn
H2 = create_fragment(coords_string="""
H 0.0 0.0 0.0
H 0.0 0.0 0.74
""", charge=0, mult=1)

#More options and more printing
res_4c = jHF(H2, "sto-3g"; maxiter=200, fock_algorithm="loop", HFtype="RHF", tei_type="4c")

res_sp = jHF(H2, "sto-3g"; maxiter=200, fock_algorithm="loop", HFtype="RHF", tei_type="sparse4c")


println("Energy (sp): $(res_sp["energy"])")
println("Energy (4c): $(res_4c["energy"])")
