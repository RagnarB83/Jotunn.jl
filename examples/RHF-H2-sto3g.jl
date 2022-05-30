using Jotunn
H2 = create_fragment(coords_string="""
H 0.0 0.0 0.0
H 0.0 0.0 0.74
""", charge=0, mult=1)

#More options and more printing
@time res_sp = jSCF(H2, "sto-3g"; maxiter=200, fock_algorithm="loop", WFtype="RHF", tei_type="sparse4c")


println("Energy (sp): $(res_sp["energy"])")
println("res_sp: $res_sp")
