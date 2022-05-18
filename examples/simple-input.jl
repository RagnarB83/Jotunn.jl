using Jotunn
H2 = create_fragment(coords_string="""
H 0.0 0.0 0.0
H 0.0 0.0 0.74
""", charge=0, mult=1)
H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)


mol=H2O
basis="def2-svp"

#Measure 1st run separately as it includes compilation
result_sp4c=jHF(mol, basis; maxiter=200, fock_algorithm="loop", HFtype="RHF", tei_type="sparse4c")
result_4c=jHF(mol, basis; maxiter=200, fock_algorithm="loop", HFtype="RHF", tei_type="4c")

println("Energy (sp): $(result_sp4c["energy"])")
println("Energy (4c): $(result_4c["energy"])")
