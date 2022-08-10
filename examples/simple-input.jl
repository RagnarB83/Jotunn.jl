using Jotunn
H2 = create_fragment(coords_string="""
H 0.0 0.0 0.0
H 0.0 0.0 0.74
""", charge=0, mult=1)
#H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)
#
mols=[H2]
basis="sto-3g"

#Measure 1st run separately as it includes compilation
for mol in mols
@time result_sp_norm=jSCF(mol, basis; maxiter=200, fock_algorithm="loop", WFtype="RHF", tei_type="sparse4c",
    calc_density=false)
println("result_sp_norm: $result_sp_norm")
end
#@time result_sp_norm=jSCF(mol, basis; maxiter=200, fock_algorithm="loop", WFtype="RHF", tei_type="sparse4c")


#println("Energy (sp norm): $(result_sp_norm["energy"])")
#println("Energy (sp_perm): $(result_sp_perm["energy"])")

