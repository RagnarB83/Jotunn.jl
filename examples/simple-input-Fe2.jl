using Jotunn
Fe2 = create_fragment(coords_string="""Fe 500.0 0.0 0.0
Fe 500.0 0.0 3.0
""", charge=0, mult=1)
mols=[Fe2]
basis="def2-svp"

#Measure 1st run separately as it includes compilation
for mol in mols
@time result_sp_norm=jSCF(mol, basis; maxiter=200, fock_algorithm="loop", WFtype="RHF", tei_type="sparse4c",
    calc_density=true)
println("result_sp_norm: $result_sp_norm")
end
