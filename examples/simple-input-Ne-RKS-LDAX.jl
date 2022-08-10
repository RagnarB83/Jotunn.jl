using Jotunn
Ne = create_fragment(coords_string="""Ne 0.0 0.0 0.0
""", charge=0, mult=1)
mols=[Ne]
basis="sto-3g"

#Measure 1st run separately as it includes compilation
for mol in mols
@time result_sp_norm=jSCF(mol, basis; maxiter=200, fock_algorithm="loop", WFtype="RKS", functional="LDA-X", tei_type="sparse4c",
    calc_density=false)
println("result_sp_norm: $result_sp_norm")
end
