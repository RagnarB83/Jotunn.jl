using Jotunn
He = create_fragment(coords_string="""He 0.0 0.0 0.0
""", charge=0, mult=1)
#H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)
mols=[He]
basis="def2-svp"

#Measure 1st run separately as it includes compilation
for mol in mols
@time result_sp_norm=jSCF(mol, basis; maxiter=200, fock_algorithm="loop", WFtype="RKS", functional="NoXC", tei_type="sparse4c",
    calc_density=false, levelshift=true)
println("result_sp_norm: $result_sp_norm")
end
