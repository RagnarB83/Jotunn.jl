using Jotunn
H2 = create_fragment(coords_string="""H 0.0 0.0 0.0
H 0.0 0.0 0.74
""", charge=0, mult=1)
mols=[H2]
basis="sto-3g"

#Measure 1st run separately as it includes compilation
for mol in mols
@time result_sp_norm=jSCF(mol, basis; maxiter=200, fock_algorithm="loop", WFtype="RKS", functional="LDA-X(libxc)", tei_type="sparse4c",
    calc_density=false, diis=true, levelshift=false, printlevel=1)
println("result_sp_norm: $result_sp_norm")
end
