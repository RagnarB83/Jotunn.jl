using Jotunn
#
He = create_fragment(coords_string="""He 0.0 0.0 0.0""", charge=0, mult=1)
Ne = create_fragment(coords_string="""Ne 0.0 0.0 0.0""", charge=0, mult=1)
Ar = create_fragment(coords_string="""Ar 0.0 0.0 0.0""", charge=0, mult=1)
Kr = create_fragment(coords_string="""Kr 0.0 0.0 0.0""", charge=0, mult=1)

mols=[He,Ne,Ar,Kr]
basis="sto-3g"

final_results=Dict()
for mol in mols
    @time result=jSCF(mol, basis; maxiter=200, fock_algorithm="loop", WFtype="RHF", tei_type="sparse4c",
    calc_density=true)
    final_results[mol.formula]=result["int_electrons"]
end

println("final_results: $final_results")
