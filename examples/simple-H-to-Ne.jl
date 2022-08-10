using Jotunn
#H2 = create_fragment(coords_string="""
#H 0.0 0.0 0.0
#H 0.0 0.0 0.74
#""", charge=0, mult=1)
#H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)
#
H = create_fragment(coords_string="""H 0.0 0.0 0.0""", charge=0, mult=2)
He = create_fragment(coords_string="""He 0.0 0.0 0.0""", charge=0, mult=1)
Li = create_fragment(coords_string="""Li 0.0 0.0 0.0""", charge=0, mult=2)
Be = create_fragment(coords_string="""Be 0.0 0.0 0.0""", charge=0, mult=1)
B = create_fragment(coords_string="""B 0.0 0.0 0.0""", charge=0, mult=2)
C = create_fragment(coords_string="""C 0.0 0.0 0.0""", charge=0, mult=3)
N = create_fragment(coords_string="""N 0.0 0.0 0.0""", charge=0, mult=4)
O = create_fragment(coords_string="""O 0.0 0.0 0.0""", charge=0, mult=3)
F = create_fragment(coords_string="""F 0.0 0.0 0.0""", charge=0, mult=2)
Ne = create_fragment(coords_string="""Ne 0.0 0.0 0.0""", charge=0, mult=1)
mols=[H,He,Li,Be,B,C,N,O,F,Ne]
basis="def2-svp"

final_results=Dict()
for mol in mols
@time result=jSCF(mol, basis; maxiter=200, fock_algorithm="loop", WFtype="RHF", tei_type="sparse4c",
    calc_density=true)
    println("result: $result")
    final_results[mol.formula]=result["int_electrons"]
end

println("final_results: $final_results")
#@time result_sp_norm=jSCF(mol, basis; maxiter=200, fock_algorithm="loop", WFtype="RHF", tei_type="sparse4c")
#println("Energy (sp norm): $(result_sp_norm["energy"])")
#println("Energy (sp_perm): $(result_sp_perm["energy"])")

