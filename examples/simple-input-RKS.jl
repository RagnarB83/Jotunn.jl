using Jotunn
LiH = create_fragment(coords_string="""
H 0.0 0.0 0.0
Li 0.0 0.0 0.74
""", charge=0, mult=1)


mol=LiH
basis="sto-3g"

jSCF(mol, basis; maxiter=200, WFtype="RKS", functional="lda")
