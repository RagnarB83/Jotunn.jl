using Jotunn

# Create molecular fragment
H2 = create_fragment(coords_string="""
H 0.0 0.0 0.0
H 0.0 0.0 0.74
""", charge=0, mult=1)

#Simple call
result= jSCF(H2, "sto-3g")
energy=result["energy"]
println("Result dictionary from Jotunn: $result")
println("Energy: $energy) Eh")

#RHF/STO-3G H2 (r=0.74 A) from ORCA
ref_energy=-1.116759307204

@assert isapprox(energy,ref_energy)
