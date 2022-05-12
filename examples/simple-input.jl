using Jotunn

# Create molecular fragment
H2 = create_fragment(coords_string="""
H 0.0 0.0 0.0
H 0.0 0.0 0.74
""", charge=0, mult=1)

#Simple call
result= jHF(H2, "sto-3g")
println("Result dictionary from Jotunn: $result")
println("Energy: $(result["energy"]) Eh")