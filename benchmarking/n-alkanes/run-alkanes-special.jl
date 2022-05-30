using Jotunn

#######################
# JOTUNN
#######################
#Settings
WFmethod="RHF"
basis="sto-3g"

Resultdict=Dict{String,Vector{Any}}()

C13=create_fragment(xyzfile="./xyzfiles/a13.xyz", charge=0, mult=1)
C14=create_fragment(xyzfile="./xyzfiles/a14.xyz", charge=0, mult=1)
C15=create_fragment(xyzfile="./xyzfiles/a15.xyz", charge=0, mult=1)
C16=create_fragment(xyzfile="./xyzfiles/a16.xyz", charge=0, mult=1)

res = jSCF(C16, basis; printlevel=1, damping=true, maxiter=200, WFtype=WFmethod)

#for i in [1,2,3,4,5,6]
#    println("DIIS startiter: $i")
#    res = jSCF(C14, basis; printlevel=1,diis=true, diis_startiter=i, damping=true, maxiter=200, WFtype=WFmethod)
#    println("For DIIS startiter: $i  iterations: $(res["finaliter"])")
#end
