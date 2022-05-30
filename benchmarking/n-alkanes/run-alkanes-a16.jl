using Jotunn
using NaturalSort
using Plots
using PrettyTables

#Reading XYZ-files into list of Jotunn fragments
xyzdir="./xyzfiles"
xyzfiles=sort(readdir(xyzdir), lt=natural) #natural sorting 

compxyzfile="a1.xyz"
xyzfile="a16.xyz"
println("xyzfile: $xyzfile")
comp_mol = create_fragment(xyzfile=xyzdir*"/"*compxyzfile, charge=0, mult=1)
mol = create_fragment(xyzfile=xyzdir*"/"*xyzfile, charge=0, mult=1)

#######################
# JOTUNN
#######################
#Settings
WFmethod="RHF"
basis="sto-3g"

Resultdict=Dict{String,Vector{Any}}()

#Compilation
println("Running compilation calc first") 
@time jSCF(comp_mol, basis; printlevel=1,diis=true, damping=true, maxiter=200, WFtype=WFmethod, nopop=true)
println("NOW RUNNING molecule :", mol.prettyformula)

    @time res = jSCF(mol, basis; printlevel=1, diis=true, damping=true, maxiter=200, WFtype=WFmethod,
                nopop=true)

println("res:", res)
