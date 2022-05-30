using Jotunn
using PrettyTables
using Crayons

#Acetic acid fragment
@time H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)


basisname="sto-3g"
maxiter=12

#Comparing DIIS with different startup times
    @time result = jSCF(H2O, basisname; maxiter=maxiter, WFtype="RHF",
        diis=false, levelshift=false, damping=false,printlevel=1)
