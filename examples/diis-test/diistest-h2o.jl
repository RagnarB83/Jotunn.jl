using Jotunn
using PrettyTables
using Crayons

#Acetic acid fragment
@time H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)


basisname="sto-3g"
maxiter=60

#Comparing DIIS with different startup times
    @time result = jSCF(H2O, basisname; maxiter=maxiter, WFtype="RHF",
        diis=true, diis_startiter=4, levelshift=false, damping=true,printlevel=1)

#Comparing DIIS with different startup times
    @time result = jSCF(H2O, basisname; maxiter=maxiter, WFtype="UHF",
        diis=true, diis_startiter=4, levelshift=false, damping=true,printlevel=1)
