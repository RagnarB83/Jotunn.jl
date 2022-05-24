using Jotunn
using PrettyTables
using Crayons

#Acetic acid fragment
@time H2O = create_fragment(xyzfile="h2o.xyz", charge=0, mult=1)


basisname="sto-3g"
maxiter=30

#Comparing DIIS with different startup times
    @time result = jHF(H2O, basisname; maxiter=maxiter, HFtype="RHF",
        diis=true, diis_startiter=4, levelshift=false, damping=false,printlevel=1)
