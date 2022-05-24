using Jotunn
using PrettyTables
using Crayons

#Acetic acid fragment
@time Acetic = create_fragment(xyzfile="acetic.xyz", charge=0, mult=1)


basisname="def2-tzvpp"
maxiter=50

#Comparing DIIS with different startup times
    @time result = jHF(Acetic, basisname; maxiter=maxiter, HFtype="RHF",
        diis=true, diis_startiter=4, levelshift=false, damping=false,printlevel=1)
