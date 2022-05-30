using Jotunn
using PrettyTables
using Crayons

#H2o cation geometry from deprince: https://www.chem.fsu.edu/~deprince/programming_projects/uks_dft/
#For comparison to their results
H2O_cation = create_fragment(coords_string="""
O       -1.328922301      0.000000000     -2.658617225
H       -1.328922301      1.430522516     -1.566094763
H       -1.328922301     -1.415897278     -1.547206005
""", charge=1, mult=2)


basisname="sto-3g"
maxiter=60

    @time result = jSCF(H2O_cation, basisname; maxiter=maxiter, WFtype="UHF",
        diis=true, diis_startiter=4, diis_size=5, levelshift=false, damping=true,printlevel=1)
