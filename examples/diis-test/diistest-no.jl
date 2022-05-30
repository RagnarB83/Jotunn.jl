using Jotunn
using PrettyTables
using Crayons

NO = create_fragment(coords_string="""
N 0.0 0.0 0.0
O 0.0 0.0 1.15
""", charge=0, mult=2)

basisname="sto-3g"
maxiter=260

#No DIIS
@time result= jSCF(NO, basisname; WFtype="UHF", levelshift=true, maxiter=200, tei_type="sparse4c")

#Comparing DIIS with different startup times
    @time result = jSCF(NO, basisname; maxiter=maxiter, WFtype="UHF",
        diis=true, diis_startiter=2, diis_size=5, levelshift=true, damping=true,printlevel=1)
