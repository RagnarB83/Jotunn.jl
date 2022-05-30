using Jotunn
using PrettyTables
using Crayons

#Acetic acid fragment
@time Acetic = create_fragment(xyzfile="acetic.xyz", charge=0, mult=1)


basisname="def2-svp"
maxiter=100

#Comparing DIIS with different startup times
#@time res_damp_diis1= jSCF(Acetic, basisname; maxiter=maxiter, WFtype="RHF",
#    diis=true, diis_startiter=1, levelshift=false, damping=false,printlevel=1)
@time res_damp_diis2= jSCF(Acetic, basisname; maxiter=maxiter, WFtype="RHF",
    diis=true, diis_startiter=2, levelshift=false, damping=false,printlevel=1)
@time res_damp_diis3= jSCF(Acetic, basisname; maxiter=maxiter, WFtype="RHF",
    diis=true, diis_startiter=3, levelshift=false, damping=false,printlevel=1)
@time res_damp_diis4= jSCF(Acetic, basisname; maxiter=maxiter, WFtype="RHF",
    diis=true, diis_startiter=4, levelshift=false, damping=false,printlevel=1)
@time res_damp_diis5= jSCF(Acetic, basisname; maxiter=maxiter, WFtype="RHF",
    diis=true, diis_startiter=5, levelshift=false, damping=false,printlevel=1)
@time res_damp_diis6= jSCF(Acetic, basisname; maxiter=maxiter, WFtype="RHF",
    diis=true, diis_startiter=6, levelshift=false, damping=false,printlevel=1)
@time res_damp_diis8= jSCF(Acetic, basisname; maxiter=maxiter, WFtype="RHF",
    diis=true, diis_startiter=8, levelshift=false, damping=false,printlevel=1)
@time res_damp_diis10= jSCF(Acetic, basisname; maxiter=maxiter, WFtype="RHF",
    diis=true, diis_startiter=10, levelshift=false, damping=false,printlevel=1)
@time res_damp_diis12= jSCF(Acetic, basisname; maxiter=maxiter, WFtype="RHF",
    diis=true, diis_startiter=12, levelshift=false, damping=false,printlevel=1)
@time res_damp_diis16= jSCF(Acetic, basisname; maxiter=maxiter, WFtype="RHF",
    diis=true, diis_startiter=16, levelshift=false, damping=false,printlevel=1)


#Compare in a nice table
println("\n")
labels=["DIIS-2", "DIIS-3", "DIIS-4", "DIIS-5","DIIS-6","DIIS-8","DIIS-10","DIIS-12","DIIS-16"]
iters=[res_damp_diis2["finaliter"],res_damp_diis3["finaliter"],res_damp_diis4["finaliter"],
    res_damp_diis5["finaliter"],res_damp_diis6["finaliter"], res_damp_diis8["finaliter"], res_damp_diis10["finaliter"],
    res_damp_diis12["finaliter"], res_damp_diis16["finaliter"]]
energies=[res_damp_diis2["energy"],res_damp_diis3["energy"],res_damp_diis4["energy"],
    res_damp_diis5["energy"],res_damp_diis6["energy"], res_damp_diis8["energy"], res_damp_diis10["energy"],
    res_damp_diis12["energy"], res_damp_diis16["energy"]]

data=hcat(labels,iters,energies)
pretty_table(data; crop=:none,  noheader = true, formatters = ft_printf("%14.8f", [3]),
        tf = tf_simple, border_crayon = crayon"bold yellow", header_crayon = crayon"bold green")

