using Jotunn
using PrettyTables
using Crayons

#Acetic acid fragment
@time Acetic = create_fragment(xyzfile="acetic.xyz", charge=0, mult=1)


basisname="def2-svp"
maxiter=100

#Comparing DIIS with different startup times
values=[2,4,6,8,10,12,14,16,18,20,22]
Resultdict_Damp_DIIS=Dict{String,Vector{Any}}()
for val in values
    @time result = jHF(Acetic, basisname; maxiter=maxiter, HFtype="RHF",
        diis=true, diis_startiter=val, levelshift=false, damping=true,printlevel=1)
        println("Resultdict_Damp_DIIS:", Resultdict_Damp_DIIS)
        println("result:", result)
        Resultdict_Damp_DIIS["DIIS"*string(val)] = [result["energy"],result["finaliter"]]
end

#Comparing DIIS with different startup times
values=[2,4,6,8,10,12,14,16,18,20,22]
Resultdict_DIIS=Dict{String,Vector{Float64}}()
for val in values
    @time result = jHF(Acetic, basisname; maxiter=maxiter, HFtype="RHF",
        diis=true, diis_startiter=val, levelshift=false, damping=false,printlevel=1)
        Resultdict_DIIS["DIIS"*string(val)] = [result["energy"],result["finaliter"]]
end

#Compare Damp-DIIS comparison in a nice table
println("\n")
labels = collect(keys(Resultdict_Damp_DIIS))
energies = [Resultdict_Damp_DIIS[i][1] for i in keys(Resultdict_Damp_DIIS)]
iters = [Resultdict_Damp_DIIS[i][2] for i in keys(Resultdict_Damp_DIIS)]
data=hcat(labels,iters,energies)
pretty_table(data; crop=:none,  header=["Setting", "Iter", "Energy"], formatters = ft_printf("%14.8f", [3]),
        tf = tf_simple, border_crayon = crayon"bold yellow", header_crayon = crayon"bold green")

#Compare DIIS comparison in a nice table
println("\n")
labels = collect(keys(Resultdict_DIIS))
energies = [Resultdict_DIIS[i][1] for i in keys(Resultdict_DIIS)]
iters = [Resultdict_DIIS[i][2] for i in keys(Resultdict_DIIS)]
data=hcat(labels,iters,energies)
pretty_table(data; crop=:none,  header=["Setting", "Iter", "Energy"], formatters = ft_printf("%14.8f", [3]),
        tf = tf_simple, border_crayon = crayon"bold yellow", header_crayon = crayon"bold green")
