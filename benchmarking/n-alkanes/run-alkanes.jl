using Jotunn
using NaturalSort
using Plots
using PrettyTables

#Reading XYZ-files into list of Jotunn fragments
xyzdir="./xyzfiles"
molecules=[]
xyzfiles=sort(readdir(xyzdir), lt=natural) #natural sorting 
for file in xyzfiles
    mol = create_fragment(xyzfile=xyzdir*"/"*file, charge=0, mult=1)
    push!(molecules,mol)
end


#######################
# JOTUNN
#######################
#Settings
WFmethod="RHF"
basis="sto-3g"

Resultdict=Dict{String,Vector{Any}}()

#Compilation
println("Running compilation calc first") 
@time jSCF(molecules[1], basis; printlevel=0,diis=true, damping=true, maxiter=200, WFtype=WFmethod)
#Looping over 
println("Now running benchmark set")
for mol in molecules
    println("NOW RUNNING molecule :", mol.prettyformula)
    @time res = jSCF(mol, basis; printlevel=0, diis=true, damping=true, maxiter=200, WFtype=WFmethod)
    Resultdict[mol.prettyformula] = [res["finaliter"],res["time"],mol.formula]
end

#Collecting data
labels = collect(keys(Resultdict))
iters = [Resultdict[i][1] for i in keys(Resultdict)]
jotunn_times = [Resultdict[i][2] for i in keys(Resultdict)]
formulas=[Resultdict[i][3] for i in keys(Resultdict)]
Catoms = [count("C",i) for i in formulas]
#Pretty tables
data=hcat(labels,iters,jotunn_times)
pretty_table(data; crop=:none,  header=["Mol", "Iter", "Time(s)"], formatters = ft_printf("%14.8f", [3]),
        tf = tf_simple, border_crayon = crayon"bold yellow", header_crayon = crayon"bold green")



#######################
#Plotting 
#######################
#Jotunn data
Plots.plot(Catoms,jotunn_times, lw=2, seriestype = :line, linecolor= :blue, label="Jotunn", grid=false)
Plots.plot!(Catoms, jotunn_times, seriestype = :scatter, markercolor= :blue, title = "Timings for set: n-alkanes. $WFmethod/$basis",
     label="Jotunn", grid=false, legend=:topleft)


xlabel!("Carbon chainlength")
ylabel!("Time (seconds)")
savefig("Benchmark-alkanes.pdf")
