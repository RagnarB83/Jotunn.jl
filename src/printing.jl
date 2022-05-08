#############################
#Printing stuff
##############################
function printdebug(string,var="")
    if debugflag == true
        println(string," : ",var)
    end
end

function print_program_header()
    print(Crayon(foreground = :white, bold = true), "="^50*"\n",Crayon(reset=true))
    print(Crayon(foreground = :blue, bold = true), "                   JOTUNN\n",Crayon(reset=true))
    print(Crayon(foreground = :blue, bold = false), "a simple quantum chemistry program in Julia\n",Crayon(reset=true))
    print(Crayon(foreground = :white, bold = true), "="^50*"\n",Crayon(reset=true))
    print(Crayon(foreground = :yellow, bold = true), "\njHF module: a RHF/UHF program\n",Crayon(reset=true))
end

function print_system(num_el,formula,E_ZZ,charge,mult,numatoms,unpaired_el)
    println("")
    print(Crayon(foreground = :green, bold = true), "SYSTEM\n",Crayon(reset=true))
    labels=["Number of atoms","Molecule formula","Charge","Spin multiplicity","No. electrons","No. unpaired electrons","Nuclear repulsion"]
    stuff=[string(numatoms),formula,string(charge),string(mult),string(num_el),string(unpaired_el),E_ZZ]
    data=hcat(labels,stuff)
    pretty_table(data; crop=:none,  noheader = true,
        formatters = ft_printf("%14.8f", [2]),
        tf = tf_simple, border_crayon = crayon"bold yellow", header_crayon = crayon"bold green")
end

"""
Print calculation setup
"""
function print_calculation_setup(HFtype,basisset,dim,guess,tei_type,fock_algorithm,lowest_S_eigenval)
    labels=["HF type", "Basis set","No. basis functions","Guess", "tei_type","Fock algorithm","S lowest eigenvalue"]
    stuff=[HFtype,basisset,dim,guess,tei_type,fock_algorithm,lowest_S_eigenval]
    data=hcat(labels,stuff)
    print(Crayon(foreground = :green, bold = true), "CALCULATION SETUP\n",Crayon(reset=true))
    pretty_table(data; crop=:none,  noheader = true,
        formatters = ft_printf("%14.8f", [2]),
        tf = tf_simple, border_crayon = crayon"bold yellow", header_crayon = crayon"bold green")
end

function print_iteration_header(iter)
    println("")
    println("="^30)
    println("SCF iteration $iter")
    println("="^30)
end

"""
Print energy contributions
"""
function print_energy_contributions(energy,Hcore,F,P,T,E_ZZ)
    #Printing final energy contributions
    one_elec_E=tr(Hcore*P)
    two_elec_E=0.5*tr((F-Hcore)*P)
    kin_energy=tr(T*P)
    pot_energy=energy-kin_energy
    virial_ratio=-pot_energy/kin_energy
    energies=[energy,E_ZZ,energy-E_ZZ,one_elec_E,two_elec_E,kin_energy,pot_energy,virial_ratio]
    labels=["Total energy","Nuclear repulsion","Electronic energy","1-electron energy","2-electron energy",
        "Kinetic energy","Potential energy","Virial ratio"]
    data=hcat(labels,energies)
    pretty_table(data; header=["Energy contributions", "E(Eh)"], crop=:none, body_hlines = [7],
        formatters = ft_printf("%14.8f", [2]))
    #@printf("Virial ratio: %.4f\n", virial_ratio)
end

"""
Print final results table
"""
function print_final_results(energy,fragment,num_el,basisset,HFtype,fock_algorithm,finaliter)
    println()
    print(Crayon(foreground = :green, bold = true), "FINAL RESULTS\n",Crayon(reset=true))
    labels=["Final HF energy","Molecule formula","Number of electrons","Basis set","HF type","Fock algorithm","SCF iterations"]
    stuff=[energy,fragment.formula,num_el,basisset,HFtype,fock_algorithm,string(finaliter)]
    data=hcat(labels,stuff)
    pretty_table(data; crop=:none,  noheader = true,
        formatters = ft_printf("%14.8f", [2]),
        tf = tf_simple, border_crayon = crayon"bold yellow", header_crayon = crayon"bold green")
        #border_crayon = crayon"green")
end

"""
Mulliken closed-shell with pretty tables
"""
function print_Mulliken(charges,elems)
    println("\n\nMulliken Population Analysis")
    data=hcat([i for i in 1:length(elems)],elems,charges)
    pretty_table(data; header=["Atom", "Element", "Charge"], crop=:none,
        formatters = ft_printf("%10.6f", [3]))
    @printf("\nSum of charges: %.4f\n", sum(charges))
end

"""
Mulliken open-shell with pretty tables
"""
function print_Mulliken(charges,elems,spinpops)
    println("\n\nMulliken Population Analysis (open-shell)")
    data=hcat([i for i in 1:length(elems)],elems,charges,spinpops)
    pretty_table(data; header=["Atom", "Element", "Charge", "Spin pop."], crop=:none,
        formatters = ft_printf("%10.6f", [3,4]))
    @printf("\nSum of charges: %.4f\n", sum(charges))
    @printf("\nSum of spin populations: %.4f\n", sum(spinpops))
end

"""
print Mayer analysis
only MBOs for now
"""
function print_Mayer_analysis(MBOs,elems; mbo_print_threshold=0.01)
    println("\nMayer bond orders (>0.01):")
    mbovals=[]
    bondpairs=[]
    for mbo in MBOs
        mbo_key_1=mbo[1][1]; mbo_key_2=mbo[1][2]
        mbo_val=mbo[2]
        elem1=elems[mbo_key_1]; elem2=elems[mbo_key_2]
        bondpair="$elem1$mbo_key_1 - $elem2$mbo_key_2"
        if mbo_val > mbo_print_threshold
            push!(bondpairs,bondpair)
            push!(mbovals,mbo_val)
        end
    end
    #println("bondpairs: $bondpairs")
    #println("mbovals: $mbovals")
    data=hcat(bondpairs,mbovals)
    #println("data: $data")
    pretty_table(data; header=["Bond", "MBO"], crop=:none, formatters = ft_printf("%6.4f", [2]))
end

"""
print_MO_energies closed-shell (pretty tables)
"""
function print_MO_energies(occ,mos)
    println("-"^30)
    println("MO Energies (closed-shell)")
    println("-"^30)
    print(Crayon(foreground = :red, bold = true), "                        ⍺",Crayon(reset=true))
    println("")
    mos_ev=mos*27.211399
    data=hcat([string(i) for i in 1:length(mos)],occ,mos,mos_ev)
    h1 = Highlighter((data, i, j) -> j in (4), bold = true, foreground = :red )
    pretty_table(data; header=["MO", "Occ.", "E(Eh)", "E(eV)"], crop=:none,
    highlighters = (h1), formatters = ft_printf("%12.4f", [3,4]))
end

"""
print_MO_energies: open-shell
"""
function print_MO_energies(occ_⍺, occ_β, mos_⍺, mos_β)
    println("-"^30)
    println("MO Energies (open-shell)")
    println("-"^30)
    print(Crayon(foreground = :red, bold = true), "                              ⍺", Crayon(foreground = :blue, bold = true),"                                    β",Crayon(reset=true))
    print("\n")
    mos_⍺_ev=mos_⍺*27.211399
    mos_β_ev=mos_β*27.211399
    data=hcat([string(i) for i in 1:length(mos_⍺)],occ_⍺,mos_⍺,mos_⍺_ev,occ_β,mos_β,mos_β_ev)
    h1 = Highlighter((data, i, j) -> j in (4), bold = true, foreground = :red )
    h2 = Highlighter((data, i, j) -> j in (7), bold = true, foreground = :blue )
    pretty_table(data; header=["MO", "Occ.", "E(Eh)", "E(eV)", "Occ.", "E(Eh)", "E(eV)"], crop=:none, 
        highlighters = (h1, h2), formatters = ft_printf("%12.4f", [3,4,6,7]))
end

"""
write multiple matrices
"""
function write_matrices(F,C,P)
    println("Writing current Fock matrix to disk: Fmatrix")
    println("Printing current MO matrix to disk: Cmatrix")
    println("Printing current P matrix to disk: Pmatrix")
    write_matrix_to_file(F,"Fmatrix")
    write_matrix_to_file(C,"Cmatrix")
    write_matrix_to_file(P,"Pmatrix")
end

"""
Write matrix to disk
https://docs.julialang.org/en/v1/stdlib/DelimitedFiles/
"""
function write_matrix_to_file(X,name)
    open(name, "a") do io
        #writedlm(io, X)
        pretty_table(io,X; header=[string(i) for i in 1:size(X)[2]], tf = tf_matrix, show_row_number = true)
    end
end