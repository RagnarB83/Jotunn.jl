#############################
#Printing stuff
##############################

"""
print_if_level2: Simple print function that only prints if printlevel is met or higher
"""
function print_if_level2(string,printlevel)
    if printlevel >= 2
        println(string)
    end
end

"""
printdebug: Simple print function that only prints if global debugflag is true
"""
function printdebug(string)
    if debugflag == true
        println("$string")
    end
end

"""
colorvalue_threshold: Return red/green color symbol depending on value smaller/larger than threshold
"""
function colorvalue_threshold(var,threshold)
    if var < threshold
        return :green
    else
        #return :red
        if var/threshold > 10000
            return :red
        elseif 1000 < var/threshold < 10000
            return :light_magenta
        else
            return :light_cyan
        end
    end
end

"""
colorvalue_bool: Return yellow color symbol depending on if boolean is true or false
"""
function colorvalue_bool(var)
    if var == true
        return :yellow
    else
        return :default
    end
end

"""
print_program_header: Print Jotunn program header
"""
function print_program_header()
    print(Crayon(foreground = :white, bold = true), "="^50*"\n",Crayon(reset=true))
    print(Crayon(foreground = :blue, bold = true), "                   JOTUNN\n",Crayon(reset=true))
    print(Crayon(foreground = :white, bold = false), "a simple quantum chemistry program in Julia\n",Crayon(reset=true))
    print(Crayon(foreground = :white, bold = true), "="^50*"\n",Crayon(reset=true))
    print(Crayon(foreground = :yellow, bold = true), "\njHF: a RHF/UHF program\n",Crayon(reset=true))
end

"""
print_system: print molecular system properties
"""
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
print_calculation_setup: Print calculation setup
"""
function print_calculation_settings(HFtype,basisset,dim,guess,tei_type,fock_algorithm,lowest_S_eigenval,
        levelshift,levelshift_val,lshift_thresh,damping,damping_val,damping_thresh,diis,diis_size,diis_startiter)
    labels=["HF type", "Basis set","No. basis functions","Guess", "2-electron type","Fock algorithm","S lowest eigenvalue", 
    "Levelshift", "Levelshift parameter", "Lshift-turnoff thresh.","Damping","Damping parameter", "Damp-turnoff thresh.",
    "DIIS","DIIS vec size","DIIS start iter."]
    stuff=[HFtype,basisset,string(dim),guess,tei_type,fock_algorithm,lowest_S_eigenval,
        string(levelshift),levelshift_val,lshift_thresh,
        string(damping),damping_val,damping_thresh,string(diis),string(diis_size),string(diis_startiter)]
    data=hcat(labels,stuff)
    print(Crayon(foreground = :green, bold = true), "CALCULATION SETTINGS\n",Crayon(reset=true))
    pretty_table(data; crop=:none,  noheader = true,
        formatters = ft_printf("%14.4f", [2]),
        tf = tf_simple, border_crayon = crayon"bold yellow", header_crayon = crayon"bold green")
end

"""
print_iteration_header: Print SCF iteration header (for printlevel > 1)
"""
function print_iteration_header(iter)
    println("")
    println(Crayon(foreground = :yellow, bold = true),"="^30,Crayon(reset=true))
    println(Crayon(foreground = :yellow, bold = true),"SCF iteration $iter",Crayon(reset=true))
    println(Crayon(foreground = :yellow, bold = true),"="^30,Crayon(reset=true))
end

"""
iteration_printing: Printing during SCF iterations.
printlevel=1 => minimal printing
printlevel >1 => more elaborate printing with SCF-iteration header
"""
function iteration_printing(iter,printlevel,energy,deltaE,energythreshold,P_RMS,
    rmsDP_threshold,P_MaxE,maxDP_threshold,levelshiftflag,damping_flag,diis_flag,FP_comm)
    if printlevel > 1
        #Fair amount of printing
        println("Damping: ", damping_flag)
        println("Levelshift: ", levelshiftflag)
        println("DIIS: ", diis_flag)
        println("Current energy: $energy")
        print("Energy change:",Crayon(foreground = colorvalue_threshold(abs(deltaE),energythreshold)), "$deltaE ",
            Crayon(reset=true)," (threshold: $energythreshold)\n")
        print("RMS-DP: ",Crayon(foreground = colorvalue_threshold(abs(P_RMS),rmsDP_threshold)), "$P_RMS ",
            Crayon(reset=true)," (threshold: $rmsDP_threshold)\n")
        print("Max-DP: ",Crayon(foreground = colorvalue_threshold(abs(P_MaxE),maxDP_threshold)), "$P_MaxE ",
            Crayon(reset=true)," (threshold: $maxDP_threshold)\n")
        #[F,P] commut
        println("Max[F,P]: ", maximum(FP_comm))
        rms_comm = sqrt(sum(x -> x*x, FP_comm) / length(FP_comm))
        println("RMS[F,P]: ", rms_comm)
    elseif printlevel == 1
        #Minimal printing
        #@printf("%6d %17.10f %17.10f %17.10f %17.10f %10s %10s\n", iter, energy, deltaE,P_RMS,P_MaxE,levelshiftflag, "no")
        #Note: Crayon output might add ~1-1.5 sec in total per 70 iterations
        @printf("%4d%15.9f", iter, energy)
        @printf("%s%16.9f%s",Crayon(foreground = colorvalue_threshold(abs(deltaE),energythreshold)),deltaE,Crayon(reset=true))
        @printf("%s%14.9f%s",Crayon(foreground = colorvalue_threshold(abs(P_RMS),rmsDP_threshold)),P_RMS,Crayon(reset=true))
        @printf("%s%14.9f%s",Crayon(foreground = colorvalue_threshold(abs(P_MaxE),maxDP_threshold)),P_MaxE,Crayon(reset=true))
        @printf("%s%8s%s",Crayon(foreground = colorvalue_bool(levelshiftflag)),levelshiftflag,Crayon(reset=true))
        @printf("%s%8s%s",Crayon(foreground = colorvalue_bool(damping_flag)) ,damping_flag,Crayon(reset=true))
        @printf("%s%8s%s",Crayon(foreground = colorvalue_bool(diis_flag)),diis_flag,Crayon(reset=true))
        @printf("\n")
    else
        #no SCF iteration printing
    end
end

"""
print_energy_contributions: Print energy contributions if SCF converged
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
print_final_results: Print final results table in the end of job
"""
function print_final_results(energy,fragment,num_el,basisset,HFtype,fock_algorithm,finaliter)
    println()
    print(Crayon(foreground = :green, bold = true), "FINAL RESULTS\n",Crayon(reset=true))
    labels=["Final HF energy","Molecule formula","Number of electrons","Basis set","HF type","Fock algorithm","SCF iterations"]
    stuff=[energy,fragment.formula,string(num_el),basisset,HFtype,fock_algorithm,string(finaliter)]
    data=hcat(labels,stuff)
    pretty_table(data; crop=:none,  noheader = true,
        formatters = ft_printf("%14.8f", [2]),
        tf = tf_simple, border_crayon = crayon"bold yellow", header_crayon = crayon"bold green")
        #border_crayon = crayon"green")
end

"""
print_Mulliken: print Mulliken analysis for closed-shell system. Using pretty tables
"""
function print_Mulliken(charges,elems)
    #println("\n\nMulliken Population Analysis")
    print(Crayon(foreground = :yellow, bold = true), "\nMulliken Population Analysis (closed-shell)\n",Crayon(reset=true))
    data=hcat([i for i in 1:length(elems)],elems,charges)
    pretty_table(data; header=["Atom", "Element", "Charge"], crop=:none,
        formatters = ft_printf("%10.6f", [3]))
    @printf("\nSum of charges: %.4f\n", sum(charges))
end

"""
print_Mulliken: print Mulliken analysis for open-shell system. Using pretty tables
"""
function print_Mulliken(charges,elems,spinpops)
    #println("\n\nMulliken Population Analysis (open-shell)")
    print(Crayon(foreground = :yellow, bold = true), "\nMulliken Population Analysis (open-shell)\n",Crayon(reset=true))
    data=hcat([i for i in 1:length(elems)],elems,charges,spinpops)
    pretty_table(data; header=["Atom", "Element", "Charge", "Spin pop."], crop=:none,
        formatters = ft_printf("%10.6f", [3,4]))
    @printf("\nSum of charges: %.4f\n", sum(charges))
    @printf("Sum of spin populations: %.4f\n", sum(spinpops))
end

"""
print_Mayer_analysis: print Mayer analysis bond orders
"""
function print_Mayer_analysis(MBOs,elems; mbo_print_threshold=0.01)
    #println("\nMayer bond orders (>0.01):")
    print(Crayon(foreground = :yellow, bold = true), "\nMayer bond orders\n",Crayon(reset=true))
    println("Threshold: $mbo_print_threshold")
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
print_MO_energies: Print MO energies for closed-shell system (using pretty tables)
"""
function print_MO_energies(occ,mos)
    println("\n")
    #println("MO Energies (closed-shell)")
    #print(Crayon(foreground = :yellow, bold = true), "-"^30*"\n",Crayon(reset=true))
    print(Crayon(foreground = :yellow, bold = true), "MO Energies (closed-shell)\n",Crayon(reset=true))
    #print(Crayon(foreground = :yellow, bold = true), "-"^30*"\n",Crayon(reset=true))
    #print(Crayon(foreground = :red, bold = true), "                        ⍺",Crayon(reset=true))
    println("")
    mos_ev=mos*27.211399
    data=hcat([string(i) for i in 1:length(mos)],occ,mos,mos_ev)
    h1 = Highlighter((data, i, j) -> j in (4), bold = true, foreground = :red )
    pretty_table(data; header=["MO", "Occ.", "E(Eh)", "E(eV)"], crop=:none,
    highlighters = (h1), formatters = ft_printf("%12.4f", [3,4]))
end

"""
print_MO_energies: Print MO energies for open-shell system (using pretty tables)
"""
function print_MO_energies(occ_⍺, occ_β, mos_⍺, mos_β)
    println("\n")
    #println("-"^30)
    #println("MO Energies (open-shell)")
    #println("-"^30)
    print(Crayon(foreground = :yellow, bold = true), "MO Energies (open-shell)\n",Crayon(reset=true))
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
write_matrices: write multiple matrices to disk
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
write_matrix_to_file: Write matrix to disk
https://docs.julialang.org/en/v1/stdlib/DelimitedFiles/
"""
function write_matrix_to_file(X,name)
    open(name, "a") do io
        #writedlm(io, X)
        pretty_table(io,X; header=[string(i) for i in 1:size(X)[2]], tf = tf_matrix, show_row_number = true)
    end
end