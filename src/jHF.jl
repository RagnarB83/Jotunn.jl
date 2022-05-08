export jHF
"""
Jotunn: jHF: a RHF/UHF program
"""
function jHF(fragment, basisset="STO-3G"; HFtype="RHF", guess="hcore", basisfile="none", maxiter=120, 
    print_final_matrices=false, rmsDP_threshold=5e-9, maxDP_threshold=1e-7, tei_type="4c",
    energythreshold=1e-8, debugprint=false, fock_algorithm="turbo", levelshift=1.0, lshift_thresh=0.001,
    printlevel=2)

    print_program_header()
    global debugflag  = debugprint

    #Removing old matrix files if present
    if debugprint == true || print_final_matrices == true
        rm("Fmatrix", force=true)
        rm("Cmatrix", force=true)
        rm("Pmatrix", force=true)
    end

    ##########################
    # BASIC SYSTEM SETUP
    ##########################
    #Num. electrons and nuc-nuc repulsion from geometry
    sum_nuccharge=sum([elem_to_nuccharge(el) for el in fragment.elems])
    num_el=floor(Int64,sum_nuccharge-fragment.charge)
    E_ZZ=nuc_nuc_repulsion(fragment.elems,fragment.coords)

    #RHF vs. UHF
    if HFtype == "RHF"
        #Doubly occupied orbitals
        numoccorbs=Int64(num_el/2)
        unpaired_electrons=0
    elseif HFtype == "UHF"
        unpaired_electrons= fragment.mult - 1
        paired_el=num_el-unpaired_electrons
        paired_el_half=paired_el/2
        numoccorbs_⍺=trunc(Int64,paired_el_half+unpaired_electrons)
        numoccorbs_β=trunc(Int64,paired_el_half)
    else
        println("Unknown HFtype! Exiting.")
        exit()
    end
    #Print basic system properties
    print_system(num_el,fragment.formula,E_ZZ,fragment.charge,fragment.mult,fragment.numatoms,unpaired_electrons)

    ##########################
    # INTEGRALS
    ##########################
    #Setting up 1-electron and 2-electron integrals
    println("Integrals provided via GaussianBasis.jl library")       
    #Basis set object creation
    bset = basis_set_create(basisset,fragment.elems,fragment.coords; basisfile=basisfile)
    dim = bset.nbas
    #Simple array of atom indices that maps onto bfs
    bset_atom_mapping = bf_atom_mapping(bset)
    println("\nCalculating 1-electron integrals")
    T = kinetic(bset)
    V = nuclear(bset)
    S = overlap(bset)

    println("Calculating 2-electron integrals")
    #Two-electron integrals: calculate and store
    #tei_type: 4c, sparse4c
    @time tei = tei_calc(bset,tei_type)

    ##########################
    # CHOOSING FOCK ALGORITHM
    ##########################
    #Choosing Fock algorithm (RHF vs. UHF, user-defined bs. best for small-system)
    Fock,fock_algorithm = choose_Fock(HFtype,fock_algorithm,dim,tei_type)
    #Fock=Fock_loop_sparse
    #fock_algorithm="sparse"


    ##########################
    # stuff
    ##########################
    #Forming Hcore matrix from kinetic and nuclear-attraction matrices
    Hcore = T + V 
    
    #Overlap diagonalization and transformation matrices
    Sval,Svec = eigen(S)
    lowest_S_eigenval=minimum(Sval)
    
    SVAL_minhalf = Diagonal(Sval)^-0.5
    Stemp = SVAL_minhalf*transpose(Svec)
    S_minhalf = Svec * Stemp

    ##########################
    # GUESS
    ##########################
    #Create initial guess for density matrix
    println("Providing guess for density matrix\n")
    if guess == "hcore"
        #Setting P to 0. Means that Fock matrix becomes F = Hcore + 0. See Fock functions.
        if HFtype=="RHF"
            P = zeros(dim,dim)
        else
            P_⍺ = zeros(dim,dim)
            P_β = zeros(dim,dim)
        end
    else
        println("unknown guess")
        exit()
    end

    print_calculation_setup(HFtype,basisset,dim,guess,tei_type,fock_algorithm,lowest_S_eigenval)

    ##########################
    # SCF
    ##########################
    println("\nBeginning SCF iterations")
    #Initializing
    energy_old=0.0; energy=0.0; P_RMS=9999; finaliter=nothing
    eps_⍺=zeros(dim); eps_β=zeros(dim); eps=zeros(dim);
    global levelshift_flag = true
    #SCF loop beginning
    if printlevel == 1
        #println("Iteration     Energy deltaE  RMS-DP Max-DP Levelshift")
        @printf("%6s %17s %17s %17s %17s %10s %10s\n", "Iter", "Energy", "deltaE", "RMS-DP", "Max-DP", "Levelshift", "Damping")
    end
    @time for iter in 1:maxiter
        #Printing per iteration
        if printlevel > 1
            #Fair amount of printing
            print_iteration_header(iter)
        end
        if HFtype == "RHF"
            F = Fock(Hcore,P,dim,tei) #Update Fock-matrix
            #if iter == 2
            #    println("iter: $iter")
            #    println("F: $F")
            #    exit()
            #end
            #Possible levelshifting
            F = levelshift_control(F,levelshift,numoccorbs,dim,P_RMS,rmsDP_threshold,iter, turnoff_threshold=lshift_thresh,printlevel)
            F′ = transpose(S_minhalf)*F*S_minhalf #Transform Fock matrix
            eps, C′ = eigen(F′) #Diagonalize transformed Fock to get eps and C'
            C = S_minhalf*C′ # Get C from C'
            P_old=deepcopy(P) #Keep copy of old P
            P = makedensity(C, dim, numoccorbs) #Calculate new P from C
            energy = E_ZZ + 0.5 * tr((Hcore+F)*P) #Calculate energy

            if debugprint == true 
                #TODO: Appending instead
                write_matrices(F,C,P) 
            end
        else
            #Solve ⍺ equations
            F_⍺ = Fock(Hcore,P_⍺,P_β,dim,tei) #Update Fock-matrix alpha
            F_⍺ = levelshift_control(F_⍺,levelshift,numoccorbs_⍺,dim,P_RMS,rmsDP_threshold,iter, turnoff_threshold=lshift_thresh)
            F′_⍺ = transpose(S_minhalf)*F_⍺*S_minhalf #Transform Fock matrix
            eps_⍺, C′_⍺ = eigen(F′_⍺) #Diagonalize transformed Fock to get eps and C'
            C_⍺ = S_minhalf*C′_⍺ # Get C from C'
            P_⍺_old=deepcopy(P_⍺) #Keep copy of old P
            P_⍺ = makedensity(C_⍺, dim, numoccorbs_⍺, 1.0) #Calculate new P from C
            #Solve β equations
            F_β = Fock(Hcore,P_β,P_⍺,dim,tei) #Update Fock-matrix beta
            F_β = levelshift_control(F_β,levelshift,numoccorbs_β,dim,P_RMS,rmsDP_threshold,iter, turnoff_threshold=lshift_thresh)
            F′_β = transpose(S_minhalf)*F_β*S_minhalf #Transform Fock matrix
            eps_β, C′_β = eigen(F′_β) #Diagonalize transformed Fock to get eps and C'
            C_β = S_minhalf*C′_β # Get C from C'
            P_β_old=deepcopy(P_β) #Keep copy of old P
            P_β = makedensity(C_β, dim, numoccorbs_β, 1.0) #Calculate new P from C
            
            #Combined density matrix
            P=P_⍺+P_β
            P_old=P_⍺_old+P_β_old
            P_⍺_β=P_⍺-P_β #spin density matrix
            energy = E_ZZ + 0.5 * tr((Hcore+F_⍺)*P_⍺) + 0.5 * tr((Hcore+F_β)*P_β)  #Calculate energy
        end



        ##########################
        # CONVERGENCE CHECK
        ##########################
        deltaE = energy-energy_old
        P_RMS, P_MaxE = deltaPcheck(P, P_old)
        energy_old=energy

        #Printing per iteration
        iteration_printing(iter,printlevel,energy,deltaE,energythreshold,P_RMS,rmsDP_threshold,
            P_MaxE,maxDP_threshold,levelshift_flag)

        if P_MaxE < maxDP_threshold && P_RMS < rmsDP_threshold && abs(deltaE) < energythreshold
            #println("\n                              SCF converged in $iter iterations! Hell yeah! 🎉")
            print(Crayon(foreground = :green, bold = true), "\n                              SCF converged in $iter iterations! Hell yeah! 🎉\n\n",Crayon(reset=true))
            finaliter=iter
            if HFtype == "RHF"
                print_energy_contributions(energy,Hcore,F,P,T,E_ZZ)
                #Printing of matrices if requested. 
                if print_final_matrices == true 
                    write_matrix_to_file(F,"Fmatrix")
                    write_matrix_to_file(C,"Cmatrix")
                    write_matrix_to_file(P,"Pmatrix")
                end
            else
                #    print_energy_contributions(energy,Hcore,F,P,T,E_ZZ)
                if print_final_matrices == true 
                    write_matrix_to_file(P,"Pmatrix")
                    write_matrix_to_file(C_⍺,"C_amatrix")
                    write_matrix_to_file(C_β,"C_bmatrix")
                    write_matrix_to_file(F_⍺,"F_amatrix")
                    write_matrix_to_file(F_β,"F_bmatrix")
                    write_matrix_to_file(P_⍺,"P_amatrix")
                    write_matrix_to_file(P_β,"P_bmatrix")
                end
            end
            break #Break from loop
        end
        if iter == maxiter
            println("Failed to converge in $maxiter iterations!")
            #TODO: Print more here?
            return nothing
        end
    end
    #####################################
    # POPULATION ANALYSIS AND PROPERTIES
    #####################################
    #SCF loop done. Calculate and print properties
    #ORBITALS AND POPULATION ANALYSIS
    if HFtype=="RHF"
        #Orbitals
        occupations=makeoccupationarray(numoccorbs,dim,2.0) #Occupation array
        print_MO_energies(occupations,eps)
        #Mulliken
        charges = mulliken(S,P,bset,fragment.elems)
        print_Mulliken(charges,fragment.elems)
        P_⍺_β=zeros(dim,dim) #dummy spin-density
        #Mayer
        MBOs = Mayer_BO(S,P,P_⍺_β, bset_atom_mapping)
        print_Mayer_analysis(MBOs,fragment.elems)
    else
        #Orbitals
        occupations_⍺=makeoccupationarray(numoccorbs_⍺,dim,1.0) #Occupation array
        occupations_β=makeoccupationarray(numoccorbs_β,dim,1.0) #Occupation array
        print_MO_energies(occupations_⍺,occupations_β,eps_⍺,eps_β)
        #Mulliken
        P_⍺_β=P_⍺_β
        charges, spinpops = mulliken(S,P,P_⍺_β,bset,fragment.elems)
        print_Mulliken(charges,fragment.elems,spinpops)
        #Mayer
        MBOs = Mayer_BO(S,P,P_⍺_β, bset_atom_mapping)
        print_Mayer_analysis(MBOs,fragment.elems)
    end

    #####################################
    # FINAL RESULTS
    #####################################
    print_final_results(energy,fragment,num_el,basisset,HFtype,fock_algorithm,finaliter)

    return energy
end