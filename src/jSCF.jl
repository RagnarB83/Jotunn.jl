export jSCF


"""
jSCF: the Jotunn SCF program (RHF,UHF,RKS,UKS)

"""
function jSCF(fragment, basisset="sto-3g"; WFtype::String="RHF", 
    functional::String="none", libxc_keyword::String="none", manual_func::String="none",
    guess::String="hcore", basisfile::String="none", maxiter::Int64=120, 
    print_final_matrices::Bool=false, rmsDP_threshold::Float64=5e-9, maxDP_threshold::Float64=1e-7, 
    tei_type::String="sparse4c", energythreshold::Float64=1e-8,
    fock_algorithm::String="loop-sparse", levelshift::Bool=false, levelshift_val::Float64=0.10, 
    lshift_thresh::Float64=0.01,damping::Bool=true, damping_val::Float64=0.4, 
    damping_thresh::Float64=0.01, diis::Bool=true, diis_size::Int64=5, diis_startiter::Int64=4, 
    DIISBfac::Float64=1.05, diis_error_conv_threshold::Float64=5e-7, calc_density::Bool=false,
    printlevel::Int64=1, nopop::Bool=false, iDMFT_kappa=0.1109136,iDMFT_b=0.18121047 )

    #Beginning of time block
    totaltime=@elapsed begin

    print_program_header(printlevel)


    #Removing old matrix files if present
    if printlevel > 2 || print_final_matrices == true
        for matrixfile in ["Fmatrix","Cmatrix","Pmatrix","C_a_matrix","C_b_matrix","F_a_matrix",
            "F_b_matrix","P_a_matrix","P_b_matrix"]
            rm(matrixfile, force=true)
        end
    end

    #############################
    # BASIS SET
    #############################
    #Setting up 1-electron and 2-electron integrals
    print_if_level("Integrals provided via GaussianBasis.jl library",1,printlevel)  
    #Basis set object creation by GaussianBasis
    bset = basis_set_create(basisset,fragment.elems,fragment.coords; basisfile=basisfile,printlevel)
    dim = bset.nbas
    #bset_atom_mapping = bf_atom_mapping(bset)
    #Array of tuples with atom,shell info for each BF
    bf_atom_shell_map=create_bf_shell_map(bset)
    #Simple array of atom indices that maps onto bfs
    bset_atom_mapping = [j[1] for j in bf_atom_shell_map]
    # Create array of l and ml values
    mlvalues=create_ml_values(bset)
    #lvalues=create_l_values(bset)
    #println("Here")
    #for i in 1:9
    #    println("set.basis[i]:", bset.basis[i])
    #    println("bset.basis[i].l:", bset.basis[i].l)
    #end
    lvalues = [bset.basis[i].l for i in 1:bset.nshells]
    #Print basis set information here
    if printlevel >0
        println(GaussianBasis.string_repr(bset))
    end

    #############################
    # BASIC SYSTEM SETUP
    #############################
    #Num. electrons and nuc-nuc repulsion from total charge and nuclear charges
    #sum_nuccharge=sum([elem_to_nuccharge(el) for el in fragment.elems])
    #num_el=floor(Int64,sum_nuccharge-fragment.charge)
    num_el=fragment.numelectrons
    E_ZZ=nuc_nuc_repulsion(fragment.elems,fragment.coords)

    #Check whether chosen multiplicity makes sense before continuing
    check_multiplicity(num_el,fragment.charge,fragment.mult)
    if fragment.mult == 2 && (WFtype=="RHF" || WFtype=="RKS")
        print_if_level("RHF/RKS and multiplicity > 1 is not possible. Switching to UHF/UKS",1,printlevel)
        if WFtype=="RHF" WFtype="UHF" end
        if WFtype=="RKS" WFtype="UKS" end
    end

    #Calculation setup for RHF vs. UHF vs. RKS vs. UKS
    if WFtype == "RHF"
        #Num occ orbitals in RHF
        numoccorbs=Int64(num_el/2)
        #Array of occupation numbers
        orb_occupations=makeoccupationarray(numoccorbs,dim,2.0)
        #
        unpaired_electrons=0
        DFTobj=nothing
        grid=false
        Intobj=Jint_HF
    elseif WFtype == "i-DMFT"
        #Num occ orbitals in RHF
        numoccorbs=Int64(num_el/2)
        #Array of occupation numbers
        orb_occupations=makeoccupationarray(numoccorbs,dim,2.0)
        println("Initial orb_occupations:", orb_occupations)
        #
        unpaired_electrons=0
        DFTobj=nothing
        grid=false
        Intobj=Jint_HF
    elseif WFtype == "RKS"
        grid=true
        #restricted Kohn-Sham
        #Num occ orbitals in RHF
        numoccorbs=Int64(num_el/2)
        #Array of occupation numbers
        orb_occupations=makeoccupationarray(numoccorbs,dim,2.0)
        unpaired_electrons=0
        println("RKS chosen.")
        println("functional: $functional and libxc_keyword: $libxc_keyword")
        if functional == "none" && libxc_keyword == "none"
            println("WFtype is RKS but no functional or libxc_keyword chosen. Exiting.")
            exit()
        elseif functional != "none" && libxc_keyword != "none"
            println("Error: Both functional and libxc_keyword can not be used. Exiting.")
            exit()
        end
        println("Calling choose_functional")
        libxc_functional,func_rung, manual_option = choose_functional(functional,libxc_keyword)
        #Non-hybrid
        if func_rung < 4
            DFTobj=JDFT(manual=manual_option, manual_func=manual_func, hybrid=false, libxc_functional=libxc_functional,
            functional_rung=func_rung)
            Intobj=Jint_KS
        #Hybrid
        else
            DFTobj=JDFT(manual=manual_option, manual_func=manual_func, hybrid=true, libxc_functional=libxc_functional,
            functional_rung=func_rung)
            Intobj=Jint_hybridKS
        end

    elseif WFtype == "UKS"
        grid=true
        unpaired_electrons= fragment.mult - 1
        paired_el=num_el-unpaired_electrons
        paired_el_half=paired_el/2
        numoccorbs_⍺=trunc(Int64,paired_el_half+unpaired_electrons)
        numoccorbs_β=trunc(Int64,paired_el_half)

        println("functional: $functional and libxc_keyword: $libxc_keyword")
        if functional == "none" && libxc_keyword == "none"
            println("WFtype is UKS but no functional or libxc_keyword chosen. Exiting.")
            exit()
        elseif functional != "none" && libxc_keyword != "none"
            println("Error: Both functional and libxc_keyword can not be used. Exiting.")
            exit()
        end

        println("Calling choose_functional")
        libxc_functional,func_rung, manual_option = choose_functional(functional,libxc_keyword, openshell=true)
        #Non-hybrid
        if func_rung < 4
            DFTobj=JDFT(manual=manual_option, manual_func=manual_func, hybrid=false, libxc_functional=libxc_functional,
            functional_rung=func_rung)
            Intobj=Jint_KS
        #Hybrid
        else
            DFTobj=JDFT(manual=manual_option, manual_func=manual_func, hybrid=true, libxc_functional=libxc_functional,
            functional_rung=func_rung)
            Intobj=Jint_hybridKS
        end
        #TODO: broken-symmetry spin-coupled case
    elseif WFtype == "UHF"
        unpaired_electrons= fragment.mult - 1
        paired_el=num_el-unpaired_electrons
        paired_el_half=paired_el/2
        numoccorbs_⍺=trunc(Int64,paired_el_half+unpaired_electrons)
        numoccorbs_β=trunc(Int64,paired_el_half)
        #TODO: broken-symmetry spin-coupled case
        DFTobj=nothing
        grid=false
        Intobj=Jint_HF
    else
        println("Unknown WFtype! Exiting.")
        exit()
    end

    #Print basic system properties
    if printlevel > 0
        print_system(num_el,fragment.prettyformula,E_ZZ,fragment.charge,fragment.mult,fragment.numatoms,unpaired_electrons)
    end
    print_geometry(fragment,printlevel)


    #############################
    # INTEGRALS
    #############################
    print_if_level("\nCalculating 1-electron integrals",1,printlevel) 
    #Calculating 1-electron integrals
    time_1el=@elapsed begin
        T = kinetic(bset)
        V = nuclear(bset)
        S = overlap(bset)
        Hcore = T + V #Combining T and V into Hcore
    end
    print_if_level("Time calculating 1-electron integrals: $time_1el",1,printlevel)
    #Overlap
    Sval,Svec = eigen(S)
    lowest_S_eigenval=minimum(Sval)
    #Transformation matrix
    SVAL_minhalf = Diagonal(Sval)^-0.5
    Stemp = SVAL_minhalf*transpose(Svec)
    S_minhalf = Svec * Stemp
    print_if_level("Calculating 2-electron integrals (tei_type: $tei_type)",1,printlevel)
    #Calculating two-electron integrals: tei_type: 4c, sparse4c
    time_tei=@elapsed tei = tei_calc(bset,tei_type,printlevel)
    print_if_level("Time calculating 2-electron integrals: $time_tei",1,printlevel)

    #Creating Integrals object (contains Hcore, TEI and basis set)
    if tei_type == "sparse4c"
        #Note: Takes 0.10-0.13 s for C16 RHF/STO-3G
        #due to list comprehension, speed up?
        integrals = Intobj(Hcore=Hcore,S=S, S_minhalf=S_minhalf, unique_indices=[i .+ 1 for i in tei[1]],
            values=tei[2], length=length(tei[2]),bset=bset,bset_atom_mapping=bset_atom_mapping,
            bf_atom_shell_map=bf_atom_shell_map,DFT=DFTobj, mlvalues=mlvalues, lvalues=lvalues)
        fock_algorithm="loop-sparse" #Only a string for printing
    else
        #To be deleted. Only kept on for RHF comparison
        integrals = Jint_4rank(Hcore=Hcore,tensor=tei,bset=bset,bset_atom_mapping=bset_atom_mapping,
        bf_atom_shell_map=bf_atom_shell_map, mlvalues=mlvalues, lvalues=lvalues)
        fock_algorithm="loop-4rank" #Only a string for printing
    end

    #############################
    # CREATING GRID
    #############################
    #grid always true for RKS/UKS
    # Also turn-on grid for RHF/UHF if density is requested
    if calc_density == true grid=true end
    if grid == true
        println("Creating grid")
        #TODO: Create gridkeywords: Grid1-7 or something?
        time_grid=@elapsed gridpoints,gridweights = numgrid_call("atomgrid",fragment,bset;
            radial_precision=1e-12,min_num_angular_points=434,max_num_angular_points=770)
        print_if_level("Number of gridpoints: $(length(gridpoints))",1,printlevel)
        print_if_level("Number of gridweights: $(length(gridweights))",1,printlevel)
        #NOTE: Gridpoint coords are in Bohrs
        #write_gridpoints_to_disk("gridpoints.xyz",gridpoints,fragment)
        #Adding gridpoints to object
        integrals.gridpoints=gridpoints
        integrals.gridweights=gridweights
        print_if_level("Time calculating grid: $time_grid",1,printlevel)

        #If grid is calculated we should also pre-calculate BFvalues at gridpoints
        BFvalues_calc(integrals) # occupies integrals.BFvalues
        #println("Size of object BFvalues:", varinfo(r"integrals.BFvalues"))

    else
        time_grid=0.0
    end
    #end

    #############################
    # CHOOSING FOCK ALGORITHM
    #############################
    #Multiple dispatch chooses Fock method based on integrals object type and number of P,F input matrices
    #TODO: Once density fitting and direct-SCF we need to revisit this

    #############################
    # GUESS
    #############################
    #Create initial guess for density matrix
    print_if_level("Providing guess: $guess",1,printlevel)
    if guess == "hcore"
        if WFtype=="RHF" || WFtype=="RKS" || WFtype=="i-DMFT"
            C = compute_core_guess(Hcore,S_minhalf)
            #P = makeP(C, dim, numoccorbs) #Calculate new P from C
            #println("P:", P)
            #Make density matrix using array of occupation numbers instead
            #NOTE: Useful for i-DMFT
            P = makeP2(C, dim, numoccorbs, orb_occupations) #Calculate new P from C
            #println("P2:", P2)
            println("Initial orb_occupations: ", orb_occupations)
            #println("Initial P: ", P)
        elseif WFtype=="UHF" || WFtype=="UKS"
            #UHF or UKS
            #A=zeros(dim,dim)
            #for i in 1:numoccorbs_⍺; A[i,i]=1 end
            #B=zeros(dim,dim)
            #for i in 1:numoccorbs_β; B[i,i]=1 end

            C = compute_core_guess(Hcore,S_minhalf)

            #C_⍺ = compute_core_guess(A*Hcore,S_minhalf)
            #C_β = compute_core_guess(B*Hcore,S_minhalf)
            #println("numoccorbs_⍺: $numoccorbs_⍺")
            #println("numoccorbs_β: $numoccorbs_β")
            #println("C_⍺: $C_⍺")
            #println("C_β: $C_β")
            P_⍺ = makeP(C, dim, numoccorbs_⍺, 1.0) #Calculate new P from C
            P_β = makeP(C, dim, numoccorbs_β, 1.0) #Calculate new P from C
            P=P_⍺+P_β
        else
            println("Unknown WFtype option.")
            exit()
        end
    else
        println("unknown guess")
        exit()
    end

    if printlevel > 0
        println("")
        print_calculation_settings(WFtype,basisset,dim,guess,tei_type,
            fock_algorithm,lowest_S_eigenval,levelshift,levelshift_val,lshift_thresh,
            damping,damping_val,damping_thresh,diis,diis_size,diis_startiter)
    end
    ##########################
    # SCF
    ##########################
    print_if_level("\nBeginning SCF iterations",1,printlevel)
    
    #Initializing some variables that will change during the iterations
    energy_old=0.0; energy=0.0; P_RMS=9999; finaliter=nothing
    Resultsdict=Dict()
    if WFtype=="RHF" || WFtype=="RKS" || WFtype=="i-DMFT"
        eps=zeros(dim)
        P_old=deepcopy(P)
    else
        eps_⍺=zeros(dim); eps_β=zeros(dim); 
        P_⍺_old=deepcopy(P_⍺); P_β_old=deepcopy(P_β)
    end

    #Initializing levelshift, damping and DIIS settings
    if levelshift == true global levelshift_flag = true else levelshift_flag = false end
    if damping == true global damping_flag = true; else damping_flag = false end
    if diis == true
        diisobj=JDIIS(active=true,diis_size=diis_size,diis_startiter=diis_startiter,DIISBfac=DIISBfac)
    else
        diisobj=JDIIS(active=false) #Creating dummy JDIIS object
    end
    #Print 1-elec matrices before SCF if printlevel > 2
    if printlevel > 2 print_1_elmatrices(T,V,Hcore,S,C,P) end

    #SCF loop beginning
    if printlevel == 1 @printf("%4s%15s%16s%14s%14s%8s%8s%8s%10s\n", "Iter", "Energy", "deltaE", "RMS-DP", "Max-DP", "Lshift", "Damp", "DIIS", "Max[F,P]") end

    #println("P_⍺")
    #print_matrix(P_⍺)
    #println("P_β")
    #print_matrix(P_β)
    #println("P")
    #print_matrix(P)

    time_scf=@elapsed for iter in 1:maxiter
        if printlevel > 1 print_iteration_header(iter) end
        if WFtype == "RHF" || WFtype=="RKS" || WFtype == "i-DMFT"

            #Optional damping of P before making Fock
            P = damping_control(P,P_old,damping,damping_val,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=damping_thresh)
            F = Fock(P,dim,integrals) #Make Fock-matrix
            #println("F: ")
            #print_matrix(F)
            #println("C: ")
            #print_matrix(C)
            #println("P: ")
            #print_matrix(P)
            FP_comm = FP_commutator(F,P,S,S_minhalf) #Calculating [F,P] commutator (for DIIS)
            #Possible levelshifting of Fock before diagonalization
            F = levelshift_control(F,levelshift,levelshift_val,numoccorbs,dim,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=lshift_thresh)
            F′ = transpose(S_minhalf)*F*S_minhalf #Transform Fock matrix
            energy =  calc_energy(E_ZZ,integrals,F,P) #Calculate RHF/RKS energy after Fock formation
            # Possible DIIS extrapolation of F′ matrix before diagonalization
            F′ = diis_control(diisobj,F′,energy,FP_comm,iter,printlevel)
            
            eps, C′ = eigen(F′) #Diagonalize transformed Fock to get eps and C'
            eps=real(eps) #Eigenvalues can be complex (tiny imaginary part) so taking real part
            C′=real(C′) #Eigenvectors can be complex (tiny imaginary part) so taking real part 
            #print_matrix(C′)

            C = S_minhalf*C′ # Get C from C'

            P_old=deepcopy(P) #Keep copy of old P

            if WFtype == "i-DMFT"
                println("i-DMFT active")
                #i-DMFT:
                #iDMFT_b=0.18121047 #fitted. Just used to shift energy
                #Correct??
                #mu=[e+iDMFT_kappa*(log(n)-log(1-n)) for (e,n) in zip(eps,orb_occupations)]
                #print("Mu:", mu)
                mu=-0.1
                #Update occupations
                println("Old orb_occupations:", orb_occupations)
                
                #sum([1/(1+exp((e-mu)/iDMFT_kappa)) for e in eps])
                println("eps:", eps)
                println("iDMFT_kappa:", iDMFT_kappa)
                println("energy:", energy)
                #for mu in [0.2921,0.2922,0.2923,0.2924]
                    orb_occupations=[1/(1+exp((e-mu)/iDMFT_kappa)) for e in eps]
                    println("mu: $mu orb_occupations: $orb_occupations  Sum: $(sum(orb_occupations))")
                #end
                print("New orb_occupations:", orb_occupations)
                print("New orb_occupations sum :", sum(orb_occupations))
                exit()
                #Update P from new C and new orb_occupations
                P = makeP2(C, dim, numoccorbs, orb_occupations)
            else
                #Switch once confirmed
                P = makeP2(C, dim, numoccorbs, orb_occupations)
                #P = makeP(C, dim, numoccorbs) #Calculate new P from C
            end

            if printlevel > 2 write_matrices(F,C,P,S) end
        else
            #Possible damping of P matrices before making Fock
            P_⍺ = damping_control(P_⍺,P_⍺_old,damping,damping_val,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=damping_thresh)
            P_β = damping_control(P_β,P_β_old,damping,damping_val,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=damping_thresh)
            
            #Solve ⍺ part
            F_⍺ = Fock(P_⍺,P_β,dim,integrals,1) #Update Fock-matrix alpha
            FP_comm_⍺ = FP_commutator(F_⍺,P_⍺,S,S_minhalf) #Calculating [F,P] commutator (for DIIS)
            F_⍺ = levelshift_control(F_⍺,levelshift,levelshift_val,numoccorbs_⍺,dim,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=lshift_thresh)
            F′_⍺ = transpose(S_minhalf)*F_⍺*S_minhalf #Transform Fock matrix
            #Solve β part
            F_β = Fock(P_β,P_⍺,dim,integrals,2) #Update Fock-matrix beta
            FP_comm_β = FP_commutator(F_β,P_β,S,S_minhalf) #Calculating [F,P] commutator (for DIIS)
            F_β = levelshift_control(F_β,levelshift,levelshift_val,numoccorbs_β,dim,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=lshift_thresh)
            F′_β = transpose(S_minhalf)*F_β*S_minhalf #Transform Fock matrix
            energy =  calc_energy(E_ZZ,integrals,F_⍺,P_⍺,F_β,P_β) #Calculate UHF/UKS energy after Fock formation
            #energy = E_ZZ + 0.5 * tr((Hcore+F_⍺)*P_⍺) + 0.5 * tr((Hcore+F_β)*P_β)  #Calculate energy

            # Possible DIIS extrapolation of F′ matrix before diagonalization
            if diisobj.active == true
                Fcomb=hcat(F′_⍺,F′_β) #Combining F′ matrices into super matrix
                FP_comm_⍺β=hcat(FP_comm_⍺,FP_comm_β) #Combining error matrices into super matrix
                Fcombextrap = diis_control(diisobj,Fcomb,energy,FP_comm_⍺β,iter,printlevel) #Do DIIS to get extrapolated super-matrix
                F′_⍺ = Fcombextrap[1:size(F′_⍺)[1],1:size(F′_⍺)[2]] #split back into F′_⍺
                F′_β = Fcombextrap[1:size(F′_β)[1],size(F′_β)[1]+1:size(Fcombextrap)[2]] #split back into F′_β
            end
            #Now continuing
            eps_⍺, C′_⍺ = eigen(F′_⍺) #Diagonalize transformed Fock to get eps and C'
            C_⍺ = S_minhalf*C′_⍺ # Get C from C'
            P_⍺_old=deepcopy(P_⍺) #Keep copy of old P
            P_⍺ = makeP(C_⍺, dim, numoccorbs_⍺, 1.0) #Calculate new P from C
            eps_β, C′_β = eigen(F′_β) #Diagonalize transformed Fock to get eps and C'
            C_β = S_minhalf*C′_β # Get C from C'
            P_β_old=deepcopy(P_β) #Keep copy of old P
            P_β = makeP(C_β, dim, numoccorbs_β, 1.0) #Calculate new P from C
            
            #Combined density and spin density matrices
            P=P_⍺+P_β
            P_old=P_⍺_old+P_β_old
            P_⍺_β=P_⍺-P_β

            #Average FP_comm for UHF
            #????
            FP_comm = (FP_comm_⍺ + FP_comm_β)/2
        end

        ##########################
        # CONVERGENCE CHECK
        ##########################
        deltaE = energy-energy_old
        #NOTE: P_RMS is still off
        P_RMS, P_MaxE = deltaPcheck(P, P_old)
        energy_old=energy

        #Printing per iteration
        #NOTE: For UHF, only checking diisobj.diis_flag 
        iteration_printing(iter,printlevel,energy,deltaE,energythreshold,P_RMS,rmsDP_threshold,
            P_MaxE,maxDP_threshold,levelshift_flag,damping_flag,diisobj.diis_flag,FP_comm)

        #NOTE: UHF, only checking convergence for diisobj_⍺ !
        #TODO: Make more general
        if check_for_convergence(deltaE,energythreshold,FP_comm,diis_error_conv_threshold,iter,printlevel) == true
            finaliter=iter
            if printlevel > 0
                if WFtype == "RHF" || WFtype == "RKS"
                    print_energy_contributions(energy,Hcore,P,T,E_ZZ,integrals.E_xc,V,num_el)
                    #Printing of matrices if requested. 
                    if print_final_matrices == true 
                        write_matrix_to_file(F,"Fmatrix")
                        write_matrix_to_file(C,"Cmatrix")
                        write_matrix_to_file(P,"Pmatrix")
                    end
                else
                    print_energy_contributions(energy,Hcore,P,T,E_ZZ,integrals.E_xc,V,num_el)
                    if print_final_matrices == true 
                        write_matrix_to_file(P,"Pmatrix")
                        write_matrix_to_file(C_⍺,"C_a_matrix")
                        write_matrix_to_file(C_β,"C_b_matrix")
                        write_matrix_to_file(F_⍺,"F_a_matrix")
                        write_matrix_to_file(F_β,"F_b_matrix")
                        write_matrix_to_file(P_⍺,"P_a_matrix")
                        write_matrix_to_file(P_β,"P_b_matrix")
                    end
                end
            end
            break
        end

        if iter == maxiter
            println("Failed to converge in $maxiter iterations!")
            #Gathering results into Dict and returning
            Resultsdict["energy"] = "None"
            Resultsdict["finaliter"]="Failed ($maxiter iters)"
            return Resultsdict
        end
    end
    print_if_level("Time calculating SCF: $time_scf",1,printlevel)

    #################################################
    # ORBITALS, POPULATION ANALYSIS AND PROPERTIES
    #################################################
    #SCF loop done. Calculate and print properties
    #DENSITY
    #println("Density matrix:")
    #pretty_table(P; header=[string(i) for i in 1:size(P)[2]], tf = tf_matrix, 
    #    show_row_number = true, formatters = ft_printf("%3.1f"))

    if nopop == false
        if grid == true
            
            #TODO: Case HF, calculate density here, already done if DFT
            #TODO: Put in table?
            #println("Now creating density")
            #println("Integrated no. of electrons: $N")
            #@time ⍴ = create_density(P,integrals)
            #@time N = integrate_density(⍴,integrals.gridweights)
            #Printing integrated no. of electrons
            println("Integrated no. of electrons: $(integrals.N_⍺+integrals.N_β)")
            println("Integrated no. of ⍺ electrons: $(integrals.N_⍺)")
            println("Integrated no. of β electrons: $(integrals.N_β)")
            Resultsdict["int_electrons"] = integrals.N_⍺ + integrals.N_β
            Resultsdict["int_electrons_⍺"] = integrals.N_⍺
            Resultsdict["int_electrons_β"] = integrals.N_β
        end

    #ORBITALS AND POPULATION ANALYSIS
        if WFtype=="RHF" || WFtype=="RKS" || WFtype == "i-DMFT"

            #Orbitals
            #TODO: Needed??
            #occupations=makeoccupationarray(numoccorbs,dim,2.0) #Occupation array
            #Now defined from beginning
            occupations=orb_occupations
            #Mulliken
            charges = mulliken(S,P,bset,fragment.elems)
            P_⍺_β=zeros(dim,dim) #dummy spin-density
            #Mayer
            MBOs = Mayer_BO(S,P,P_⍺_β, bset_atom_mapping)
            #Print
            if printlevel > 0
                print_MO_energies(occupations,eps)
                print_Mulliken(charges,fragment.elems)
                print_Mayer_analysis(MBOs,fragment.elems)
            end
        else
            #Orbitals
            occupations_⍺=makeoccupationarray(numoccorbs_⍺,dim,1.0) #Occupation array
            occupations_β=makeoccupationarray(numoccorbs_β,dim,1.0) #Occupation array
            #Mulliken
            P_⍺_β=P_⍺_β
            charges, spinpops = mulliken(S,P,P_⍺_β,bset,fragment.elems)
            #Mayer
            MBOs = Mayer_BO(S,P,P_⍺_β, bset_atom_mapping)
            
            if printlevel > 0
                print_MO_energies(occupations_⍺,occupations_β,eps_⍺,eps_β)
                print_Mulliken(charges,fragment.elems,spinpops)
                print_Mayer_analysis(MBOs,fragment.elems)
            end
        end
    end

    #####################################
    # PRINT FINAL RESULTS
    #####################################
    if printlevel > 0
        print_final_results(energy,fragment,num_el,basisset,WFtype,fock_algorithm,finaliter)
        if WFtype == "i-DMFT"
            println("iDMFT_kappa: ",iDMFT_kappa)
            println("iDMFT_b: ", iDMFT_b)
            println("Final occupations: ", orb_occupations)
            println("Final sum of occupations:", sum(orb_occupations))
            S_term = -sum([n*log(n)+(1-n)*log(1-n) for n in orb_occupations])
            println("S_term: ", S_term)
            E_cum = -iDMFT_kappa*S_term-iDMFT_b
            println("E_cum: ", E_cum)
            DMFT_energy= energy + E_cum
            println("Final i-DMFT energy: ", DMFT_energy)
        end
    end
    #Gathering results into Dict and returning.
    #TODO: Add more here. 
    #Mulliken charges/spinpops, Mayer MBOs, dipole, all energy contributions,
    #orbital energies, occupation numbers, WFtype, numelectrons
    
    Resultsdict["energy"]=energy
    Resultsdict["finaliter"]=finaliter

    #TIMINGS: TODO properly
    if printlevel >0
        println("Timings:")
        print_if_level("Time calculating 1-electron integrals: $time_1el",1,printlevel)
        print_if_level("Time calculating 2-electron integrals: $time_tei",1,printlevel)
        print_if_level("Time calculating SCF: $time_scf",1,printlevel)
        print_if_level("Time calculating grid: $time_grid",1,printlevel)
        all_times=[time_1el,time_tei,time_scf,time_grid]
        sum_of_all_times=sum(all_times)
        println("Sum of all recorded times:", sum_of_all_times)
    end
end #end of timing block
Resultsdict["time"]=totaltime #addint total time
    return Resultsdict
end



####################
# Some structs
###################

"""
JDFT: Jotunn DFT functional object
"""
Base.@kwdef mutable struct JDFT
    manual::Bool
    hybrid::Bool
    libxc_functional::Any
    functional_rung::Int64
    alpha::Float64=0.0 #HF exchange
    #For manual=True, manual_func distinguishes between options
    manual_func::String="LDA"
end

"""
Jint_DMFT: Jotunn DMFT object
"""
Base.@kwdef mutable struct Jint_DMFT

end

"""
Jint_HF: Jotunn integral struct in sparse form for RHF/UHF
"""
Base.@kwdef mutable struct Jint_HF
    #1-electron integrals
    Hcore::Matrix{Float64}
    S::Matrix{Float64}
    S_minhalf::Matrix{Float64}
    #TWO-electron integrals
    #Original unique list of integral indices from GaussianBasis but with 1-indexing
    unique_indices::Vector{NTuple{4, Int16}}=[]
    #Integral values
    values::Vector{Float64}=[]
    length::Int64=0
    #Basis set object and data
    bset::BasisSet
    bset_atom_mapping::Vector{Int64}
    bf_atom_shell_map::Vector{Tuple{Int64, Int64}}
    mlvalues::Vector{Int64}
    lvalues::Vector{Int64}
    #Gridpoints (used for density generation in HF)
    gridpoints::Vector{Tuple{Float64, Float64, Float64}}=[]
    gridweights::Vector{Float64}=[]
    #Pre-calculated basis function values at gridpoints
    BFvalues::Array{Float64}=[]
    ∇BFvalues::Array{Float64}=[] #derivative of BF values at gridpoints
    #
    DFT::Nothing
    E_xc::Float64=0.0
end

"""
Jint_KS: Jotunn integral struct in sparse form
"""
Base.@kwdef mutable struct Jint_KS
    #1-electron integrals
    Hcore::Matrix{Float64}
    S::Matrix{Float64}
    S_minhalf::Matrix{Float64}
    #TWO-electron integrals
    #Original unique list of integral indices from GaussianBasis but with 1-indexing
    unique_indices::Vector{NTuple{4, Int16}}=[]
    #Integral values
    values::Vector{Float64}=[]
    length::Int64=0
    #Basis set object and data
    bset::BasisSet
    bset_atom_mapping::Vector{Int64}
    bf_atom_shell_map::Vector{Tuple{Int64, Int64}}
    mlvalues::Vector{Int64}
    lvalues::Vector{Int64}
    #Gridpoints (only used by DFT or density generation)
    gridpoints::Vector{Tuple{Float64, Float64, Float64}}=[]
    gridweights::Vector{Float64}=[]
    #Pre-calculated basis function values at gridpoints
    BFvalues::Array{Float64}=[]
    ∇BFvalues::Array{Float64}=[] #derivative of BF values at gridpoints
    #DFT functional
    DFT::JDFT
    E_xc::Float64=0.0
    E_xc_⍺::Float64=0.0 #only used in UKS
    E_xc_β::Float64=0.0 #only used in UKS
    J::Array{Float64}=[]
    N_⍺::Float64=0.0
    N_β::Float64=0.0
end

"""
Jint_hybridKS: Jotunn integral struct in sparse form for hybrid Kohn-Sham
"""
Base.@kwdef mutable struct Jint_hybridKS
    #1-electron integrals
    Hcore::Matrix{Float64}
    S::Matrix{Float64}
    S_minhalf::Matrix{Float64}
    #TWO-electron integrals
    #Original unique list of integral indices from GaussianBasis but with 1-indexing
    unique_indices::Vector{NTuple{4, Int16}}=[]
    #Integral values
    values::Vector{Float64}=[]
    length::Int64=0
    #Basis set object and data
    bset::BasisSet
    bset_atom_mapping::Vector{Int64}
    bf_atom_shell_map::Vector{Tuple{Int64, Int64}}
    mlvalues::Vector{Int64}
    lvalues::Vector{Int64}
    #Gridpoints (only used by DFT or density generation)
    gridpoints::Vector{Tuple{Float64, Float64, Float64}}=[]
    gridweights::Vector{Float64}=[]
    #Pre-calculated basis function values at gridpoints
    BFvalues::Array{Float64}=[]
    ∇BFvalues::Array{Float64}=[] #derivative of BF values at gridpoints
    #DFT functional
    DFT::JDFT
    E_xc::Float64=0.0
    N_⍺::Float64=0.0
    N_β::Float64=0.0
end

"""
Jint: Jotunn integral struct in 4-rank form for RHF/UHF. Rarely used. 
Only kept for RHF speed-comparisons for now
To be deleted
"""
Base.@kwdef mutable struct Jint_4rank
    #1-electron integrals
    Hcore::Matrix{Float64}
    #TWO-electron integrals
    tensor::Array{Float64, 4}
    #Basis set object and data
    bset::BasisSet
    bset_atom_mapping::Vector{Int64}
    bf_atom_shell_map::Vector{Tuple{Int64, Int64}}
    mlvalues::Vector{Int64}
    lvalues::Vector{Int64}
    #Gridpoints (only used by DFT or density generation)
    gridpoints::Vector{Tuple{Float64, Float64, Float64}}=[]
    gridweights::Vector{Float64}=[]
    #Pre-calculated basis function values at gridpoints
    BFvalues::Array{Float64}=[]
    ∇BFvalues::Array{Float64}=[] #derivative of BF values at gridpoints
end




#Calculate energy function with multiple methods (multiple dispatch)
#Dispatches based on type of integrals object and the number of P,F matrices provided (R vs. U)
#RHF/RKS versions
##################
function calc_energy(E_ZZ::Float64,integrals::Jint_HF,F::Matrix{Float64},P::Matrix{Float64})
    return E_ZZ + 0.5 * tr((integrals.Hcore+F)*P)
end
#RHF-4c version (to be deleted)
function calc_energy(E_ZZ::Float64,integrals::Jint_4rank,F::Matrix{Float64},P::Matrix{Float64})
    return E_ZZ + 0.5 * tr((integrals.Hcore+F)*P)
end
#RKS version
function calc_energy(E_ZZ::Float64,integrals::Jint_KS,F::Matrix{Float64},P::Matrix{Float64})
    println("integrals.E_xc:", integrals.E_xc)
    #return E_ZZ + 0.5 * tr((integrals.Hcore+integrals.Hcore+integrals.J)*P) + integrals.E_xc
    # Nuc-repulsion + Hcore + J + Exc
    return E_ZZ + tr(integrals.Hcore*P) + 0.5 * tr(integrals.J*P) + integrals.E_xc
end
#RKS-hybrid version
#function calc_energy(E_ZZ::Float64,integrals::Jint_hybridKS,F::Matrix{Float64},P::Matrix{Float64})
#    return E_ZZ + 0.5 * tr((integrals.Hcore+integrals.Hcore+integrals.J)*P) + integrals.E_xc
#end



#UHF/UKS versions
##################
#UHF
function calc_energy(E_ZZ::Float64,integrals::Jint_HF,F_⍺::Matrix{Float64},
    P_⍺::Matrix{Float64},F_β::Matrix{Float64},P_β::Matrix{Float64})
    return E_ZZ + 0.5 * tr((integrals.Hcore+F_⍺)*P_⍺) + 0.5 * tr((integrals.Hcore+F_β)*P_β)  #Calculate energy
end
#UHF-4c version (to be deleted)
function calc_energy(E_ZZ::Float64,integrals::Jint_4rank,F_⍺::Matrix{Float64},
    P_⍺::Matrix{Float64},F_β::Matrix{Float64},P_β::Matrix{Float64})
    return E_ZZ + 0.5 * tr((integrals.Hcore+F_⍺)*P_⍺) + 0.5 * tr((integrals.Hcore+F_β)*P_β)  #Calculate energy
end
#UKS
function calc_energy(E_ZZ::Float64,integrals::Jint_KS,F_⍺::Matrix{Float64},
    P_⍺::Matrix{Float64},F_β::Matrix{Float64},P_β::Matrix{Float64})
    P=P_⍺+P_β
    return E_ZZ + tr(integrals.Hcore*P) + 0.5 * tr(integrals.J*P) + integrals.E_xc
    #return E_ZZ + 0.5 * tr((integrals.Hcore+F_⍺)*P_⍺) + 0.5 * tr((integrals.Hcore+F_β)*P_β) #Calculate energy
end
#UKS-hybrid version
#function calc_energy(E_ZZ::Float64,integrals::Jint_hybridKS,F_⍺::Matrix{Float64},
#    P_⍺::Matrix{Float64},F_β::Matrix{Float64},P_β::Matrix{Float64})
#    return E_ZZ + 0.5 * tr((integrals.Hcore+F_⍺)*P_⍺) + 0.5 * tr((integrals.Hcore+F_β)*P_β)  #Calculate energy
#end