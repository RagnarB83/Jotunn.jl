export jHF

"""
Jint: Jotunn integral struct in sparse form
"""
struct Jint_sparse
    #Original unique list of integral indices from GaussianBasis but with 1-indexing
    unique_indices::Vector{NTuple{4, Int16}}
    #Integral values
    values::Vector{Float64}
    length::Int64
end


"""
Jotunn: jHF: a RHF/UHF program
"""

function jHF(fragment, basisset="sto-3g"; HFtype::String="RHF", guess::String="hcore", basisfile::String="none", maxiter::Int64=120, 
    print_final_matrices::Bool=false, rmsDP_threshold::Float64=5e-9, maxDP_threshold::Float64=1e-7, tei_type::String="sparse4c",
    energythreshold::Float64=1e-8, debugprint::Bool=false, fock_algorithm::String="loop", 
    levelshift::Bool=false, levelshift_val::Float64=0.10, lshift_thresh::Float64=0.01,
    damping::Bool=true, damping_val::Float64=0.4, damping_thresh::Float64=0.01,
    diis::Bool=false, diis_size::Int64=7, diis_thresh::Float64=0.01,
    printlevel::Int64=1, fock4c_speedup::String="simd")

    print_program_header()
    global debugflag  = debugprint


    #Removing old matrix files if present
    if debugprint == true || print_final_matrices == true
        for matrixfile in ["Fmatrix","Cmatrix","Pmatrix","C_a_matrix","C_b_matrix","F_a_matrix",
            "F_b_matrix","P_a_matrix","P_b_matrix"]
            rm(matrixfile, force=true)
        end
    end

    ##########################
    # BASIC SYSTEM SETUP
    ##########################
    #Num. electrons and nuc-nuc repulsion from total charge and nuclear charges
    sum_nuccharge=sum([elem_to_nuccharge(el) for el in fragment.elems])
    num_el=floor(Int64,sum_nuccharge-fragment.charge)
    E_ZZ=nuc_nuc_repulsion(fragment.elems,fragment.coords)

    #Check whether chosen multiplicity makes sense before continuing
    check_multiplicity(num_el,fragment.charge,fragment.mult)
    if fragment.mult != 1
        println("RHF and multiplicity > 1 is not possible. Switching to UHF")
        HFtype="UHF"
    end
    #RHF vs. UHF
    if HFtype == "RHF"
        #Num occ orbitals in RHF
        numoccorbs=Int64(num_el/2)
        unpaired_electrons=0
    elseif HFtype == "UHF"
        unpaired_electrons= fragment.mult - 1
        paired_el=num_el-unpaired_electrons
        paired_el_half=paired_el/2
        numoccorbs_⍺=trunc(Int64,paired_el_half+unpaired_electrons)
        numoccorbs_β=trunc(Int64,paired_el_half)
        #TODO: broken-symmetry spin-coupled case
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
    #Basis set object creation by GaussianBasis
    bset = basis_set_create(basisset,fragment.elems,fragment.coords; basisfile=basisfile)
    dim = bset.nbas
    #Simple array of atom indices that maps onto bfs
    bset_atom_mapping = bf_atom_mapping(bset)
    println("\nCalculating 1-electron integrals")
    #Calculating 1-electron integrals
    T = kinetic(bset)
    V = nuclear(bset)
    S = overlap(bset)
    Hcore = T + V #Combining T and V into Hcore
    
    #Overlap diagonalization
    Sval,Svec = eigen(S)
    lowest_S_eigenval=minimum(Sval)
    #Transformation matrix
    SVAL_minhalf = Diagonal(Sval)^-0.5
    Stemp = SVAL_minhalf*transpose(Svec)
    S_minhalf = Svec * Stemp

    println("Calculating 2-electron integrals (tei_type: $tei_type)")
    #Calculating two-electron integrals: tei_type: 4c, sparse4c
    @time tei = tei_calc(bset,tei_type)

    ##########################
    # CHOOSING FOCK ALGORITHM
    ##########################
    #Choosing Fock algorithm (based on RHF vs. UHF, 2el-int-type etc.)
    Fock,fock_algorithm = choose_Fock(HFtype,fock_algorithm,dim,tei_type)

    ##########################
    # GUESS
    ##########################
    #Create initial guess for density matrix
    println("Providing guess for density matrix")
    if guess == "hcore"
        println("Hcore guess used")
        #Setting P to 0. Means that Fock matrix becomes F = Hcore + 0. See Fock functions.
        if HFtype=="RHF"
            P = zeros(dim,dim)
            guess_energy= E_ZZ + 0.5 * tr((Hcore+Hcore)*P) #Calculate guess energy
            println("Energy of guess: $guess_energy Eh\n")
        else
            P_⍺ = zeros(dim,dim)
            P_β = zeros(dim,dim)
        end
    else
        println("unknown guess")
        exit()
    end
    
    print_calculation_settings(HFtype,basisset,dim,guess,tei_type,
        fock_algorithm,lowest_S_eigenval,levelshift,levelshift_val,lshift_thresh,
        damping,damping_val,damping_thresh,diis,diis_size,diis_thresh)

    ##########################
    # SCF
    ##########################
    println("\nBeginning SCF iterations")
    #Initializing some variables tha will change during the iterations
    energy_old=0.0; energy=0.0; P_RMS=9999; finaliter=nothing
    if HFtype=="RHF" 
        eps=zeros(dim)
        P_old=deepcopy(P)
    else
        eps_⍺=zeros(dim); eps_β=zeros(dim); 
        P_⍺_old=deepcopy(P_⍺); P_β_old=deepcopy(P_β)
    end
    #Initializing levelshift, damping and DIIS flags
    if levelshift == true global levelshift_flag = true else levelshift_flag = false end
    if damping == true global damping_flag = true; else damping_flag = false end
    if diis == true global diis_flag = true; diis_error_matrices=nothing; Fockmatrices=nothing; 
    else global diis_flag = false;diis_error_matrices=nothing; Fockmatrices=nothing end
    
    #SCF loop beginning
    if printlevel == 1 @printf("%4s%15s%16s%14s%14s%8s%8s%8s\n", "Iter", "Energy", "deltaE", "RMS-DP", "Max-DP", "Lshift", "Damp", "DIIS") end
    
    @time for iter in 1:maxiter
        if printlevel > 1 print_iteration_header(iter) end
        if HFtype == "RHF"
            #Possible damping of P before making Fock
            P = damping_control(P,P_old,damping,damping_val,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=damping_thresh)
            @time F = Fock(Hcore,P,dim,tei) #Make Fock-matrix

            #Possible levelshifting of Fock before diagonalization
            F = levelshift_control(F,levelshift,levelshift_val,numoccorbs,dim,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=lshift_thresh)
            F′ = transpose(S_minhalf)*F*S_minhalf #Transform Fock matrix

            # Possible DIIS extrapolation of F′ matrix before diagonalization
            F′ = diis_control(F′,F,diis_error_matrices,Fockmatrices,diis,diis_size,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=diis_thresh)

            eps, C′ = eigen(F′) #Diagonalize transformed Fock to get eps and C'
            C = S_minhalf*C′ # Get C from C'
            P_old=deepcopy(P) #Keep copy of old P
            P = makedensity(C, dim, numoccorbs) #Calculate new P from C
            energy = E_ZZ + 0.5 * tr((Hcore+F)*P) #Calculate energy

            if debugprint == true write_matrices(F,C,P) end
        else
            #Possible damping of P matrices before making Fock
            P_⍺ = damping_control(P_⍺,P_⍺_old,damping,damping_val,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=damping_thresh)
            P_β = damping_control(P_β,P_β_old,damping,damping_val,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=damping_thresh)
            
            #Solve ⍺ part
            F_⍺ = Fock(Hcore,P_⍺,P_β,dim,tei) #Update Fock-matrix alpha
            F_⍺ = levelshift_control(F_⍺,levelshift,levelshift_val,numoccorbs_⍺,dim,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=lshift_thresh)
            F′_⍺ = transpose(S_minhalf)*F_⍺*S_minhalf #Transform Fock matrix
            eps_⍺, C′_⍺ = eigen(F′_⍺) #Diagonalize transformed Fock to get eps and C'
            C_⍺ = S_minhalf*C′_⍺ # Get C from C'
            P_⍺_old=deepcopy(P_⍺) #Keep copy of old P
            P_⍺ = makedensity(C_⍺, dim, numoccorbs_⍺, 1.0) #Calculate new P from C
            #Solve β part
            F_β = Fock(Hcore,P_β,P_⍺,dim,tei) #Update Fock-matrix beta
            F_β = levelshift_control(F_β,levelshift,levelshift_val,numoccorbs_β,dim,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=lshift_thresh)
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

            #if debugprint == true write_matrices(F,C,P) end
        end

        ##########################
        # CONVERGENCE CHECK
        ##########################
        deltaE = energy-energy_old
        P_RMS, P_MaxE = deltaPcheck(P, P_old)
        energy_old=energy

        #Printing per iteration
        iteration_printing(iter,printlevel,energy,deltaE,energythreshold,P_RMS,rmsDP_threshold,
            P_MaxE,maxDP_threshold,levelshift_flag,damping_flag,diis_flag)

        if P_MaxE < maxDP_threshold && P_RMS < rmsDP_threshold && abs(deltaE) < energythreshold
            print(Crayon(foreground = :green, bold = true), 
                "\n                              SCF converged in $iter iterations! Hell yeah! 🎉\n\n",Crayon(reset=true))
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
                    write_matrix_to_file(C_⍺,"C_a_matrix")
                    write_matrix_to_file(C_β,"C_b_matrix")
                    write_matrix_to_file(F_⍺,"F_a_matrix")
                    write_matrix_to_file(F_β,"F_b_matrix")
                    write_matrix_to_file(P_⍺,"P_a_matrix")
                    write_matrix_to_file(P_β,"P_b_matrix")
                end
            end
            break #Break from loop
        end
        if iter == maxiter
            println("Failed to converge in $maxiter iterations!")
            #Gathering results into Dict and returning
            Resultsdict=Dict("energy"=>"None","finaliter"=>"Failed ($maxiter iters)")
            return Resultsdict
        end
    end
    #################################################
    # ORBITALS, POPULATION ANALYSIS AND PROPERTIES
    #################################################
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
    # PRINT FINAL RESULTS
    #####################################
    print_final_results(energy,fragment,num_el,basisset,HFtype,fock_algorithm,finaliter)

    #Gathering results into Dict and returning.
    #TODO: Add more here. 
    #Mulliken charges/spinpops, Mayer MBOs, dipole, all energy contributions,
    #orbital energies, occupation numbers, HFtype, numelectrons
    Resultsdict=Dict("energy"=>energy,"finaliter"=>finaliter)

    #TIMINGS: TODO

    return Resultsdict


end