export jSCF


"""
jSCF: the Jotunn SCF program (RHF,UHF,RKS,UKS)
"""

function jSCF(fragment, basisset="sto-3g"; WFtype::String="RHF", functional::String="none",
    guess::String="hcore", 
    basisfile::String="none", maxiter::Int64=120, print_final_matrices::Bool=false, 
    rmsDP_threshold::Float64=5e-9, maxDP_threshold::Float64=1e-7, tei_type::String="sparse4c",
    energythreshold::Float64=1e-8, debugprint::Bool=false, fock_algorithm::String="loop", 
    levelshift::Bool=false, levelshift_val::Float64=0.10, lshift_thresh::Float64=0.01,
    damping::Bool=true, damping_val::Float64=0.4, damping_thresh::Float64=0.01,
    diis::Bool=true, diis_size::Int64=5, diis_startiter::Int64=4, DIISBfac::Float64=1.05,
    diis_error_conv_threshold::Float64=5e-7, calc_density::Bool=false,
    printlevel::Int64=1, fock4c_speedup::String="simd")

    #Timing whole function
    totaltime=@elapsed begin

    print_program_header(printlevel)
    global debugflag  = debugprint


    #Removing old matrix files if present
    if debugprint == true || print_final_matrices == true
        for matrixfile in ["Fmatrix","Cmatrix","Pmatrix","C_a_matrix","C_b_matrix","F_a_matrix",
            "F_b_matrix","P_a_matrix","P_b_matrix"]
            rm(matrixfile, force=true)
        end
    end

    #############################
    # BASIC SYSTEM SETUP
    #############################
    #Num. electrons and nuc-nuc repulsion from total charge and nuclear charges
    sum_nuccharge=sum([elem_to_nuccharge(el) for el in fragment.elems])
    num_el=floor(Int64,sum_nuccharge-fragment.charge)
    E_ZZ=nuc_nuc_repulsion(fragment.elems,fragment.coords)

    #Check whether chosen multiplicity makes sense before continuing
    check_multiplicity(num_el,fragment.charge,fragment.mult)
    if fragment.mult != 1 && WFtype=="RHF"
        print_if_level("RHF and multiplicity > 1 is not possible. Switching to UHF",1,printlevel)
        WFtype="UHF"
    end
    #RHF vs. UHF
    if WFtype == "RHF"
        #Num occ orbitals in RHF
        numoccorbs=Int64(num_el/2)
        unpaired_electrons=0
        DFTobj=nothing
        grid=false
    elseif WFtype == "RKS"
        #restricted Kohn-Sham
        #Num occ orbitals in RHF
        numoccorbs=Int64(num_el/2)
        unpaired_electrons=0
        if functional == "none"
            println("WFtype is RKS but no functional chosen. Exiting.")
            exit()
        end
        grid=true
        #TODO: Choose and create DFT functional object here
        DFTobj=JDFT(functional,false)
    elseif WFtype == "UKS"
        println("UKS not yet ready")
        exit()
        DFTobj=JDFT(functional,false)
        grid=true
        unpaired_electrons= fragment.mult - 1
        paired_el=num_el-unpaired_electrons
        paired_el_half=paired_el/2
        numoccorbs_⍺=trunc(Int64,paired_el_half+unpaired_electrons)
        numoccorbs_β=trunc(Int64,paired_el_half)
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
    else
        println("Unknown WFtype! Exiting.")
        exit()
    end

    #Print basic system properties
    if printlevel > 0
        print_system(num_el,fragment.prettyformula,E_ZZ,fragment.charge,fragment.mult,fragment.numatoms,unpaired_electrons)
    end
    #############################
    # INTEGRALS
    #############################
    #Setting up 1-electron and 2-electron integrals
    print_if_level("Integrals provided via GaussianBasis.jl library",1,printlevel)  
    #Basis set object creation by GaussianBasis
    bset = basis_set_create(basisset,fragment.elems,fragment.coords; basisfile=basisfile,printlevel)
    dim = bset.nbas
    #Simple array of atom indices that maps onto bfs
    bset_atom_mapping = bf_atom_mapping(bset)
    #Create array of tuples with atom,shell info for each BF
    bf_atom_shell_map=create_bf_shell_map(bset)

    print_if_level("\nCalculating 1-electron integrals",1,printlevel) 
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
    print_if_level("Calculating 2-electron integrals (tei_type: $tei_type)",1,printlevel)
    #Calculating two-electron integrals: tei_type: 4c, sparse4c
    #
    time_tei=@elapsed tei = tei_calc(bset,tei_type,printlevel)
    print_if_level("Time calculating 2-electron integrals: $time_tei",1,printlevel)

    #Creating Integrals object (contains Hcore, TEI and basis set)
    if tei_type == "sparse4c"
        integrals = Jint(Hcore=Hcore,unique_indices=[i .+ 1 for i in tei[1]],
        values=tei[2], length=length(tei[2]),bset=bset,bset_atom_mapping=bset_atom_mapping,
        bf_atom_shell_map=bf_atom_shell_map)
    else
        #println("tei_type other than sparse4c disabled for the moment")
        #exit()
        integrals = Jint_4rank(Hcore=Hcore,tensor=tei,bset=bset,bset_atom_mapping=bset_atom_mapping,
        bf_atom_shell_map=bf_atom_shell_map)
    end

    #############################
    # CREATING GRID
    #############################
    #grid always true for RKS/UKS
    # Also turn-on grid for RHF/UHF if density is requested
    if calc_density == true grid=true 
    if grid == true
        println("Creating grid")
        gridpoints,gridweights = numgrid_call("atomgrid",fragment,bset)
        println("Number of gridpoints: $(length(gridpoints))")
        println("Number of gridweights: $(length(gridweights))")
        #NOTE: Gridpoint coords are in Bohrs
        #Adding gridpoints to object
        integrals.gridpoints=gridpoints
        integrals.gridweights=gridweights
    end
    #end

    #############################
    # CHOOSING FOCK ALGORITHM
    #############################
    #Choosing Fock/Kohn-Sham algorithm (based on RHF vs. RKS vs. UKS vs. UHF, 2el-int-type etc.)
    Fock,fock_algorithm = choose_Fock(WFtype,DFTobj,fock_algorithm,tei_type,printlevel)

    #############################
    # GUESS
    #############################
    #Create initial guess for density matrix
    print_if_level("Providing guess: $guess",1,printlevel)
    if guess == "hcore"
        if WFtype=="RHF" || WFtype=="RKS"
            C = compute_core_guess(Hcore,S_minhalf)
            P = makeP(C, dim, numoccorbs) #Calculate new P from C
        else
            #UHF
            C = compute_core_guess(Hcore,S_minhalf)
            P_⍺ = makeP(C, dim, numoccorbs_⍺) #Calculate new P from C
            P_β = makeP(C, dim, numoccorbs_β) #Calculate new P from C
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
    if WFtype=="RHF" || WFtype=="RKS"
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

    #SCF loop beginning
    if printlevel == 1 @printf("%4s%15s%16s%14s%14s%8s%8s%8s%10s\n", "Iter", "Energy", "deltaE", "RMS-DP", "Max-DP", "Lshift", "Damp", "DIIS", "Max[F,P]") end
    
    time_scf=@elapsed for iter in 1:maxiter
        if printlevel > 1 print_iteration_header(iter) end
        if WFtype == "RHF" || WFtype=="RKS"
            #Optional damping of P before making Fock
            P = damping_control(P,P_old,damping,damping_val,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=damping_thresh)
            F = Fock(P,dim,integrals) #Make Fock-matrix
            FP_comm = FP_commutator(F,P,S,S_minhalf) #Calculating [F,P] commutator (for DIIS)
            #Possible levelshifting of Fock before diagonalization
            F = levelshift_control(F,levelshift,levelshift_val,numoccorbs,dim,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=lshift_thresh)
            F′ = transpose(S_minhalf)*F*S_minhalf #Transform Fock matrix
            energy = E_ZZ + 0.5 * tr((Hcore+F)*P) #Calculate energy after Fock formation
            # Possible DIIS extrapolation of F′ matrix before diagonalization
            F′ = diis_control(diisobj,F′,energy,FP_comm,iter,printlevel)
            eps, C′ = eigen(F′) #Diagonalize transformed Fock to get eps and C'
            C = S_minhalf*C′ # Get C from C'
            P_old=deepcopy(P) #Keep copy of old P
            P = makeP(C, dim, numoccorbs) #Calculate new P from C
            if debugprint == true write_matrices(F,C,P) end
        else
            #Possible damping of P matrices before making Fock
            P_⍺ = damping_control(P_⍺,P_⍺_old,damping,damping_val,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=damping_thresh)
            P_β = damping_control(P_β,P_β_old,damping,damping_val,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=damping_thresh)
            
            #Solve ⍺ part
            F_⍺ = Fock(P_⍺,P_β,dim,integrals) #Update Fock-matrix alpha
            FP_comm_⍺ = FP_commutator(F_⍺,P_⍺,S,S_minhalf) #Calculating [F,P] commutator (for DIIS)
            F_⍺ = levelshift_control(F_⍺,levelshift,levelshift_val,numoccorbs_⍺,dim,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=lshift_thresh)
            F′_⍺ = transpose(S_minhalf)*F_⍺*S_minhalf #Transform Fock matrix
            #Solve β part
            F_β = Fock(P_β,P_⍺,dim,integrals) #Update Fock-matrix beta
            FP_comm_β = FP_commutator(F_β,P_β,S,S_minhalf) #Calculating [F,P] commutator (for DIIS)
            F_β = levelshift_control(F_β,levelshift,levelshift_val,numoccorbs_β,dim,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=lshift_thresh)
            F′_β = transpose(S_minhalf)*F_β*S_minhalf #Transform Fock matrix

            energy = E_ZZ + 0.5 * tr((Hcore+F_⍺)*P_⍺) + 0.5 * tr((Hcore+F_β)*P_β)  #Calculate energy

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
            
            #if debugprint == true write_matrices(F,C,P) end

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
                if WFtype == "RHF"
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
            end
            break
        end

        if iter == maxiter
            println("Failed to converge in $maxiter iterations!")
            #Gathering results into Dict and returning
            Resultsdict=Dict("energy"=>"None","finaliter"=>"Failed ($maxiter iters)")
            return Resultsdict
        end
    end
    print_if_level("Time calculating SCF: $time_scf",1,printlevel)
    #################################################
    # ORBITALS, POPULATION ANALYSIS AND PROPERTIES
    #################################################
    #SCF loop done. Calculate and print properties
    #ORBITALS AND POPULATION ANALYSIS
    if printlevel > 0
        if WFtype=="RHF" || WFtype=="RKS"
            if grid == true 
                println("Now creating density")
                @time ⍴ = create_density(P,integrals)
                #println("⍴: $⍴")
                @time N = integrate_density(⍴,integrals.gridweights)
                println("Integrated no. of electrons: $N")
            end
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
    end

    #####################################
    # PRINT FINAL RESULTS
    #####################################
    if printlevel > 0
        print_final_results(energy,fragment,num_el,basisset,WFtype,fock_algorithm,finaliter)
    end
    #Gathering results into Dict and returning.
    #TODO: Add more here. 
    #Mulliken charges/spinpops, Mayer MBOs, dipole, all energy contributions,
    #orbital energies, occupation numbers, WFtype, numelectrons
    Resultsdict=Dict("energy"=>energy,"finaliter"=>finaliter)

    #TIMINGS: TODO
end
Resultsdict["time"]=totaltime
    return Resultsdict
end



####################
# Some structs
###################

"""
Jint: Jotunn integral struct in sparse form
NOTE: Separate object for 4c and DF ??
"""
Base.@kwdef mutable struct Jint
    #1-electron integrals
    Hcore::Matrix{Float64}
    #TWO-electron integrals
    #Original unique list of integral indices from GaussianBasis but with 1-indexing
    unique_indices::Vector{NTuple{4, Int16}}=[]
    #Integral values
    values::Vector{Float64}=[]
    length::Int64=0
    #Basis set object and data
    bset::BasisSet
    bset_atom_mapping::Vector{Int64}
    bf_atom_shell_map::Vector{Any}
    #Gridpoints (only used by DFT or density generation)
    gridpoints::Vector{Any}=[]
    gridweights::Vector{Any}=[]
end

"""
Jint: Jotunn integral struct in 4-rank form. Rarely used.
"""
Base.@kwdef mutable struct Jint_4rank
    #1-electron integrals
    Hcore::Matrix{Float64}
    #TWO-electron integrals
    tensor::Array{Float64, 4}
    #Basis set object and data
    bset::BasisSet
    bset_atom_mapping::Vector{Int64}
    bf_atom_shell_map::Vector{Any}
    #Gridpoints (only used by DFT or density generation)
    gridpoints::Vector{Any}=[]
    gridweights::Vector{Any}=[]
end

"""
JDFT: Jotunn DFT functional object
"""
Base.@kwdef mutable struct JDFT
    name::String
    hybrid::Bool
end
