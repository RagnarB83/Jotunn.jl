"""
JDIIS: Jotunn DIIS struct
"""
Base.@kwdef mutable struct JDIIS
    #Settings
    active::Bool #whether DIIS method is used in this job
    diis_size::Int64=5
    diis_startiter::Int64=2
    DIISBfac::Float64=1.05
    #Data
    diis_flag::Bool=false #whether DIIS is active in this iteration
    errorvectors::Vector{Matrix{Float64}}=[]
    maxerrors::Vector{Float64}=[]
    rmserrors::Vector{Float64}=[]
    max_diis_error::Float64=999999
    Fockmatrices::Vector{Matrix{Float64}}=[]
    energies::Vector{Float64}=[]
end

"""
compute_core_guess: Compute MO coefficients from F=Hcore
"""
function compute_core_guess(Hcore,S_minhalf)
    F=Hcore #initial Fock matrix
    F′ = transpose(S_minhalf)*F*S_minhalf #Transform Fock matrix
    eps, C′ = eigen(F′) #Diagonalize transformed Fock to get eps and C'
    C = S_minhalf*C′ # Get C from C'
    return C
end

"""
deltaPcheck: Compare density matrices by RMSE and MaxE
"""
function deltaPcheck(P, Pold)
    diff = [abs(i) for i in (P-Pold)]
    #rms_e=sqrt(sum([i^2 for i in diff])/length(diff))
    #rms_e = sqrt(sum(x -> x*x, diff) / (length(diff)))
    #Only unique values used in RMS calc.
    rms_e = sqrt(sum([abs(i)^2 for i in LowerTriangular(diff)])/(size(diff)[1]*(size(diff)[1]+1)))
    max_e=maximum([abs(i) for i in diff])
    return rms_e, max_e
end





"""
diis_control: Do DIIS extrapolation of F′

TODO: We should start collecting error vectors immediately if DIIS method is true.
Only when DIIS-startiter or other threshold met do we start extrapolation though
"""
function diis_control(diisobj,F′,energy,FP_comm,iter,printlevel)
    if diisobj.active == true
        print_if_level("DIIS method active. Iteration: $iter  (DIIS startiter: $(diisobj.diis_startiter))",2,printlevel)
        
        #DIIS error vector
        diis_err = FP_comm
        max_diis_err=maximum(diis_err)
        rms_diis_err= sqrt(sum(x -> x*x, diis_err) / length(diis_err))
        #Setting Max DIIS error in object for convergence check
        diisobj.max_diis_error=max_diis_err
        print_if_level("Max DIIS error: $max_diis_err",2,printlevel)
        print_if_level("RMS DIIS error: $rms_diis_err",2,printlevel)
        #Add error vector to array and F′ to Fockmatrices
        print_if_level("DIIS is active. Storing Fock and error vector for iteration.",2,printlevel)

        #Throwing out option
        #TODO: Throw out vector from high-energy iteration instead?
        throwout_option="maxerror" #Options: maxerror or oldest vector

        #Throwing out DIIS vector if we are at diis_size
        #NOTE: Throw out high-energy vector instead? or largest error vector
        if length(diisobj.errorvectors) == diisobj.diis_size
            #println("Throwing out before adding")
            if throwout_option == "maxerror"
                #println("diisobj.maxerrors: ", diisobj.maxerrors)
                #println("diisobj.rmserrors: ", diisobj.rmserrors)
                #Not sure about abs
                delindex = findmax([abs(i) for i in diisobj.maxerrors])[2]
                #println("Deleting index:", delindex)
                deleteat!(diisobj.errorvectors,delindex)
                deleteat!(diisobj.Fockmatrices,delindex)
                deleteat!(diisobj.energies,delindex)
                deleteat!(diisobj.maxerrors,delindex)
                deleteat!(diisobj.rmserrors,delindex)
                #println("After DIIS deletion")
                #println("diisobj.maxerrors: ", diisobj.maxerrors)
                #println("diisobj.rmserrors: ", diisobj.rmserrors)
            elseif throwout_option == "oldest"
                #Deleting first vector (i.e. oldest)
                diisobj.errorvectors = diisobj.errorvectors[2:end]
                diisobj.Fockmatrices = diisobj.Fockmatrices[2:end]
                diisobj.energies = diisobj.energies[2:end]
                diisobj.maxerrors = diisobj.maxerrors[2:end]
                diisobj.rmserrors = diisobj.rmserrors[2:end]
            end
            
        end

        #Add current error vector,Fock matrix and energy to lists
        if iter > 1 #Skipping bogus data from first iteration (at least bad when guess=Hcore)
            push!(diisobj.errorvectors,diis_err)
            push!(diisobj.Fockmatrices,F′)
            push!(diisobj.energies,energy)
            push!(diisobj.maxerrors,max_diis_err)
            push!(diisobj.rmserrors,rms_diis_err)
        end
        println("DIIS current size:", length(diisobj.errorvectors))
        #else
        #    print_if_level("DIIS has not yet reached diis_startiter=$(diisobj.diis_startiter). Not storing data.",2,printlevel)
        #end

        #DO DIIS extrapolation if active and as soon as there are enough vectors
        if iter >= diisobj.diis_startiter
            #println("DIIS startiter reached")
            if length(diisobj.errorvectors) > 1
                print_if_level("DIIS extrapolation will be performed",2,printlevel)
                diisobj.diis_flag=true #Now setting DIIS flag to true
                global damping_flag = false #Makes sure  damping is off when DIIS extrapolation is on
                #TODO: Revisit above

                #Bias DIIS to the Fock matrix with lowest energy by multiplying diagonal elements of
                #the other error vectors with 1.05 or so
                lowestE_index=findmin(diisobj.energies)[2]
                for diis_m_index in eachindex(diisobj.errorvectors)
                    if diis_m_index != lowestE_index
                        diis_m=diisobj.errorvectors[diis_m_index]
                        for i in 1:size(diis_m)[1]
                            diis_m[i,i] = diisobj.DIISBfac*diis_m[i,i]
                        end
                    end
                end

                #Defining B
                B=zeros(length(diisobj.errorvectors)+1,length(diisobj.errorvectors)+1)
                B[end, :].=-1
                B[:, end].=-1
                B[end, end]=0.0
                for i in eachindex(diisobj.errorvectors)
                    for j in eachindex(diisobj.errorvectors)
                        B[i,j] = dot(diisobj.errorvectors[i],diisobj.errorvectors[j])
                    end
                end
                #Z, residual vector
                Z = zeros(length(diisobj.errorvectors)+1)
                Z[end] = -1 #Last value is -1

                #Solve linear equation to get ci coefficeints
                try
                    coeffs = B\Z
                catch problem
                    if isa(problem,LinearAlgebra.SingularException)
                        print_if_level("Singular matrix during DIIS solve. Skipping extrapolation in this step",2,printlevel)
                        return F′
                    end
                    print_if_level("Unknown DIIS problem. Skipping extrapolation in this step.",2,printlevel)
                    return F′
                end
                
                coeffs = coeffs[1:end-1] #Removing last value (lambda)
                print_if_level("DIIS coefficients: $coeffs",2,printlevel)

                #Extrapolate F′ from old F′ and ci coefficients
                newF′=zeros(size(F′))
                for i in 1:length(coeffs)
                    newF′ += coeffs[i]*diisobj.Fockmatrices[i]
                end
                return newF′
            else
                #Not enough matrices to do extrapolation. Returning  F′
                print_if_level("Not enough matrices ($(length(diisobj.errorvectors))) to do DIIS extrapolation",2,printlevel)
                return F′
            end
        else
            println("DIIS active but startiter not yet reached. No extrapolation yet (only storing).")
            return F′
        end
    else
        #No DIIS used in job: unchanged F′
        return F′
    end
end


"""
levelshift_control: Control and do levelshifting of Fock matrix during SCF
"""
function levelshift_control(F,levelshift,levelshift_val,numorbs,dim,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=0.001)
    global levelshift_flag
    print_if_level("levelshift option: $levelshift",2,printlevel)
    if levelshift == true
        print_if_level("Levelshift option is on",2,printlevel)
        if printlevel > 1 println("Levelshift flag: $levelshift_flag") end
        if levelshift_flag == true
            #println("Iter: $iter P_RMS: $P_RMS (threshold to turn off levelshift:", turnoff_threshold)
            if P_RMS > turnoff_threshold
                if printlevel > 1 println("Level shift active. Adding levelshift: $levelshift_val in iteration $iter") end
                F = levelshift_Fock(F,levelshift_val,numorbs,dim)
            else
                if printlevel > 1
                    println("P_RMS=$P_RMS less than threshold!")
                    println("Levelshifting is off!")
                end
                #println("Removing levelshift in iteration $iter")
                #F = levelshift_Fock(F,-levelshift,numorbs,dim)
                levelshift_flag=false
            end
        end
    else
        print_if_level("Levelshift option is off",2,printlevel)
    end
    return F
end

"""levelshift_Fock: levleshift Fock matrix
"""
function levelshift_Fock(F,parameter, occ, dim)
    #Shift the diagonal elements of the virtual block of the Fock matrix
    for i in occ+1:dim
        F[i,i] = F[i,i] + parameter
    end
    return F
end


"""
damping_control: Control and do damping during SCF
"""
function damping_control(P,Pold,damping,dampingpar,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=0.001)
    global damping_flag
    #if printlevel > 1 println("Damping: $damping_flag") end
    if isnothing(dampingpar) == false
        if damping_flag == true
            #println("Iter: $iter P_RMS: $P_RMS (threshold to turn off damping:", turnoff_threshold)
            if P_RMS > turnoff_threshold
                if printlevel > 1 println("Damping active. Adding damping $dampingpar in iteration $iter") end
                #println("P before damping: $P")
                P = damping_P(P,Pold;mixpar=dampingpar)
                #println("P after damping: $P")
            else
                if printlevel > 1
                    println("P_RMS=$P_RMS less than threshold!")
                    println("Damping is off!")
                end
                #println("Removing damping in iteration $iter")
                damping_flag=false
            end
        end
    end
    return P
end



"""
damping: mixing of density
Currently static
"""
function damping_P(Pnew,Pold;mixpar=0.4)
    Pnew_damped=(1-mixpar)Pnew+mixpar*Pold
end

"""
[F,P] commutator
"""
function FP_commutator(F,P,S,Sminhalf)
    # regular
    #commut = F*P*S-S*P*F
    #orthonormalized basis
    commut= transpose(Sminhalf)*(F*P*S-S*P*F)*Sminhalf
    return commut
end


function check_for_convergence(deltaE,energythreshold,FP_comm,diis_error_conv_threshold,iter,printlevel)

    #Current behaviour: if either deltaE or MaxDIISerror condition is fulfilled
    #with reasonable thresholds on other, we signal convergence
    
    max_diis_error=maximum(FP_comm)
    if abs(deltaE) < energythreshold && max_diis_error < diis_error_conv_threshold*10
        print_if_level("Energy convergence threshold reached: $(abs(deltaE)) < $energythreshold",1,printlevel)
        #println("diisobj.max_diis_error: $diisobj.max_diis_error and diis_error_conv_threshold: $diis_error_conv_threshold ")
        if printlevel >= 1
            print(Crayon(foreground = :green, bold = true), 
                "\n                              SCF converged in $iter iterations! Hell yeah! 🎉\n\n",
                Crayon(reset=true))
        end
            return true
    #If DIIS error is converged and energy threshold off by 10
    elseif max_diis_error < diis_error_conv_threshold && abs(deltaE) < energythreshold*10
        print_if_level("DIIS convergence threshold reached: $(max_diis_error) < $diis_error_conv_threshold",1,printlevel)
        if printlevel >= 1
            print(Crayon(foreground = :green, bold = true), 
                "\n                              SCF converged in $iter iterations! Hell yeah! 🎉\n\n",
                Crayon(reset=true))
        end
        return true
    end
    return false
end