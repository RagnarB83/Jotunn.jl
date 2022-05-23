"""
deltaPcheck: Compare density matrices by RMSE and MaxE
"""
function deltaPcheck(P, Pold)
    diff = Pold - P
    mae=abs(mean(diff))
    #rms_e=sqrt(sum([i^2 for i in diff])/length(diff))
    rms_e = sqrt(sum(x -> x*x, diff) / length(diff))
    max_e=maximum(diff)
    return rms_e, max_e
end


"""
diis_control: Do DIIS extrapolation or not
"""
function diis_control(F′,F,P,S,S_minhalf,energy,diis_error_matrices,Fockmatrices,energies,diis,diis_size,DIISBfac,
        P_RMS,rmsDP_threshold,iter,printlevel, FP_comm,diis_startiter)
        global diis_flag
    if diis == true
        println("DIIS method active. Iteration: $iter  (DIIS startiter: $diis_startiter)")
        
        #Calculate current-iteration DIIS error vector based on [F,P]
        #diis_err=FP_comm #Direct
        diis_err = diis_error_vector(FP_comm,S_minhalf) # orthonormal basis
        #diis_err = 2*diis_error_vector(FP_comm,S_minhalf) #orthonormal basis scaled by X
        
        println("Max DIIS error: $(maximum(diis_err))")
        rms_diis_err= sqrt(sum(x -> x*x, diis_err) / length(diis_err))
        println("RMS DIIS error: $rms_diis_err")

        #Add error vector to array and F′ to Fockmatrices
        if iter >= diis_startiter
            println("DIIS is active")
            #Throwing out old DIIS vector if we are at diis_size
            #NOTE: Throw out high-energy vector instead?
            if length(diis_error_matrices) == diis_size
                diis_error_matrices = diis_error_matrices[2:end]
                Fockmatrices = Fockmatrices[2:end]
                energies = energies[2:end]
            end
            #Add current error vector,Fock matrix and energy to lists
            push!(diis_error_matrices,diis_err)
            push!(Fockmatrices,F′)
            push!(energies,energy)
        else
            println("DIIS has not yet reached diis_startiter=$diis_startiter")
        end

        #DO DIIS extrapolation if enough data
        if length(diis_error_matrices) > 1
            println("DIIS error matrix size: $(length(diis_error_matrices))")
            diis_flag=true #Now setting DIIS flag to true

            #Bias DIIS to Fock matrix with lowest energy
            #relenergies=[i-minimum(energies) for i in energies]
            lowestE_index=findmin(energies)[2]
            for diis_m_index in eachindex(diis_error_matrices)
                if diis_m_index != lowestE_index
                    diis_m=diis_error_matrices[diis_m_index]
                    for i in 1:size(diis_m)[1]
                        diis_m[i,i] = DIISBfac*diis_m[i,i]
                    end
                end
            end

            #Defining B
            B=zeros(length(diis_error_matrices)+1,length(diis_error_matrices)+1)
            B[end, :].=-1
            B[:, end].=-1
            B[end, end]=0.0
            for i in eachindex(diis_error_matrices)
                for j in eachindex(diis_error_matrices)
                    B[i,j] = dot(diis_error_matrices[i],diis_error_matrices[j])
                end
            end
            #Z, residual vector
            Z = zeros(length(diis_error_matrices)+1)
            Z[end] = -1 #Last value is -1

            #Solve linear equation to get ci coefficeints
            coeffs = B\Z
            coeffs = coeffs[1:end-1] #Removing last value (lambda)
            println("DIIS coefficients: ", coeffs)

            #Extrapolate F′ from old F′ and ci coefficients
            newF′=zeros(size(F′))
            for i in 1:length(coeffs)
                newF′ += coeffs[i]*Fockmatrices[i]
            end
            return newF′
        else
            #Not enough matrices to do extrapolation. Returning  F′
            println("Not enough matrices ($(length(diis_error_matrices))) to do DIIS extrapolation")
            return F′
        end
    else
        #No DIIS: unchanged F′
        return F′
    end
end


"""
levelshift_control: Control and do levelshifting of Fock matrix during SCF
"""
function levelshift_control(F,levelshift,levelshift_val,numorbs,dim,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=0.001)
    global levelshift_flag
    printdebug("levelshift option: $levelshift")
    if levelshift == true
        printdebug("Levelshift option is on")
        if printlevel > 1 println("Levelshift flag: $levelshift_flag") end
        if levelshift_flag == true
            #println("Iter: $iter P_RMS: $P_RMS (threshold to turn off levelshift:", turnoff_threshold)
            if P_RMS > turnoff_threshold
                if printlevel > 1 println("Level shift active. Adding levelshift: $levelshift in iteration $iter") end
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
        printdebug("Levelshift option is off")
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
function FP_commutator(F,P,S)
    commut = F*P*S-S*P*F
    return commut
end

"DIIS error vector"
function diis_error_vector(FP_comm,Sminhalf)
    #commut= FP_commutator(F,P,S)
    error_vec = transpose(Sminhalf)*FP_comm*Sminhalf
    return error_vec
end