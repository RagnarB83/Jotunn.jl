"""
deltaPcheck: Compare density matrices by RMSE and MaxE
"""
function deltaPcheck(P, Pold)
    diff = Pold - P
    mae=abs(mean(diff))
    rms_e=sqrt(sum([i^2 for i in diff])/length(diff))
    max_e=maximum(diff)
    return rms_e, max_e
end


"""
diis_control: Do DIIS extrapolation or not
"""
function diis_control(F′,F,diis_error_matrices,Fockmatrices,diis,diis_size,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=diis_thresh)

    if diis == true
        println("DIIS true")
        #Calculate new DIIS error vector based on current F and P
        diis_err = diis_error_vector(F,P,S,S_minhalf)

        #Add error vector to array and F to 
        push!(diis_error_matrices,diis_err)
        push!(Fockmatrices,F)

        if length(diis_error_matrices) > 1
            println("length(diis_error_matrices) :", length(diis_error_matrices))

            #B = zeros()
            #for i in length(diis_error_matrices)
            #    B
            #end

            #Determine ci coefficeints

            #SVD
            U,Sval, Vt = svd(a)

            #Do extrapolation based on F′ and ci coefficients
            newF=zeros(size(F′))
            for i in 1:size(ci_s)
                newF += ci_s[i]*Fockmatrices[i]
            end
            newF′=F′
            return newF′
        else
            #Note enough matrices to do extrapolation.
            #Returning unchanged matrix
            println("Not enough matrices ($(length(diis_error_matrices))) to do DIIS extrapolation")
            return F′
        end
    else
        #No DIIS: unchanged F′
        #println("No DIIS")
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
    if printlevel > 1 println("Damping: $damping_flag") end
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
function diis_error_vector(F,P,S,Sminhalf)
    commut= FP_commutator(F,P,S)
    error_vec = transpose(Sminhalf)*commut*Sminhalf
    return error_vec
end