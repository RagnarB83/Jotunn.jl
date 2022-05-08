"""
deltaPcheck: Compare density matrices
"""
function deltaPcheck(P, Pold)
    diff = Pold - P
    mae=abs(mean(diff))
    rms_e=sqrt(sum([i^2 for i in diff])/length(diff))
    max_e=maximum(diff)
    return rms_e, max_e
end

"""
levelshift_control
"""
function levelshift_control(F,levelshift,numorbs,dim,P_RMS,rmsDP_threshold,iter,printlevel; turnoff_threshold=0.001)
    global levelshift_flag
    if printlevel > 1 println("Levelshift: $levelshift_flag") end
    if isnothing(levelshift) == false
        if levelshift_flag == true
            #println("Iter: $iter P_RMS: $P_RMS (threshold to turn off levelshift:", turnoff_threshold)
            if P_RMS > turnoff_threshold
                if printlevel > 1 println("Level shift active. Adding levelshift: $levelshift in iteration $iter") end
                F = levelshift_Fock(F,levelshift,numorbs,dim)
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
    end
    return F
end

"""levelshift_Fock:
"""
function levelshift_Fock(F,parameter, occ, dim)
    #Shift the diagonal elements of the virtual block of the Fock matrix
    for i in occ+1:dim
        F[i,i] = F[i,i] + parameter
    end
    return F
end