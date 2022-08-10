###################################################################
#Fock algorithms based on pre-computed 2-electron integrals as 
#either full rank-4 tensor or condensed sparse version
###################################################################



# Fock functions that use full sparse form of 2-electron integrals
# via Jint struct
###################################################################

"""
Fock_loop_sparse_RHF: Multiple matrix-element update per unique set of indices
For RHF
"""
function Fock(P,dim,integrals::Jint_HF)
    #println("Inside RFock for Jint_HF")
    #- Next combine assignments into simpler ones
    G = zeros(dim,dim)
    #Looping over unique sets of indices and values from sparse2e4c
    
    @inbounds @fastmath for i in eachindex(integrals.unique_indices)
        #This is the unique set of indices provided by GaussianBasis/libcint (sparse)
        µ,ν,λ,σ=integrals.unique_indices[i]
        value=integrals.values[i] #integral value
        #Logic to handle Fock matrix contributions for all permutations
        if µ > ν
            µν = µ*(µ+1)/2 + ν
        else
            µν = ν*(ν+1)/2 + µ
        end
        if λ > σ
            λσ = λ*(λ+1)/2 + σ
        else
            λσ = σ*(σ+1)/2 + λ
        end
        γµν = µ !== ν
        γλσ = λ !== σ
        γxy = µν !== λσ

        #Now adding value to correct matrix element
        if γµν && γλσ && γxy
            #println("Case1: e.g. µ!=ν!=c!=d  and 1232 ")
            G[µ,ν] += 2*P[λ,σ]*value # J µνλσ
            G[ν,µ] += 2*P[λ,σ]*value # J νµλσ
            #G[µ,ν] += P[σ,λ]*value # J µνσλ
            G[λ,σ] += 2*P[µ,ν]*value # J λσµν
            #G[λ,σ] += P[ν,µ]*value # J λσνµ
            #G[ν,µ] += P[σ,λ]*value # J νµσλ
            G[σ,λ] += 2*P[ν,µ]*value # J σλνµ
            #G[σ,λ] += P[µ,ν]*value # J σλµν

            G[µ,λ] -= 0.5*P[ν,σ]*value # K   µνλσ
            G[ν,λ] -= 0.5*P[µ,σ]*value # K   νµλσ
            G[µ,σ] -= 0.5*P[ν,λ]*value # K   µνσλ
            G[λ,µ] -= 0.5*P[σ,ν]*value # K   λσµν
            G[λ,ν] -= 0.5*P[σ,µ]*value # K   λσνµ
            G[ν,σ] -= 0.5*P[µ,λ]*value # K   νµσλ
            G[σ,ν] -= 0.5*P[λ,µ]*value # K   σλνµ
            G[σ,µ] -= 0.5*P[λ,ν]*value # K   σλµν

        elseif γλσ && γxy
            #println("Case2")
            G[µ,ν] += 2*P[λ,σ]*value # J µνλσ
            #G[µ,ν] += P[σ,λ]*value # J µνσλ
            #Strange diis+damp+levelshift fail here
            G[λ,σ] += P[µ,ν]*value # J λσµν
            G[σ,λ] += P[ν,µ]*value # J σλνµ

            G[µ,λ] -= 0.5*P[ν,σ]*value # K   µνλσ
            G[µ,σ] -= 0.5*P[ν,λ]*value # K   µνσλ
            G[λ,µ] -= 0.5*P[σ,ν]*value # K   λσµν
            G[σ,ν] -= 0.5*P[λ,µ]*value # K   σλνµ
        elseif γµν && γxy
            #println("Case3 (incomplete)")
            G[µ,ν] += P[λ,σ]*value # J µνλσ
            G[ν,µ] += P[λ,σ]*value # J νµλσ
            G[λ,σ] += 2*P[µ,ν]*value # J λσµν
            #G[λ,σ] += P[ν,µ]*value # J λσνµ

            G[µ,λ] -= 0.5*P[ν,σ]*value # K   µνλσ
            G[ν,λ] -= 0.5*P[µ,σ]*value # K   νµλσ
            G[λ,µ] -= 0.5*P[σ,ν]*value # K   λσµν
            G[λ,ν] -= 0.5*P[σ,µ]*value # K   λσνµ

        elseif γµν && γλσ
            #println("Case4")
            G[µ,ν] += 2*P[λ,σ]*value # J µνλσ
            G[ν,µ] += 2*P[λ,σ]*value # J νµλσ
            #G[µ,ν] += P[σ,λ]*value # J µνσλ
            #G[ν,µ] += P[σ,λ]*value # J νµσλ

            G[µ,λ] -= 0.5*P[ν,σ]*value # K   µνλσ
            G[ν,λ] -= 0.5*P[µ,σ]*value # K   νµλσ
            G[µ,σ] -= 0.5*P[ν,λ]*value # K   µνσλ
            G[ν,σ] -= 0.5*P[µ,λ]*value # K   νµσλ
        elseif γxy #µ=ν; λ=σ µν != λσ
            #println("Case5: example: 1122 and 2211")
            G[µ,ν] += P[λ,σ]*value # J µνλσ
            G[λ,σ] += P[µ,ν]*value # J λσµν

            G[µ,λ] -= 0.5*P[ν,σ]*value # K   µνλσ
            G[λ,µ] -= 0.5*P[σ,ν]*value # K   λσµν
        else
            #println("Case6, Example: 1111. ")
            G[µ,ν] += P[λ,σ]*value # J µνλσ
            G[µ,λ] -= 0.5*P[ν,σ]*value # K µνλσ
        end
    end


    F = integrals.Hcore + G
    return F
end

"""
Fock_loop_sparse_UHF: Multiple matrix-element update per unique set of indices
"""
function Fock(Pi,Pj,dim,integrals::Jint_HF,manifold)
    #println("Inside UFock for Jint_HF")
    G = zeros(dim,dim)
    #Looping over unique sets of indices and values from sparse2e4c
    @inbounds @fastmath for i in eachindex(integrals.unique_indices)
        #This is the unique set of indices provided by GaussianBasis/libcint (sparse)
        µ,ν,λ,σ=integrals.unique_indices[i]
        value=integrals.values[i] #integral value
        #Logic to handle Fock matrix contributions for all permutations
        if µ > ν
            µν = µ*(µ+1)/2 + ν
        else
            µν = ν*(ν+1)/2 + µ
        end
        if λ > σ
            λσ = λ*(λ+1)/2 + σ
        else
            λσ = σ*(σ+1)/2 + λ
        end
        γµν = µ !== ν
        γλσ = λ !== σ
        γxy = µν !== λσ

        #Now adding value to correct matrix element
        if γµν && γλσ && γxy
            #println("Case1: e.g. µ!=ν!=c!=d  and 1232 ")
            G[µ,ν] += Pi[λ,σ]*value+Pj[λ,σ]*value # J µνλσ
            G[ν,µ] += Pi[λ,σ]*value+Pj[λ,σ]*value # J νµλσ
            G[µ,ν] += Pi[σ,λ]*value+Pj[σ,λ]*value # J µνσλ
            G[λ,σ] += Pi[µ,ν]*value+Pj[µ,ν]*value # J λσµν
            G[λ,σ] += Pi[ν,µ]*value+Pj[ν,µ]*value # J λσνµ
            G[ν,µ] += Pi[σ,λ]*value+Pj[σ,λ]*value # J νµσλ
            G[σ,λ] += Pi[ν,µ]*value+Pj[ν,µ]*value # J σλνµ
            G[σ,λ] += Pi[µ,ν]*value+Pj[µ,ν]*value # J σλµν

            G[µ,λ] -= Pi[ν,σ]*value # K   µνλσ
            G[ν,λ] -= Pi[µ,σ]*value # K   νµλσ
            G[µ,σ] -= Pi[ν,λ]*value # K   µνσλ
            G[λ,µ] -= Pi[σ,ν]*value # K   λσµν
            G[λ,ν] -= Pi[σ,µ]*value # K   λσνµ
            G[ν,σ] -= Pi[µ,λ]*value # K   νµσλ
            G[σ,ν] -= Pi[λ,µ]*value # K   σλνµ
            G[σ,µ] -= Pi[λ,ν]*value # K   σλµν

        elseif γλσ && γxy
            #println("Case2")
            G[µ,ν] += Pi[λ,σ]*value+Pj[λ,σ]*value # J µνλσ
            G[µ,ν] += Pi[σ,λ]*value+Pj[σ,λ]*value # J µνσλ
            G[λ,σ] += Pi[µ,ν]*value+Pj[µ,ν]*value # J λσµν
            G[σ,λ] += Pi[ν,µ]*value+Pj[ν,µ]*value # J σλνµ

            G[µ,λ] -= Pi[ν,σ]*value # K   µνλσ
            G[µ,σ] -= Pi[ν,λ]*value # K   µνσλ
            G[λ,µ] -= Pi[σ,ν]*value # K   λσµν
            G[σ,ν] -= Pi[λ,µ]*value # K   σλνµ
        elseif γµν && γxy
            #println("Case3 (incomplete)")
            G[µ,ν] += Pi[λ,σ]*value+Pj[λ,σ]*value # J µνλσ
            G[ν,µ] += Pi[λ,σ]*value+Pj[λ,σ]*value # J νµλσ
            G[λ,σ] += Pi[µ,ν]*value+Pj[µ,ν]*value # J λσµν
            G[λ,σ] += Pi[ν,µ]*value+Pj[ν,µ]*value # J λσνµ

            G[µ,λ] -= Pi[ν,σ]*value # K   µνλσ
            G[ν,λ] -= Pi[µ,σ]*value # K   νµλσ
            G[λ,µ] -= Pi[σ,ν]*value # K   λσµν
            G[λ,ν] -= Pi[σ,µ]*value # K   λσνµ

        elseif γµν && γλσ
            #println("Case4")
            G[µ,ν] += Pi[λ,σ]*value+Pj[λ,σ]*value # J µνλσ
            G[ν,µ] += Pi[λ,σ]*value+Pj[λ,σ]*value # J νµλσ
            G[µ,ν] += Pi[σ,λ]*value+Pj[σ,λ]*value # J µνσλ
            G[ν,µ] += Pi[σ,λ]*value+Pj[σ,λ]*value # J νµσλ

            G[µ,λ] -= Pi[ν,σ]*value # K   µνλσ
            G[ν,λ] -= Pi[µ,σ]*value # K   νµλσ
            G[µ,σ] -= Pi[ν,λ]*value # K   µνσλ
            G[ν,σ] -= Pi[µ,λ]*value # K   νµσλ
        elseif γxy #µ=ν; λ=σ µν != λσ
            #println("Case5: example: 1122 and 2211")
            G[µ,ν] += Pi[λ,σ]*value+Pj[λ,σ]*value # J µνλσ
            G[λ,σ] += Pi[µ,ν]*value+Pj[µ,ν]*value # J λσµν

            G[µ,λ] -= Pi[ν,σ]*value # K   µνλσ
            G[λ,µ] -= Pi[σ,ν]*value # K   λσµν
        else
            #println("Case6, Example: 1111. ")
            G[µ,ν] += Pi[λ,σ]*value+Pj[λ,σ]*value # J µνλσ
            G[µ,λ] -= Pi[ν,σ]*value # K µνλσ
        end
    end
    F = integrals.Hcore + G
    return F
end


###################################################################
# Fock functions that use full 4-rank tensor
#TODO: Enable macro that can switch from @simd to @turbo
###################################################################
"""
Fock_loop_RHF: Fock-matrix RHF case loop-version
Benchmark: RHF/def2-QZVPP on H2O with MBA: 
243.604541  seconds total, 3.481820 seconds per Fock  (1 THREAD)
"""
function Fock(P,dim,integrals::Jint_4rank)
    #println("Inside RFock for Jint_4rank")
    JK = zeros(dim,dim)
    @simd for µ in 1:dim
        for ν in 1:dim
            for λ in 1:dim
                for σ in 1:dim
                    @inbounds JK[µ,ν] += P[λ,σ]*(integrals.tensor[µ,ν,λ,σ]-0.5*integrals.tensor[µ,λ,ν,σ])
                end
            end
        end
    end
    F = integrals.Hcore + JK
    return F
end


"""
Fock_loop_UHF: Fock-matrix UHF case using simple loop
"""
function Fock(Pi,Pj,dim,integrals::Jint_4rank,manifold)
    #println("Inside UFock for Jint_4rank")
    JK = zeros(dim,dim)
    @simd for µ in 1:dim
        for ν in 1:dim
            for λ in 1:dim
                for σ in 1:dim
                    @inbounds JK[µ,ν] += Pi[λ,σ]*integrals.tensor[ν,µ,λ,σ]+Pj[λ,σ]*integrals.tensor[ν,µ,λ,σ]-Pi[λ,σ]*integrals.tensor[ν,λ,µ,σ]
                end
            end
        end
    end
    F = integrals.Hcore + JK
    return F
end


###################################################################
# Fock "Kohn-Sham" functions
###################################################################

"""
Fock for DFT and RKS: J and XC
"""
function Fock(P,dim,integrals::Jint_KS)
    println("Inside RFock for Jint_KS")
    J = zeros(dim,dim)
    #Looping over unique sets of TEI indices and values from sparse2e4c
    @inbounds @fastmath for i in eachindex(integrals.unique_indices)
        #This is the unique set of TEI indices provided by GaussianBasis/libcint (sparse)
        µ,ν,λ,σ=integrals.unique_indices[i]
        value=integrals.values[i] #integral value
        #Logic to handle Fock matrix contributions for all permutations
        if µ > ν
            µν = µ*(µ+1)/2 + ν
        else
            µν = ν*(ν+1)/2 + µ
        end
        if λ > σ
            λσ = λ*(λ+1)/2 + σ
        else
            λσ = σ*(σ+1)/2 + λ
        end
        γµν = µ !== ν
        γλσ = λ !== σ
        γxy = µν !== λσ

        if γµν && γλσ && γxy
            J[µ,ν] += 2*P[λ,σ]*value # J µνλσ
            J[ν,µ] += 2*P[λ,σ]*value # J νµλσ
            J[λ,σ] += 2*P[µ,ν]*value # J λσµν
            J[σ,λ] += 2*P[ν,µ]*value # J σλνµ
        elseif γλσ && γxy
            J[µ,ν] += 2*P[λ,σ]*value # J µνλσ
            J[λ,σ] += P[µ,ν]*value # J λσµν
            J[σ,λ] += P[ν,µ]*value # J σλνµ
        elseif γµν && γxy
            J[µ,ν] += P[λ,σ]*value # J µνλσ
            J[ν,µ] += P[λ,σ]*value # J νµλσ
            J[λ,σ] += 2*P[µ,ν]*value # J λσµν
        elseif γµν && γλσ
            J[µ,ν] += 2*P[λ,σ]*value # J µνλσ
            J[ν,µ] += 2*P[λ,σ]*value # J νµλσ
        elseif γxy #µ=ν; λ=σ µν != λσ
            J[µ,ν] += P[λ,σ]*value # J µνλσ
            J[λ,σ] += P[µ,ν]*value # J λσµν
        else
            J[µ,ν] += P[λ,σ]*value # J µνλσ
        end
    end

    #Calculate Vxc contribution from P, bfs on gridpoints
    Fxc,Exc = Fxc_form(integrals,P)
    #println("Exc:", Exc)
    integrals.E_xc=Exc
    integrals.J=J
    #println("Inside Fock")
    #println("----------------------------")
    #println("Hcore:")
    #print_matrix(integrals.Hcore)
    #println("J:")
    #print_matrix(J)
    #println("integrals.Hcore + J:")
    #print_matrix(integrals.Hcore + J)
    #println("Fxc:")
    #print_matrix(Fxc)
    #println("----------------------------")
    F = integrals.Hcore + J + Fxc
    return F
end

"""
Fock for DFT and UKS (two sets of P): J, and XC

TODO: Truncate assignments down by symmetry
"""
function Fock(Pi,Pj,dim,integrals::Jint_KS,manifold)
    #println("Inside UFock for Jint_KS")
    J = zeros(dim,dim)
    #Looping over unique sets of TEI indices and values from sparse2e4c
    @inbounds @fastmath for i in eachindex(integrals.unique_indices)
        #This is the unique set of TEI indices provided by GaussianBasis/libcint (sparse)
        µ,ν,λ,σ=integrals.unique_indices[i]
        value=integrals.values[i] #integral value
        #Logic to handle Fock matrix contributions for all permutations
        if µ > ν
            µν = µ*(µ+1)/2 + ν
        else
            µν = ν*(ν+1)/2 + µ
        end
        if λ > σ
            λσ = λ*(λ+1)/2 + σ
        else
            λσ = σ*(σ+1)/2 + λ
        end
        γµν = µ !== ν
        γλσ = λ !== σ
        γxy = µν !== λσ

            #Now adding value to correct matrix element
            if γµν && γλσ && γxy
                #println("Case1: e.g. µ!=ν!=c!=d  and 1232 ")
                J[µ,ν] += Pi[λ,σ]*value+Pj[λ,σ]*value # J µνλσ
                J[ν,µ] += Pi[λ,σ]*value+Pj[λ,σ]*value # J νµλσ
                J[µ,ν] += Pi[σ,λ]*value+Pj[σ,λ]*value # J µνσλ
                J[λ,σ] += Pi[µ,ν]*value+Pj[µ,ν]*value # J λσµν
                J[λ,σ] += Pi[ν,µ]*value+Pj[ν,µ]*value # J λσνµ
                J[ν,µ] += Pi[σ,λ]*value+Pj[σ,λ]*value # J νµσλ
                J[σ,λ] += Pi[ν,µ]*value+Pj[ν,µ]*value # J σλνµ
                J[σ,λ] += Pi[µ,ν]*value+Pj[µ,ν]*value # J σλµν
    
            elseif γλσ && γxy
                #println("Case2")
                J[µ,ν] += Pi[λ,σ]*value+Pj[λ,σ]*value # J µνλσ
                J[µ,ν] += Pi[σ,λ]*value+Pj[σ,λ]*value # J µνσλ
                J[λ,σ] += Pi[µ,ν]*value+Pj[µ,ν]*value # J λσµν
                J[σ,λ] += Pi[ν,µ]*value+Pj[ν,µ]*value # J σλνµ
    
            elseif γµν && γxy
                #println("Case3 (incomplete)")
                J[µ,ν] += Pi[λ,σ]*value+Pj[λ,σ]*value # J µνλσ
                J[ν,µ] += Pi[λ,σ]*value+Pj[λ,σ]*value # J νµλσ
                J[λ,σ] += Pi[µ,ν]*value+Pj[µ,ν]*value # J λσµν
                J[λ,σ] += Pi[ν,µ]*value+Pj[ν,µ]*value # J λσνµ
    
    
            elseif γµν && γλσ
                #println("Case4")
                J[µ,ν] += Pi[λ,σ]*value+Pj[λ,σ]*value # J µνλσ
                J[ν,µ] += Pi[λ,σ]*value+Pj[λ,σ]*value # J νµλσ
                J[µ,ν] += Pi[σ,λ]*value+Pj[σ,λ]*value # J µνσλ
                J[ν,µ] += Pi[σ,λ]*value+Pj[σ,λ]*value # J νµσλ
    
            elseif γxy #µ=ν; λ=σ µν != λσ
                #println("Case5: example: 1122 and 2211")
                J[µ,ν] += Pi[λ,σ]*value+Pj[λ,σ]*value # J µνλσ
                J[λ,σ] += Pi[µ,ν]*value+Pj[µ,ν]*value # J λσµν
            else
                #println("Case6, Example: 1111. ")
                J[µ,ν] += Pi[λ,σ]*value+Pj[λ,σ]*value # J µνλσ
            end
        end

    #Calculate Vxc contribution from P, bfs on gridpoints
    Fxc,Exc = Fxc_form_UKS(integrals,Pi,Pj,manifold)
    #println("Exc: $Exc")
    integrals.E_xc=Exc
    integrals.J=J
    #println("Hcore+J:")
    #print_matrix(integrals.Hcore + J)
    F = integrals.Hcore + J + Fxc
    #println("J:")
    #print_matrix(J)
    #println("Fxc:")
    #print_matrix(Fxc)
    #println("F:")
    #print_matrix(F)
    #println("Pi:")
    #print_matrix(Pi)
    #println("Pj:")
    #print_matrix(Pj)
    return F

end


"""
Fock for HybridDFT and RKS: J,K and XC
"""
function Fock(P,dim,integrals::Jint_hybridKS)
    println("Inside RFock for Jint_hybridKS")
    println("Disabled")
    exit()
    G = zeros(dim,dim)
    #HF exchange scaling
    alpha=Jint_hybridKS.DFT.alpha
    rest=1-alpha
    #Looping over unique sets of TEI indices and values from sparse2e4c
    @inbounds @fastmath for i in eachindex(integrals.unique_indices)
        #This is the unique set of TEI indices provided by GaussianBasis/libcint (sparse)
        µ,ν,λ,σ=integrals.unique_indices[i]
        value=integrals.values[i] #integral value
        #Logic to handle Fock matrix contributions for all permutations
        if µ > ν
            µν = µ*(µ+1)/2 + ν
        else
            µν = ν*(ν+1)/2 + µ
        end
        if λ > σ
            λσ = λ*(λ+1)/2 + σ
        else
            λσ = σ*(σ+1)/2 + λ
        end
        γµν = µ !== ν
        γλσ = λ !== σ
        γxy = µν !== λσ

        if γµν && γλσ && γxy
            G[µ,ν] += 2*P[λ,σ]*value # J µνλσ
            G[ν,µ] += 2*P[λ,σ]*value # J νµλσ
            G[λ,σ] += 2*P[µ,ν]*value # J λσµν
            G[σ,λ] += 2*P[ν,µ]*value # J σλνµ

            G[µ,λ] -= 0.5*P[ν,σ]*value*alpha # K   µνλσ
            G[ν,λ] -= 0.5*P[µ,σ]*value*alpha # K   νµλσ
            G[µ,σ] -= 0.5*P[ν,λ]*value*alpha # K   µνσλ
            G[λ,µ] -= 0.5*P[σ,ν]*value*alpha # K   λσµν
            G[λ,ν] -= 0.5*P[σ,µ]*value*alpha # K   λσνµ
            G[ν,σ] -= 0.5*P[µ,λ]*value*alpha # K   νµσλ
            G[σ,ν] -= 0.5*P[λ,µ]*value*alpha # K   σλνµ
            G[σ,µ] -= 0.5*P[λ,ν]*value*alpha # K   σλµν

        elseif γλσ && γxy
            G[µ,ν] += 2*P[λ,σ]*value # J µνλσ
            G[λ,σ] += P[µ,ν]*value # J λσµν
            G[σ,λ] += P[ν,µ]*value # J σλνµ

            G[µ,λ] -= 0.5*P[ν,σ]*value*alpha # K   µνλσ
            G[µ,σ] -= 0.5*P[ν,λ]*value*alpha # K   µνσλ
            G[λ,µ] -= 0.5*P[σ,ν]*value*alpha # K   λσµν
            G[σ,ν] -= 0.5*P[λ,µ]*value*alpha # K   σλνµ
        elseif γµν && γxy
            G[µ,ν] += P[λ,σ]*value # J µνλσ
            G[ν,µ] += P[λ,σ]*value # J νµλσ
            G[λ,σ] += 2*P[µ,ν]*value # J λσµν

            G[µ,λ] -= 0.5*P[ν,σ]*value*alpha # K   µνλσ
            G[ν,λ] -= 0.5*P[µ,σ]*value*alpha # K   νµλσ
            G[λ,µ] -= 0.5*P[σ,ν]*value*alpha # K   λσµν
            G[λ,ν] -= 0.5*P[σ,µ]*value*alpha # K   λσνµ

        elseif γµν && γλσ
            G[µ,ν] += 2*P[λ,σ]*value # J µνλσ
            G[ν,µ] += 2*P[λ,σ]*value # J νµλσ

            G[µ,λ] -= 0.5*P[ν,σ]*value*alpha # K   µνλσ
            G[ν,λ] -= 0.5*P[µ,σ]*value*alpha # K   νµλσ
            G[µ,σ] -= 0.5*P[ν,λ]*value*alpha # K   µνσλ
            G[ν,σ] -= 0.5*P[µ,λ]*value*alpha # K   νµσλ
        elseif γxy #µ=ν; λ=σ µν != λσ
            G[µ,ν] += P[λ,σ]*value # J µνλσ
            G[λ,σ] += P[µ,ν]*value # J λσµν

            G[µ,λ] -= 0.5*P[ν,σ]*value*alpha # K   µνλσ
            G[λ,µ] -= 0.5*P[σ,ν]*value*alpha # K   λσµν
        else
            G[µ,ν] += P[λ,σ]*value # J µνλσ
            G[µ,λ] -= 0.5*P[ν,σ]*value*alpha # K µνλσ
        end
    end

    #Calculate Fxc contribution
    #Fxc = Fxc_form(integrals,P)
    #F = integrals.Hcore + G - Fxc
    #Need to scale K separately and Fx 
    return F
end

"""
Fock for HybridDFT and UKS: J,K and XC
"""
function Fock(Pi,Pj,dim,integrals::Jint_hybridKS)
    println("Inside UFock for Jint_hybridKS")
    println("This one is not ready")
    exit()
end


