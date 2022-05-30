
#Fock algorithms based on pre-computed 2-electron integrals as 
#either full rank-4 tensor or condensed sparse version
###################################################################

"""
choose_Fock: Choosing Fock algorithm based on user-chosen WFtype,fock_algorithm variables etc.
"""
function choose_Fock(WFtype,DFTobj,fock_algorithm,tei_type,printlevel)
    print_if_level("Choosing Fock algorithm.",1,printlevel)
    #DFTobj only used if RKS/UKS
    if WFtype == "RHF"
        if tei_type == "4c"
            if fock_algorithm == "loop" #fast for small systems
                print_if_level("Using Fock_loop",1,printlevel)
                Fock=Fock_loop_RHF
            elseif fock_algorithm == "turbo" #faster for larger systems (but long compilation).
                print_if_level("Using Fock_turbo (requires LoopVectorization to be loaded",1,printlevel)
                Fock=Fock_turbo_RHF
            end
        elseif tei_type == "sparse4c"
            if fock_algorithm == "loop" || fock_algorithm == "loop_sparse" #generally recommended
                print_if_level("Using Fock_loop_sparse",1,printlevel)
                Fock=Fock_loop_sparse_RHF
                fock_algorithm="loop_sparse"
            else
                println("unknown fock_algorithm choice")
                exit()
            end
        end
    elseif WFtype == "RKS"
        if tei_type == "sparse4c"
            #Non-hybrid DFT
            if DFTobj.hybrid == false
                Fock=NonHybridDFT_RKS
                fock_algorithm="NonHybridDFT_RKS"
                print_if_level("Using $fock_algorithm",1,printlevel)
            #Hybrid-DFT
            else
                Fock=HybridDFT_RKS
                fock_algorithm="HybridDFT_RKS"
                print_if_level("Using $fock_algorithm",1,printlevel)
            end
        else
            println("Only tei_type=sparse4c allowed for RKS")
            exit()
        end
    elseif WFtype == "UKS"
        println("UKS not ready")
        exit()
    else #UHF
        if tei_type == "4c"
            if fock_algorithm == "loop"
                Fock=Fock_loop_UHF
            elseif fock_algorithm == "turbo"
                print_if_level("Using Fock_turbo_UHF (requires LoopVectorization to be loaded",1,printlevel)
                Fock=Fock_UHF_turbo_UHF
            end
        elseif tei_type == "sparse4c"
            if fock_algorithm == "loop"
                print_if_level("Using Fock_loop_sparse_UHF",1,printlevel)
                Fock=Fock_loop_sparse_UHF
                fock_algorithm="loop_sparse(UHF)"
            else
                println("unknown fock_algorithm choice")
                exit()
            end
        end
    end
    return Fock,fock_algorithm
end


###################################################################
# Fock functions that use full sparse form of 2-electron integrals
# via Jint struct
###################################################################

"""
Fock_loop_sparse: Multiple matrix-element update per unique set of indices
"""
function Fock_loop_sparse_RHF(P,dim,integrals::Jint)
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
function Fock_loop_sparse_UHF(Pi,Pj,dim,integrals::Jint)
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
Fock_loop: Fock-matrix RHF case loop-version
Benchmark: RHF/def2-QZVPP on H2O with MBA: 
243.604541  seconds total, 3.481820 seconds per Fock  (1 THREAD)
"""
function Fock_loop_RHF(P,dim,integrals::Jint_4rank)
    printdebug("This is regular Fock loop")
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
Fock_UHF_loop: Fock-matrix UHF case using simple loop
"""
function Fock_loop_UHF(Pi,Pj,dim,integrals::Jint_4rank)
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

"""
Fock_turbo: Fock-matrix RHF case with @turbo macro
Benchmark: RHF/def2-QZVPP on H2O with MBA: 
63.693465  seconds total, 0.325363 seconds per Fock (3-4 threads)
77.170879 seconds total, 0.320401 per Fock (1 thread)
Warning: Uses threading by default. Turn off by export OMP_NUM_THREADS=1
"""
#   function Fock_turbo_RHF(Hcore,P,dim,tei::Array{Float64,4})
#      #println("Fock_turbo disabled.")
#      #exit()
#       JK = zeros(dim,dim)
#       @turbo for µ in 1:dim
#           for ν in 1:dim
#               for λ in 1:dim
#                   for σ in 1:dim
#                       JK[µ,ν] += P[λ,σ]*(tei[ν,µ,λ,σ]-0.5*tei[ν,λ,µ,σ])
#                   end
#               end
#           end
#       end
#       F = Hcore + JK
#       return F
#   end

 """
 Fock_UHF_turbo: Fock-matrix UHF-case with @turbo macro
 """
#  function Fock_turbo_UHF(Hcore,Pi,Pj,dim,tei::Array{Float64,4})
#     JK = zeros(dim,dim)
#     @turbo for µ in 1:dim
#         for ν in 1:dim
#             for λ in 1:dim
#                 for σ in 1:dim
#                     JK[µ,ν] += Pi[λ,σ]*tei[ν,µ,λ,σ]+Pj[λ,σ]*tei[ν,µ,λ,σ]-Pi[λ,σ]*tei[ν,λ,µ,σ]
#                 end
#             end
#         end
#     end
#     F = Hcore + JK
#     return F
# end

#"""
#Fock_tullio: Fock with Tullio and Loopvectorization library (n)
#Benchmark: RHF/def2-QZVPP on H2O with MBA: 
#268.934002 seconds total, 3.69 seconds per Fock  (ONLY TULLIO LIBRARY LOADED)
#63.543055 seconds total, 0.354601 seconds per Fock  (LOOPVECT AND TULLIO LIBRARY LOADED, THREADED)
#72.814342 seconds total, 0.398519 seconds per Fock  (LOOPVECT AND TULLIO LIBRARY LOADED, 1 THREAD)
#"""
#function Fock_tullio(Hcore,P,dim,tei::Array{Float64,4})
#    JK = zeros(dim,dim)
#    @tullio JK[µ,ν] += P[λ,σ]*(tei[ν,µ,λ,σ]-0.5*tei[ν,λ,µ,σ])
#    F = Hcore + JK
#    return F
#end


#"""
#Fock matrix using tensor macro (REQUIRES TensorOperations)
#Benchmark: RHF/def2-QZVPP on H2O with MBA: 
#109.18 seconds total, 1.13 seconds per Fock (1 thread)
#"""
#function Fock_tensor(Hcore,P,dim,tei::Array{Float64,4})
#    JK = zeros(dim,dim)
#    @tensor JK[µ,ν] += P[λ,σ]*(tei[ν,µ,λ,σ]-0.5*tei[ν,λ,µ,σ])
#    F = Hcore + JK
#    return F
#end


###################################################################
# Kohn-Sham matrix fomrations
###################################################################

"""
Fock_loop_sparse: Multiple matrix-element update per unique set of indices
"""
function NonHybridDFT_RKS(P,dim,integrals::Jint)
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
    println("J: $J")
    println("")

    #Calculate Vxc contribution from P, bfs on gridpoints
    Fxc = Fxc_form(integrals,P)

    F = integrals.Hcore + J - Fxc
    return F
end

"""
HybridDFT_RKS: 
"""
function HybridDFT_RKS(P,dim,integrals::Jint)
    G = zeros(dim,dim)
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

            G[µ,λ] -= 0.5*P[ν,σ]*value # K   µνλσ
            G[ν,λ] -= 0.5*P[µ,σ]*value # K   νµλσ
            G[µ,σ] -= 0.5*P[ν,λ]*value # K   µνσλ
            G[λ,µ] -= 0.5*P[σ,ν]*value # K   λσµν
            G[λ,ν] -= 0.5*P[σ,µ]*value # K   λσνµ
            G[ν,σ] -= 0.5*P[µ,λ]*value # K   νµσλ
            G[σ,ν] -= 0.5*P[λ,µ]*value # K   σλνµ
            G[σ,µ] -= 0.5*P[λ,ν]*value # K   σλµν

        elseif γλσ && γxy
            G[µ,ν] += 2*P[λ,σ]*value # J µνλσ
            G[λ,σ] += P[µ,ν]*value # J λσµν
            G[σ,λ] += P[ν,µ]*value # J σλνµ

            G[µ,λ] -= 0.5*P[ν,σ]*value # K   µνλσ
            G[µ,σ] -= 0.5*P[ν,λ]*value # K   µνσλ
            G[λ,µ] -= 0.5*P[σ,ν]*value # K   λσµν
            G[σ,ν] -= 0.5*P[λ,µ]*value # K   σλνµ
        elseif γµν && γxy
            G[µ,ν] += P[λ,σ]*value # J µνλσ
            G[ν,µ] += P[λ,σ]*value # J νµλσ
            G[λ,σ] += 2*P[µ,ν]*value # J λσµν

            G[µ,λ] -= 0.5*P[ν,σ]*value # K   µνλσ
            G[ν,λ] -= 0.5*P[µ,σ]*value # K   νµλσ
            G[λ,µ] -= 0.5*P[σ,ν]*value # K   λσµν
            G[λ,ν] -= 0.5*P[σ,µ]*value # K   λσνµ

        elseif γµν && γλσ
            G[µ,ν] += 2*P[λ,σ]*value # J µνλσ
            G[ν,µ] += 2*P[λ,σ]*value # J νµλσ

            G[µ,λ] -= 0.5*P[ν,σ]*value # K   µνλσ
            G[ν,λ] -= 0.5*P[µ,σ]*value # K   νµλσ
            G[µ,σ] -= 0.5*P[ν,λ]*value # K   µνσλ
            G[ν,σ] -= 0.5*P[µ,λ]*value # K   νµσλ
        elseif γxy #µ=ν; λ=σ µν != λσ
            G[µ,ν] += P[λ,σ]*value # J µνλσ
            G[λ,σ] += P[µ,ν]*value # J λσµν

            G[µ,λ] -= 0.5*P[ν,σ]*value # K   µνλσ
            G[λ,µ] -= 0.5*P[σ,ν]*value # K   λσµν
        else
            G[µ,ν] += P[λ,σ]*value # J µνλσ
            G[µ,λ] -= 0.5*P[ν,σ]*value # K µνλσ
        end
    end

    #Calculate Fxc contribution
    #Fxc = Fxc_form(integrals,P)
    #F = integrals.Hcore + G - Fxc
    #Need to scale K separately and Fx 
    exit()
    return F
end



"""
Fxc formation
"""
function Fxc_form(integrals,P)
    println("Inside Fxc_form")
    dim=integrals.bset.nbas
    #Gridpoints inside integrals object
    gridpoints=integrals.gridpoints
    gridweights=integrals.gridweights
    println("gridpoints: $gridpoints")
    println("gridweights: $gridweights")
    #Create density array at selected gridpoints
    ⍴ = create_density(P,integrals,gridpoints)

    println("⍴:", ⍴)
    println("gridpoints: $gridpoints")


    #Calculate xc contribution from ⍴ at each gridpoint
    #Call manual LDA exchange functional

    if xc_method == "manual"
        e = LDA_Exchange(⍴,gridweights)
    elseif xc_method == "libxc"
        #Call LibXC
    else
        Fxc=zeros(integrals.bset.nbas,integrals.bset.nbas)
    end
    exit()

    return Fxc
end


function LDA_Exchange(⍴,gridweights,dim)
    E_xc=0.0
    V_xc=zeros(dim,dim)
    prefactor_E=-(3/4)*(3/π)^(1/3)
    prefactor_V=-(6/π)^(1/3)
    #for ⍴_val,weight in zip(⍴,gridweights)
    #    E_xc+= prefactor*⍴_val^(4/3)*weight
    #    #V_xc+= prefactor*⍴_val^(1/3)*weight
    #end
    return E_xc
end

function basisfunc_sum(shell,gp,atom_xyz)
    #println("Inside basisfunc_sum")    
    # r from gp and atom_xyz
    #println("gp :$gp and atom_xyz:$atom_xyz")
    #Gridpoint XYZ position relative to atom 
    x = gp[1] - (atom_xyz[1]*1.88973)
    y = gp[2] - (atom_xyz[2]*1.88973)
    z = gp[3] - (atom_xyz[3]*1.88973)

    #
    #println("Shell is: $shell with angmom l: $(shell.l)")
    #println("shell.coef: ", shell.coef)
    #println("shell.exp: ", shell.exp)
    #println("shell.l: ", shell.l)

    spherharm=1 #?? TODO
    val=0.0
    for (coef,⍺) in zip(shell.coef,shell.exp)
        #println("coef: $coef , ⍺: $⍺")
        val+= (x^2+y^2+z)^shell.l * coef*spherharm*exp(-1*⍺*(x^2+y^2+z^2))
    end
    return val
end



"""
create_density: Create density from P (density matrix) and basis functions 
for input gridpoints
"""
function create_density(P,integrals)
    gridpoints=integrals.gridpoints
    ⍴=zeros(length(gridpoints)) # array of rho-values for each gridpoint
    bf_atom_shell_map=integrals.bf_atom_shell_map
    dim=integrals.bset.nbas
    #println("bf_atom_shell_map: $bf_atom_shell_map")
    for (gp_i,gp) in enumerate(gridpoints)
        #println("=====================================")
        #println("gp: $gp")
        #println("⍴ (current)", ⍴)
        x,y,z = gp
        #println("x,y,z: $x, $y, $z")
        #Looping over basis function indices
        for µ in 1:dim
            for ν in 1:dim
                #println("Basis functions µ: $µ and ν: $ν")
                #Get atom and shell tuple for each bf
                atom_shell_µ_index=bf_atom_shell_map[µ] #(atomindex,shell) for bf µ
                atom_shell_ν_index=bf_atom_shell_map[ν] #(atomindex,shell) for bf ν
                atom_µ = atom_shell_µ_index[1] #atom-index that µ belongs to
                atom_ν = atom_shell_ν_index[1] #atom-index that ν belongs to
                #Atom XYZ coords
                atom_µ_xyz=integrals.bset.atoms[atom_µ].xyz
                atom_ν_xyz=integrals.bset.atoms[atom_ν].xyz
                #For each bf µ we get the shell 
                #TODO: Multiply by degeneracy instead ??
                shell_µ = integrals.bset[atom_µ][atom_shell_µ_index[2]]
                shell_ν = integrals.bset[atom_ν][atom_shell_ν_index[2]]

                #Now get the basis function value at position
                ϕ_µ_bfval = basisfunc_sum(shell_µ,gp,atom_µ_xyz)
                ϕ_ν_bfval = basisfunc_sum(shell_ν,gp,atom_ν_xyz)
               
                #Calculating density at point
                densval=P[µ,ν] * ϕ_µ_bfval * ϕ_ν_bfval
                #println("densval: $densval")
                ⍴[gp_i] += densval
            end
        end
    end
   return ⍴
end

"""
make_grid()
"""
function make_grid()
    gridpoints=[[0.0,0.0,0.0],[0.0,0.0,0.5],[0.0,0.5,0.0],[0.5,0.0,0.0],
    [0.0,0.5,0.5],[0.0,0.5,0.5],[0.5,0.0,0.5],[0.5,0.5,0.5],[10.0,10.0,10.0]]
    return gridpoints
end

