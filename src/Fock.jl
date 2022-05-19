
#Fock algorithms based on pre-computed 2-electron integrals as full
# rank-4 tensor or condensed sparse version
###################################################################


"""
choose_Fock: Choosing Fock algorithm based on user-chosen HFtype,fock_algorithm variables etc.
"""
function choose_Fock(HFtype,fock_algorithm,dim,tei_type)
    println("Choosing Fock algorithm.")
    #Choosing loop algo for small system (Loopvectorization compilation takes a bit long)
    #Disabling for now since we are not using LoopVectorization
    #H2O def2-TZVPP and smaller benefits
    #if dim < 50
    #    println("Small system (Basis dim: $dim). Choosing loop Fock.")
    #    fock_algorithm="loop"
    #end
    if HFtype == "RHF"
        if tei_type == "4c"

            # Change macro 
            #if fock4c_speedup =="turbo"
            #    println("fock4c_speedup requires LoopVectorization library")
            #    var"@speedup" = var"@turbo"
            #elseif fock4c_speedup =="simd"
            #    var"@speedup" = var"@simd"
            #else
            #    #Have nothing
            #    var"@speedup" = var"@show"
            #end

            if fock_algorithm == "loop" #fast for small systems
                println("Using Fock_loop")
                Fock=Fock_loop_RHF
            elseif fock_algorithm == "turbo" #faster for larger systems (but long compilation).
                println("Using Fock_turbo (requires LoopVectorization to be loaded")
                Fock=Fock_turbo
            end
        elseif tei_type == "sparse4c"
            if fock_algorithm == "loop" #generally recommended
                println("Using Fock_loop_sparse")
                Fock=Fock_loop_sparse_RHF
                fock_algorithm="loop_sparse"
            elseif fock_algorithm == "loop2" #elegant but slow version
                #SLOWER (to be removed)
                println("Using Fock_loop_sparse_perm")
                Fock=Fock_loop_sparse_perm_RHF
                fock_algorithm="loop_sparse_perm"
            else
                println("unknown choice")
                exit()
            end
        end
    else #UHF
        if tei_type == "4c"
            if fock_algorithm == "loop"
                Fock=Fock_loop_UHF
            elseif fock_algorithm == "turbo"
                println("Using Fock_turbo_UHF (requires LoopVectorization to be loaded")
                Fock=Fock_UHF_turbo_UHF
            end
        elseif tei_type == "sparse4c"
            if fock_algorithm == "loop"
                println("Using Fock_loop_sparse_UHF")
                Fock=Fock_loop_sparse_UHF
                fock_algorithm="loop_sparse(UHF)"
            else
                println("unknown choice")
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
function Fock_loop_sparse_RHF(Hcore,P,dim,tei::Jint_sparse)

    #- STILL NEED TO LOOK INTO CORRECTNESS
    #- Are exchange terms definitely correct ??? Look at Xik things
    #- Next combine assignments into simpler ones
    G = zeros(dim,dim)
    #Looping over unique sets of indices and values from sparse2e4c
    for i in eachindex(tei.unique_indices)
        #This is the unique set of indices provided by GaussianBasis/libcint (sparse)
        µ,ν,λ,σ=tei.unique_indices[i]
        value=tei.values[i] #integral value
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
            G[µ,ν] += P[λ,σ]*value # J µνλσ
            G[ν,µ] += P[λ,σ]*value # J νµλσ
            G[µ,ν] += P[σ,λ]*value # J µνσλ
            G[λ,σ] += P[µ,ν]*value # J λσµν
            G[λ,σ] += P[ν,µ]*value # J λσνµ
            G[ν,µ] += P[σ,λ]*value # J νµσλ
            G[σ,λ] += P[ν,µ]*value # J σλνµ
            G[σ,λ] += P[µ,ν]*value # J σλµν

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
            G[µ,ν] += P[λ,σ]*value # J µνλσ
            G[µ,ν] += P[σ,λ]*value # J µνσλ
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
            G[λ,σ] += P[µ,ν]*value # J λσµν
            G[λ,σ] += P[ν,µ]*value # J λσνµ

            G[µ,λ] -= 0.5*P[ν,σ]*value # K   µνλσ
            G[ν,λ] -= 0.5*P[µ,σ]*value # K   νµλσ
            G[λ,µ] -= 0.5*P[σ,ν]*value # K   λσµν
            G[λ,ν] -= 0.5*P[σ,µ]*value # K   λσνµ

        elseif γµν && γλσ
            #println("Case4")
            G[µ,ν] += P[λ,σ]*value # J µνλσ
            G[ν,µ] += P[λ,σ]*value # J νµλσ
            G[µ,ν] += P[σ,λ]*value # J µνσλ
            G[ν,µ] += P[σ,λ]*value # J νµσλ

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
    F = Hcore + G
    return F
end

"""
Fock_loop_sparse_UHF: Multiple matrix-element update per unique set of indices
"""
function Fock_loop_sparse_UHF(Hcore,Pi,Pj,dim,tei::Jint_sparse)
    G = zeros(dim,dim)
    #Looping over unique sets of indices and values from sparse2e4c
    for i in eachindex(tei.unique_indices)
        #This is the unique set of indices provided by GaussianBasis/libcint (sparse)
        µ,ν,λ,σ=tei.unique_indices[i]
        value=tei.values[i] #integral value
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
    F = Hcore + G
    return F
end
        #@inbounds JK[µ,ν] += Pi[λ,σ]*tei[ν,µ,λ,σ]+Pj[λ,σ]*tei[ν,µ,λ,σ]-Pi[λ,σ]*tei[ν,λ,µ,σ]


"""
Fock_loop_sparse_perm: RHF Sparse-integral loop-version RHF of Fock-matrix with permutations
Allocates too much
"""
#Fock_loop(Hcore,P,dim,indices,tei)
function Fock_loop_sparse_perm_RHF(Hcore,P,dim,tei::Jint_sparse)
    #println("THIS IS Fock_loop_sparse")
    G = zeros(dim,dim)
    #Looping over tuples of indices from sparse2e4c
    for i in eachindex(tei.unique_indices)
        #This is the unique set of indices provided by GaussianBasis/libcint (sparse)
        #uniq_tuple=tei.unique_indices[i]
        perms=permutations(tei.unique_indices[i]...)
        #println("perms: $perms type: $(typeof(perms))")
        # This is the associated integral value
        value=tei.values[i] 
        #println("here") 
        #Looping over all symmetry-related sets of indices of uniq_tuple
        for p in eachindex(perms)
        #for tuple in [(1,1,1,1),(1,1,1,2)]
            µ,ν,λ,σ=(perms[p])
            @inbounds G[µ,ν] += P[λ,σ]*value # J
            @inbounds G[µ,λ] -= 0.5*P[ν,σ]*value # K
        end
    end
    F = Hcore + G
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
function Fock_loop_RHF(Hcore,P,dim,tei::Array{Float64,4})
    printdebug("This is regular Fock loop")
    JK = zeros(dim,dim)
    @simd for µ in 1:dim
        for ν in 1:dim
            for λ in 1:dim
                for σ in 1:dim
                    @inbounds JK[µ,ν] += P[λ,σ]*(tei[µ,ν,λ,σ]-0.5*tei[µ,λ,ν,σ])
                end
            end
        end
    end
    F = Hcore + JK
    return F
end


"""
Fock_UHF_loop: Fock-matrix UHF case using simple loop
"""
function Fock_loop_UHF(Hcore,Pi,Pj,dim,tei::Array{Float64,4})
    JK = zeros(dim,dim)
    @simd for µ in 1:dim
        for ν in 1:dim
            for λ in 1:dim
                for σ in 1:dim
                    @inbounds JK[µ,ν] += Pi[λ,σ]*tei[ν,µ,λ,σ]+Pj[λ,σ]*tei[ν,µ,λ,σ]-Pi[λ,σ]*tei[ν,λ,µ,σ]
                end
            end
        end
    end
    F = Hcore + JK
    return F
end

"""
Fock_turbo: Fock-matrix RHF case with @turbo macro
Benchmark: RHF/def2-QZVPP on H2O with MBA: 
63.693465  seconds total, 0.325363 seconds per Fock (3-4 threads)
77.170879 seconds total, 0.320401 per Fock (1 thread)
Warning: Uses threading by default. Turn off by export OMP_NUM_THREADS=1
"""
#  function Fock_turbo_RHF(Hcore,P,dim,tei::Array{Float64,4})
#     println("Fock_turbo disabled.")
#     exit()
#      JK = zeros(dim,dim)
#      @turbo for µ in 1:dim
#          for ν in 1:dim
#              for λ in 1:dim
#                  for σ in 1:dim
#                      JK[µ,ν] += P[λ,σ]*(tei[ν,µ,λ,σ]-0.5*tei[ν,λ,µ,σ])
#                  end
#              end
#          end
#      end
#      F = Hcore + JK
#      return F
#  end

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

