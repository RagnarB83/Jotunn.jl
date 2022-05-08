#Basic Fock algorithms based on pre-computed 2-electron integrals as full rank-4 tensor

"""
Choosing Fock algorithm based on HFtype, user-input and system-size
"""
function choose_Fock(HFtype,fock_algorithm,dim,tei_type)
    println("Choosing Fock algorithm.")
    #Choosing loop algo for small system (Loopvectorization compilation takes a bit long)
    #H2O def2-TZVPP and smaller benefits
    if dim < 70
        println("Small system (Basis dim: $dim). Choosing loop Fock.")
        fock_algorithm="loop"
    end
    if HFtype == "RHF"
        if tei_type == "4c"
            if fock_algorithm == "loop" #fast for small systems
                Fock=Fock_loop
            elseif fock_algorithm == "turbo" #fastest for larger systems. LoopVectorization
                Fock=Fock_turbo
            end
        elseif tei_type == "sparse4c"
            if fock_algorithm == "loop" #fast for small systems
                Fock=Fock_loop_sparse
                fock_algorithm="loop_sparse"
            else
                println("unknown choice")
                exit()
            end
        end
    else #UHF
        if fock_algorithm == "loop" #fast for small systems
            Fock=Fock_UHF_loop
        elseif fock_algorithm == "turbo" #fastest for larger systems. LoopVectorization
            Fock=Fock_UHF_turbo
        end
    end
    return Fock,fock_algorithm
end

"""
Basic loop-version of the Fock-matrix
Benchmark: RHF/def2-QZVPP on H2O with MBA: 
243.604541  seconds total, 3.481820 seconds per Fock  (1 THREAD)
"""
function Fock_loop(Hcore,P,dim,tei)
    println("This is regular Fock loop")
    JK = zeros(dim,dim)
    @simd for µ in 1:dim
        for ν in 1:dim
            for λ in 1:dim
                for σ in 1:dim
                    #println("µ: $µ, ν: $ν, λ: $λ, σ: $σ  Value:", P[λ,σ]*(tei[ν,µ,λ,σ]-0.5*tei[ν,λ,µ,σ]))
                    #println("J tei[ν,µ,λ,σ]:", tei[ν,µ,λ,σ])
                    #println("K tei[ν,λ,µ,σ]:", tei[ν,λ,µ,σ])
                    @inbounds JK[µ,ν] += P[λ,σ]*(tei[µ,ν,λ,σ]-0.5*tei[µ,λ,ν,σ])
                    #JK[µ,ν] += 2*P[λ,σ]*tei[µ,ν,λ,σ] wrong
                    #JK[µ,ν] -= P[λ,σ]*tei[µ,λ,ν,σ] wrong
                end
            end
        end
    end
    F = Hcore + JK
    return F
end

"""
Sparse-integral loop-version of the Fock-matrix. 
"""
#Fock_loop(Hcore,P,dim,indices,tei)
function Fock_loop_sparse(Hcore,P,dim,tei)
    println("THIS IS Fock_loop_sparse")
    JK = zeros(dim,dim)
    sparse_length=length(tei[1])
    #Looping over tuples of indices known not to be zero
    println("tei: $tei")
    println("sparse_length: $sparse_length")
    for i in 1:sparse_length
        println("\ni: $i")
        indices=tei[1][i]
        val=tei[2][i]
        println("indices: $indices")
        µ=indices[1]+1; ν=indices[2]+1; λ=indices[3]+1; σ=indices[4]+1 
        println("µ: $µ, ν: $ν, λ: $λ, σ: $σ")
        #JK[µ,ν] += P[λ,σ]*(tei[2][ν,µ,λ,σ]-0.5*tei[ν,λ,µ,σ])
        if µ == ν == λ == σ
            println("case 1")
            println("val: $val")
            J=P[λ,σ]*val
            K=-1*0.5*P[λ,σ]*val
            JK[µ,ν] += J+K
            println("JKcont added:", J+K)
        elseif µ != ν && λ != σ && µν != λσ #Example: [1,2,3,1] => [1,2,3,1],[2,1,3,1],[2,1,1,3],[1,2,1,3]
            
            println("val: $cal")
            J=P[λ,σ]*val
            K=-1*0.5*P[λ,σ]*val
            #Coulomb
            JK[µ,ν] += 2*P[λ,σ]*val
            JK[ν,µ] += 2*P[λ,σ]*val #could be done later?
            JK[λ,σ] += 2*P[µ,ν]*val
            JK[σ,λ] += 2*P[ν,µ]*val #could be done later?


        else
            println("case 3")
            JK[µ,ν] += P[λ,σ]*(val-0.5*val)
            JK[ν,µ] += P[σ,λ]*(val-0.5*val)
            println("Val added:", P[λ,σ]*(val-0.5*val))
            println("Val added:", P[λ,σ]*(val-0.5*val))
        end
        print("---------------------------------------")
    end
    F = Hcore + JK
    return F
end


"""
Unrestricted Fock matrix.
"""
function Fock_UHF_loop(Hcore,Pi,Pj,dim,tei)
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
Using LoopVectorization turbo macro
Benchmark: RHF/def2-QZVPP on H2O with MBA: 
63.693465  seconds total, 0.325363 seconds per Fock (3-4 threads)
77.170879 seconds total, 0.320401 per Fock (1 thread)
Warning: Uses threading by default. Turn off by export OMP_NUM_THREADS=1
"""
 function Fock_turbo(Hcore,P,dim,tei)
     JK = zeros(dim,dim)
     @turbo for µ in 1:dim
         for ν in 1:dim
             for λ in 1:dim
                 for σ in 1:dim
                     JK[µ,ν] += P[λ,σ]*(tei[ν,µ,λ,σ]-0.5*tei[ν,λ,µ,σ])
                 end
             end
         end
     end
     F = Hcore + JK
     return F
 end

 function Fock_UHF_turbo(Hcore,Pi,Pj,dim,tei)
    JK = zeros(dim,dim)
    @turbo for µ in 1:dim
        for ν in 1:dim
            for λ in 1:dim
                for σ in 1:dim
                    JK[µ,ν] += Pi[λ,σ]*tei[ν,µ,λ,σ]+Pj[λ,σ]*tei[ν,µ,λ,σ]-Pi[λ,σ]*tei[ν,λ,µ,σ]
                end
            end
        end
    end
    F = Hcore + JK
    return F
end

"""
Fock_tullio: Fock with Tullio and Loopvectorization library (n)
Benchmark: RHF/def2-QZVPP on H2O with MBA: 
268.934002 seconds total, 3.69 seconds per Fock  (ONLY TULLIO LIBRARY LOADED)
63.543055 seconds total, 0.354601 seconds per Fock  (LOOPVECT AND TULLIO LIBRARY LOADED, THREADED)
72.814342 seconds total, 0.398519 seconds per Fock  (LOOPVECT AND TULLIO LIBRARY LOADED, 1 THREAD)
"""
#function Fock_tullio(Hcore,P,dim,tei)
#    JK = zeros(dim,dim)
#    @tullio JK[µ,ν] += P[λ,σ]*(tei[ν,µ,λ,σ]-0.5*tei[ν,λ,µ,σ])
#    F = Hcore + JK
#    return F
#end


"""
Fock matrix using tensor macro (REQUIRES TensorOperations)
Benchmark: RHF/def2-QZVPP on H2O with MBA: 
109.18 seconds total, 1.13 seconds per Fock (1 thread)
"""
#function Fock_tensor(Hcore,P,dim,tei)
#    JK = zeros(dim,dim)
#    @tensor JK[µ,ν] += P[λ,σ]*(tei[ν,µ,λ,σ]-0.5*tei[ν,λ,µ,σ])
#    F = Hcore + JK
#    return F
#end