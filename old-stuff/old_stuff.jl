
function get_energy(P,Hcore,F,dim,E_ZZ)
    E=E_ZZ #E containing nuc-repuls
    for mu in 1:dim
        for nu in 1:dim
            E += 0.5*P[mu,nu]*(Hcore[mu,nu]+F[mu,nu])
        end
    end
    return E
end

function newget_energy(P,Hcore,F,dim,E_ZZ)
    E=E_ZZ #E containing nuc-repuls
    for mu in 1:dim
        for nu in 1:dim
            E += 2*Hcore[mu,nu]*P[mu,nu]
            for r in 1:dim
                for s in 1:dim
                    E += G[μ,ν,ρ,σ]*D[σ,ρ]*D[μ,ν]
                end
            end
        end
    end
    return E
end

function convert_to_matrix(c)
    d=zeros(length(c),3)
    for i in 1:length(c)
        d[i,1]=c[i][1]
        d[i,2]=c[i][2]
        d[i,3]=c[i][3]
    end
    return d
end


function print_Mulliken(charges,elems)
    println("\n\nMulliken Population Analysis")
    println("-"^40)
    println("             Charge")
    println("-"^40)
    for i in 1:length(elems)
            charge=charges[i]
            spin=""
            el=elems[i]
            @printf("%4i %4s %10.5f\n", i, el, charge)
    end
    println("-"^40)
    @printf("\nSum of charges: %.3f\n", sum(charges))
end

function print_Mulliken(charges,elems,spinpops)
    println("\n\nMulliken Population Analysis (open-shell)")
    println("-"^40)
    println("             Charge     Spin.pop")
    println("-"^40)
    for i in 1:length(elems)
        charge=charges[i]
        spin=spinpops[i]
        el=elems[i]
        @printf("%4i %4s %10.5f %10.5f\n", i, el, charge, spin)
    end
    println("-"^40)
    @printf("\nSum of charges: %.3f\n", sum(charges))
    @printf("Sum of spin populations: %.3f\n", sum(spinpops))
end

function print_Mayer_analysis(MBOs,elems)
    println("\nMayer bond orders:")
    for mbo in MBOs
        mbo_key_1=mbo[1][1]; mbo_key_2=mbo[1][2]
        mbo_val=mbo[2]
        elem1=elems[mbo_key_1]; elem2=elems[mbo_key_2]
        println("$elem1$mbo_key_1 - $elem2$mbo_key_2 : $mbo_val")
    end
end



"""
print_MO_energies
"""
function print_MO_energies(occ,mos)
    println("-"^30)
    println("MO Energies (closed-shell y)")
    println("-"^30)
    println("                          ⍺")
    println("-"^50)
    @printf "%5s%12s%12s%12s" "MO" "Occ." "E(Eh)" "E(eV)\n"
    #for (i,orb) in enumerate(zip(eps)
    for (i,(n,orb)) in enumerate(zip(occ,mos))
        orbev=orb*27.211399
        @printf "%5i%12.4f%12.6f%12.4f\n" i n orb orbev
    end
    println("-"^50)
end

"""
print_MO_energies: open-shell
"""
function print_MO_energies(occ_⍺, occ_β, mos_⍺, mos_β)
    println("-"^30)
    println("MO Energies (open-shell)")
    println("-"^30)
    println("                          ⍺                                         β")
    println("-"^100)
    @printf "%5s%12s%12s%12s      %12s%12s%12s" "MO" "Occ." "E(Eh)" "E(eV)" "Occ." "E(Eh)" "E(eV)\n"
    for (i,(n_⍺,n_β,orb_⍺,orb_β)) in enumerate(zip(occ_⍺,occ_β, mos_⍺, mos_β))
        @printf "%5i%12.4f%12.6f%12.4f      %12.4f%12.6f%12.4f\n" i n_⍺ orb_⍺ orb_⍺*27.211399 n_β orb_β orb_β*27.211399
    end
    println("-"^100)
end




"""
Fock_loop_sparse: Sparse-integral loop-version RHF case of the Fock-matrix. 
Not ready
"""
#Fock_loop(Hcore,P,dim,indices,tei)
function oldFock_loop_sparse(Hcore,P,dim,tei)
    println("THIS IS Fock_loop_sparse (does not work yet")
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
Yoshimine sort: gives compound index from four indices.
From: https://adambaskerville.github.io/posts/HartreeFockGuide/
"""
function yoshimine_sort(a,b,c,d)
    if a > b
        ab = a*(a+1)/2 + b
    else
        ab = b*(b+1)/2 + a
    end
    if c > d
        cd = c*(c+1)/2 + d
    else
        cd = d*(d+1)/2 + c
    end
    if ab > cd
        abcd = ab*(ab+1)/2 + cd
    else
        abcd = cd*(cd+1)/2 + ab
    end
    return floor(Int64,abcd)
end

#Extending getindex function so that object of type Jint support 4-index 
#function getindex(x::Jint,i::Int64,j::Int64,k::Int64,l::Int64)
    #println("x: $x")
    #println("i: $i j:$j k:$k l:$l")
    #bla=yoshimine_sort(i,j,k,l)
    #println("bla: $bla")

    #NOTE: if lookupindex in dict fails that means the integral was never put in dict, meaning it must be 0.0
    #Hence returning 0.0
    #THIS IS VERY SLOW
    #try
    #    return x.integraldict[yoshimine_sort(i,j,k,l)]
    #catch
    #    return 0.0
    #end
    #return x.integraldict[yoshimine_sort(i,j,k,l)]
#end


#NOTE: Idea, to integral cutoffs here?
function yosh_lib(integrals,dim)
    intdict=Dict{Int64,Array{Tuple}}()
    #Looping over all permutations
    for µ in 1:dim
        for ν in 1:dim
            for λ in 1:dim
                for σ in 1:dim
                    res = yoshimine_sort(µ,ν,λ,σ)
                    if haskey(intdict, res)
                        push!(intdict[res],(µ,ν,λ,σ))
                    else
                        intdict[res]=[]
                        push!(intdict[res],(µ,ν,λ,σ))
                    end
                end
            end
        end
    end
    return intdict
end

"""
Reads unique indices-tuples and integral values of sparse GaussianBasis and creates 
a Jint object with a list of all indices-permutations and values in same order
"""
function ordered_int_table_sparse(integrals)
    #Basis dimension
    dim=integrals[1][end][1] + 1
    #
    println("Starting yosh_lib")
    @time intdict = yosh_lib(integrals,dim)
    println("Dond with yosh_lib")
    #
    integral_all_indices=Vector{NTuple{4, Int16}}[]
    println("Starting loop")
    #Looping over ordered unique tuples and getting all permutations in order
    @time for i in 1:length(integrals[1])
        uniq_tuple=integrals[1][i] .+ 1
        #println("uniq_tuple: $uniq_tuple")
        res = yoshimine_sort(uniq_tuple...)
        val = intdict[res]
        #println("val: $val")
        #exit()
        #for (key,val) in intdict
        #    println("key: $key val: $val")
        #    if uniq_tuple in val
        #        push!(integral_all_indices,val)
        #    end
        #end
        push!(integral_all_indices,val)
    end
    println("loop done")
    #Creating Jint object: original unique indices (but with 1-based indexing),
    #All index permutations, all value
    #println("integral_all_indices: $integral_all_indices")
    #println("typeof $(typeof(integral_all_indices))")
    println("Creating Jint object")
    @time jintobject = Jint_sparse([i .+ 1 for i in integrals[1]],integral_all_indices,integrals[2])
    #println("jintobject: $jintobject")
    #println(jintobject.values)
    #println(jintobject.all_indices)
    #println(jintobject.unique_indices)
    #exit()
    return jintobject
end


"""
Takes sparseERI_2e4c object from GaussianBasis/libcint and creates a dict of compound indices as 
keys and integral values as values
"""
function lookup_table_sparse_create(integrals)
    integraldict=Dict{Int64, Float64}()
    for (count,i) in enumerate(integrals[1])
        z = i .+ 1 #shifted by 1
        index = yoshimine_sort(z...)
        val=integrals[2][count]
        integraldict[index] = val
    end
    #Create integral object
    jintobject = Jint(integraldict)
    return jintobject
end

"""
Grabs 2-el integral value from sparse integral using all 4 indices
"""
function get_tei(i,j,k,l)
    return sp[yoshimine_sort(i,j,k,l)]
end



# Fock_loop_sparse2: Version without creating permutations.
#Instead we do multiple F matrix-element assignments for each unique set

#NOTE: Keeping copy of this function before modifications. Works well
function Fock_loop_sparse2(Hcore,P,dim,tei)

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
Fock_loop_sparse_perm: RHF Sparse-integral loop-version RHF of Fock-matrix with permutations
Allocates too much and hence slower
"""
#Fock_loop(Hcore,P,dim,indices,tei)
function Fock_loop_sparse_perm_RHF(Hcore,P,dim,tei::Jint_sparse)
    #println("THIS IS Fock_loop_sparse")
    G = zeros(dim,dim)
    perms= zeros(Int8,32) #8*4=32 is the max number of indices
    #perms = SVector{32}(zeros(32))
    #Looping over tuples of indices from sparse2e4c
    for i in eachindex(tei.unique_indices)
        #This is the unique set of indices provided by GaussianBasis/libcint (sparse)
        #uniq_tuple=tei.unique_indices[i]
        #println("Calling permutations")
        permutations!(perms,tei.unique_indices[i]...)
        #println("After perms:", perms)
        #println("type perms :", typeof(perms))
        #exit()
        #println("perms: $perms type: $(typeof(perms))")
        # This is the associated integral value
        #println("val: ")
        value=tei.values[i] 
        #println("here") 
        #Looping over all symmetry-related sets of indices of uniq_tuple
        #for p in perms
        #println("part loop begin")
        for (µ,ν,λ,σ) in partition(perms, 4)
            if µ == 0
                break
            end
            #println("type of perms[p]", typeof(perms[p]))
            #µ,ν,λ,σ=p
            @inbounds G[µ,ν] += P[λ,σ]*value # J
            @inbounds G[µ,λ] -= 0.5*P[ν,σ]*value # K
        end
    end
    F = Hcore + G
    return F
end