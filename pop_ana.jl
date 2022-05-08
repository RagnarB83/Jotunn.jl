"""
Mulliken for closed-shells
"""
function mulliken(S,P,bset,elems)
    PS=P*S
    #Basis function charges
    #diagDS #array of diagonal elements of DS
    #Reduced orbital charges?
    #Overlap charges?

    atom_populations=zeros(length(bset.atoms))
    charges=zeros(length(bset.atoms))
    num_bf_per_atoms=zeros(length(bset.atoms))
    j=1
    #Iterating over atom and summing over bfs belonging to atom
    for atom in 1:length(bset.atoms)
        #Number of basis functions for this atom
        atom_bfs=sum([bset[atom][i].l for i in 1:length(bset[atom])].*2 .+1)
        num_bf_per_atoms[atom] = atom_bfs
        #Sum of diagonal elements for bfs belonging to atom
        total_pop=sum([PS[i,i] for i in j:j+atom_bfs-1])
        atom_populations[atom] = total_pop
        charges[atom] = bset.atoms[atom].Z - total_pop
        j+=atom_bfs
    end
    return charges
end

"""
Mulliken for open-shells
"""
function mulliken(S,P,R,bset,elems)
    PS=P*S; RS=R*S
    charges=zeros(length(bset.atoms))
    spinpops=zeros(length(bset.atoms))
    num_bf_per_atoms=zeros(length(bset.atoms))
    j=1
    #Iterating over atom and summing over bfs belonging to atom
    for atom in 1:length(bset.atoms)
        atom_bfs=sum([bset[atom][i].l for i in 1:length(bset[atom])].*2 .+1)
        num_bf_per_atoms[atom] = atom_bfs
        total_pop_P=sum([PS[i,i] for i in j:j+atom_bfs-1])
        total_pop_R=sum([RS[i,i] for i in j:j+atom_bfs-1])
        charges[atom] = bset.atoms[atom].Z - total_pop_P
        spinpops[atom] = total_pop_R
        j+=atom_bfs
    end
    return charges, spinpops
end


"""
Mayer bond order function. Returns dictionary of all atom pairs
"""
function Mayer_BO(S,P,R,bset_atom_mapping)
    numatoms=bset_atom_mapping[end]
    PS=P*S; RS=R*S
    MBOdict=Dict{Tuple,Float64}() #dictionary: (A,B) = B_AB
    for A_atom in 1:numatoms
        #atom A bf indices
        A_indices = findall( x -> x == A_atom, bset_atom_mapping)
        for B_atom in A_atom+1:numatoms
            B_AB=0.0
            B_indices = findall( x -> x == B_atom, bset_atom_mapping)
            for µ in A_indices
                for ν in B_indices
                    B_AB += PS[µ,ν]*PS[ν,µ]+RS[µ,ν]*RS[ν,µ]
                end
            end
            MBOdict[(A_atom,B_atom)] = B_AB
        end
    end
    return MBOdict
end
