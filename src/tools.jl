
"""
elem_to_nuccharge: Get nuclear charge for element-string
"""
function elem_to_nuccharge(elem)
    elem_nuccharge_dict=Dict("H"=>1.0, "He"=>2.0, "Li"=>3.0, "Be"=>4.0, "B"=>5.0, "C"=>6.0, "N"=>7.0, "O"=>8.0, "F"=>9.0, "Ne"=>10.0 )
    if haskey(elem_nuccharge_dict,elem) == false
        println("Error: elem_to_nuccharge does not have this element defined. Needs to be fixed!")
        exit()
    end
    charge=elem_nuccharge_dict[elem]
    return charge
end

"""
Calculate nuclear-nuclear repulsion from elements and and coordinates
"""
function nuc_nuc_repulsion(elems,coords)
    ang2bohr = 1.88972612546
    num=length(elems)
    charges=[elem_to_nuccharge(el) for el in elems]
    coords_b=coords*ang2bohr
    VC=0.0
    @inbounds for j in 1:num
        for i in j+1:num
            @inbounds rij_x = coords_b[i,1] - coords_b[j,1]
            @inbounds rij_y = coords_b[i,2] - coords_b[j,2]
            @inbounds rij_z = coords_b[i,3] - coords_b[j,3]
            @fastmath r = rij_x*rij_x+rij_y*rij_y+rij_z*rij_z
            @fastmath d = sqrt(r)
            @inbounds @fastmath VC += charges[i] * charges[j] / (d)
        end
    end
    return VC
end

"""
Calculate density matrix P from MO coefficient matrix C
Scaling: 2.0 for RHF, 1.0 for UHF
"""
function makedensity(C, dim, Norb,scaling=2.0)
    P=zeros(dim,dim)
    for µ in 1:dim
        for ν in 1:dim
            for o in 1:Norb
                P[µ,ν] += scaling*C[µ,o]*C[ν,o]
            end
        end
    end
    return P
end

"""
Create array of occupation numbers
"""
function makeoccupationarray(norb,dim,val)
    occ=zeros(dim)
    for i in 1:norb
        occ[i]=val
    end
    return occ
end

