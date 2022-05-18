
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
nuc_nuc_repulsion: Calculate nuclear-nuclear repulsion from elements and and coordinates
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
makedensity: Calculate density matrix P from MO coefficient matrix C
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
makeoccupationarray: Create array of occupation numbers
"""
function makeoccupationarray(norb,dim,val)
    occ=zeros(dim)
    for i in 1:norb
        occ[i]=val
    end
    return occ
end

#Check whether spin multiplicity is consistent with the nuclear charge and total charge
"""
check_multiplicity: Checks whether spin multiplicity is compatible with number of electrons
"""
function check_multiplicity(num_el,charge,mult)
    unpaired_electrons=mult-1
    result = map(iseven, (num_el,unpaired_electrons))
    if result[1] != result[2]
        println("The spin multiplicity $mult ($unpaired_electrons unpaired electrons) is incompatible with the total number of electrons $num_el")
        println("Exiting!")
        exit()
    end
end


#plot_MO and plot_density. Plots from arrays themselves.
#Write wrapper tools that are standalone from basisfile or add basifile option below

"""
Plot MO with index MOindex
"""
function plot_MO(MOindex,bset,C)

    #Read basis functions from bset and coordinates

    #MO coefficents from C

    #x_range=range(-5,5, step=0.1)
    #x_range=range(-5,5, length=100)

end

"""
Plot density
"""
function plot_density(P,bset,C; box_size=[-5,5], grid_size=100)

    isovalue=0.05
    #Ranges for x,y and z
    x_min = y_min =z_min = box_size[1]
    x_max =y_max = z_max = box_size[2]

    x_space=range(x_min,x_max, length=grid_size)
    y_space=range(y_min,y_max, length=grid_size)
    z_space=range(z_min,z_max, length=grid_size)
    dim=bset.nbas

    #rho=zeros(min_lim)

    #bset.basis. Basis functions ordered onto atoms in atom-ordering
    #Length .basis is numatoms. Vector{Vector{BasisFunction}}

    for x in range(x_lim_min,x_lim_max)
        for µ in 1:dim
            for ν in 1:dim
                ϕ_µ=bset[µ]
                ϕ_ν=bset[ν]
                rho[x] = P[µ,ν] * ϕ_µ * ϕ_ν
            end
        end
    end


end