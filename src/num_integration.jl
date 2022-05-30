using PyCall

#Interface to Numgrid: a Rust library with Python API
#TODO: Call Rust library directly via C-Interface
#Write separate Julia Interface



"""
numgrid_call: Function that calls the numgrid py interface 
"""
function numgrid_call(gridtype,fragment, bset, radial_precision=1.0e-12,min_num_angular_points=86,max_num_angular_points=302,hardness=3)

    #Defining numgrid Python part
    py"""
    import numgrid
    def numgrid_atomgrid(coordinates_bohr, alpha_max=None, alpha_min = None, proton_charges=None, 
        radial_precision=1.0e-12, min_num_angular_points=86,max_num_angular_points=302, hardness=3):
        print("radial_precision: ", radial_precision)
        print("min_num_angular_points: ", min_num_angular_points)
        print("max_num_angular_points: ", max_num_angular_points)
        print("Center coordinates(au): ", coordinates_bohr)
        print("proton_charges:", proton_charges)
        print("alpha_max:", alpha_max)
    
        for center_index in range(len(coordinates_bohr)):
            # atom grid using explicit basis set parameters
            coordinates, weights = numgrid.atom_grid(
                alpha_min[center_index],
                alpha_max[center_index],
                radial_precision,
                min_num_angular_points,
                max_num_angular_points,
                proton_charges,
                center_index,
                coordinates_bohr,
                hardness=hardness,
            )
        return coordinates, weights
    
    def LMG_radial_grid(nuc_charge=None, alpha_min=None, alpha_max=None, radial_precision=1.0e-12):
        
        # radial grid (LMG) using explicit basis set parameters
        radii, weights = numgrid.radial_grid_lmg(
            alpha_min=alpha_min,
            alpha_max=alpha_max,
            radial_precision=radial_precision,
            proton_charge=nuc_charge,
        )
        return radii, weights
    def KK_radial_grid(num_points):
        # radial grid with 100 points using Krack-Koster approach
        radii, weights = numgrid.radial_grid_kk(num_points=num_points)
    
    def angular_grid(num_points):
        # angular grid with 14 points
        coordinates, weights = numgrid.angular_grid(num_points=num_points)
    """

    #Converting coordinates to Bohr and into vector of Tuples
    coordinates_bohr = [Tuple(fragment.coords[i,:]*1.8897161646321) for i in 1:size(fragment.coords,1)]
    proton_charges = fragment.nuccharges #array of nuclear charges

    #Getting basis set highest exponents for each atom
    atoms_largest_expons_per_atom=[]
    atoms_smallest_expons_per_atom_and_l=[]
    for atom in eachindex(fragment.elems)
        smallest_exp_for_angmom=Dict()
        #Collecting exponents for all primitives for all bfs
        allexponents_atom=[bset[atom][i].exp for i in 1:length(bset[atom])]
        #Collecting angular momentums for each shell
        allangmoms_atom=[bset[atom][i].l for i in 1:length(bset[atom])]
        unique_angmoms=unique(allangmoms_atom)
        #Finding largest (regardless of angmom)
        largest_expon = maximum(collect(Iterators.flatten(allexponents_atom)))
        #Looping over unique angmoms and finding smallest exponent for each angmom
        for angmom in unique_angmoms
            smallest = minimum(collect(Iterators.flatten([bset[atom][i].exp for i in 1:length(bset[atom]) if bset[atom][i].l == angmom])))
            #Adding to dict
            smallest_exp_for_angmom[angmom]=smallest
        end

        push!(atoms_largest_expons_per_atom,largest_expon)
        push!(atoms_smallest_expons_per_atom_and_l,smallest_exp_for_angmom)
    end
    if gridtype == "atomgrid"
        #alpha_min and alpha_max gives:
        #coordinates, weights = py"numgrid_atomgrid(30,min_num_angular_points=86,max_num_angular_points=434)"
        
        #Getting grid for all atoms
        grid_coords_bohr, gridweights = py"numgrid_atomgrid"(coordinates_bohr, alpha_min=atoms_smallest_expons_per_atom_and_l, alpha_max=atoms_largest_expons_per_atom, proton_charges=proton_charges, radial_precision=radial_precision, 
            min_num_angular_points=min_num_angular_points,max_num_angular_points=max_num_angular_points,hardness=hardness)
    elseif gridtype == "LMG+Leb"
        radii, gridweights = py"LMG_radial_grid"(basis_name=basis_name)
        # An angular grid with e.g. 14 points
        grid_coords_bohr, gridweights = py"angular_grid"(14)
        #gridcoords=grid_coords_bohr*0.529177249
    elseif gridtype == "KK+Leb"
        # Krack-Koster radial grid. e.g. 100 points
        radii, gridweights = py"KK_radial_grid"(100)
        # An angular grid with e.g. 14 points
        grid_coords_bohr, gridweights = py"angular_grid"(14)
        #gridcoords=grid_coords_bohr*0.529177249

    end
    return grid_coords_bohr, gridweights
end

"""
integrate_density
"""
function integrate_density(⍴,gridweights)
    N=0.0
    for (⍴_val,weight) in zip(⍴,gridweights)
        N+= ⍴_val*weight
    end
    return N
end