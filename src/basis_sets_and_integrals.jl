#Calculate 2-electron integrals via GaussianBasis.jl
"""
tei_calc: Calculate 2-electron integrals in 4 different ways
"""
function tei_calc(bset,tei_type,printlevel)
    if tei_type == "4c"
        #Full 4-rank  tensor (1-based indexing already used here)
        #TODO: Create Jint object instead ??
        return ERI_2e4c(bset)
    elseif tei_type == "sparse4c"
        #Sparse form. Only unique sets of indices and integral values
        sparseintegrals = sparseERI_2e4c(bset)
        print_if_level("Sparse integral calculation done.",1,printlevel)
        return sparseintegrals
    elseif tei_type == "3c"
        #3c version. Requires auxiliary basis set
        #TODO: Create Jint object instead ??
        return ERI_2e3c(bset, aux)
    elseif tei_type == "2c"
        #2c version
        #TODO: Create Jint object instead ??
        return ERI_2e2c(bset)
    end
    return integrals
end

"""
Calculate all 2-el index permutations for a given set of indices. Called by Fock.
"""
function permutations!(perms,a,b,c,d)
    #println("Inside permutations, 1perms: $perms")
    #a,b,c,d=tup
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
    γab = a !== b
    γcd = c !== d
    γxy = ab !== cd
    #Setting everything to zero

    perms .= 0 #Setting vector elements to zero
    
    #println("Inside permutations, 2perms: $perms")
    #println("Here: perms: $perms")
    if γab && γcd && γxy
        #println("Case1: e.g. a!=b!=c!=d  and 1232 ")
        #return [(a,b,c,d),(b,a,c,d),(a,b,d,c),(c,d,a,b),(c,d,b,a),(b,a,d,c),(d,c,b,a),(d,c,a,b)]
        perms[1:32] .= a,b,c,d,b,a,c,d,a,b,d,c,c,d,a,b,c,d,b,a,b,a,d,c,d,c,b,a,d,c,a,b
        return nothing
        #return SA[(a,b,c,d),(b,a,c,d),(a,b,d,c),(c,d,a,b),(c,d,b,a),(b,a,d,c),(d,c,b,a),(d,c,a,b)]
    elseif γcd && γxy
        #println("Case2")
        perms[1:16] .= a,b,c,d,a,b,d,c,c,d,a,b,d,c,b,a
        return nothing
        #return [(a,b,c,d),(a,b,d,c),(c,d,a,b),(d,c,b,a)]
        #return SA[(a,b,c,d),(a,b,d,c),(c,d,a,b),(d,c,b,a)]
    elseif γab && γxy
        #println("Case3 (incomplete)")
        perms[1:16] .= a,b,c,d,b,a,c,d,c,d,a,b,c,d,b,a
        return nothing
        #return [(a,b,c,d),(b,a,c,d),(c,d,a,b),(c,d,b,a)]
        #return SA[(a,b,c,d),(b,a,c,d),(c,d,a,b),(c,d,b,a)]
    elseif γab && γcd
        #println("Case4")
        perms[1:16] .= a,b,c,d,b,a,c,d,a,b,d,c,b,a,d,c
        return nothing
        #return [(a,b,c,d),(b,a,c,d),(a,b,d,c),(b,a,d,c)]
        #return SA[(a,b,c,d),(b,a,c,d),(a,b,d,c),(b,a,d,c)]
    elseif γxy #a=b; c=d ab != cd
        #println("Case5: example: 1122")
        perms[1:8] .= a,b,c,d,c,d,a,b
        return nothing
        #return [(a,b,c,d),(c,d,a,b)]
        #return SA[(a,b,c,d),(c,d,a,b)]
    else
        #println("Case6, else. a=a=a=a. Adding nothing(already there)")
        perms[1:4] .= a,b,c,d
        #println("Inside permutations, 3perms: $perms")
        return nothing
        #return [(a,b,c,d)]
        #return SA[(a,b,c,d)]
    end
end


"""
read_basis_file: Read ORCA-style basis set file and return dictionary of element:basis-info-dict
"""
function read_basis_file(filename,elems; format="orca")
    print_if_level("Reading basis-set file: $filename",1,printlevel)
    #Element names for ORCA-style basis set files
    element_names=Dict("HYDROGEN" => "H", "CARBON" => "C", "HELIUM" => "He", 
    "HELIUM" => "He", "CARBON" => "C", "NITROGEN" => "N", "OXYGEN" => "O", "FLUORINE" => "F")

    #Dictionary of system elements with empty list to be filled
    elem_basis_dict=Dict()
    for el in elems
        elem_basis_dict[el]=[]
    end

    grab_basis=false
    grab_bf=false
    num_prims=0
    el=nothing
    ang_mom=nothing
    #Found elements
    found_elements=[]
    for (count,line) in enumerate(eachline(filename))
        if grab_basis == true
            #println("grab_basis is true")
            linesplit=split(line)
            if grab_bf == true
                #println("yes grabbing bf")
                gaussian = split(line)
                num_prim = parse(Int,gaussian[1])
                exponent= parse(Float64,gaussian[2])
                coeff = parse(Float64,gaussian[3])
                push!(elem_basis_dict[el][end],[num_prim,exponent,coeff])
                if num_prim == num_prims
                    #println("last primitive. setting grab_bf to false")
                    grab_bf=false
                end
            end
            if length(linesplit) == 2
                #println("NEW BASIS FUNCTION STARTING")
                grab_bf=true
                ang_mom=linesplit[1]
                #Adding new basis function: e.g. [S,3,[]]
                num_prims=parse(Int,String(linesplit[2]))
                #println("elem_basis_dict:", elem_basis_dict)
                push!(elem_basis_dict[el],[ang_mom,num_prims])
            end
        end

        #Trim line and then check if in element_names
        line_trimmed=line
        #print
        try
            el = element_names[line_trimmed]
            #println("FOUND ELEMENT $(line_trimmed) in BASIS FILE: $el")
            grab_basis=false
            if el in elems
                #Only grab basis set for element if system contains that element
                println("Found system element $el in basis-set-file")
                push!(found_elements,el)
                grab_basis=true
            end
        catch E
        end

    end

    #Check if basis-set was found for all system elements in file
    for el in elems
        if (el in found_elements) == false
            println("Error:System element $el was not found in basis-set file")
            println("Exiting!")
            exit()
        end
    end

    return elem_basis_dict
end

"""
basis_set_create: From basis-set-name, elements and coordinates, create a BasisSet object (GaussianBasis.jl)
Option: built-in basis, manual definition or read from ORCA-format basis file
"""
function basis_set_create(basis,elems,coordinates; basisfile="none",printlevel)
    print_if_level("Creating basis set object",1,printlevel)
    #get coordinates as a multi-linestring (for Molecules.jl):
    coords_string = array_to_string(elems,coordinates)
    if basis == "manual"
        #Coords as Molecules atoms object
        atoms = Molecules.parse_string(coords_string)

        #Basis functions
        #Example: Manual STO-3G
        #Angular moment number, expansion coeff, exponents
        # 1 s basis function from 3 primitive Gaussians 
        s_bf = BasisFunction(0, [0.1543289673,0.5353281423,0.4446345422],
        [.425250914, 0.6239137298, 0.1688554040])
        print("s_bf : $s_bf")

        #Map  between basis functions and atoms
        #Here single H-atom only
        bmap = [[s_bf], [s_bf]]
        #bmap = [[s_bf], [s_bf]]   # S bf on each H-atom

        bset = BasisSet(basis, atoms, bmap)
    elseif basis == "read_file"
        angmom_dict = Dict("S" => 0, "P" => 1, "D" => 2, "F" => 3, "G" => 4, "H" => 5, "I" => 6)
        basisname="basis-from-file"
        #Coords as Molecules atoms object
        atoms = Molecules.parse_string(coords_string)
        #0. Read basis set from file into dictionary
        basis_dict = read_basis_file(basisfile,elems,format="orca")
        #Initialize Basis function element dictionary
        #Dict of element: [GaussianBasis,GaussianBasis]
        elems_BasisFunctions_dict=Dict()
        for el in elems
            elems_BasisFunctions_dict[el] = []
        end
        for elbas in basis_dict
            list_of_bfs=[]
            if length(elbas[2]) > 0
                for bs in elbas[2]
                    element=elbas[1]
                    angmom=bs[1]
                    numprims=bs[2]
                    allprims=[bs[i] for i in 3:3+numprims-1]
                    exponents=hcat(allprims...)[2,:]
                    coeffs=hcat(allprims...)[3,:]
                    bf = BasisFunction(angmom_dict[angmom], coeffs, exponents)
                    push!(elems_BasisFunctions_dict[element],bf)
                end
            end
        end
        #Now create vector of basis functions in order of atoms
        basismap=Vector{BasisFunction}[]
        for el in elems
            push!(basismap,[])
            for basisfunc in elems_BasisFunctions_dict[el]
                push!(basismap[end],basisfunc)
            end
        end
        bset = BasisSet(basisname, atoms, basismap)

    #Use internal basis sets in GaussianBasis.jl
    else
        #Using basis variable as basis-set name
        #TODO: coords_string fails if scientific notation present in number
        bset = BasisSet(basis, coords_string)
        try
            #Basis set via internal basis-set name and coordinates-string
            bset = BasisSet(basis, coords_string)
        catch e
            println("Unknown basis set name. Not found internally in GaussianBasis.jl")
            exit()
        end

    end

    return bset
end

"""
array_to_string: Converting coordinate-array to multi-line string.
Used in BasisSet creation
"""
function array_to_string(elems,coords)
    coords_string=""
    for i in 1:length(elems)
        el = elems[i]
        c_row = coords[i,:]
        c1=@sprintf "%.14f" c_row[1];c2=@sprintf "%.14f" c_row[2];c3=@sprintf "%.14f" c_row[3]
        c_line=c1*" "*c2*" "*c3
        coords_string = coords_string * el * " " * c_line * "\n"
    end
    return coords_string
end


"""
bf_atom_mapping: Created simply array (size of basis-set dimension) of atom indices in basis-function order.
i.e. map of which basis-function belongs to which atom
TODO: Get rid of ?? Since create_bf_shell_map contains same info??
"""
#function bf_atom_mapping(bset)
#    println("yy")
#    mapping=Int64[]
#    for atom in 1:length(bset.atoms)
#        atom_bfs=sum([bset[atom][i].l for i in 1:length(bset[atom])].*2 .+1)
#        #atom_bfs=sum([bset[atom][i].l for i in 1:bset.shells_per_atom[atom]].*2 .+1)
#        for j in 1:atom_bfs
#            push!(mapping,atom)
#        end
#    end
#    return mapping
#end


"""
create_bf_shell_map: Create array (dim:basis dim) of (atom, shell) tuples,
i.e. what atom and shell each BF belong to
"""
function old_create_bf_shell_map(bset)
    bf_atom_shell_mapping=Tuple{Int64, Int64}[]
    for atom_ind in 1:bset.natoms
        shells = bset[atom_ind]
        #numshells = bset.shells_per_atom[atom_ind]
        ind=0
        for shell in shells
            angmom=shell.l
            #degen=angmom_degen[angmom]
            degen=2*angmom+1
            ind+=1
            for i in 1:degen
                tup=(atom_ind,ind)
                push!(bf_atom_shell_mapping,tup)
            end
            
        end
    end
    return bf_atom_shell_mapping
end



"""
create_bf_shell_map: Create array (dim:basis dim) of (atom, shell) tuples,
i.e. what atom and shell each BF belong to
"""
function create_bf_shell_map(bset)
    bf_atom_shell_mapping=Tuple{Int64, Int64}[]
    #Looping over shells
    atom_ind=1 #atom-index counter
    shellcounter=0 #shell-index counter that is reset
    for (i,shell) in enumerate(bset.basis)
        shellcounter+=1
        #if shell index larger than shells per atom
        if shellcounter > bset.shells_per_atom[atom_ind]
            atom_ind+=1 #new atom
            shellcounter=1 #resetting shell counter
        end
        for j in 1:2*shell.l+1 #looping over degenerate bfs in shell
            push!(bf_atom_shell_mapping,(atom_ind,i))
        end
    end
    return bf_atom_shell_mapping
end

"""
get_basis_composition: could be used to print primitive composition
"""
function get_basis_composition(bset)
    #Basis set composition dicts per atom index
    BFs_contracted_dict=Dict{Int64,Vector{Any}}()
    BFs_primitives_dict=Dict{Int64,Vector{Any}}()
    basis_func_angmoms=[]; basis_func_prims=[]
    #Looping over shells
    atom_ind=1 #atom-index counter
    shellcounter=0 #shell-index counter that is reset
    for shell in bset.basis
        shellcounter+=1
        #if shell index larger than shells per atom
        if shellcounter > bset.shells_per_atom[atom_ind]
            BFs_contracted_dict[atom_ind] = basis_func_angmoms
            BFs_primitives_dict[atom_ind] = basis_func_prims
            atom_ind+=1 #new atom
            shellcounter=1 #resetting shell counter
            basis_func_prims=[]
            basis_func_angmoms=[]
        end
        push!(basis_func_angmoms,shell.l)
        for i in range(1,length(shell.exp))
            push!(basis_func_prims,shell.l)
        end
    end
    #Last item
    BFs_contracted_dict[atom_ind] = basis_func_angmoms
    BFs_primitives_dict[atom_ind] = basis_func_prims
    return BFs_contracted_dict,BFs_primitives_dict
end


"""
create_ml_values: Create lists of all m_l values for each basis function.
"""
function create_ml_values(bset)
    mlvalues=[]
    for shell in bset.basis
        for i in shell.l:-1:-shell.l
            push!(mlvalues,i) 
        end
    end
    return mlvalues
end

"""
create_l_values: Create lists of all angmom l values for each basis function.
"""
function create_l_values(bset)
    lvalues=[]
    for shell in bset.basis
        for i in shell.l:-1:-shell.l
            push!(lvalues,shell.l)
        end
    end
    return lvalues
end

"""
get_bf: For basis function i, get GB shell object, atom-position of that bf and ml value.
"""
function get_bf(integrals,i)
    #Get atom and shell indices for bf i
    atom_shell_i_index=integrals.bf_atom_shell_map[i] #(atomindex,shell) for bf i
    atom_i = atom_shell_i_index[1] #atom-index that i belongs to
    atom_i_xyz=integrals.bset.atoms[atom_i].xyz
    shell_i = integrals.bset[atom_i][atom_shell_i_index[2]] #Get shell for i
    m_i=integrals.mlvalues[i]
    return shell_i, atom_i_xyz, m_i
end

"""
Evaluate basis function at gridpoint in spherical harmonics.
NOTE: Incomplete, have not figured out angular part
See: https://pyscf.org/_modules/pyscf/gto/mole.html#cart2sph
"""
function basisfunc_eval_sph(shell,gp,atom_xyz,m)
    #Gridpoint XYZ position relative to atom in Bohrs
    x = gp[1] - (atom_xyz[1]*1.88973)
    y = gp[2] - (atom_xyz[2]*1.88973)
    z = gp[3] - (atom_xyz[3]*1.88973)
    #println("gridpoint xyz in Bohr: $x  $y  $z")
    #Simple normalization from GaussianBasis (angular momenta already incoporated into coefficients)
    #N_sp=sqrt(1/(4pi))
    l=shell.l
    if shell.l == 0
        #N_sp=1.0
        N_sp=sqrt(1/(4*pi)) #0.282094791773878143
    elseif shell.l == 1
        N_sp=0.488602511902919921*m
        #if m == 1
        #    N_sp=0.488602511902919921
        #    #N_sp=sqrt(1/(4*pi))
        #elseif m== 0
        #    N_sp=sqrt(1/(4*pi)) #0.282094791773878143
        #else
        #    N_sp=sqrt(1/(4*pi)) #0.282094791773878143
        #end
            #N_sp=sqrt((2*l+1)/(4*pi)*sqrt((factorial(l-abs(m)))/factorial(l+abs(m))))
        #println("m: $m")
        #N_sp=1.0
        #N_sp=0.2102
        #N_sp=0.28
        #N_sp=sqrt(1/(4*pi))
        #N_sp=((-1)^m)*sqrt(((2*l+1)*factorial((l-abs(m))))/(4*pi*factorial((l+abs(m)))))
        #N_sp=sqrt(1/(4*pi))
        #N_sp=(-1)^(m)*sqrt(((1))/(4*pi))
        #N_sp=sqrt((2*l+1)/(4*pi)*sqrt((factorial(l-m))/factorial(l+m)))
        #println("N_sp:", N_sp)
    elseif shell.l == 2
        #N_sp=sqrt(5/(16*pi))
        #N_sp=0.21157109383040862 #sqrt(1/(4*pi))*3/4
        N_sp=0.2102 #between this and 0.2103
        #N_sp=sqrt(1/(4*pi))
    #    l=shell.l
        #https://www2.atmos.umd.edu/~dkleist/docs/shtns/doc/html/spec.html
        #N_sp=sqrt((2*l+1)/(4*pi)*sqrt((factorial(l-m))/factorial(l+m)))
        #N_sp=(-1)^(m)*sqrt(((1))/(4*pi))
        #N_sp=0.745*sqrt(1/(4pi))
        #N_sp=sqrt((2*l+1)/(4*pi)*sqrt((factorial(l-abs(m)))/factorial(l+abs(m))))
    end

    #Actual normalization used
    #sqrt((2^(2*l+3)*factorial(l+1)*(2*a)^(l+1.5))/(factorial(2*l+2)*sqrt(pi)))

    #N_sp=sqrt((2*l+1)/(4*pi)*sqrt((factorial(l-m))/factorial(l+m)))
    #N_sp=((-1)^m)*sqrt(((2*l+1)*factorial((l-m)))/(4*pi*factorial((l+m))))
    #N_sp=sqrt((2*l+1)/(4*pi)*sqrt((factorial(l-abs(m)))/factorial(l+abs(m))))
    #if shell.l == 2
    #    blox=3
    #else
    blox=shell.l
    #end
    #Looping over primitives
    val=0.0
    for (coef,⍺) in zip(shell.coef,shell.exp)
        val+= ((x+y+z)^blox) * coef * N_sp * exp(-1*⍺*(x^2+y^2+z^2))
    end
    return val
end

"""
basisfunc_eval2: Doing BF evaluation differently
"""
function basisfunc_eval2(shell,gp,atom_xyz,m)
    #Gridpoint XYZ position relative to atom in Bohrs
    x = gp[1] - (atom_xyz[1]*1.8897161646321)
    y = gp[2] - (atom_xyz[2]*1.8897161646321)
    z = gp[3] - (atom_xyz[3]*1.8897161646321)

    r2=(x^2+y^2+z^2)

    val=0.0
    #S-function
    if shell.l == 0
        i=0;j=0;k=0
        N=0.282094791773878143
        #Looping over primitives
        for (coef,⍺) in zip(shell.coef,shell.exp)
            val += N * coef * exp(-⍺*r2)
        end
        return val
    #P-function
    elseif shell.l == 1
        N=0.488602511902919921
        if m == 1
            #i=0;j=0;k=0
            i=0;j=1;k=0
        elseif m == 0
            #i=0;j=0;k=0
            i=0;j=0;k=1
        elseif m == -1
            #i=0;j=0;k=0
            i=1;j=0;k=0
        end
        #if m == 1
        #    i=1;j=0;k=0
        #elseif m==0
        #    i=0;j=1;k=0
        #elseif m==-1
        #    i=0;j=0;k=1
        #end
        #Looping over primitives
        for (coef,⍺) in zip(shell.coef,shell.exp)
            val += N * x^i * y^j * z^k * coef * exp(-⍺*r2)
            #val += N * (x+y+z)^shell.l * coef * exp(-⍺*r2)
        end
        return val
    #D-function
    elseif shell.l == 2
        N=0.6307831305050401 # sqrt(5/(4pi))
        #https://www.theochem.ru.nl/~pwormer/Knowino/knowino.org/wiki/Gaussian_type_orbitals.html
        #N=1.145*0.488602511902919921
        if m == 2 #xy
            M=1.7320508075688772 #sqrt(3)
            i=1;j=1;k=0
            #Looping over primitives
            for (coef,⍺) in zip(shell.coef,shell.exp)
                val += N * M * x*y * coef * exp(-⍺*r2)
                #val += N * M * x^i*y^j*z^k* coef * exp(-⍺*r2)
            end
        return val
        elseif m==1 #yz
            M=1.7320508075688772 #sqrt(3)
            i=0;j=1;k=1
            #Looping over primitives
            for (coef,⍺) in zip(shell.coef,shell.exp)
                val += N * M * y*z * coef * exp(-⍺*r2)
                #val += N * M * x^i*y^j*z^k* coef * exp(-⍺*r2)
            end
        elseif m==0 #xz
            M=1.7320508075688772 #sqrt(3)
            i=1;j=0;k=1
            #Looping over primitives
            for (coef,⍺) in zip(shell.coef,shell.exp)
                val += N * M * x*z * coef * exp(-⍺*r2)
                #val += N * M * x^i*y^j*z^k* coef * exp(-⍺*r2)
            end
        elseif m==-1 #3*z2-r^2
            M=0.5
            #Looping over primitives
            for (coef,⍺) in zip(shell.coef,shell.exp)
                val+= N * M * (3*z^2-r2) * coef * exp(-⍺*r2)
            end
            return val
        elseif m==-2 # x2-y2
            M=0.8660254037844386 #0.5*sqrt(3)
            #Looping over primitives
            for (coef,⍺) in zip(shell.coef,shell.exp)
                val+= N * M * (x^2-y^2) * coef * exp(-⍺*r2)
            end
            return val
        else
            println("WTF!")
            exit()
        end
    #F-function
    #http://www.rsc.org/suppdata/c7/cp/c7cp00194k/c7cp00194k1.pdf
    elseif shell.l == 3
        #TODO: Incomplete. returning 0.0 for now
        #println("Warning: f-basis functions barely supported. Result will be inaccurate")
        return 0.0
    #G-function
    elseif shell.l == 4
        #TODO: Incomplete. returning 0.0 for now
        #println("Warning: g-basis functions barely supported. Result will be inaccurate")
        return 0.0
    else
        println("unsupported angmom")
        exit()
        
    end

    return val
end


"""
BFvalues_calc:Precalculate basisfunction values at gridpoints
Creates array of size numgridpoints x numbfs

TODO: Derivative

TODO: Add exclamation mark
"""
function BFvalues_calc(integrals)
    dim=integrals.bset.nbas

    BF=zeros(length(integrals.gridpoints),dim)
    #∇BF=zeros(length(integrals.gridpoints),dim)
    for (i,gp) in enumerate(integrals.gridpoints)
        for µ in 1:dim
            #Get basis functions µ
            shell_µ,atom_µ_xyz,m_µ = get_bf(integrals,µ)

            #Now get the basis function value at gridpoint position
            #ϕ_µ = basisfunc_eval_sph(shell_µ,gp,atom_µ_xyz,m_µ)  DOES NOT WORK
            ϕ_µ = basisfunc_eval2(shell_µ,gp,atom_µ_xyz,m_µ) #DOES WORK
            #ϕ_µ = basisfunc_eval_libcint(shell_µ,gp,atom_µ_xyz,m_µ) #TODO: Call libcint routine directly

            #Add to array
            BF[i,µ] = ϕ_µ

            #Derivative of BF         TODO
            #∇BF[i,µ] =

        end
    end
    integrals.BFvalues=BF
    #integrals.∇BFvalues=∇BF
end
