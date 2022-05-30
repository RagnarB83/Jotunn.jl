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
function bf_atom_mapping(bset)
    mapping=Int64[]
    for atom in 1:length(bset.atoms)
        atom_bfs=sum([bset[atom][i].l for i in 1:length(bset[atom])].*2 .+1)
        for j in 1:atom_bfs
            push!(mapping,atom)
        end
    end
    return mapping
end


"""
create_bf_shell_map: Create array (dim:basis dim) of (atom, shell) tuples,
i.e. what atom and shell each BF belong to
"""
function create_bf_shell_map(bset)
    bf_atom_shell_mapping=Tuple{Int64, Int64}[]
    for atom_ind in 1:bset.natoms
        shells = bset[atom_ind]
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




############################
# MANUAL INTEGRALS
############################

#print("Reading integrals from file")
#Set basis dimension manually
#dim=2
#Read 1-electron integrals
#println("Reading file: tint.dat")
#T = read_integrals_from_file("tint.dat",dim)
#println("Reading file: vint.dat")
#V = read_integrals_from_file("vint.dat",dim)
#Read and diagonalize overlap matrix S. Get transformation matrix rij_x
#println("Reading file: Sint.dat")
#S = read_integrals_from_file("Sint.dat",dim)
#Read 2-electron integrals
#println("Reading file: twoelecint.dat")
#tei = read_twoel_integrals_from_file("twoelecint.dat", dim) #TODO: Not ready

"""
read_integrals_from_file: Reading 1-electron integrals from file, e.g. vint.dat
"""
function read_integrals_from_file(file,dim)
    M=zeros(dim,dim)
    for (count,line) in enumerate(eachline(file))
        splitline=split(line)
        index1=parse(Int64,splitline[1])
        index2=parse(Int64,splitline[2])
        val=parse(Float64,splitline[3])
        M[index1,index2] = val
        if index1 != index2
            M[index2,index1] = val
        end
    end
    return M
end

"""
read_twoel_integrals_from_file: Reading 2-electron integrals from file
"""
#TODO: UNFINISHED
function read_twoel_integrals_from_file(file,dim)
    M=zeros(dim*2,dim*2)
    #X=zeros(dim,dim,dim,dim)
    for (count,line) in enumerate(eachline(file))
        println("line:", line)
        splitline=split(line)
        numcols=length(splitline)
        val=parse(Float64,splitline[numcols])
        indices=[parse(Int64,i) for i in splitline[1:end-1]]
        println("indices:", indices)
        println("val:", val)
        
        exit()
        #M[index1,index2] = 
    end
    return M
end
