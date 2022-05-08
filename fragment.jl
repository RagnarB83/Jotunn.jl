export Fragment, create_fragment, read_xyzfile

#Python-function usage in this file:
"""
py_functions_coords.totmasslist(elems)
py_functions_coords.nucchargelist(elems)
"""


#Keyword arguments with structs (without Parameters.jl)
#Field without defaults are required keyword arguments at object creation
#Future Julia will likely not need @kwdef
Base.@kwdef mutable struct Fragment
    elems::Array{String,1}
    coords::Array{Float64,2}
    energy::Float64 = 0.0
    connectivity::Array{Array{Int64,1}} = []
    atomcharges::Array{Float64,1} = []
    atomtypes::Array{String,1} = []
    Centralmainfrag::Array{Float64,1} = []
    fragmenttype_labels::Array{String,1} = []
    formula::String = ""
    prettyformula::String = ""
    label::String = ""
    charge::Any=missing
    mult::Any=missing
    nuccharge::Int64=0
    mass::Float64 = 0.0
    numatoms::Int64=0
    atomlist::Array{Float64,1} = []
    allatoms::Array{Float64,1} = []
    list_of_masses::Array{Float64,1} = []
end


eldict_covrad=Dict(["H" => 0.31, "He" => 0.28, "Li" => 1.28, "Be" => 0.96, "B" => 0.84, "C" => 0.76, "N" => 0.71, "O" => 0.66, "F" => 0.57, "Ne" => 0.58, "Na" => 1.66, "Mg" => 1.41, "Al" => 1.21, "Si" => 1.11, "P" => 1.07, "S" => 1.05, "Cl" => 1.02, "Ar" => 1.06, "K" => 2.03, "Ca" => 1.76, "Sc" => 1.70, "Ti" => 1.6, "V" => 1.53, "Cr" => 1.39, "Mn" => 1.61, "Fe" => 1.52, "Co" => 1.50, "Ni" => 1.24, "Cu" => 1.32, "Zn" => 1.22, "Ga" => 1.22, "Ge" => 1.20, "As" => 1.19, "Se" => 1.20, "Br" => 1.20, "Kr" => 1.16, "Rb" => 2.2, "Sr" => 1.95, "Y" => 1.9, "Zr" => 1.75, "Nb" => 1.64, "Mo" => 1.54, "Tc" => 1.47, "Ru" => 1.46, "Rh" => 1.42, "Pd" => 1.39, "Ag" => 1.45, "Cd" => 1.44, "In" => 1.42, "Sn" => 1.39, "Sb" => 1.39, "Te" => 1.38, "I" => 1.39, "Xe" => 1.40])


#Function to create fragment object of type Fragment. Better for user than constructor or object-creation statement
function create_fragment(;coords_string=nothing,xyzfile=nothing,pdbfile=nothing,fragfile=nothing, coords=nothing,
    elems=nothing, calc_connectivity=false, label=nothing, charge=nothing, mult=nothing)

    #If file argument passed. Other
    if xyzfile != nothing
        elems,coords=read_xyzfile(xyzfile)
    elseif fragfile != nothing
        elems,coords=read_fragfile(fragfile)
    elseif pdbfile != nothing
        elems,coords=read_pdbfile(pdbfile)
    elseif coords_string != nothing
        elems,coords=read_coordsstring(coords_string)
    elseif coords == nothing
        println("No file argument or coords array passed to create_fragment")
        error()
    end

    #First creation of fragment
    println("Creating fragment")
    fragment=Fragment(elems=elems,coords=coords, charge=charge, mult=mult)

    #Set attributes like formula and mass
    #fragment.mass=py_functions_coords.totmasslist(elems)
    #fragment.list_of_masses=py_functions_coords.list_of_masses(elems)
    #fragment.nuccharge=py_functions_coords.nucchargelist(elems)
    fragment.numatoms=length(elems)
    fragment.formula=join(elems)

    if label != nothing
        fragment.label=label
    end
    #Connectivity if desired
    if calc_connectivity == true
        println("Calculating connectivity")
        conndepth=40
        scale=1.0
        tol=0.1
        eldict_covrad
        fragment.connectivity=Coordinates.calc_connectivity(coords,elems,conndepth,scale, tol,eldict_covrad)
    end
    return fragment
end


#Read multi-line string and return elements and coordinates
function read_coordsstring(mlstring)
    elems=[]
    numlines=count("\n",mlstring)
    #println("Reading $numlines lines")
    coords=zeros(numlines,3)
    atomcount=0
    for (index,line) in enumerate(split(mlstring,"\n"))
        if length(line) >0
           arr=split(line)
           el=string(arr[1])
           push!(elems,el)
           x=parse(Float64,arr[2])
           y=parse(Float64,arr[3])
           z=parse(Float64,arr[4])

           coords[index,1]=x
           coords[index,2]=y
           coords[index,3]=z
        end
    end 
    return elems,coords
end

#Read XYZ file
function read_xyzfile(filename)
    #Will accept atom-numbers as well as symbols
    elements = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K",
            "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
            "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs",
            "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
            "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
            "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"]
    println("Reading coordinates from XYZfile $filename ")
    elems=[]
    file = open(filename, "r")
    numatoms=nothing
    coords=nothing
    for (count,line) in enumerate(eachline(file))

        if count == 1
            numatoms=parse(Int64,line)
            coords=zeros(numatoms,3)
        elseif count == 2
            #title-line
        elseif count > 2
            atomnumber=count-2
            splitline=split(line)
            col1=splitline[1]
            #Check if element-symbol or atomic-numbers
            if tryparse(Int64,col1) == nothing
                push!(elems,col1)
            else
                el=parse(Int64,col1)
                elem=elements[el]
                push!(elems,elem)
            end

            x=parse(Float64,splitline[2])
            y=parse(Float64,splitline[3])
            z=parse(Float64,splitline[4])
            coords[atomnumber,1]=x
            coords[atomnumber,2]=y
            coords[atomnumber,3]=z
        end
    end
    close(file)
    @assert size(coords)[1] == numatoms "Number of coordinates does not match header line"
    @assert size(coords)[1] == length(elems) "Number of coordinates does not match elements."
    println("Read $numatoms atoms")
    return elems,coords
end