"""
Integrate XC energy from e_xc values on a grid
"""
function integrate_XC_energy(e_xc_grid,⍴,integrals)
    E_xc=0.0
    #Looping over gridpoints and calculate Exc
    for i in 1:length(e_xc_grid)
        gw=integrals.gridweights[i]
        e_xc=e_xc_grid[i]
        E_xc += gw * e_xc * ⍴[i]
    end
    return E_xc
end

"""
Integrate XC potential from v_xc values on a grid
"""
function integrate_XC_potential(vrho_grid,⍴,integrals)
    dim=integrals.bset.nbas
    V=zeros(dim,dim)
    for µ in 1:dim
        for ν in 1:dim
            for i in 1:length(integrals.gridpoints)
                gw=integrals.gridweights[i]
                v_xc=vrho_grid[i]
                ϕ_µ = integrals.BFvalues[i,µ]
                ϕ_ν = integrals.BFvalues[i,ν]
                if integrals.lvalues[µ] == integrals.lvalues[ν]
                    #println("X1. shell_µ.l: $shell_µ.l and shell_ν.l:$shell_ν.l")
                    if integrals.mlvalues[µ] == integrals.mlvalues[ν]
                        #println("integrals.mlvalues[µ]: ", integrals.mlvalues[µ])
                        #println("integrals.mlvalues[ν]: ", integrals.mlvalues[ν])
                        V[µ,ν] +=  v_xc * ϕ_µ * ϕ_ν * gw
                    else
                        #println("x2. else. skip")
                    end
                else
                    #println("x3 else")
                    #V[µ,ν] +=  v_xc * ϕ_µ * ϕ_ν * gw
                end
            end
        end
    end
    return V
end


function E_x_LDA(integrals,Cx,⍴)
    E=0.0
    #Looping over gridpoints and calculate Exc
    for i in 1:length(integrals.gridpoints)
        gw=integrals.gridweights[i]
        E += gw * Cx * ⍴[i]^(1/3)*⍴[i]
    end
    return E
end

####
function V_x_LDA(integrals,CxLDA,⍴)
    println("V_x_LDA:", V_x_LDA)
    dim=integrals.bset.nbas
    V=zeros(dim,dim)
    for µ in 1:dim
        for ν in 1:dim
            #println("-------------------------------------------")
            #println("µ: $µ and ν:$ν")
            #Get basis functions µ and ν
            shell_µ,atom_µ_xyz,m_µ = get_bf(integrals,µ)
            shell_ν,atom_ν_xyz,m_ν = get_bf(integrals,ν)
            for i in 1:length(integrals.gridpoints)
                gw=integrals.gridweights[i]

                if shell_µ.l == shell_ν.l
                    #println("X1. shell_µ.l: $shell_µ.l and shell_ν.l:$shell_ν.l")
                    if integrals.mlvalues[µ] == integrals.mlvalues[ν]
                        #println("integrals.mlvalues[µ]: ", integrals.mlvalues[µ])
                        #println("integrals.mlvalues[ν]: ", integrals.mlvalues[ν])
                        #Using pre-calculated values
                        ϕ_µ = integrals.BFvalues[i,µ]
                        ϕ_ν = integrals.BFvalues[i,ν]
                        V[µ,ν] += 4/3 * CxLDA * ⍴[i]^(1/3) * ϕ_µ * ϕ_ν * gw
                    else
                        #println("x2. else")
                        # e.g. px and py : nothing
                    end
                else
                    #println("x3 else")
                    ϕ_µ = integrals.BFvalues[i,µ]
                    ϕ_ν = integrals.BFvalues[i,ν]
                    V[µ,ν] += 4/3 * CxLDA * ⍴[i]^(1/3) * ϕ_µ * ϕ_ν * gw
                end
            end
        end
    end
    return V
end



"""
create_density: Create density from P (density matrix) and basis functions 
for provided gridpoints
"""
function create_density(P,integrals)
    gridpoints=integrals.gridpoints
    ⍴=zeros(length(gridpoints)) # array of rho-values for each gridpoint
    dim=integrals.bset.nbas

    for gp_i in eachindex(gridpoints)
        for µ in 1:dim
            for ν in 1:dim
                ϕ_µ = integrals.BFvalues[gp_i,µ]
                ϕ_ν = integrals.BFvalues[gp_i,ν]
                ⍴[gp_i] += P[µ,ν] * ϕ_µ * ϕ_ν
            end
        end
    end
    #Making sure gridpoints dens values are never negative
    for gp_i in eachindex(gridpoints)
        if ⍴[gp_i] < 0 
            println("negative:")
            println("⍴[gp_i]:", ⍴[gp_i])
            ⍴[gp_i] = max(0.0,⍴[gp_i])
            println("⍴[gp_i]:", ⍴[gp_i])
            exit()
        end
        ⍴[gp_i] = max(0.0,⍴[gp_i])
    end
   return ⍴
end

function create_∇⍴(P,integrals)
    gridpoints=integrals.gridpoints
    ∇⍴=zeros(length(gridpoints)) # array of drho-values for each gridpoint
    dim=integrals.bset.nbas
    for gp_i in eachindex(gridpoints)
        for µ in 1:dim
            for ν in 1:dim
                ϕ_µ = integrals.BFvalues[gp_i,µ]
                ϕ_ν = integrals.BFvalues[gp_i,ν]
                ∇ϕ_µ = integrals.∇BFvalues[xyz,µ,gp_i]
                ∇ϕ_ν = integrals.∇BFvalues[xyz,ν,gp_i]
                ∇⍴ += P[µ,ν]*(∇ϕ_µ*ϕ_ν + ϕ_µ*∇ϕ_ν)
            end
        end
    end
    return ∇⍴
end



"""
Fxc formation in each Fock matrix formation: RKS
"""
function Fxc_form(integrals,P)
    println("Inside Fxc_form")
    dim=integrals.bset.nbas

    #Initializing Fxc matrix and Exc energy
    Fxc=zeros(dim,dim)
    Exc=0.0

    #Create density on grid from current P (needed by all rungs)
    ⍴ = create_density(P,integrals)

    #Integrate ⍴
    N_⍴ = integrate_density(⍴,integrals.gridweights)
    println("Integrated no. electrons: $N_⍴")

    integrals.N_⍺=N_⍴/2
    integrals.N_β=N_⍴/2


    #Manual functional specification
    #Activated by functional keyword="manual" and manual_func
    if integrals.DFT.manual == true
        println("DFT manual functional chosen with manual_func:$(integrals.DFT.manual_func)")
        if integrals.DFT.manual_func=="LDA"
            println("manual_func: LDA chosen")
            CxLDA=-0.7385587663820223
            #Exchange-correlation energy
            Exc = E_x_LDA(integrals,CxLDA,⍴)
            #Exchange-correlation potential
            Fxc = V_x_LDA(integrals,CxLDA,⍴)
        else
            println("Unknown manual_func")
            exit()
        end
        return Fxc,Exc
    end

    if integrals.DFT.functional_rung == 0
        println("Functional of rung 0 chosen, i.e. no XC. Pure Hartree-hell.")
        return Fxc,Exc
    elseif integrals.DFT.functional_rung == 1
        println("LibXC functional of rung 1 chosen")
        #Evaluate using input rho
        vrho_grid,e_xc_grid = Libxc.evaluate(integrals.DFT.libxc_functional, rho=⍴)
        #Integrate values on the grid
        Exc = integrate_XC_energy(e_xc_grid,⍴,integrals)
        Fxc = integrate_XC_potential(vrho_grid,⍴,integrals)
        return Fxc,Exc
    elseif integrals.DFT.functional_rung == 2
        println("LibXC functional of rung 2 chosen. Now we need sigma")
        println("Not ready. Exiting")
        exit()
        #sigma = functionalderivative(P,integrals)

        #Evaluate using input rho and sigma
        vrho_grid,e_xc_grid,vsigma_grid = Libxc.evaluate(integrals.DFT.libxc_functional, rho=⍴, sigma=sigma)
        #Integrate values on the grid
        Exc = integrate_XC_energy(e_xc_grid,⍴,integrals)
        Fxc = integrate_XC_potential(vrho_grid,⍴,integrals)
        #TODO: sigma
        println("here in gga")
        exit()
        return Fxc,Exc
    elseif integrals.DFT.functional_rung == 3
        println("LibXC functional 3 rung chosen. Now we need tau")
        println("Not ready. Exiting")
        exit()
        #Evaluate using input rho, sigma and tau
        vrho_grid,e_xc_grid,vsigma_grid,vtau_grid = Libxc.evaluate(integrals.DFT.libxc_functional, rho=⍴, sigma=sigma, tau=tau)

        #TODO: sigma and tau
        println("here in meta-gga")
        println("Not ready. Exiting")
        exit()
        return Fxc,Exc
    elseif integrals.DFT.functional_rung == 4
        println("LibXC functional 4 rung chosen")
        #Need to differentiate between hybrid-GGA and hybrid-meta-GGA
        #TODO: Use different functional rung 4a vs. 4b or 4.0 vs. 4.5 ??
        println("not yet hybrid")
        exit()
        return Fxc,Exc
    else
        println("functional_rung not supported")
        exit()
    end
end



"""
UKS Fxc formation: UKS

TODO: P_i and P_j refer to either P_⍺ or P_β
"""
function Fxc_form_UKS(integrals,P_i,P_j,manifold)
    dim=integrals.bset.nbas

    #Initializing Fxc matrix and Exc energy
    Fxc=zeros(dim,dim)
    Exc=0.0
    
    println("P_i:", P_i)
    println("P_j:", P_j)
    #Create density on grid from current P_⍺ and P_β (needed by all rungs)
    ⍴_i = create_density(P_i,integrals)
    ⍴_j = create_density(P_j,integrals)

    #Full electron density
    ⍴=⍴_i+⍴_j
    #Combined density (for libxc)
    ⍴_comb=zeros(length(⍴_i)*2)
    ⍴_comb[1:2:end] = ⍴_i
    ⍴_comb[2:2:end] = ⍴_j

    #Integrate ⍴ spin manifold and store
    N_⍴_i = integrate_density(⍴_i,integrals.gridweights)
    #N_⍴_j = integrate_density(⍴_j,integrals.gridweights)   No need to do here
    if manifold == 1
        integrals.N_⍺=N_⍴_i
        println("Integrated no. ⍺ electrons: $N_⍴_i")
    elseif manifold == 2
        integrals.N_β=N_⍴_i
        println("Integrated no. β electrons: $N_⍴_i")
    end

    #Manual functional specification
    #Activated by functional keyword="manual" and manual_func
    if integrals.DFT.manual == true
        println("DFT manual functional chosen with manual_func:$(integrals.DFT.manual_func)")
        if integrals.DFT.manual_func=="LDA"
            #println("manual_func: LDA chosen")
            #TODO: This needs to be tested

            CxLDA = -0.7385587663820223 # -(3/4)*(3/pi)^(1/3)
            println("CxLDA: $CxLDA")
            #Exchange-correlation energy
            Exc = E_x_LDA(integrals,CxLDA,⍴)
            #Exchange-correlation potential
            Fxc = V_x_LDA(integrals,CxLDA,⍴)
        #Not sure where this alternative definition comes from and how much it is used
        elseif integrals.DFT.manual_func=="LSDA"
            CxLSDA = -0.9305257363491 #-(3/2)*(3/(4*pi))^(1/3)
            println("CxLSDA: $CxLSDA")
            #Exchange-correlation energy
            Exc = E_x_LDA(integrals,CxLSDA,⍴)
            #Exchange-correlation potential
            Fxc = V_x_LDA(integrals,CxLSDA,⍴)
        else
            println("Unknown manual_func")
            exit()
        end
        return Fxc,Exc
    end

    #If not manual then NoXC or LibXC
    if integrals.DFT.functional_rung == 0
        #println("Functional of rung 0 chosen, i.e. no XC. Pure Hartree-hell.")
        return Fxc,Exc
    elseif integrals.DFT.functional_rung == 1
        println("LibXC functional of rung 1 chosen")
        #Evaluate using input rho combined
        vrho_grid,e_xc_grid = Libxc.evaluate(integrals.DFT.libxc_functional, rho=⍴_comb)
        numgridpoints=length(integrals.gridpoints)
        vrho_grid_i=vrho_grid[1:numgridpoints]

        #Integrate values on the grid
        Exc = integrate_XC_energy(e_xc_grid,⍴,integrals)
        Fxc = integrate_XC_potential(vrho_grid_i,⍴_i,integrals)

        return Fxc,Exc
    elseif integrals.DFT.functional_rung == 2
        println("LibXC functional of rung 2 chosen. Now we need sigma")
        println("Not ready. Exiting")
        exit()
        #sigma = functionalderivative(P,integrals)
        #Evaluate using input rho and sigma
        #vrho_grid,e_xc_grid,vsigma_grid = Libxc.evaluate(integrals.DFT.libxc_functional, rho=⍴, sigma=sigma)
        #Integrate values on the grid
        #Exc = integrate_XC_energy(e_xc_grid,⍴,integrals)
        #Fxc = integrate_XC_potential(vrho_grid,⍴,integrals)
        #TODO: sigma
        #println("here in gga")
        exit()
        return Fxc,Exc
    elseif integrals.DFT.functional_rung == 3
        println("LibXC functional 3 rung chosen. Now we need tau")
        println("Not ready. Exiting")
        exit()
        #Evaluate using input rho, sigma and tau
        vrho_grid,e_xc_grid,vsigma_grid,vtau_grid = Libxc.evaluate(integrals.DFT.libxc_functional, rho=⍴, sigma=sigma, tau=tau)

        #TODO: sigma and tau
        println("here in meta-gga")
        println("Not ready. Exiting")
        exit()
        return Fxc,Exc
    elseif integrals.DFT.functional_rung == 4
        println("LibXC functional 4 rung chosen")
        #Need to differentiate between hybrid-GGA and hybrid-meta-GGA
        #TODO: Use different functional rung 4a vs. 4b or 4.0 vs. 4.5 ??
        println("not yet hybrid")
        exit()
        return Fxc,Exc
    else
        println("functional_rung not supported")
        exit()
    end
end

"""
Choose functional definition based on functional string and libxc_keyword string.
Returns Libxc functional and functional rung index unless manual or NoXC option
"""
function choose_functional(functional_string,libxc_keyword;openshell=false)
    if openshell == true
        n_spin=2 #UKS: n_spin=2
    else
        n_spin=1 #RKS: n_spin=1, 
    end
    
    #TODO: Manual definition here also???
    #NoXC is Hartree only, i.e. rung 0
    if functional_string =="NoXC"
        return nothing, 0, false
    #Manual option: could add more manual options here
    elseif functional_string=="manual"
        return nothing, -1, true
    #Choosing functional based on Jotunn-defined functional keyword or libxc manual keyword
    elseif functional_string == "LDA-X" || functional_string == "LDA-X(libxc)"
        println("Using LDA-X(libxc)")
        #Creating functional
        libxc_functional = Libxc.Functional(:lda_x, n_spin=n_spin)
    #Defining functional 
    elseif functional_string == "PBE-X"     
        libxc_functional = Libxc.Functional(:gga_x_pbe, n_spin=n_spin)
    #No functional keyword string. Defining functional by libxc name.
    elseif functional_string == "none"
        println("Using manual libxc option (libxc_keyword: $libxc_keyword)")
        #Manual libxc option chosen. Converting to Symbol
        libxc_symbol=Symbol(libxc_keyword)
        println("libxc_symbol:", libxc_symbol)
        libxc_functional = Libxc.Functional(libxc_symbol, n_spin=n_spin)
        println("libxc_functional: $libxc_functional")
    else
        println("choose_functional else option")
        exit()
    end

    functional_rung = determine_rung(libxc_functional)
    return libxc_functional,functional_rung, false
end

"""
determine_rung: Determine functional rung of libxc-created functional object
functional: Libxc.jl Functional type
"""
function determine_rung(functional)

    #Checking if VV10 correlation is present and exit
    if Libxc.is_vv10(functional) == true
        println("VV10 correlation component not yet implemented")
        exit()
    end

    #Check for rung before continuing
    if Libxc.is_lda(functional) == true
        return 1
    elseif Libxc.is_gga(functional) == true
        if Libxc.is_hybrid(functional) == true
            println("libxc_functional.exx_coefficient:", libxc_functional.exx_coefficient)
            return 4
            exit()
            if Libxc.global_hybrid(functional) == true
                println("global hybrid")
                #Use 4.0 ?
                exit()
                #TODO
            elseif Libxc.is_range_separated(functional) == true
                println("Range-separated hybrids not yet implemented in Jotunn!")
                exit()
                #TODO
            end
        else
            return 2
        end
    elseif Libxc.is_mgga(functional) == true
        #Note: For now assuming tau-dependent functionals only (no Laplacian)
        if Libxc.is_hybrid(functional) == true
            return 4
            println("libxc_functional.exx_coefficient:", libxc_functional.exx_coefficient)
            exit()
            if Libxc.global_hybrid(functional) == true
                println("global hybrid")
            elseif Libxc.is_range_separated(functional) ==true
                println("Range-separated hybrids not yet implemented in Jotunn!")
                exit()
            end
        else
            return 3
        end
    else
        println("Unknown functional-rung")
        exit()
    end
end