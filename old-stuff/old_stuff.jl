
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