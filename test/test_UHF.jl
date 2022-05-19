
H2 = create_fragment(coords_string="""
    H 0.0 0.0 0.0
    H 0.0 0.0 0.74
    """, charge=0, mult=1)
    
NO = create_fragment(coords_string="""
N 0.0 0.0 0.0
O 0.0 0.0 1.15
""", charge=0, mult=2)

@testset "UHF/STO-3g on H2 (4c)" begin
    #UHF on a closed-shell system. Same result as RHF
    result= jHF(H2, "sto-3g"; HFtype="UHF", tei_type="4c")
    #RHF/STO-3G H2 (r=0.74 A) from ORCA 5.0.3
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "UHF/STO-3g on H2 (sparse4c)" begin
    #UHF on a closed-shell system. Same result as RHF
    result= jHF(H2, "sto-3g"; HFtype="UHF", tei_type="sparse4c")
    #RHF/STO-3G H2 (r=0.74 A) from ORCA 5.0.3
    ref_energy=-1.116759307204
    @test isapprox(result["energy"],ref_energy)
end

@testset "UHF/STO-3g on NO (4c)" begin
    #UHF on a doublet radical
    result= jHF(NO, "sto-3g"; HFtype="UHF", levelshift=true, maxiter=200, tei_type="4c")
    #UHF/STO-3G NO (r=0.1.15 A) from ORCA 5.0.3
    ref_energy=-127.530134374969
    @test isapprox(result["energy"],ref_energy)
end


@testset "UHF/STO-3g on NO (sparse4c)" begin
    #UHF on a doublet radical
    result= jHF(NO, "sto-3g"; HFtype="UHF", levelshift=true, maxiter=200, tei_type="sparse4c")
    #UHF/STO-3G NO (r=0.1.15 A) from ORCA 5.0.3
    ref_energy=-127.530134374969
    @test isapprox(result["energy"],ref_energy)
end
