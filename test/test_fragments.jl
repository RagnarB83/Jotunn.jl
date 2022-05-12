#Some basic tests on fragment creation


#Testing various formatted strings
H2 = create_fragment(coords_string="""
    H 0.0 0.0 0.0
    H 0.0 0.0 0.74
    """, charge=0, mult=1)
H2_2 = create_fragment(coords_string="""H 0.0 0.0 0.0
H 0.0 0.0 0.74
""", charge=0, mult=1)
H2_3 = create_fragment(coords_string="""
H 0.0 0.0 0.0
H 0.0 0.0 0.74""", charge=0, mult=1)
H2_4 = create_fragment(coords_string="""H 0.0 0 0.0
H 0 0 0.74""", charge=0, mult=1)

@testset "fragments tests" begin
    @test H2.numatoms == 2
    @test length(H2.elems) == 2
    @test size(H2.coords) == (2,3)

    @test H2_2.numatoms == 2
    @test H2_3.numatoms == 2
    @test H2_4.numatoms == 2

    @test size(H2_2.coords) == (2,3)
    @test size(H2_3.coords) == (2,3)
    @test size(H2_4.coords) == (2,3)
end