@testset verbose=true "Issues" begin
    # This is just to collect the hodgepodge of random tests inspired by issues that have
    # crept up.  The numbers refer to github issues.

    @testset "Issue #15" begin
        t = collect(LinRange(-10, 10, 201))
        @test t .* imz == imz .* t
        @test t * imz == imz * t
        @test t * imz == imz .* t
        @test t .* imz == imz * t
    end

    @testset "Issue #70" begin
        p = Rotor(3,-1,2, 1.2)
        vr = [Rotor(-2,1,1,3), Rotor(2,0,-1,3)]
        @test typeof(p * vr) === typeof(vr)
        @test typeof(vr * p) === typeof(vr)
    end
end
