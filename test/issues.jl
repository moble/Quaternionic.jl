# Tests inspired by GitHub issues; the numbers refer to the issues.

@testitem "Issues: #15" tags=[:unit, :fast] begin
    t = collect(LinRange(-10, 10, 201))
    @test t .* imz == imz .* t
    @test t * imz == imz * t
    @test t * imz == imz .* t
    @test t .* imz == imz * t
end

@testitem "Issues: #70" tags=[:unit, :fast] begin
    p = Rotor(3,-1,2, 1.2)
    vr = [Rotor(-2,1,1,3), Rotor(2,0,-1,3)]
    @test typeof(p * vr) === typeof(vr)
    @test typeof(vr * p) === typeof(vr)
end
