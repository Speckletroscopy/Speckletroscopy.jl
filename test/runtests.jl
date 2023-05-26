using Speckletroscopy
using Test

@testset "Speckletroscopy.jl" begin
    @testset "LightSource.jl" begin
        @test σTemp(1.0,1.0)^2 ≈ 9.178e-14
    end
end;
