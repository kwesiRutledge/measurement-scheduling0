using Test
using LinearAlgebra

@testset "compute_Sx0() Function Tests" begin
    #Test The function
    A = [1 0; 0 2]
    @test compute_Sx0(A,3) == [1 0; 0 1; 1 0 ; 0 2; 1 0 ; 0 4; 1 0 ; 0 8]
end;

@testset "compute_Sw() Function Tests" begin
    A = [1 0; 0 2]
    #Test what happens for time horizon 1
    @test compute_Sw(A,1) == [0.0 0; 0 0; 1 0 ; 0 1]
    @test compute_Sw(A,2) == [0.0 0 0 0; 0 0 0 0; 1 0 0 0; 0 1 0 0; 1 0 1 0; 0 2 0 1 ]
end;

@testset "compute_C_M() Function Tests" begin
    C = [1 0]
    M = [0 1]
    T = 2

    #Test What happens for multiple time horizons
    @test compute_C_M( C , M , T) == [1.0 0 0 0 0 0; 0 0 1 0 0 0]
end