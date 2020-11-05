using Test
using LinearAlgebra

@testset "LinearSystem Class Struct" begin
    #Test The Constructor
    A = [1 0; 0 2]
    B = [1;1]
    C = [1 0]
    eta_w = 1
    eta_v = 1

    #Create System 1
    system1 = LinearSystem(A,B,C,eta_w,eta_v)
    @test system1.A == A
end;

@testset "check() tests" begin
    #Test1: Bad A Matrix
    A = [1 0; 0 2; 0 3]
    B = [1;1]
    C = [1 0]
    eta_w = 1
    eta_v = 1

    #Create System 1
    system1 = LinearSystem(A,B,C,eta_w,eta_v)
    @test !check(system1)

    #Test2
    A = [1 0; 0 3]
    system2 = LinearSystem(A,B,C,eta_w,eta_v)
    @test check(system2)

    #Test3
    C = [1 0 0]
    system3 = LinearSystem(A,B,C,eta_w,eta_v)
    @test !check(system3)
end

@testset "n_x() Tests" begin
     #Test1: Bad A Matrix
     A = [1 0; 0 2; 0 3]
     B = [1;1]
     C = [1 0]
     eta_w = 1
     eta_v = 1
 
     #Create System 1
     system1 = LinearSystem(A,B,C,eta_w,eta_v)
     @test_throws ArgumentError n_x(system1)

     #Test2
     A = [1 0; 0 3]
     system2 = LinearSystem(A,B,C,eta_w,eta_v)
     @test n_x(system2) == 2
end

@testset "n_y() Tests" begin
     #Test1: Bad A Matrix
     A = [1 0; 0 2; 0 3]
     B = [1;1]
     C = [1 0]
     eta_w = 1
     eta_v = 1
 
     #Create System 1
     system1 = LinearSystem(A,B,C,eta_w,eta_v)
     @test_throws ArgumentError n_y(system1)

     #Test2
     A = [1 0; 0 3]
     system2 = LinearSystem(A,B,C,eta_w,eta_v)
     @test n_y(system2) == 1

    #Test3
    C = [1 0 0]
    system3 = LinearSystem(A,B,C,eta_w,eta_v)
    @test_throws ArgumentError n_y(system3)

    #Test4
    C = [1 0;2 2;3 1]
    system4 = LinearSystem(A,B,C,eta_w,eta_v)
    @test n_y(system4) == 3
end

@testset "find_ALAP_time_from_X0_to_bound Tests" begin
    #Test the error handling
    @test_throws MethodError find_ALAP_time_from_X0_to_bound(1)
end