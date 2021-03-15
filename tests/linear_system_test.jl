using Test
using LinearAlgebra
using Polyhedra

@testset "LinearSystem Class Struct" begin
    #Test The Constructor
    A = [1 0; 0 2]
    B = [1;1]
    C = [1 0]
    D = [1 0]
    d = 0.0
    # eta_w = polyhedron(MixedMatHRep(HalfSpace([1.0, 0.2], 1.0) ∩ HalfSpace([-1.0, -0.2], 1.0) ))
    Pw = polyhedron(hrep( [1.0 0.2 ; -1.0 -0.2] , [1.0,1.0] ))
    # eta_v = polyhedron(MixedMatHRep(HalfSpace([1.0], 1.0) ∩ HalfSpace([-1.0], 1.0) ))
    Pv = polyhedron(hrep( [HalfSpace([1.0], 1.0), HalfSpace([-1.0], 1.0)] ))
    Pu = Pv

    #Create System 1
    system1 = LinearSystem(A,B,C,D,d,Pw,Pv,Pu)
    @test system1.A == A
end;

@testset "check() tests" begin
    #Test1: Bad A Matrix
    A = [1 0; 0 2; 0 3]
    B = [1;1]
    C = [1 0]
    D = [1 0]
    d = 0.0

    Pw = polyhedron(hrep( [1.0 0.2 ; -1.0 -0.2] , [1.0,1.0] ))
    Pv = polyhedron(hrep( [HalfSpace([1.0], 1.0), HalfSpace([-1.0], 1.0)] ))
    Pu = Pv

    #Create System 1
    system1 = LinearSystem(A,B,C,D,d,Pw,Pv,Pu)
    @test !check(system1)

    #Test2
    A = [1 0; 0 3]
    system2 = LinearSystem(A,B,C,D,d,Pw,Pv,Pu)
    @test check(system2)

    #Test3
    C = [1 0 0]
    system3 = LinearSystem(A,B,C,D,d,Pw,Pv,Pu)
    @test !check(system3)
end

@testset "x_dim() Tests" begin
    #Test1: Bad A Matrix
    A = [1 0; 0 2; 0 3]
    B = [1;1]
    C = [1 0]
    D = [1 0]
    d = 0.0

    Pw = polyhedron(hrep( [1.0 0.2 ; -1.0 -0.2] , [1.0,1.0] ))
    Pv = polyhedron(hrep( [HalfSpace([1.0], 1.0), HalfSpace([-1.0], 1.0)] ))
    Pu = Pv
    
    #Create System 1
    system1 = LinearSystem(A,B,C,D,d,Pw,Pv,Pu)
    @test_throws ArgumentError x_dim(system1)

    #Test2
    A = [1 0; 0 3]
    system2 = LinearSystem(A,B,C,D,d,Pw,Pv,Pu)
    @test x_dim(system2) == 2

    # Test 3: One Dimensional System
    A = 7
    B = 5
    C = 3
    D = 1.0
    d = 0.0

    Pw = polyhedron(hrep( [HalfSpace([1.0], 1.0), HalfSpace([-1.0], 1.0)] ))
    Pv = Pw
    Pu = Pv

    system3 = LinearSystem(A,B,C,D,d,Pw,Pv,Pu)
    @test x_dim(system3) == 1
end

@testset "y_dim() Tests" begin
    #Test1: Bad A Matrix
    A = [1 0; 0 2; 0 3]
    B = [1;1]
    C = [1 0]
    D = [1 0]
    d = 0.0

    Pw = polyhedron(hrep( [1.0 0.2 ; -1.0 -0.2] , [1.0,1.0] ))
    Pv = polyhedron(hrep( [HalfSpace([1.0], 1.0), HalfSpace([-1.0], 1.0)] ))
    Pu = Pv
    
    #Create System 1
    system1 = LinearSystem(A,B,C,D,d,Pw,Pv,Pu)
    @test_throws ArgumentError y_dim(system1)

    #Test2
    A = [1 0; 0 3]
    system2 = LinearSystem(A,B,C,D,d,Pw,Pv,Pu)
    @test y_dim(system2) == 1

    #Test3
    C = [1 0 0]
    system3 = LinearSystem(A,B,C,D,d,Pw,Pv,Pu)
    @test_throws ArgumentError y_dim(system3)

    #Test4
    C = [1 0;2 2;3 1]
    Pv = polyhedron(hrep([HalfSpace([1.0,1.5,2.0], 1.0) , HalfSpace([-1.0,-1.5,-2.0], 1.0)] ))
    system4 = LinearSystem(A,B,C,D,d,Pw,Pv,Pu)
    @test y_dim(system4) == 3
end

@testset "define_simple_eta_HPolytope() Tests" begin
    #Test 1: Simple 1d system 
    A = 3
    B = 7
    C = 1.0
    D = 1.0
    d = 0.0

    Pw = polyhedron(hrep([HalfSpace([1.0], 0.1) , HalfSpace([-1.0], 0.1)] ))
    Pv = polyhedron(hrep([HalfSpace([1.0], 0.5) , HalfSpace([-1.0], 0.5)] ))
    Pu = Pv

    #Create system1
    system1 = LinearSystem(A,B,C,D,d,Pw,Pv,Pu)

    eta_x0 = 0.2
    P_x0 = polyhedron(hrep( [HalfSpace([1.0],eta_x0),HalfSpace([-1.0],eta_x0)] ))
    T1 = 1
    H_eta1, h_eta1 = define_simple_eta_HPolytope(system1,P_x0,T1)
    @test H_eta1 == [1.0 0 0;-1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1 ]
    @test transpose(h_eta1) == transpose([0.1; 0.1 ; 0.5; 0.5; 0.2; 0.2])

end

@testset "find_est_error_bound_from_time_0_to_T Tests" begin
    #Test with simple scalar system
    A = 1.0
    B = 1
    C = 1
    eta_w = 0.1
    P_w = MixedMatHRep(HalfSpace([1.0], eta_w) ∩ HalfSpace([-1.0], eta_w) )
    P_v = MixedMatHRep(HalfSpace([1.0], 0.5) ∩ HalfSpace([-1.0], 0.5) )

    #Create system1
    system1 = LinearSystem(A,B,C,P_w,P_v)
    eta_x0 = 0.2
    T1 = 1

    out1 = find_est_error_bound_from_time_0_to_T(system1, T1 , eta_x0)
    feasibilityflag1 = out1[1]
    opt_val1 = out1[2]
    @test opt_val1 == eta_x0+eta_w 
end

@testset "evaluate_schedule_wrt Tests" begin
    #Test with simple scalar system
    A = 1.0
    B = 1
    C = 1
    eta_w = 0.1
    eta_v = 0.5

    #Create system1
    system1 = LinearSystem(A,B,C,eta_w,eta_v)
    eta_x0 = 0.2
    T1 = 1

end