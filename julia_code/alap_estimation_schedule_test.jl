@testset "norm_of_est_error_after Tests" begin
    #Test1: Bad A Matrix
    A = [1 0; 0 2; 0 3]
    B = [1;1]
    C = [1 0]
    eta_w = 1
    eta_v = 1

    #Create System 1
    system1 = LinearSystem(A,B,C,eta_w,eta_v)
    @test_throws ArgumentError norm_of_est_error_after( 2 , system1 , 0.1 )

    #Test 2: Scalar System
    A = 1
    B = 1
    C = 1
    eta_w = 0.1
    eta_v = 0.2

    system2 = LinearSystem(A,B,C,eta_w,eta_v) #Create system2

    eta_x0 = 0.3
    estErrorBound2, feasibFlag2 = norm_of_est_error_after( 2 , system2 , eta_x0 )
    @test estErrorBound2 == 0.5


end

@testset "find_time_from_X0_to_bound Tests" begin
    #Test the error handling
    @test_throws MethodError find_time_from_X0_to_bound(1,2,3)
end

@testset "alap_estimation_schedule_alg1 Tests" begin
    # Test with nice Time Horizon and MeasurementBudget values.
    T = 12
    B = 3

    # Create Dummy System

    schedule1 = alap_estimation_schedule_alg1(T,B) #This should be the values [3,6,9]
    @test [3,6,9] == schedule1
end