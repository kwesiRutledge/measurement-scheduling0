"""
    constraint_generator_test.jl
    Description:
        
"""

using Test
using LinearAlgebra

@testset "GetPolytopeInclusionConstraints() Tests" begin
    # ================================ 
    # Test1: Normal Hx and Hy Matrices
     Hx1 = [1 0; 0 1; -1 0;0 -1]
     hx1 = [1;1;1;1]
     Hx2 = Hx1
     hx2 = [2;2;2;2]
 
     #Create System 1
     #cg0 = ConstraintGenerator()
     GetPolytopeInclusionConstraints(Hx1,hx1,Hx2,hx2)


     #Test2
     Hx1 = [1 0; 0 1; -1 0;0 -1]
     hx1 = [1;1;1;1]
     Hx2 = Hx1
     hx2 = [2;2;2;2]

     A = [1 0; 0 3]

end