% Class-based unit testing for SolveLinearEquations
%
%   run these tests with two command-line commands:
%   - testCase = testSolveLinearEquations;
%   - res = run(testCase);
%
classdef testSolveLinearEquations < matlab.unittest.TestCase
    properties
        localFunctions
    end
    
    methods (TestClassSetup)
        function getFunctionHandles(testCase)
             [handles]=SolveLinearEquations('-test',{});
            testCase.localFunctions = handles;
         end
    end
    
    methods (Test)
       
        function testLinearEquationSolver (testCase)
            kTol = 1e-10;
            
            aa = [2, -3; 1, 2];
            bb = [5, 13];
            [act] = SolveLinearEquations(aa, length(aa), bb);
            exp = [7,3];
            testCase.verifyEqual(act,exp,'AbsTol',kTol);
            
            aa = [1,-2,3;2,1,-1;3,-3,2];
            bb = [-5,19,8];
            [act] = SolveLinearEquations(aa, length(aa), bb);
            exp = [7, 3, -2];
            testCase.verifyEqual(act,exp,'AbsTol',kTol);
           
            aa = [1, -2, 3, 4; ...
                  2, 1, -1, -2; ...
                  3, -3, 2, 1; ...
                  -1, -2, 1, -4];
              
             bb =  [ 15, 9, 13, -35];
            [act] = SolveLinearEquations(aa, length(aa), bb);
            exp = [7, 3, -2, 5];
            testCase.verifyEqual(act,exp,'AbsTol',kTol);
            
        end
        
    end
end
        
    
    
    