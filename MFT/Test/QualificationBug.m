classdef QualificationBug < matlab.unittest.TestCase
    
    methods (Test)        
        function testQualificationBug (testCase)
            kTol = 0.003;            
            act = 5.022;
            exp = 5.026;
            testCase.verifyEqual(act,exp,'AbsTol',kTol)
        end
    end
end
            