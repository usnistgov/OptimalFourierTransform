% Class-based unit testing for MFT.m
%
%   run these tests with two command-line commands:
%   - testCase = testMFT;
%   - res = run(testCase);
%
classdef testMFT < matlab.unittest.TestCase
     properties
        localFunctions
    end
    
    methods (TestClassSetup)
        function getFunctionHandles(testCase)
            [f,~]=MFT('-test',{},{},{});
            testCase.localFunctions = f;
        end
    end
    
    methods (Test)
        
        function testConsolidateFreqs (testCase)
            handle=testCase.localFunctions(1);
            ConsolidateFreqs = handle{1};
            kNuConsolidate = 0.1;
            nu_MFT = [1, 5, 5.05, 6, 16, 16.025, 25, 25];
            [act] = ConsolidateFreqs(nu_MFT,kNuConsolidate);
            exp = [1, 5.0250, 6, 16.0125, 25];
            testCase.verifyEqual(act,exp);
        end
        
        function testConsolidateFreqsByBracket (testCase)
            handle=testCase.localFunctions(2);
            ConsolidateFreqsByBracket = handle{1};
            kNuConsolidate = 0.1;
            nu_MFT = [1, 5, 5.05, 6, 16, 16.025, 25, 25, 30, 32, 32.075, 45, 59];
            bracket = [4,5,4];
            [act] = ConsolidateFreqsByBracket(nu_MFT,bracket,kNuConsolidate);
            exp = [1, 5.025, 6. 16.0125, 25, 30, 32.0375, 45, 59];
            testCase.verifyEqual(act,exp);
        end
    end
end
        
        