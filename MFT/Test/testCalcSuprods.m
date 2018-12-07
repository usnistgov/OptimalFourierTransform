% Class-based unit testing for CalcSuprods
%
%   run these tests with two command-line commands:
%   - testCase = testCalcSuprods;
%   - res = run(testCase);
%
classdef testCalcSuprods < matlab.unittest.TestCase
    properties
        localFunctions
    end
    
    methods (TestClassSetup)
        function getFunctionHandles(testCase)
            [cc,~,~,~]=CalcRegSuprods('-test',{},{});
            testCase.localFunctions = cc;
        end
    end
    
    methods (Test)
        
        function testFactorize(testCase)
            handle=testCase.localFunctions(1);
            Factorize = handle{1};
            [lastN_SU,extensionsSU,factorsSU]=Factorize(92);
            act = {lastN_SU,extensionsSU,factorsSU};
            exp = {92,92,[2 2 23]};
            testCase.verifyEqual(act,exp); 
        end
        
        function testSmallestPrimeFactorOf(testCase)
            handle=testCase.localFunctions(2);
            SmallestPrimeFactorOf = handle{1};
            act = SmallestPrimeFactorOf(27);
            exp = 3;
            testCase.verifyEqual(act,exp);
            act = SmallestPrimeFactorOf(29);
            exp = 29;
            testCase.verifyEqual(act,exp);
        end
        
        function testCalcRegSuprods(testcase)
            
            % integral arguments
            testOneSuprodComparison(testcase,3,4,20); 
            testOneSuprodComparison(testcase,3,0,19);
            testOneSuprodComparison(testcase,0,4,37); 
            testOneSuprodComparison(testcase,3,10,20); 
            testOneSuprodComparison(testcase,10,10,20);
            testOneSuprodComparison(testcase,4,4,20); 
            testOneSuprodComparison(testcase,4,4,19); 
            testOneSuprodComparison(testcase,0,0,20); 
            testOneSuprodComparison(testcase,0,0,20); 
            
            %Arguments outside [0,N/2] 
            testOneSuprodComparison(testcase,23,34,20); 
            testOneSuprodComparison(testcase,17,-3,20); 
            testOneSuprodComparison(testcase,0,10,20); 
            testOneSuprodComparison(testcase,0.001,10.001,20); 
            testOneSuprodComparison(testcase,123.01,147.2,19); 
            testOneSuprodComparison(testcase,-123.01,-147.2,20); 
           
            % non-integral arguments
            testOneSuprodComparison(testcase,3.141,4.252,20);            
            testOneSuprodComparison(testcase,0.141,4.252,20);
            testOneSuprodComparison(testcase,9.941,4.252,20);
            testOneSuprodComparison(testcase,3.141,4.252,19);
            testOneSuprodComparison(testcase,0.141,4.252,19);
            testOneSuprodComparison(testcase,9.941,4.252,19); 
            
            %N as factors
            testOneSuprodComparison(testcase,3.141,4.252,256); 
            testOneSuprodComparison(testcase,3.141,4.252,39); 
            testOneSuprodComparison(testcase,3.141,4.252,600); 
            testOneSuprodComparison(testcase,3.141,4.252,97);  % 97 is prime
            
            
        end
         
    end
    
    methods
        function testOneSuprodComparison(testCase, nu, mu, n)
            handle=testCase.localFunctions(3);
            CalcRegSuprods_Naive = handle{1};
            
            kTol = 0.0000000001;
            
            [cc, cs, sc, ss] = CalcRegSuprods(nu, mu, n);
            [ccRef, csRef, scRef, ssRef] = CalcRegSuprods_Naive(nu, mu, n);
            testCase.verifyEqual(cc,ccRef,'AbsTol',1^-10);
            testCase.verifyEqual(cs,csRef,'AbsTol',1^-10);
            testCase.verifyEqual(sc,scRef,'AbsTol',1^-10);
            testCase.verifyEqual(ss,ssRef,'AbsTol',1^-10);   
        end
    end
    
end
    
