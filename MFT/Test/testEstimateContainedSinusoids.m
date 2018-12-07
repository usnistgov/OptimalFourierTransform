% Class-based unit testing for EstimateContainedSinusoids
%
%   run these tests with two command-line commands:
%   - testCase = testEstimateContainedSinusoids;
%   - res = run(testCase);
%
classdef testEstimateContainedSinusoids < matlab.unittest.TestCase
    properties
        localFunctions
        TS
    end
    
    methods (TestClassSetup)
        function getFunctionHandles(testCase)
            [handles,~]=EstimateContainedSinusoids('-test',{});
            testCase.localFunctions = handles;
        end
    end
    
    methods (TestMethodSetup)
        function setTsDefaults(testCase)
            nSinusoidsMax = 20;
            Fs = 600/2000;   % Sample rate   
            duration = 2000;
            testCase.TS = ArtificialTS;
            testCase.TS.Name = 'testTS';
            testCase.TS.Description = ' Time Series for Testing';
            testCase.TS.T0 = 0;
            testCase.TS.Extent = duration;
            testCase.TS.nSamples = uint32(duration * Fs);
            testCase.TS.Freqs = [1/300];
            testCase.TS.Amps = [1];
            testCase.TS.Phases = [0];
            for i=1:nSinusoidsMax  
                testCase.TS.Amps(i) = 20 + i;
                testCase.TS.Phases(i) = (117 + 10*i)*pi/180; 
            end
        end
    end
    
    methods (Test)
%         function testClassifyRegFreqsAsEdgeOrNormal (testCase)
%             handle=testCase.localFunctions(1);
%             ClassifyRegFreqsAsEdgeOrNormal = handle{1};
%         end
        
        function testComputeTheAAMatrix (testCase)
            handle=testCase.localFunctions(2);
            ComputeTheAAMatrix = handle{1};
            n = 16;
            nu = [2 3.5 4.1];
            isEdgeNu = [false false false];
            aa = ComputeTheAAMatrix(n,nu,isEdgeNu);
        end
        
        function testFailedCase (testCase)
           testCase.TS.Extent = 2000;
           testCase.TS.Freqs = [1/400, 1/402, 1/404];
           testEstimateContainedSinusoidOnce (testCase);            
        end
        
        function testContainsSinusoidMultiple (testCase)
            kNuEgdeWidth = 0.00005;
            
            testCase.TS.Freqs = [1/300];
            testEstimateContainedSinusoidOnce (testCase);

            % Degenerates
            testCase.TS.Freqs = [0];
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.Freqs = [0.5 * testCase.TS.nSamples/testCase.TS.Extent];
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.Freqs = [kNuEgdeWidth + 1e-7];
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.Freqs = [kNuEgdeWidth - 1e-7];
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.nSamples = uint32(7 * 7 * 11);  % odd
            testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth - 1e-7) / testCase.TS.Extent];
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth + 1e-7) / testCase.TS.Extent];
            testEstimateContainedSinusoidOnce (testCase);
            
            % Tests on 2 Sinusoids
            testCase.TS.Freqs = [1/300, 1/49];
            testCase.TS.Phases = [67 * pi/180, 27 * pi/180];
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.Phases = [117 * pi/180, 27 * pi/180];
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.Phases = [217 * pi/180, 327 * pi/180];
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.Phases = [297 * pi/180, 127 * pi/180];
            testEstimateContainedSinusoidOnce (testCase);
            
            % Random
            testCase.TS.Freqs = [1/1000, 1/12];
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.Freqs = [1/1000, 1/2000];
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.Freqs = [1/20, 1/10];
            testEstimateContainedSinusoidOnce (testCase);
            
           % Degeneerates
           testCase.TS.Extent = 2000;
           testCase.TS.Freqs = [0, 1/120];
           testEstimateContainedSinusoidOnce (testCase);
           testCase.TS.Freqs = [1/1000, 0];
           testEstimateContainedSinusoidOnce (testCase);
           testCase.TS.Freqs = [0.5 * testCase.TS.nSamples/testCase.TS.Extent, 1/23];
           testEstimateContainedSinusoidOnce (testCase);

           % Just outside the Edge Zone
           testCase.TS.Freqs = [(kNuEgdeWidth + 1e-9)/testCase.TS.Extent, 0.5 * testCase.TS.nSamples/testCase.TS.Extent];
           testEstimateContainedSinusoidOnce (testCase);
           
           % Just inside the Edge Zone
           testCase.TS.Freqs = [(kNuEgdeWidth - 1e-9)/testCase.TS.Extent, 0.5 * testCase.TS.nSamples/testCase.TS.Extent];
           testEstimateContainedSinusoidOnce (testCase);
           
           % Just outside the Edge Zone
           testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth - 1e-9)/testCase.TS.Extent, 0];
           testEstimateContainedSinusoidOnce (testCase);
           
           % Just inside the Edge Zone
           testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth + 1e-9)/testCase.TS.Extent, 0];
           testEstimateContainedSinusoidOnce (testCase);
           
           testCase.TS.Freqs = [0.5 * testCase.TS.nSamples / testCase.TS.Extent, 1/23];
           testEstimateContainedSinusoidOnce (testCase);
           testCase.TS.Freqs = [0.5 * testCase.TS.nSamples / testCase.TS.Extent, 0];
           testEstimateContainedSinusoidOnce (testCase);
           testCase.TS.Freqs = [0, 0.5 * testCase.TS.nSamples / testCase.TS.Extent];
           testEstimateContainedSinusoidOnce (testCase);
           
           % Tests on 3 sinusoids
           setTsDefaults(testCase)
           testCase.TS.Freqs = [1/400, 1/402, 1/404];
           testEstimateContainedSinusoidOnce (testCase);
           testCase.TS.Freqs = [1/1000, 1/12, 1/129];
           testEstimateContainedSinusoidOnce (testCase);
           testCase.TS.Freqs = [0, 1/120, 1/200];
           testEstimateContainedSinusoidOnce (testCase);
           testCase.TS.Freqs = [1/120, 0, 1/130];
           testEstimateContainedSinusoidOnce (testCase);
           testCase.TS.Freqs = [1/120, 1/130, 0];
           testEstimateContainedSinusoidOnce (testCase);
           
           % Tests on 4 sinusoids
           
           testCase.TS.Freqs = [1/1000, 1/12, 1/129, 1/239];
           testEstimateContainedSinusoidOnce (testCase);
           testCase.TS.Freqs = [1/380, 1/403, 1/407, 1/416];
           testEstimateContainedSinusoidOnce (testCase);
           
           % Tests on 5 sinusoids
           testCase.TS.Freqs = [1/1000, 1/12, 1/129, 1/239, 1/339];
           testEstimateContainedSinusoidOnce (testCase);
           testCase.TS.Freqs = [1/380, 1/402, 1/410, 1/435, 1/440];
           testEstimateContainedSinusoidOnce (testCase);
           testCase.TS.Freqs = [1/606, 1/405, 1/201, 1/102, 1/49];
           testCase.TS.Amps = [8,20,20,9,20];
           testCase.TS.Phases = [45,0,17,0,217]*pi/180;
           testEstimateContainedSinusoidOnce (testCase);
                      
       end
        
    end
    
    methods (Access=private)
        function testEstimateContainedSinusoidOnce (testCase)
           % This runs one iteration of test ArtificialTS parameters for
           % each test run have been set up by TestMethodSetup (automatically 
           % set by matlab before each test) and by
           % TestMultipleContainedSinusoids.
           kTol = 0.003;
           kNuEgdeWidth = 0.00005;     % <- make sure this is the same as declared in EstimateContainedSinusoids
           
           % for each test, show a figure of the time series
           testCase.TS = testCase.TS.makeTime;
           testCase.TS = testCase.TS.makeTS;
           Figure = figure;
           plot(testCase.TS.time, testCase.TS.Ts)
           title(testCase.TS.Freqs);
           pause;
           close(Figure);
           
           
           for i = 1:length(testCase.TS.Freqs)
                nu(i) = testCase.TS.Freqs(i) * testCase.TS.Extent;
            end
            [actCosPart,actSinPart] = EstimateContainedSinusoids(testCase.TS.Ts,nu);
            
            % create sine and cosine parts from the ArtificialTS parameters
            for i =1:length(testCase.TS.Freqs)
                expCosPart(i) = testCase.TS.Amps(i) * cos (testCase.TS.Phases(i));
                expSinPart(i) = testCase.TS.Amps(i) * sin (testCase.TS.Phases(i));
                if nu(i) <= kNuEgdeWidth || nu(i) >= length(testCase.TS.Ts)/2 - kNuEgdeWidth
                    expSinPart(i) = 0;
                end
            end
            testCase.verifyEqual(actCosPart,expCosPart,'AbsTol',kTol)
            testCase.verifyEqual(actSinPart,expSinPart,'AbsTol',kTol)
        end 
        
    end
    
end
            
            