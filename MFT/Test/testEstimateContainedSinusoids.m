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
        FLAGS   % Object that holds information about the possibility of round-off corruption or absurd values
    end
    
    methods (TestClassSetup)
        function getFunctionHandles(testCase)
            global flags
            [handles,~]=EstimateContainedSinusoids('-test',{});
            testCase.localFunctions = handles;
            flags = Flags();
        end
    end
    
    methods (TestMethodSetup)
        function setTsDefaults(testCase)
            nSinusoidsMax = 20;
            Fs = 600/2000;   % Sample rate   
            duration = 2000;
            for i=1:nSinusoidsMax  
                Amps(i) = 20 + i;
                Phases(i) = (117 + 10*i)*pi/180; 
            end
            
            testCase.TS = ArtificialTS(...
                'Name', 'testEstimateContainedSinusoids', ...
                'Description',' Time Series for Testing EstimateContainedSinusoids', ...
                'T0', 0, ...
                'Extent', duration, ...
                'nSamples', uint32(duration * Fs), ...
                'Freqs', [1/300], ...
                'Amps', Amps, ...
                'Phases', Phases, ...
                'NoiseUniformLow', 0, ...
                'NoiseUniformHi', 0, ...
                'NoiseGaussMean', 0, ...
                'NoiseGaussSD', 0, ...'=
                'RngState', -1);                                 
        end
    end
    
    
    
    methods (Test)
        function regressionTests (testCase)
            tuneNuEdgeWidth (testCase)
%           testComputeTheAAMatrix (testCase); 
%           testContainsSinusoidMultiple (testCase)            
        end
    
    end
    
    
    
    methods (Access = private)
        
        function tuneNuEdgeWidth (testCase)
            % test to help find the best kNuEgdeWidth
           kNuEgdeWidth = 0.00005;     % <- make sure this is the same as declared in EstimateContainedSinusoids

            % 'Choosing the edgewidth
            %----------------------
            %- As a nu approaches an edge (0 or N/2):
            %   1. The sine average at that nu approaches zero.
            %   2. The row of the "aa" matrix for the sine average in that nu approaches all zeroes.
            %   3. The col of the "aa" matrix for the sine part for that nu approaches all zeroes.
            %- In this case, can get absurdly high values of the sine part at that freq, as it makes little difference to
            %  anything in the equations described by aa. The amplitudes of the sinusoids thus found appear to "blow up",
            %  becoming absurdly high.
            %- The equations are correct, but it seems the suprod values that populate aa are wrong -- simply
            %  because they are very small (typically E-10 and E-20) and thus dominated by roundoff error. So the values in
            %  the aa matrix become random, so the computed values of the sine part at this nu becomes absurdly wrong.
            %- The solution is to declare the zone close to the edge "the edge zone", where the edge zone covers all
            %  the nu's where the roundoff error in suprod values is significant. In the edge zone, simply treat the
            %  nu as at the edge, and thus omit the equation from the "aa" matrix.
            %- Set the edge zone by trial and error, guided by testEstimateContainedSinusoids.tuneNuEdgeWidth.
            %  Obvously want the edge width to be as small as possible to get more accurate results, but big enough
            %  to avoid the effect of roundoff in computing suprod values.
            %- The value of aa that is the smallest, and thus most vulnerable to suprod roundoff, is at the
            %  intersection of the row for the sine average of nu and the col for the sine part for nu, the ss suprod
            %  value. At the nu near the edge, cc ~ 1, cs = sc ~ 10^-d, ss ~ 10^-2d. The roundoff error of a suprod
            %  value using double arithemtic (which is good for 15 decimal places) for time series of 1000 data points or so
            %  creates a noise floor around 10^-12 or so. So d is about 6, i.e. set edgewidth so ss(edgewidth) is around
            %  10^-12.
            %- Set the tol in TestEstContainedSinusoidOnce very low (maybe 10^-7) then watch the failure messages
            %  for the degenerate cases as adjust kNuEdgeWidth. To stop the misses blowing up, need kNuEdgeWidth >= 10^-6,
            %  or 10^-5 to be surer. Make it too large however, and get other errors ramping up due to artificial movement of
            %  frequencies to the edge. Best compromise seems to be about 0.00002.
            %- Note that the nu values for EstimateContainedSinusoids must be sufficiently distinct that the rows and
            %  cols of aa are not close to being linearly dependent. So consolidate the frequencies first.
            %-----------------------------------------------------------------------------------------------
            testCase.TS.Freqs = [0];
            testCase.TS.Name = 'DC';
            testEstimateContainedSinusoidOnce (testCase);
           testCase.TS.Freqs = [kNuEgdeWidth + 1e-6];
            testCase.TS.Name = sprintf('DC + %f',kNuEgdeWidth + 1e-6);
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.Freqs = [kNuEgdeWidth - 1e-6];
            testCase.TS.Name = sprintf('DC + %f',kNuEgdeWidth - 1e-6);            
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.Freqs = [0.5 * testCase.TS.nSamples/testCase.TS.Extent];
            testCase.TS.Name = 'Nyquist';
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth - 1e-6) / testCase.TS.Extent];
            testCase.TS.Name = sprintf('NyQuist - %f',kNuEgdeWidth - 1e-6);                        
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth + 1e-6) / testCase.TS.Extent];
            testCase.TS.Name = sprintf('NyQuist - %f',kNuEgdeWidth + 1e-6);                        
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.nSamples = uint32(7 * 7 * 11);  % odd
            testCase.TS.Freqs = [0.5 * testCase.TS.nSamples/testCase.TS.Extent];
            testCase.TS.Name = 'Nyquist';
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth - 1e-6) / testCase.TS.Extent];
            testCase.TS.Name = sprintf('NyQuist - %f',kNuEgdeWidth - 1e-6);                        
            testEstimateContainedSinusoidOnce (testCase);
            testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth + 1e-6) / testCase.TS.Extent];
            testCase.TS.Name = sprintf('NyQuist - %f',kNuEgdeWidth + 1e-6);                        
            testEstimateContainedSinusoidOnce (testCase);
          
        end
        
        
        function testComputeTheAAMatrix (testCase)
            handle=testCase.localFunctions(2);
            ComputeTheAAMatrix = handle{1};
            n = 16;
            nu = [2 3.5 4.1];
            isEdgeNu = [false false false];
            aa = ComputeTheAAMatrix(n,nu,isEdgeNu);
        end
        
        function testContainsSinusoidMultiple (testCase)
            kNuEgdeWidth = 0.00005;
            
            testCase.TS.Freqs = [1/300];
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
           global flags
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
           title(testCase.TS.Name);
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
            disp(sprintf('flags = %s',flags.flagString))
        end 
        
    end
    
end
            
            