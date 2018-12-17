% Class-based unit testing for OFT
%
%   run these tests with two command-line commands:
%   - testCase = testOFT;
%   - res = run(testCase);
%
classdef testOFT < matlab.unittest.TestCase
    properties
        localFunctions
        TS
    end
    
    methods (TestClassSetup)
        function getFunctionHandles(testCase)
            [handles,~,~]=OFT('-test',{},{},{}, {});
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
        
%         function testWaitBar (testCase)
%             handle=testCase.localFunctions(1);
%             waitBarOFT = handle{1};
%             
%             waitBarOFT(0,'Wait Bar Test',true);
%             pause(1)
%             waitBarOFT(.33,'are you waiting?',true)
%             pause(1)
%             waitBarOFT(.66,'what are you waiting for?',true)
%             pause(1)
%             waitBarOFT(1,'Passed Wait Bar Test!',true)
%             pause(1)
%             waitBarOFT(2,'you should not be seeing this',true)                      
%         end
%         
        function testOFTCalc (testCase)
            global Figure
            
            
            % load the test data and print it
            load('tsData.mat');
            t = double(tsData.t);
            t0 = t(1);
            deltaT = mean(t(2:end)-t(1:end-1));
            for tau=0:length(t)-1
                regT(tau+1) = t0+deltaT*tau;
            end
            ts = double(tsData.g);
            
             Figure = figure;
%             plot(regT,ts,'b')
%             hold on
%               
            [actFreqs, actOFT, actFracErr] = OFT(ts, regT, 0.000001, true, 6, false); 
            
        end
    end
    
    methods (Access=private)
        function testOftOnce (testCase)
           testCase.TS = testCase.TS.makeTime;
            testCase.TS = testCase.TS.makeTS;
            Figure = figure;
            plot(testCase.TS.time,testCase.TS.Ts)
            title(testCase.TS.Freqs);
            pause;
            close(Figure);
            
            [~, ~, ~] = OFT(testCase.TS.Ts,testCase.TS.time, .1, true, 5, true);
        end
    end
            
    
    
end