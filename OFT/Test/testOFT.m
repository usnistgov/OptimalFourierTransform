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
            testCase.TS.Name = 'Test #';
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
        
%         function testOFTCalc_1 (testCase)
%         % This loads and analyses the same time series data used in the 
%         % original EXCEL "climate.xlsm by David Evans.
%                        
%             load('tsData.mat');  % load the test data
%             t = double(tsData.t);
%             t0 = t(1);
%             deltaT = mean(t(2:end)-t(1:end-1));
%             for tau=0:length(t)-1
%                 regT(tau+1) = t0+deltaT*tau;
%             end
%             ts = double(tsData.g);
%                           
%             [actFreqs, actOFT, actFracErr] = OFT(ts, regT, 0.000001, true, 6, true); 
%             [act_TS] = testCase.SynthesizeOFT(actFreqs, t, actOFT);
%             
%             close all
%             figure(2)
%             subplot(2,1,1)
%             hold on
%             plot (t,ts,'b')
%             plot(t,act_TS,'r')
%             hold off
%             subplot (2,1,2)
%             ts_diff = act_TS' - ts;
%             plot (t, ts_diff,'g')
%             pause
%             close all
%             
% % TODO:  set up some expected values to verify                        
% %           actAmp = abs(actOFT)                        
% %           actAngle = radtodeg(angle(actOFT));
% %           for i = 1:length(actAngle)
% %               if actAngle(i) < 0; actAngle(i) = 360 + actAngle(i); end
% %           end
%                      
%         end
         
        function testMultipleOFT (testCase)
            kNuEgdeWidth = 0.00005;   
            
%             testCase.TS.Name = 'Test #1-1 (1/300)';
%             testCase.TS.Freqs = [1/300];
%             testOftOnce (testCase);
%             
%             testCase.TS.Name = 'Test #1-2 (0)';                     
%             testCase.TS.Freqs = [0];
%             testOftOnce (testCase);            
            
            % Close to Nyquist frequency does not work so need to revisit and
            % find out how to handle it or throw an error.
            testCase.TS.Name = 'Test #1-3 (Nyquist)';            
            testCase.TS.Freqs = [0.5 * testCase.TS.nSamples/testCase.TS.Extent];
            testOftOnce (testCase);           

%             testCase.TS.Name = 'Test #1-4 (Nyquist - kNuEgdeWidth - 1e-7)';            
%             testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth - 1e-7) / testCase.TS.Extent];
%             testOftOnce (testCase);
% 
%             testCase.TS.Name = 'Test #1-5 (Nyquist - kNuEgdeWidth + 1e-7)';            
%             testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth + 1e-7) / testCase.TS.Extent];
%             testOftOnce (testCase);
% 
%             testCase.TS.Name = 'Test #1-6 (Nyquist + 1e-7)';            
%             testCase.TS.Freqs = [kNuEgdeWidth + 1e-7];
%             testOftOnce (testCase);
% 
%             testCase.TS.Name = 'Test #1-7 (Nyquist - 1e-7)';            
%             testCase.TS.Freqs = [kNuEgdeWidth - 1e-7];
%             testOftOnce (testCase);
% 
%             % an odd number of samples
%             testCase.TS.nSamples = uint32(7 * 7 * 11);  % odd
% 
%             testCase.TS.Name = 'Test #1-8 (1/300) (odd)';            
%             testCase.TS.Freqs = [1/300];
%             testOftOnce (testCase);                       
%              
%             % odd samples near nyquist.  Still has a hard time but looks
%             % different than with even samples.
%             testCase.TS.Name = 'Test #1-9 (Nyquist) (odd)';            
%             testCase.TS.Freqs = [0.5 * testCase.TS.nSamples/testCase.TS.Extent];
%             testOftOnce (testCase);                       
% 
%             testCase.TS.Name = 'Test #1-10 (Nyquist - kNuEgdeWidth - 1e-7) (odd)';            
%             testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth - 1e-7) / testCase.TS.Extent];
%             testOftOnce (testCase);
% 
%             testCase.TS.Name = 'Test #1-11 (Nyquist - kNuEgdeWidth + 1e-7) (odd)';            
%             testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth + 1e-7) / testCase.TS.Extent];
%             testOftOnce (testCase);        
%             
%             % Tests on 2 Sinusoids
%             testCase.TS.Name = 'Test #2-1 (1/300(67), 1/49(27)) (odd)';            
%             testCase.TS.Freqs = [1/300, 1/49];
%             testCase.TS.Phases = [67 * pi/180, 27 * pi/180];
%             testOftOnce (testCase);
% 
%             testCase.TS.Name = 'Test #2-2 (1/300(117), 1/49(27)) (odd)';            
%             testCase.TS.Phases = [117 * pi/180, 27 * pi/180];
%             testOftOnce (testCase);
% 
%             testCase.TS.Name = 'Test #2-3 (1/300(217), 1/49(327)) (odd)';            
%             testCase.TS.Phases = [217 * pi/180, 327 * pi/180];
%             testOftOnce (testCase);
% 
%             testCase.TS.Name = 'Test #2-4 (1/300(297), 1/49(127)) (odd)';            
%             testCase.TS.Phases = [297 * pi/180, 127 * pi/180];
%             testOftOnce (testCase);
%             
% 
%             % Random
%             testCase.TS.Name = 'Test #2-5 (1/1000(297), 1/12(127)) (odd)';            
%             testCase.TS.Freqs = [1/1000, 1/12];
%             testOftOnce (testCase);
% 
%             testCase.TS.Name = 'Test #2-6 (1/1000(297), 1/2000(127)) (odd)';            
%             testCase.TS.Freqs = [1/1000, 1/2000];
%             testOftOnce (testCase);
% 
%             testCase.TS.Name = 'Test #2-7 (1/1000(297), 1/2000(127)) (odd)';            
%             testCase.TS.Freqs = [1/20, 1/10];
%             testOftOnce (testCase);
%             
%            % Degeneerates
%            testCase.TS.Name = 'Test #2-8 (0, 1/120(127)) (odd)';            
%            testCase.TS.Extent = 2000;
%            testCase.TS.Freqs = [0, 1/120];
%            testOftOnce (testCase);
% 
%            testCase.TS.Name = 'Test #2-9 (1/1000(297), 0) (odd)';            
%            testCase.TS.Freqs = [1/1000, 0];
%            testOftOnce (testCase);
% 
%            testCase.TS.Name = 'Test #2-10 (Nyquist, 1/23) (odd)';            
%            testCase.TS.Freqs = [0.5 * testCase.TS.nSamples/testCase.TS.Extent, 1/23];
%            testOftOnce (testCase);
% 
%            % Just outside the Edge Zone
%            testCase.TS.Name = 'Test #2-11 (kNuEgdeWidth + 1e-9, Nyquist) (odd)';            
%            testCase.TS.Freqs = [(kNuEgdeWidth + 1e-9)/testCase.TS.Extent, 0.5 * testCase.TS.nSamples/testCase.TS.Extent];
%            testOftOnce (testCase);
%            
%            % Just inside the Edge Zone
%            testCase.TS.Name = 'Test #2-12 (kNuEgdeWidth - 1e-9, Nyquist) (odd)';            
%            testCase.TS.Freqs = [(kNuEgdeWidth - 1e-9)/testCase.TS.Extent, 0.5 * testCase.TS.nSamples/testCase.TS.Extent];
%            testOftOnce (testCase);
%            
%            % Just outside the Edge Zone
%            testCase.TS.Name = 'Test #2-13 (Nyquist- kNuEgdeWidth - 1e-9, 0) (odd)';            
%            testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth - 1e-9)/testCase.TS.Extent, 0];
%            testOftOnce (testCase);
%            
%            % Just inside the Edge Zone
%            testCase.TS.Name = 'Test #2-14 (Nyquist- kNuEgdeWidth + 1e-9, 0) (odd)';            
%            testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth + 1e-9)/testCase.TS.Extent, 0];
%            testOftOnce (testCase);
%            
%            testCase.TS.Name = 'Test #2-15 (Nyquist), 1/23 (odd)';            
%            testCase.TS.Freqs = [0.5 * testCase.TS.nSamples / testCase.TS.Extent, 1/23];
%            testOftOnce (testCase);
% 
%            testCase.TS.Name = 'Test #2-16 (Nyquist, 0) (odd)';            
%            testCase.TS.Freqs = [0.5 * testCase.TS.nSamples / testCase.TS.Extent, 0];
%            testOftOnce (testCase);
% 
%            testCase.TS.Name = 'Test #2-17 (0, Nyquist) (odd)';                       
%            testCase.TS.Freqs = [0, 0.5 * testCase.TS.nSamples / testCase.TS.Extent];
%            testOftOnce (testCase);
%            
%            % Tests on 3 sinusoids
%            setTsDefaults(testCase)
%            
%            testCase.TS.Name = 'Test #3-1 (1/400, 1/402, 1/404) (odd)';                       
%            testCase.TS.Freqs = [1/400, 1/402, 1/404];
%            testOftOnce (testCase);
% 
%            testCase.TS.Name = 'Test #3-2 (1/1000, 1/12, 1/129) (odd)';                       
%            testCase.TS.Freqs = [1/1000, 1/12, 1/129];
%            testOftOnce (testCase);
% 
%            testCase.TS.Name = 'Test #3-3 (0, 1/120, 1/200) (odd)';                       
%            testCase.TS.Freqs = [0, 1/120, 1/200];
%            testOftOnce (testCase);
% 
%            testCase.TS.Name = 'Test #3-4 (1/120, 0, 1/130) (odd)';                       
%            testCase.TS.Freqs = [1/120, 0, 1/130];
%            testOftOnce (testCase);
%            
%            testCase.TS.Name = 'Test #3-5 (1/120, 1/130, 0) (odd)';                       
%            testCase.TS.Freqs = [1/120, 1/130, 0];
%            testOftOnce (testCase);
%            
%            % Tests on 4 sinusoids
%            
%            testCase.TS.Name = 'Test #4-1 (1/1000, 1/12, 1/129, 1/239) (odd)';                       
%            testCase.TS.Freqs = [1/1000, 1/12, 1/129, 1/239];
%            testOftOnce (testCase);
% 
%            testCase.TS.Name = 'Test #4-2 (1/380, 1/403, 1/407, 1/416) (odd)';                       
%            testCase.TS.Freqs = [1/380, 1/403, 1/407, 1/416];
%            testOftOnce (testCase);
%            
%            % Tests on 5 sinusoids
%            testCase.TS.Name = 'Test #5-1 (1/1000, 1/12, 1/129, 1/239, 1/339) (odd)';                       
%            testCase.TS.Freqs = [1/1000, 1/12, 1/129, 1/239, 1/339];
%            testOftOnce (testCase);
% 
%            testCase.TS.Name = 'Test #5-2 (1/380, 1/402, 1/410, 1/435, 1/440) (odd)';                       
%            testCase.TS.Freqs = [1/380, 1/402, 1/410, 1/435, 1/440];
%            testOftOnce (testCase);
% 
%            testCase.TS.Name = 'Test #5-3 (8;1/606(45), 20;1/405(0);, 20;1/201(17), 9;1/102(0), 20;1/49(217)) (odd)';                       
%            testCase.TS.Freqs = [1/606, 1/405, 1/201, 1/102, 1/49];
%            testCase.TS.Amps = [8,20,20,9,20];
%            testCase.TS.Phases = [45,0,17,0,217]*pi/180;
%            testOftOnce (testCase);
                                              
        end
        
        
        
    end
    
    methods (Access=private)

    end
            
    methods (Access = private)
        
        function testOftOnce (testCase)
           % This runs one iteration of test ArtificialTS parameters for
           % each test run have been set up by TestMethodSetup (automatically 
           % set by matlab before each test) and by testMultipleOFT
           % TestMultipleContainedSinusoids.            
            disp(testCase.TS.Name)
            testCase.TS = testCase.TS.makeTime;
            testCase.TS = testCase.TS.makeTS;
            
%             Figure = figure;
%             plot(testCase.TS.time,testCase.TS.Ts)
%             title(testCase.TS.Freqs);
%             pause;
%             close(Figure);
            
             t = testCase.TS.time;
            [actFreqs, actOFT, fracErr] = OFT(testCase.TS.Ts, t, .1, true, 5, true);  % recon, plot progress plot p
%            [actFreqs, actOFT, fracErr] = OFT(testCase.TS.Ts, t, .1, false, 5, true);  % no recon, plot progress plot p
%            [actFreqs, actOFT, fracErr] = OFT(testCase.TS.Ts, t, .1, true, 5, false); % recon, no plot progress plot p
            [act_TS] = testCase.SynthesizeOFT(actFreqs, t, actOFT);
            orig_TS = testCase.TS.Ts;

            
            close all
            figure(2)
            s(1)=subplot(2,1,1);
            hold on
            plot (t,orig_TS,'b')
            plot(t,act_TS,'r')
            hold off

            s(2)=subplot (2,1,2);
            ts_diff = act_TS - orig_TS;

            plot (t, ts_diff,'g')
            title(s(1),testCase.TS.Name);
            title(s(2),'Residual');
            pause
            close all
        end
        
        function [ts_OFT] = SynthesizeOFT(testCase, Freqs, t, actOFT)
           n = length(t);
           extent = t(end)+mean(diff(t))-t(1);
           twoPiExtentON = 2* pi * extent / n;
           ts_OFT = zeros(1,length(t));
           
           for i = 1:length(Freqs)
              mag = abs(actOFT(i));
              ang = angle(actOFT(i));
              %ang = mod(angle(actOFT(i)),pi);
              twoPiTime = zeros(1,n);
              for tau = n-1:-1:0
                  twoPiTime(tau+1) = twoPiExtentON * tau;
              end
              radians = twoPiTime * Freqs(i);
              ts_OFT = ts_OFT + mag * cos(radians-ang);              
           end
        end
    end
    
end