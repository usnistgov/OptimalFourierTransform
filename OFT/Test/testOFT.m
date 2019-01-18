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
            [handles,~,~]=OFT_fn('-test',{},{},{}, {});
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
            testCase.TS.NoiseUniformLow = 0;
            testCase.TS.NoiseUniformHi = 0;
            testCase.TS.NoiseGaussMean = 0;
            testCase.TS.NoiseGaussSD = 0;      
        end
    end
    
    methods (Test)
        function regressionTests (testCase)
%            testWaitBar (testCase)
%            testOFTCalc_1 (testCase)
%            testMultipleOFT (testCase)
%             testNyquist (testCase)
%             testNearDC (testCase)
%             testLab (testCase)
            testNoise (testCase)
        end
    end
%--------------------------------------------------------------------------                
    methods (Access = private)
    % these private methods are called by the regression tests
    
        function testWaitBar (testCase)
            handle=testCase.localFunctions(1);
            waitBarOFT = handle{1};
            
            waitBarOFT(0,'Wait Bar Test',true);
            pause(1)
            waitBarOFT(.33,'1/3 done',true)
            pause(1)
            waitBarOFT(.66,'2/3 done',true)
            pause(1)
            waitBarOFT(1,'Passed Wait Bar Test!',true)
            pause(1)
            waitBarOFT(2,'you should not be seeing this',true)                      
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
    
    
    methods (Access=private)
    % these private methods are groups of artificial time series that are
    % called by the regression tests
    
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
            [actFreqs, actOFT, fracErr] = OFT_fn(testCase.TS.Ts, t, .000001, true, 5, true);  % recon, plot progress plot p
            %[actFreqs, actOFT, fracErr] = OFT_fn(testCase.TS.Ts, t, .000001, false, 5, true);  % no recon, plot progress plot p
            %[actFreqs, actOFT, fracErr] = OFT_fn(testCase.TS.Ts, t, .000001, true, 5, false); % recon, no plot progress plot p
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
            
            disp(mat2str(actFreqs));
            disp(mat2str(actOFT));
            
            pause
            close all
        end
      
    
        function testOFTCalc_1 (testCase)
        % This loads and analyses the same time series data used in the 
        % original EXCEL "climate.xlsm by David Evans.
                       
            load('tsData.mat');  % load the test data
            t = double(tsData.t);
            t0 = t(1);
            deltaT = mean(t(2:end)-t(1:end-1));
            for tau=0:length(t)-1
                regT(tau+1) = t0+deltaT*tau;
            end
            ts = double(tsData.g);
                          
            [actFreqs, actOFT, actFracErr] = OFT_fn(ts, regT, 0.000001, true, 6, true); 
            [act_TS] = testCase.SynthesizeOFT(actFreqs, t, actOFT);
            
            close all
            figure(2)
            subplot(2,1,1)
            hold on
            plot (t,ts,'b')
            plot(t,act_TS,'r')
            hold off
            subplot (2,1,2)
            ts_diff = act_TS' - ts;
            plot (t, ts_diff,'g')
            pause
            close all
%             
% % TODO:  set up some expected values to verify                        
% %           actAmp = abs(actOFT)                        
% %           actAngle = radtodeg(angle(actOFT));
% %           for i = 1:length(actAngle)
% %               if actAngle(i) < 0; actAngle(i) = 360 + actAngle(i); end
% %           end
%                      
        end
         
        function testNyquist (testCase)
            kNuEgdeWidth = 0.00005;
            
            testCase.TS.Name = 'Test Nyquist (600)';
            testCase.TS.nSamples = uint32(600);
            testCase.TS.Freqs = [0.5 * testCase.TS.nSamples/testCase.TS.Extent];
            testOftOnce (testCase);
            
            testCase.TS.Name = 'Test Nyquist - kNuEgdeWidth - 1e-7 (600)';            
            testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth - 1e-7) / testCase.TS.Extent];
            testOftOnce (testCase);

            testCase.TS.Name = 'Test Nyquist - kNuEgdeWidth + 1e-7 (600)';            
            testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth + 1e-7) / testCase.TS.Extent];
            testOftOnce (testCase);

            % odd samples near nyquist.  
            testCase.TS.Name = 'Test Nyquist (539)';            
            testCase.TS.nSamples = uint32(7 * 7 * 11);  % odd
            testCase.TS.Freqs = [0.5 * testCase.TS.nSamples/testCase.TS.Extent];
            testOftOnce (testCase);                            
            
        end
        
        function testNearDC (testCase)
            kNuEgdeWidth = 0.00005;            
            testCase.TS.nSamples = uint32(600);    
            
            testCase.TS.Name = 'Test DC  (600)';            
            testCase.TS.Freqs = [0];
            testOftOnce (testCase);            
            
            testCase.TS.Name = 'Test DC + 1e-7 (600)';            
            testCase.TS.Freqs = [kNuEgdeWidth + 1e-7];
            testOftOnce (testCase);

            testCase.TS.Name = 'Test DC- 1e-7 (600)';            
            testCase.TS.Freqs = [kNuEgdeWidth - 1e-7];
            testOftOnce (testCase);
            
            testCase.TS.nSamples = uint32(7 * 7 * 11);  % odd  
            testCase.TS.Name = 'Test DC  (539)';            
            testCase.TS.Freqs = [0];
            testOftOnce (testCase);                        
            
            testCase.TS.Name = 'Test DC + 1e-7 (539)';            
            testCase.TS.Freqs = [kNuEgdeWidth + 1e-7];
            testOftOnce (testCase);

            testCase.TS.Name = 'Test DC - 1e-7 (539)';            
            testCase.TS.Freqs = [kNuEgdeWidth - 1e-7];
            testOftOnce (testCase);
            
                        
        end

    
       function testMultipleOFT (testCase)
            kNuEgdeWidth = 0.00005;  
            testCase.TS.nSamples = uint32(600);                
            
            testCase.TS.Name = 'Test #1-1 (1/300)';
            testCase.TS.Freqs = [1/300];
            testOftOnce (testCase);
            
             % an odd number of samples
            testCase.TS.nSamples = uint32(7 * 7 * 11);  % odd

            testCase.TS.Name = 'Test #1-2 (1/300) (odd)';            
            testCase.TS.Freqs = [1/300];
            testOftOnce (testCase);                       
             
            % Tests on 2 Sinusoids
            testCase.TS.Name = 'Test #2-1 (1/300(67), 1/49(27)) (odd)';            
            testCase.TS.Freqs = [1/300, 1/49];
            testCase.TS.Phases = [67 * pi/180, 27 * pi/180];
            testOftOnce (testCase);

            testCase.TS.Name = 'Test #2-2 (1/300(117), 1/49(27)) (odd)';            
            testCase.TS.Phases = [117 * pi/180, 27 * pi/180];
            testOftOnce (testCase);

            testCase.TS.Name = 'Test #2-3 (1/300(217), 1/49(327)) (odd)';            
            testCase.TS.Phases = [217 * pi/180, 327 * pi/180];
            testOftOnce (testCase);

            testCase.TS.Name = 'Test #2-4 (1/300(297), 1/49(127)) (odd)';            
            testCase.TS.Phases = [297 * pi/180, 127 * pi/180];
            testOftOnce (testCase);
            

            % Random
            testCase.TS.Name = 'Test #2-5 (1/1000(297), 1/12(127)) (odd)';            
            testCase.TS.Freqs = [1/1000, 1/12];
            testOftOnce (testCase);

            testCase.TS.Name = 'Test #2-6 (1/1000(297), 1/2000(127)) (odd)';            
            testCase.TS.Freqs = [1/1000, 1/2000];
            testOftOnce (testCase);

            testCase.TS.Name = 'Test #2-7 (1/1000(297), 1/2000(127)) (odd)';            
            testCase.TS.Freqs = [1/20, 1/10];
            testOftOnce (testCase);
            
           % Degeneerates
           testCase.TS.Name = 'Test #2-8 (0, 1/120(127)) (odd)';            
           testCase.TS.Extent = 2000;
           testCase.TS.Freqs = [0, 1/120];
           testOftOnce (testCase);

           testCase.TS.Name = 'Test #2-9 (1/1000(297), 0) (odd)';            
           testCase.TS.Freqs = [1/1000, 0];
           testOftOnce (testCase);

           testCase.TS.Name = 'Test #2-10 (Nyquist, 1/23) (odd)';            
           testCase.TS.Freqs = [0.5 * testCase.TS.nSamples/testCase.TS.Extent, 1/23];
           testOftOnce (testCase);

           % Just outside the Edge Zone
           testCase.TS.Name = 'Test #2-11 (kNuEgdeWidth + 1e-9, Nyquist) (odd)';            
           testCase.TS.Freqs = [(kNuEgdeWidth + 1e-9)/testCase.TS.Extent, 0.5 * testCase.TS.nSamples/testCase.TS.Extent];
           testOftOnce (testCase);
           
           % Just inside the Edge Zone
           testCase.TS.Name = 'Test #2-12 (kNuEgdeWidth - 1e-9, Nyquist) (odd)';            
           testCase.TS.Freqs = [(kNuEgdeWidth - 1e-9)/testCase.TS.Extent, 0.5 * testCase.TS.nSamples/testCase.TS.Extent];
           testOftOnce (testCase);
           
           % Just outside the Edge Zone
           testCase.TS.Name = 'Test #2-13 (Nyquist- kNuEgdeWidth - 1e-9, 0) (odd)';            
           testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth - 1e-9)/testCase.TS.Extent, 0];
           testOftOnce (testCase);
           
           % Just inside the Edge Zone
           testCase.TS.Name = 'Test #2-14 (Nyquist- kNuEgdeWidth + 1e-9, 0) (odd)';            
           testCase.TS.Freqs = [(0.5 * testCase.TS.nSamples - kNuEgdeWidth + 1e-9)/testCase.TS.Extent, 0];
           testOftOnce (testCase);
           
           testCase.TS.Name = 'Test #2-15 (Nyquist), 1/23 (odd)';            
           testCase.TS.Freqs = [0.5 * testCase.TS.nSamples / testCase.TS.Extent, 1/23];
           testOftOnce (testCase);

           testCase.TS.Name = 'Test #2-16 (Nyquist, 0) (odd)';            
           testCase.TS.Freqs = [0.5 * testCase.TS.nSamples / testCase.TS.Extent, 0];
           testOftOnce (testCase);

           testCase.TS.Name = 'Test #2-17 (0, Nyquist) (odd)';                       
           testCase.TS.Freqs = [0, 0.5 * testCase.TS.nSamples / testCase.TS.Extent];
           testOftOnce (testCase);
           
           % Tests on 3 sinusoids
           setTsDefaults(testCase)
           
           testCase.TS.Name = 'Test #3-1 (1/400, 1/402, 1/404) (odd)';                       
           testCase.TS.Freqs = [1/400, 1/402, 1/404];
           testOftOnce (testCase);

           testCase.TS.Name = 'Test #3-2 (1/1000, 1/12, 1/129) (odd)';                       
           testCase.TS.Freqs = [1/1000, 1/12, 1/129];
           testOftOnce (testCase);

           testCase.TS.Name = 'Test #3-3 (0, 1/120, 1/200) (odd)';                       
           testCase.TS.Freqs = [0, 1/120, 1/200];
           testOftOnce (testCase);

           testCase.TS.Name = 'Test #3-4 (1/120, 0, 1/130) (odd)';                       
           testCase.TS.Freqs = [1/120, 0, 1/130];
           testOftOnce (testCase);
           
           testCase.TS.Name = 'Test #3-5 (1/120, 1/130, 0) (odd)';                       
           testCase.TS.Freqs = [1/120, 1/130, 0];
           testOftOnce (testCase);
           
           % Tests on 4 sinusoids
           
           testCase.TS.Name = 'Test #4-1 (1/1000, 1/12, 1/129, 1/239) (odd)';                       
           testCase.TS.Freqs = [1/1000, 1/12, 1/129, 1/239];
           testOftOnce (testCase);

           testCase.TS.Name = 'Test #4-2 (1/380, 1/403, 1/407, 1/416) (odd)';                       
           testCase.TS.Freqs = [1/380, 1/403, 1/407, 1/416];
           testOftOnce (testCase);
           
           % Tests on 5 sinusoids
           testCase.TS.Name = 'Test #5-1 (1/1000, 1/12, 1/129, 1/239, 1/339) (odd)';                       
           testCase.TS.Freqs = [1/1000, 1/12, 1/129, 1/239, 1/339];
           testOftOnce (testCase);

           testCase.TS.Name = 'Test #5-2 (1/380, 1/402, 1/410, 1/435, 1/440) (odd)';                       
           testCase.TS.Freqs = [1/380, 1/402, 1/410, 1/435, 1/440];
           testOftOnce (testCase);

           testCase.TS.Name = 'Test #5-3 (8;1/606(45), 20;1/405(0);, 20;1/201(17), 9;1/102(0), 20;1/49(217)) (odd)';                       
           testCase.TS.Freqs = [1/606, 1/405, 1/201, 1/102, 1/49];
           testCase.TS.Amps = [8,20,20,9,20];
           testCase.TS.Phases = [45,0,17,0,217]*pi/180;
           testOftOnce (testCase);
                                              
       end    
    
       function testLab (testCase)
           % set up the artificial time series
           testCase.TS.T0 = 0;
           testCase.TS.Extent = 2000;
           testCase.TS.nSamples = uint32(300);
           testCase.TS.Freqs = [0, ...
                                0.00123762376237624, ...
                                0.00165016501650165, ...
                                0.00247524752475248, ...
                                0.0033003300330033, ...
                                0.0048780487804878, ...
                                0.0099009900990099, ...
                                0.0204081632653061, ...
                                0.0434782608695652, ...
                                0.0714285714285714...
                                ];
           testCase.TS.Amps = [0, 0, 0, 0, 0, 9, 0, 0, 0, 0];
           testCase.TS.Phases = [0, 90, 45, 230, 0, 20, 15, 0, 215, 0] * pi/180;
           
           testCase.TS.Name = 'One Sinusoid';
           testOftOnce (testCase);
           
           testCase.TS.Name = 'Two Sinusoids';           
           testCase.TS.Amps(5) = 8;
           testOftOnce (testCase);
           
           testCase.TS.Name = 'Three Sinusoids'; 
           testCase.TS.Amps(5) = 0;
           testCase.TS.Amps(3) = 11;
           testCase.TS.Freqs(6) = 0.00495049504950495;
           testCase.TS.Amps(9) = 13;
           testOftOnce (testCase);
           
           
           testCase.TS.Name = 'Four Sinusoids (2 sec)'; 
           testCase.TS.Phases(3) = 235*pi/180;
           testCase.TS.Amps(5) = 8;
           testOftOnce (testCase);
           
           testCase.TS.Name = 'Five Sinusoids (2 sec)'; 
           testCase.TS.Freqs(3) = 0.00165837479270315;
           testCase.TS.Phases(3) = 45*pi/180;
           testCase.TS.Phases(4) = 235*pi/180;
           testCase.TS.Freqs(5) = 0.00332225913621262;
           testCase.TS.Amps(8) = 12;
           testOftOnce (testCase);
           
           testCase.TS.Name = 'Six Sinusoids (9 sec)'; 
           testCase.TS.Freqs(3) = 0.00165016501650165;
           testCase.TS.Phases(4) = 230*pi/180;           
           testCase.TS.Freqs(5) = 0.0033003300330033;
           testCase.TS.Amps(4) = 10;
           testCase.TS.Amps(5) = 0;
           testCase.TS.Amps(8) = 0;
           testCase.TS.Amps(10) = 7;
           testOftOnce (testCase);
           
           testCase.TS.Name = 'Seven Sinusoids (10 sec)'; 
           testCase.TS.Freqs(3) = 0.00166389351081531;
           testCase.TS.Phases(4) = 235*pi/180;
           testCase.TS.Freqs(5) = 0.00332225913621262;
           testCase.TS.Amps(8) = 12;
           testOftOnce (testCase);
           
           testCase.TS.Name = 'Eight Sinusoids (49 sec)'; 
           testCase.TS.Freqs(3) = 0.00165016501650165;
           testCase.TS.Phases(4) = 230*pi/180;
           testCase.TS.Freqs(5) = 0.0033003300330033;
           testCase.TS.Amps(5) = 9;
           testOftOnce (testCase);
           
           testCase.TS.Name = 'Nine Sinusoids (39 sec)'; 
           testCase.TS.Amps(2) = 8;
           testCase.TS.Freqs(2) = 0.0133333333333333;
           testCase.TS.Phases(2) = 300*pi/180;
           testCase.TS.Phases(4) = 0;
           testCase.TS.Amps(8) = 4;
           testCase.TS.Phases(8) = 340*pi/180;
           testCase.TS.Amps(9) = 8;
           testCase.TS.Phases(10) = 40*pi/180;
           testOftOnce (testCase);
          
           testCase.TS.Name = 'Ten Sinusoids (33 sec)'; 
           testCase.TS.Freqs(1) = 0.0065359477124183;
           testCase.TS.Amps(1) = 9;
           testCase.TS.Phases(1) = 150*pi/180;  
           testOftOnce (testCase);
           
           testCase.TS.Name = 'Two Sinusoids close but resolvable'; 
           testCase.TS.Freqs(1) = 0;
           testCase.TS.Freqs(2) = 0.00123762376237624;
           testCase.TS.Freqs(6) = 0.00346020761245675;
           testCase.TS.Amps = [0, 0, 0, 0, 8, 9, 0, 0, 0, 0];
           testcase.TS.Phases = [0,90,45,230,0,20,15,0,215,0] * pi/180;
           testOftOnce (testCase);
           
           testCase.TS.Name = 'Two Sinusoids too close'; 
           testCase.TS.Freqs(6) = 0.00343642611683849;
           testOftOnce (testCase);
           
           testCase.TS.Name = 'Two Sinusoids consolidated'; 
           testCase.TS.Freqs(5) = 0.00332225913621262;
           testCase.TS.Freqs(6) = 0.00334448160535117;
           testOftOnce (testCase);
                                                                                                                      
       end
       
       function testNoise (testCase)
           
           %set up a 60 fps 5 second system
           Fs = 60;   % Sample rate
           duration = 5;
           testCase.TS = ArtificialTS;
           testCase.TS.Name = 'No-Noise Test';
           testCase.TS.T0 = 0;
           testCase.TS.Extent = duration;
           testCase.TS.nSamples = uint32(duration * Fs);
           testCase.TS.Freqs = [5];
           testCase.TS.Amps = [70];
           testCase.TS.Phases = [0];
           
           testCase.TS.NoiseUniformLow = 0;
           testCase.TS.NoiseUniformHi = 0;
           testCase.TS.NoiseGaussMean = 0;
           testCase.TS.NoiseGaussSD = 0;
         
           testOftOnce (testCase);
           
           %---
           
           testCase.TS.Name = 'Noise Test';
            
           testCase.TS.NoiseUniformLow = -.01;
           testCase.TS.NoiseUniformHi = .01;
           testCase.TS.NoiseGaussMean = 0;
           testCase.TS.NoiseGaussSD = 0;
          
           testOftOnce (testCase);
                                           
       end
    end
end