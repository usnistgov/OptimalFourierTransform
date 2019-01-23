classdef testOFT_class < matlab.unittest.TestCase
% Class-based unit testing for class-based OFT
%
%   run these tests with two command-line commands:
%   - testCase = testOFT_class;
%   - res = run(testCase);
%
    
    properties
        TS
        bWaitBar
        bShowResult
    end
    
    methods (Test)
        function regressionTests (testCase)
            testCase.bWaitBar = true;      % true:  watch internal stage progress
            testCase.bShowResult = true;   % true, show final result and pause after each test
            setTsDefaults(testCase)
            
            % comment out any line below to skip those tests
%             testNyquist (testCase)
%             testNearDC (testCase)
            testLab (testCase)
        end
    end
    
    %----------------------------------------------------------------------
    methods (Access = private)
        
        function testConstructor (testCase)
            oft = OFT();            
        end
                            
        function testOftOnce (testCase)
            disp(testCase.TS.Name)
            testCase.TS = testCase.TS.makeTime;
            testCase.TS = testCase.TS.makeTS;

            oft = OFT();
            oft.bWaitBar = testCase.bWaitBar;
            [actFreqs, actOFT, actFracErr] = oft.OFT_fn(testCase.TS.Ts, testCase.TS.time);
            
            if testCase.bShowResult
                t = testCase.TS.time;
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
            end
            
            disp(mat2str(actFracErr));
            disp(mat2str(actFreqs));
            disp(mat2str(actOFT));
            
            if testCase.bShowResult; pause; end
            close all
            
        end
        %-------------------
        
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

        
        %-------------------
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
            testCase.TS.Freqs = 1/300;
            testCase.TS.Amps = 1;
            testCase.TS.Phases = 0;
            for i=1:nSinusoidsMax
                testCase.TS.Amps(i) = 20 + i;
                testCase.TS.Phases(i) = (117 + 10*i)*pi/180;
            end
            testCase.TS.NoiseUniformLow = 0;
            testCase.TS.NoiseUniformHi = 0;
            testCase.TS.NoiseGaussMean = 0;
            testCase.TS.NoiseGaussSD = 0;
        end
        %-------------------
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
           testCase.TS.Amps(7) = 8; 
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
           testCase.TS.Phases = [0,90,45,230,0,20,15,0,215,0] * pi/180;
           testOftOnce (testCase);
           
           testCase.TS.Name = 'Two Sinusoids too close'; 
           testCase.TS.Freqs(6) = 0.00343642611683849;
           testOftOnce (testCase);
           
           testCase.TS.Name = 'Two Sinusoids consolidated'; 
           testCase.TS.Freqs(5) = 0.00332225913621262;
           testCase.TS.Freqs(6) = 0.00334448160535117;
           testOftOnce (testCase);
                                                                                                                      
       end            
                   
        %-------------------
        function testNyquist (testCase)
            kNuEgdeWidth = 0.00005;
            
            testCase.TS.Name = 'Test Nyquist (600)';
            testCase.TS.nSamples = uint32(600);
            testCase.TS.Freqs = 0.5 * testCase.TS.nSamples/testCase.TS.Extent;
            testOftOnce (testCase);
            
            testCase.TS.Name = 'Test Nyquist - kNuEgdeWidth - 1e-7 (600)';
            testCase.TS.Freqs = (0.5 * testCase.TS.nSamples - kNuEgdeWidth - 1e-7) / testCase.TS.Extent;
            testOftOnce (testCase);
            
            testCase.TS.Name = 'Test Nyquist - kNuEgdeWidth + 1e-7 (600)';
            testCase.TS.Freqs = (0.5 * testCase.TS.nSamples - kNuEgdeWidth + 1e-7) / testCase.TS.Extent;
            testOftOnce (testCase);
            
            % odd samples near nyquist.
            testCase.TS.Name = 'Test Nyquist (539)';
            testCase.TS.nSamples = uint32(7 * 7 * 11);  % odd
            testCase.TS.Freqs = 0.5 * testCase.TS.nSamples/testCase.TS.Extent;
            testOftOnce (testCase);
            
        end
        %-------------------
        function testNearDC (testCase)
            kNuEgdeWidth = 0.00005;            
            testCase.TS.nSamples = uint32(600);    
            
            testCase.TS.Name = 'Test DC  (600)';            
            testCase.TS.Freqs = 0;
            testOftOnce (testCase);            
            
            testCase.TS.Name = 'Test DC + 1e-7 (600)';            
            testCase.TS.Freqs = kNuEgdeWidth + 1e-7;
            testOftOnce (testCase);

            testCase.TS.Name = 'Test DC- 1e-7 (600)';            
            testCase.TS.Freqs = kNuEgdeWidth - 1e-7;
            testOftOnce (testCase);
            
            testCase.TS.nSamples = uint32(7 * 7 * 11);  % odd  
            testCase.TS.Name = 'Test DC  (539)';            
            testCase.TS.Freqs = 0;
            testOftOnce (testCase);                        
            
            testCase.TS.Name = 'Test DC + 1e-7 (539)';            
            testCase.TS.Freqs = kNuEgdeWidth + 1e-7;
            testOftOnce (testCase);

            testCase.TS.Name = 'Test DC - 1e-7 (539)';            
            testCase.TS.Freqs = kNuEgdeWidth - 1e-7;
            testOftOnce (testCase);
            
                        
        end       
        %-------------------
        
        
        
    end
    

end

