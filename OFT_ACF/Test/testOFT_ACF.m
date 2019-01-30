classdef testOFT_ACF < matlab.unittest.TestCase
% Class-based unit testing for Optimal Fourier Transform - Autocorrelation
%
%   run these tests with two command-line commands:
%   - testCase = testOFT_ACF;
%   - res = run(testCase);
%    
    properties
        TS
        bWaitBar
        bDoRecon
        bDoAcf
        bShowResult
        kAcfThreshold
    end
    
    methods (Test)
        function regressionTests (testCase)
            testCase.bShowResult = true;   % true, show final result and pause after each test
            testCase.bWaitBar = true;      % Show the progress bar and internal progress
            testCase.bDoRecon = true;
            testCase.bDoAcf = true;
            % comment out any line below to skip those tests
%             testNyquist (testCase)
%             testNearDC (testCase)
%             testLab (testCase)
%             testACF (testCase)
            test50_45 (testCase)
        end
    end
    
    %----------------------------------------------------------------------
    methods (Access = private)
        function testOftOnce (testCase)
            disp(testCase.TS.Name)
            testCase.TS = testCase.TS.makeTime;
            testCase.TS = testCase.TS.makeTS;

            oft = OFT_ACF();
            oft.bWaitBar = testCase.bWaitBar;
            oft.bDoRecon = testCase.bDoRecon;
            oft.bDoAcf = testCase.bDoAcf;
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
        function setRegressionDefaults(testCase)
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
        function testACF(testCase)
            setACFDefaults(testCase)
            testCase.TS.NoiseGaussMean = 0;
            testCase.TS.NoiseGaussSD = 0.01;
            %testCase.TS.NoiseUniformLow = -0.005;
            %testCase.TS.NoiseUniformHi = 0.005;            
            %testAcfMinimize(testCase)
            testOftOnce(testCase)
        end       
        %-------------------              
        function testAcfMinimize(testCase)
           testCase.TS = testCase.TS.makeTS;                                 
           oft = OFT_ACF();
           
%            oft.bWaitBar = true;         
%            [freqs,MFT,fracErr] = oft.OFT_fn(testCase.TS.Ts,testCase.TS.time);
%            disp(mat2str(freqs))
%            disp(mat2str(MFT))

            %DEBUG figure
            oft.Fig1 = figure(1);
            clf(oft.Fig1);
            figure(oft.Fig1);
            %hold on

            ts = testCase.TS.Ts;
            [nuGuessW] = oft.MinimizeAcfByVaryingMultipleFreqs(ts, [4.5,15.5], 151);
            disp(mat2str(nuGuessW));
                        
%             % DEBUG ----------------------------------
%             [cosPart, sinPart] = EstimateContainedSinusoids(ts, nuGuessW);
%             [resid,~] = SubtractMultipleSinusoidsFromTS (ts,cosPart, sinPart, nuGuessW);
%             sumAcfResid = oft.SumAbsACF(resid);
%             nuGuessW = 0;
%             [cosPart, sinPart] = EstimateContainedSinusoids(ts, nuGuessW);
%             [resid,~] = SubtractMultipleSinusoidsFromTS (ts,cosPart, sinPart, nuGuessW);            
%             plot(ts-resid)
%             sumAcfOrig = oft.SumAbsACF(resid);
%             msg = sprintf('sumAcfResid = %f, sumAcfOrig = %f',sumAcfResid, sumAcfOrig);
%             disp(msg)
%             % DEBUG --------------------------------------
           
            figure(oft.Fig1);
            close(oft.Fig1);
            %hold off

%            [ax,bx,cx,foundMinimum,xMin,fxMin] = oft.AcfBracketGuess1D([14.98,5.02],[1,0],151,testCase.TS.Ts);
%            msg = sprintf('ax = %f, bx = %f, cx = %f, foundMinimum = %i, xMin = %f, fxMin = %f',ax,bx,cx,foundMinimum,xMin,fxMin);
%            disp(msg)
%            pause
           
           
           %[ax,bx,cx,foundMaximum,xMax,fxMax] = oft.BracketMaxShapiroWilk([15,5],[1,0],151,testCase.TS.Ts);
%            plot(testCase.TS.time,testCase.TS.Ts)
%            pause
           close all
        end             
        %-------------------              
        function  setACFDefaults(testCase)
            Fs = 60;
            duration = 5;
            testCase.TS = ArtificialTS;
            testCase.TS.T0 = 0;
            testCase.TS.Extent = duration;
            testCase.TS.nSamples = uint32(duration * Fs);
            testCase.TS.Freqs = [1,3];
            testCase.TS.Amps = [1,.5];
            testCase.TS.Phases = [0,45]*pi/180;
            
            testCase.TS.NoiseUniformLow = 0;
            testCase.TS.NoiseUniformHi = 0;
            testCase.TS.NoiseGaussMean = 0;
            testCase.TS.NoiseGaussSD = 0;
            
            testCase.TS = testCase.TS.makeTime;
        end
        %-------------------               
        function testLab (testCase)
            setRegressionDefaults(testCase)
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
            
            testCase.TS.NoiseGaussMean = 0;
            testCase.TS.NoiseGaussSD = .01;            
            
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
       function test50_45 (testCase)
           
           % load the time series
            load('5045.mat')                        
            
            testCase.bDoRecon = true;
            testCase.bDoAcf = true;
            testCase.kAcfThreshold = 1e-12;
            
            testCase.TS = ArtificialTS;
            Fs = 50;
            duration = 5;
            testCase.TS.Name = '50F045F';
            testCase.TS.T0 = 0;
            testCase.TS.Extent = duration;
            testCase.TS.nSamples = uint32(duration * Fs);            
            testCase.TS.Ts = TS.MagErr(:,6);
            
            testCase.TS = testCase.TS.makeTime;
            
%             % subtract the best fit line
%             P = polyfit(testCase.TS.time',testCase.TS.Ts,1);
%             yfit = P(1)* testCase.TS.time'+P(2);
%             testCase.TS.Ts = testCase.TS.Ts - yfit;
            
            
            oft = OFT_ACF();
            oft.bWaitBar = testCase.bWaitBar;
            oft.bDoRecon = testCase.bDoRecon;
            oft.bDoAcf = testCase.bDoAcf;
            oft.kAcfThreshold = testCase.kAcfThreshold;
            [actFreqs, actOFT, actFracErr] = oft.OFT_fn(testCase.TS.Ts, testCase.TS.time);
            
            if testCase.bShowResult
                t = testCase.TS.time;
                [act_TS] = testCase.SynthesizeOFT(actFreqs, t, actOFT);
                orig_TS = testCase.TS.Ts';
                
                
                close all
                figure(2)
                s(1)=subplot(2,1,1);
                hold on
                plot (t,orig_TS,'b')
                plot(t,act_TS,'r')
                hold off
                
                s(2)=subplot (2,1,2);
                ts_diff = act_TS - orig_TS;
                
                plot (t, ts_diff,'g.')
                title(s(1),testCase.TS.Name);
                title(s(2),'Residual');
                figure(3)
                hist(ts_diff);
                %hold on
                %hist((randn(length(ts_diff),1)*std(ts_diff))+mean(ts_diff));

            end
            
            disp(mat2str(actFracErr));
            disp(mat2str(actFreqs));
            disp(mat2str(abs(actOFT)));
            
            if testCase.bShowResult; pause; end
            
            
%             for i = 1:length(actFreqs)
%                 [act_TS] = testCase.SynthesizeOFT(actFreqs(1:i), t, actOFT(1:i));
%                 close all
%                 figure(2)
%                 s(1)=subplot(2,1,1);
%                 hold on
%                 plot (t,orig_TS,'b')
%                 plot(t,act_TS,'r')
%                 hold off
%                 
%                 s(2)=subplot (2,1,2);
%                 ts_diff = act_TS - orig_TS;
%                 
%                 plot (t, ts_diff,'g')
%                 title(s(1),testCase.TS.Name);
%                 title(s(2),'Residual');
%                 drawnow
%             end
                
            
            
            close all            
            
        end
        
    end
    
end

