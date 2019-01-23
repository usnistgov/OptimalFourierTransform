classdef testSW < matlab.unittest.TestCase
% Class-based unit testing for Shapiro Wilk
%
%   run these tests with two command-line commands:
%   - testCase = testSW;
%   - res = run(testCase);
%    
    properties
        TS
    end
    
    methods (Test)
        function regressionTests (testCase)
           testNoise (testCase)
        end
    end
    
    %----------------------------------------------------------------------
    methods (Access = private)
        function  setTSDefaults(testCase)
            Fs = 60;
            duration = 5;
            testCase.TS = ArtificialTS;
            testCase.TS.T0 = 0;
            testCase.TS.Extent = duration;
            testCase.TS.nSamples = uint32(duration * Fs);
            testCase.TS.Freqs = [1,3];
            testCase.TS.Amps = [1,1];
            testCase.TS.Phases = [0,45]*pi/180;
            
            testCase.TS.NoiseUniformLow = 0;
            testCase.TS.NoiseUniformHi = 0;
            testCase.TS.NoiseGaussMean = 0;
            testCase.TS.NoiseGaussSD = 0;
            
            testCase.TS = testCase.TS.makeTime;
        end
        
        function oneTest(testCase)
           testCase.TS = testCase.TS.makeTS; 
           
%            plot(testCase.TS.time,testCase.TS.Ts)
%            pause
                                
           oft = OFT_ShapiroWilk();
           %oft = OFT();
           
%            oft.bWaitBar = true;         
%            [freqs,MFT,fracErr] = oft.OFT_fn(testCase.TS.Ts,testCase.TS.time);
%            disp(mat2str(freqs))
%            disp(mat2str(MFT))
           
           
           [ax,bx,cx,foundMinimum,xMin,fxMin] = oft.AcfBracketGuess1D([14.98,5.02],[1,0],151,testCase.TS.Ts);
           msg = sprintf('ax = %f, bx = %f, cx = %f, foundMinimum = %i, xMin = %f, fxMin = %f',ax,bx,cx,foundMinimum,xMin,fxMin);
           disp(msg)
%            pause
%            [sumAbs, maxCorr] = oft.AcfBracketGuess1D([15,5.02],[1,0],151,testCase.TS.Ts);
%            msg = sprintf('sumAbs = %f, maxCorr = %f',sumAbs,maxCorr);
%            disp(msg)
%            pause
%            [sumAbs, maxCorr] = oft.AcfBracketGuess1D([15,5.02],[1,0],151,testCase.TS.Ts);
%            msg = sprintf('sumAbs = %f, maxCorr = %f',sumAbs,maxCorr);
%            disp(msg)
%            pause
%            [sumAbs, maxCorr] = oft.AcfBracketGuess1D([15,5.02],[1,0],151,testCase.TS.Ts);
%            msg = sprintf('sumAbs = %f, maxCorr = %f',sumAbs,maxCorr);
%            disp(msg)
%            pause
           

           
           %[ax,bx,cx,foundMaximum,xMax,fxMax] = oft.BracketMaxShapiroWilk([15,5],[1,0],151,testCase.TS.Ts);
%            plot(testCase.TS.time,testCase.TS.Ts)
%            pause
           close all
        end
    
        function testNoise(testCase)
            setTSDefaults(testCase)
            testCase.TS.NoiseGaussMean = 0;
            testCase.TS.NoiseGaussSD = 0;
            oneTest(testCase)
        end
    end
    
end

