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
            [f,~,~]=MFT('-test',{},{},{});
            testCase.localFunctions = f;
        end
    end
    
    methods (Test)
        
        function testConsolidateFreqs (testCase)
%             handle=testCase.localFunctions(1);
%             ConsolidateFreqs = handle{1};
            kNuConsolidate = 0.1;
            nu_MFT = [1, 5, 5.05, 6, 16, 16.025, 25, 25];
            [act] = ConsolidateFreqs(nu_MFT,kNuConsolidate);
            exp = [1, 5.0250, 6, 16.0125, 25];
            testCase.verifyEqual(act,exp);
        end
        
        function testConsolidateFreqsByBracket (testCase)
            handle=testCase.localFunctions(1);
            ConsolidateFreqsByBracket = handle{1};
            kNuConsolidate = 0.1;
            nu_MFT = [1, 5, 5.05, 6, 16, 16.025, 25, 25, 30, 32, 32.075, 45, 59];
            bracket = [4,5,4];
            [act] = ConsolidateFreqsByBracket(nu_MFT,bracket,kNuConsolidate);
            exp = [1, 5.025, 6. 16.0125, 25, 30, 32.0375, 45, 59];
            testCase.verifyEqual(act,exp);
        end
        
        function testMFT_Fn (testCase)
            load('tsData.mat');
            t = double(tsData.t); 
            t0 = t(1);
            deltaT = mean(t(2:end)-t(1:end-1));
            for tau=0:length(t)-1
                regT(tau+1) = t0+deltaT*tau;
            end
            ts = double(tsData.g);
            
            Figure = figure;
            plot(regT,ts,'b')
            hold on
            
            [actFreqs, actMFT, actFracErr] = MFT(ts, regT, Freqs, []);
            expMFT = [0.840472028559792*exp((337.811918578934i*pi/180)),...
                      0.819839127174502*exp((197.267337528315i*pi/180)),...
                      0.30879463509027*exp((351.14365854658i*pi/180)),...
                      0.0863292731948021*exp((356.72102222165i*pi/180)),...
                      0.142350768427654*exp((242.759379140154i*pi/180)),...
                      0.0545294288353872*exp((200.587845563089i*pi/180)),...
                      0.06172570381516*exp((281.21433966186i*pi/180)),...
                      0.0890214288150472*exp((243.659819109798i*pi/180))];
              testCase.verifyEqual(actMFT,expMFT,'absTol',0.00001+0.00001);
              
              expFracErr = 0.68673537826427;
              testCase.verifyEqual(actFracErr,expFracErr,'absTol',1e-6);
                                    
             ts_MFT = testCase.SynthesizeMFT(actFreqs, regT, actMFT);
             plot(regT,real(ts_MFT),'r')
                 
            hold off
            pause
            close(Figure)

            
        end
    end

    methods (Access = private)
        
        function [ts_MFT] =  SynthesizeMFT(testCase, Freqs, t, MFT)
           
            n = length(t);
            extent = t(end)-t(1);
            twoPiExtentON = 2 * pi * extent / n;
            ts_MFT = zeros(1,length(t));
            for i = 1:length(Freqs)
               mag = abs(MFT(i));
               ang = angle(MFT(i));
               twoPiTime = zeros(1,length(t));
               for tau = n-1:-1:0
                   twoPiTime(tau+1) = twoPiExtentON * tau;
               end
               radians = twoPiTime * Freqs(i);
               ts_MFT = ts_MFT + mag * cos(radians-ang);
            end
                                   
        end
        
    end
    
end


        
        