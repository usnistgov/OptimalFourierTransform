function [fracErr] = FractionalError (ts, extent, freqs, MFT)
%- Returns (the absolute deviation of the ts synthesizd from cosPart and sinPart) divided by
%  (the absolute deviation of ts).
%- Provides a figure of merit for how good a transform is. The smaller the better, less than 1% is good.
%- In:  ts     [0..N-1]      Time series of extent
%       freq   [1..nFreqs]
%       cosPart[1..nFreqs]
%       sinPart[1..nFreqs]

cosPart = real(MFT);
sinPart = imag(MFT);
n = length(ts);
twoPiExtentON = 2 * pi * extent / n;
absDevOfTS = 0;
absDevOfError = 0;

for tau = n-1:-1:0
   
    twoPiExtentTauOn = twoPiExtentON * tau;
    synthesizedValue = 0;
    for i = 1:length(freqs)
       radians =  twoPiExtentTauOn * freqs(i);
       synthesizedValue = synthesizedValue + cosPart(i) * cos(radians) + sinPart(i) * sin(radians);
    end
    
    orig = ts(tau+1);
    err = abs(orig - synthesizedValue);
    
    absDevOfError = absDevOfError + err;
    
    orig = abs(orig);
    absDevOfTS = absDevOfTS + orig;             
end
%absDevOfOriginalTS = AbsDevOfTS_Fn(ts, n)
fracErr = absDevOfError / absDevOfTS;
end


