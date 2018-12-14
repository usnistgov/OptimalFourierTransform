function [freqs, MFT] = SortSinusoidsByAmplitude (freqs, MFT, removeFraction)
%- Orders sinusoids by amplitude, and removes any whose amplitude is less than removeFraction of
%  the second largest max amplitude (there is often one huge amplitude of a near-zero frequency).

nFreqs = length(freqs);
ampSq = zeros(1,nFreqs);
maxAmpSq = 0;
secondMaxAmpSq = 0;

for i = 1:nFreqs
   c = real(MFT(i));
   s = imag(MFT(i));
   x = c^2+s^2;
   if x > secondMaxAmpSq
      if x > maxAmpSq
         secondMaxAmpSq = maxAmpSq;
         maxAmpSq = x;   
      else
       secondMaxAmpSq = x;
      end
   end
end

% count frequencies to be removed
threshold = secondMaxAmpSq * removeFraction;
nRemove = 0;
for i = 1:nFreqs
   if ampSq(i) < threshold; nRemove = nRemove+1; end 
end

% sort freqs in order of amplitude
[indexArr]=SortArrayOfIndexesToDoubles(indexArr, valueArr, 1, nFreqs);
tempF = freqs;
tempMFT = MFT;
for i = 1:nFreqs
   freq(i) = tempF(indexArr(i));
   MFT(i) = tempMFT(indexArr(i));
end

% remove frequencies whos amplitudes are too small
nFreqs = nFreqs - nRemove;
freqs = freqs(1:nFreqs);
MFT = MFT(1:nFreqs);
end
