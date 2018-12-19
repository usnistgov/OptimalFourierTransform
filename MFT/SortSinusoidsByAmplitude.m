function [freqs, Mft] = SortSinusoidsByAmplitude (freqs, Mft, removeFraction)
%- Orders sinusoids by amplitude, and removes any whose amplitude is less than removeFraction of
%  the second largest max amplitude (there is often one huge amplitude of a near-zero frequency).

nFreqs = length(freqs);
%ampSq = zeros(1,nFreqs);
maxAmpSq = 0;
secondMaxAmpSq = 0;

% for i = 1:nFreqs
%    c = real(Mft(i));
%    s = imag(Mft(i));
%    x = c^2+s^2;
%    ampSq(i) = x;
%    if x > secondMaxAmpSq
%       if x > maxAmpSq
%          secondMaxAmpSq = maxAmpSq;
%          maxAmpSq = x;   
%       else
%        secondMaxAmpSq = x;
%       end
%    end
% end

ampSq = abs(Mft).^2;
for i = 1:nFreqs
   if ampSq(i) > secondMaxAmpSq 
       secondMaxAmpSq = ampSq(i);
       if ampSq(i) > maxAmpSq
           secondMaxAmpSq = maxAmpSq;
           maxAmpSq = ampSq(i);
       end
   end
end

% count frequencies to be removed
threshold = secondMaxAmpSq * removeFraction;
nRemove = 0;
for i = 1:nFreqs
   if ampSq(i) < threshold; nRemove = nRemove+1; end 
end

indexArr = zeros(1,nFreqs);
for i = 1:nFreqs
   indexArr(i) = i;
   ampSq(i) = -ampSq(i);
end

% sort freqs in order of amplitude
[indexArr]=SortArrayOfIndexesToDoubles(indexArr, ampSq, 1, nFreqs);
tempF = freqs;
tempMFT = Mft;
for i = 1:nFreqs
   freqs(i) = tempF(indexArr(i));
   Mft(i) = tempMFT(indexArr(i));
end

% remove frequencies whos amplitudes are too small
nFreqs = nFreqs - nRemove;
freqs = freqs(1:nFreqs);
Mft = Mft(1:nFreqs);
end
