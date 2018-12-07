function [freq] = SortSinusoidsByAmplitude (freq, cosPart, sinPart, removefraction)
%- Orders sinusoids by amplitude, and removes any whose amplitude is less than removeFraction of
%  the second largest max amplitude (there is often one huge amplitude of a near-zero frequency).
