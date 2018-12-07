function [ts, absDev] = SubtractMultipleSinusoidsFromTS (ts, cosPart, sinPart, nu)
%- Subtracts sinusoids with cosPart(i), sinPart(i) at frequency index nu() from time series ts.
%  Computes absolute deviation of ts after subtraction.
%- In:  ts     [0..N-1]  Time series
%       cosPart[1..nNu]  The sinusoids to subtract [tau = 0,...,N-1] are
%       sinPart[1..nNu]  cosPart(i) * Cos(2 pi nu(i) tau / N) + sinPart(i) * Sin(2 pi nu(i) tau / N)
%       nu     [1..nNu]
%- Out: ts     [0..N-1]  Time series
%       absDev           Absolute deviation of output ts
n = length(ts);
twoPiON = 2*pi/n;
twoPiNuONArr = nu * twoPiON;

absDev = 0;
for tau = n-1:-1:0
    x = ts(tau);
    for i = 1:length(nu)
        radians = twoPiNuONArr(i) * tau;
        x = x - (cosPart(i) * cos(radians) + sinPart(i) * sin(radians)
    end
    ts(tau) = x;
    x = abs(x);
    absDev = absDev + x
end
end
