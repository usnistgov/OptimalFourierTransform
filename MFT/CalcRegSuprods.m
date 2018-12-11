function [cc,cs,sc,ss]=CalcRegSuprods(nu,mu,n)
% Ported (with permission) by Allen Goldstein, NIST from:
% http://jonova.s3.amazonaws.com/cfa/climate.xlsm/Suprods(vba)
% written by: Dr David Evans
%             david.evans@sciencespeak.com
%=========================================================================%
% The below code is not part of the function call.  It returns handles to
% all the local subfunctions for the purpose of unit testing if the value
% of the "nu" input (normally a double) is the string '-test'
if ischar(nu) && strcmp(nu,'-test')
    cs = 0;
    sc = 0;
    ss = 0;
    cc = localfunctions;
    return
end
%=========================================================================%
% The actual CalcRegSuprods call begins here
%- nu and mu can be outside [0, N/2].
%- Tries a divide and conquer approach if N is divisible by small primes.
%- Compared to CalcRegSuprods_Naive on a 5-sinusoid fit, reduced time from 
%  15 secs to 6 secs.

persistent lastN_SU;
persistent extensionsSU;
persistent factorsSU;
if isempty(lastN_SU) || lastN_SU ~= n
    [lastN_SU,extensionsSU,factorsSU] = Factorize(n);
end

% if nu and mu are integers, then there is a quicker way to calculate the Suprods
if floor(nu)==nu && floor(mu)==mu
    [cc,cs,sc,ss] = CalcRegSuprodsForIntegralArguments(nu,mu,n);
    return
end

extensions = extensionsSU;
[cc,cs,sc,ss] = CalcRegSuprods_Naive(nu/extensions, mu/extensions, n/extensions);

i = 1;
while extensions > 1
    extensions = extensions / factorsSU(i);
    if factorsSU(i) == 2
        [cc,cs,sc,ss]=ExtendRegSuprodBy2 (nu/extensions,mu/extensions,cc,cs,sc,ss);
    else
        [cc,cs,sc,ss]=ExtendRegSuprodByP (factorsSU(i),nu/extensions,mu/extensions,cc,cs,sc,ss);
    end
    i = i + 1;
    if i > length(factorsSU)
        break;
    end
end


end


function [lastN_SU,extensionsSU,factorsSU] = Factorize(n)
lastN_SU = n;
maxSmallestFactorSU = floor(n^2);
nFactors_SU = 0;
extensionsSU = 1;
while n > 1
    [p] = SmallestPrimeFactorOf(n);
    nFactors_SU = nFactors_SU+1;
    factorsSU(nFactors_SU)=p;
    n = n/p;
    if p < lastN_SU && p <= 50
        extensionsSU = extensionsSU * p;
    end
end
end

function [p] = SmallestPrimeFactorOf(n)
if 2 * floor(n/2) == n
    p = 2;
    return
end
for i = 3:2:floor(n^2)
    div = n/i;
    if floor(div) == div
        p = i;
        return
    end
end
p = n;   
end

function [cc,cs,sc,ss] = CalcRegSuprods_Naive(nu, mu, n)
% - nu and mu can be outside [0, N/2].
% - Straightforward method, as per defining formulae.

if floor(nu)==nu && floor(mu)==mu
    [cc,cs,sc,ss] = CalcRegSuprodsForIntegralArguments(nu,mu,n);
    return
end

oneOnN = 1/n;
twoPiON =  2*pi*oneOnN;
twoPiNuON = twoPiON*nu;
cc = 1;
cs = 0;
sc = 0;
ss = 0;

if nu == mu
   for tau = n-1:-1:1
       radiansNu = twoPiNuON * tau;
       c = cos(radiansNu);
       s = sin(radiansNu);
       cc = cc + c * c;
       cs = cs + c * s;
       ss = ss + s * s;
   end
    sc = cs;
else
   twoPiMuON = twoPiON*mu;
   for tau = n-1:-1:1
       radiansNu = twoPiNuON * tau;
       radiansMu = twoPiMuON * tau;
       CNu = cos(radiansNu);
       SNu = sin(radiansNu);
       CMu = cos(radiansMu);
       SMu = sin(radiansMu);
       cc = cc + CNu * CMu;
       cs = cs + CNu * SMu;
       sc = sc + SNu * CMu;
       ss = ss + SNu * SMu;
   end
end
cc = cc * oneOnN;
cs = cs * oneOnN;
sc = sc * oneOnN;
ss = ss * oneOnN;
end

function [cc, cs, sc, ss]=ExtendRegSuprodBy2(nu, mu, cc, cs, sc, ss)
radiansNu = pi*nu;
radiansMu = pi*mu;
CNu = cos(radiansNu);
SNu = sin(radiansNu);
CMu = cos(radiansMu);
SMu = sin(radiansMu);
ccHat = CNu * CMu;
csHat = CNu * SMu;
scHat = SNu * CMu;
ssHat = SNu * SMu;
ccx = cc + ccHat * cc - csHat * cs - scHat * sc + ssHat * ss;
csx = cs + csHat * cc + ccHat * cs - ssHat * sc - scHat * ss;
scx = sc + scHat * cc - ssHat * cs + ccHat * sc - csHat * ss;
ssx = ss + ssHat * cc + scHat * cs + csHat * sc + ccHat * ss;

cc = 0.5 * ccx;
cs = 0.5 * csx;
sc = 0.5 * scx;
ss = 0.5 * ssx;

end

function[cc, cs, sc, ss]=ExtendRegSuprodByP(p, nu, mu, cc, cs, sc, ss)
%- In:   p                 p >= 2
%        cc, cs, sc, ss    Suprod values for (nu/p, mu/p, N/p)
%- Out:  cc, cs, sc, ss    Suprod values for (nu, mu, N)
radiansNuNoi = 2*pi * nu / p;
radiansMuNoi = 2*pi * mu / p;
ccx = cc;
csx = cs;
scx = sc;
ssx = ss;

for i = 1:p-1
    radiansNu = radiansNuNoi * i;
    radiansMu = radiansMuNoi * i;
    CNu = cos(radiansNu);
    SNu = sin(radiansNu);
    CMu = cos(radiansMu);
    SMu = sin(radiansMu);
    ccHat = CNu * CMu;
    csHat = CNu * SMu;
    scHat = SNu * CMu;
    ssHat = SNu * SMu;
    ccx = ccx + ccHat * cc - csHat * cs - scHat * sc + ssHat * ss;
    csx = csx + csHat * cc + ccHat * cs - ssHat * sc - scHat * ss;
    scx = scx + scHat * cc - ssHat * cs + ccHat * sc - csHat * ss;
    ssx = ssx + ssHat * cc + scHat * cs + csHat * sc + ccHat * ss;
end
cc = ccx / p;
cs = csx / p;
sc = scx / p;
ss = ssx / p;
end

function [cc,cs,sc,ss]=CalcRegSuprodsForIntegralArguments(nu,mu,n)
% Ported (with permission) by Allen Goldstein, NIST from:
% http://jonova.s3.amazonaws.com/cfa/climate.xlsm/Suprods(vba)
% written by: Dr David Evans
%             david.evans@sciencespeak.com
%=============================================================

[nu,multNu]=MoveFreqIxInto0ToHalfN(n,nu);
[mu,multMu]=MoveFreqIxInto0ToHalfN(n,mu);
cs = 0;
sc = 0;
if nu == mu
    if nu == 0 || nu == 0.5 * n
        cc = 1;
        ss = 0;
    else
        cc = 0.5;
        ss = 0.5 * multNu * multMu;
    end
else
    cc = 0;
    ss = 0;
end
end

