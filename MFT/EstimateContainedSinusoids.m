function [cosPart, sinPart] = EstimateContainedSinusoids(g, nu)
% estimates sinusoids listed in nu()
% Ported (with some changes listed below) from VBA by Allen Goldstein, NIST from:
% http://jonova.s3.amazonaws.com/cfa/climate.xlsm
% written by: Dr David Evans
%             david.evans@sciencespeak.com
%
% changes:  
%   - This code wors only with "regular" time series data "g".
%     Regular means that the sampling rate of the data is constant.
%   - Rather than returning separate sine and cosine parts, this returns 
%     complex "MFT" values
%-------------------------------------------------------------------------%
% - input:
%       g [0..N-1] Time Series Data
%       t [0..N-1] time of the data points
%       nu [1..M]  distinct, real-valued frequency indices.  
%                   Algorithm fails if the nu values are extremely close 
%                   together so consolidate first.
% - output:
%       MFT [1..M] complex sine and cosine parts
%

%=========================================================================%
% The below code is not part of the function call.  It returns handles to
% all the local subfunctions for the purpose of unit testing if the value
% of the "nu" input (normally a double) is the string '-test'
if ischar(g) && strcmp(g,'-test')
    cosPart = localfunctions;
    sinPart = 0;
    return
end
%================================================================+========%


[isEdgeNu] = ClassifyRegFreqsAsEdgeOrNormal(length(g),nu);
[aa] = ComputeTheAAMatrix(length(g), nu, isEdgeNu);
[cosAvg, sinAvg] = ComputeRegCosAndSinAverages (g, nu);

nEdgeNu = 0;
rowN = 1;
bb = zeros(1,length(nu));
for i = 1:length(nu)
    bb(rowN) = cosAvg(i);
    if isEdgeNu(i)
        rowN = rowN + 1;
        nEdgeNu = nEdgeNu + 1;
    else
        bb(rowN+1) = sinAvg(i);
        rowN = rowN + 2;
    end
end

[bb] = SolveLinearEquations (aa, 2 * length(nu) - nEdgeNu, bb);

rowN = 1;
    cosPart = zeros(1,length(nu));
    sinPart = zeros(1,length(nu));
for i = 1:length(nu)
    cosPart(i) = bb(rowN);
    if isEdgeNu(i)
        sinPart(i) = 0;
        rowN = rowN + 1;
    else
        sinPart(i) = bb(rowN +1);
        rowN = rowN+2;
    end
end

end

%============================== LOCAL FUNCTIONS ==========================%
function [isEdgeNu] = ClassifyRegFreqsAsEdgeOrNormal(n,nu)
% 
%-----------------------------------------------------------------------------------------------
kNuEgdeWidth = 0.00005;
% 'Choosing the edgewidth
%----------------------
%- As a nu approaches an edge (0 or N/2):
%   1. The sine average at that nu approaches zero.
%   2. The row of the "aa" matrix for the sine average in that nu approaches all zeroes.
%   3. The col of the "aa" matrix for the sine part for that nu approaches all zeroes.
%- In this case, can get absurdly high values of the sine part at that freq, as it makes little difference to
%  anything in the equations described by aa. The amplitudes of the sinusoids thus found appear to "blow up",
%  becoming absurdly high.
%- The equations are correct, but it seems the suprod values that populate aa are wrong -- simply
%  because they are very small (typically E-10 and E-20) and thus dominated by roundoff error. So the values in
%  the aa matrix become random, so the computed values of the sine part at this nu becomes absurdly wrong.
%- The solution is to declare the zone close to the edge "the edge zone", where the edge zone covers all
%  the nu's where the roundoff error in suprod values is significant. In the edge zone, simply treat the
%  nu as at the edge, and thus omit the equation from the "aa" matrix.
%- Set the edge zone by trial and error, guided by assorted_TestEstimatingMultipleContainedSinusoids.
%  Obvously want the edge width to be as small as possible to get more accurate results, but big enough
%  to avoid the effect of roundoff in computing suprod values.
%- The value of aa that is the smallest, and thus most vulnerable to suprod roundoff, is at the
%  intersection of the row for the sine average of nu and the col for the sine part for nu, the ss suprod
%  value. At the nu near the edge, cc ~ 1, cs = sc ~ 10^-d, ss ~ 10^-2d. The roundoff error of a suprod
%  value using double arithemtic (which is good for 15 decimal places) for time series of 1000 data points or so
%  creates a noise floor around 10^-12 or so. So d is about 6, i.e. set edgewidth so ss(edgewidth) is around
%  10^-12.
%- Set the tol in TestEstContainedSinusoidOnce very low (maybe 10^-7) then watch the failure messages
%  for the degenerate cases as adjust kNuEdgeWidth. To stop the misses blowing up, need kNuEdgeWidth >= 10^-6,
%  or 10^-5 to be surer. Make it too large however, and get other errors ramping up due to artificial movement of
%  frequencies to the edge. Best compromise seems to be about 0.00002.
%- Note that the nu values for EstimateContainedSinusoids must be sufficiently distinct that the rows and
%  cols of aa are not close to being linearly dependent. So consolidate the frequencies first.
%-----------------------------------------------------------------------------------------------

halfNMinus = (n/2) - kNuEgdeWidth;
halfNPlus = (n/2) + kNuEgdeWidth;
isEdgeNu = false(1,length(nu));
for i = 1:length(nu)
    nuVal=nu(i);
    if nuVal <= kNuEgdeWidth
        if nuVal >= -kNuEgdeWidth
            isEdgeNu(i)=true;
        else
            [nuVal] = MoveFreqIxInto0ToHalfN(n,nuVal);
            isEdgeNu(i)= (nuVal < kNuEgdeWidth) | (nuVal > halfNMinus);
        end
    else
        if nuVal >= halfNMinus
            if nuVal <= halfNPlus
                isEdgeNu(i)= true;
            else
                [nuVal] = MoveFreqIxInto0ToHalfN(n,nuVal);
                isEdgeNu(i)= (nuVal < kNuEgdeWidth) | (nuVal > halfNMinus);
            end
            
        else
            isEdgeNu(i)=false;
        end
    end

end
end

function [aa]=ComputeTheAAMatrix(n,nu,isEdgeNu)
global flags
% Ported (with some changes listed below) from VBA by Allen Goldstein, NIST from:
% http://jonova.s3.amazonaws.com/cfa/climate.xlsm/EstimatingContainedSinusoids(vba)
% written by: Dr David Evans
%             david.evans@sciencespeak.com
%
% changes:  
%   - This code works only with "regular" time series data "g".  Regular means that the sampling rate of the data is constant.
%===============================================================================================%

%- Computes aa array.
%- Degenerate case of a zero row in the "aa" matrix arises iff row is for sine part of an edge nu (0 or N/2).
%- "aa" is a square symmetric matrix: 2 rows and cols for each normal nu, one row and col for each edge nu.

rowN = 1;
for i=1:length(nu)
    colN = 1;
    if isEdgeNu(i)          % row for edge nu
        for j = 1:i         % col
            [cc,cs,~,~] = CalcRegSuprods (nu(i), nu(j), n);
            if isEdgeNu(j)
                aa(rowN,colN)=cc;
                if j < i
                    aa(colN, rowN)=cc;
                end
                colN = colN+1;
            else
                aa(rowN,colN)=cc;
                aa(rowN,colN+1)=cs;
                if j < i
                    aa(colN,rowN)=cc;
                    aa(colN+1,rowN)=cs;
                end
                colN = colN+2;
            end
        end
        rowN = rowN+1;
    else
        for j = 1:i
            [cc,cs,sc,ss] = CalcRegSuprods (nu(i), nu(j), n);
            if isEdgeNu(j)
                aa(rowN,colN) = cc;
                aa(rowN+1,colN) = sc;
                if j < i 
                    aa(colN,rowN) = cc;
                    aa(colN,rowN+1) = sc;
                end
                colN = colN+1;
            else
                aa(rowN,colN)=cc;
                aa(rowN,colN+1)=cs;
                aa(rowN+1,colN)=sc;
                aa(rowN+1,colN+1)=ss;
                if j < i
                    aa(colN,rowN)=cc;
                    aa(colN+1,rowN)=cs;
                    aa(colN,rowN+1)=sc;
                    aa(colN+1,rowN+1)=ss;
                end
                 colN=colN+2;
            end 
        end
        rowN = rowN+2;
    end
end
    if isobject(flags)
        flags=flags.FlagLowElementsInMatrix(aa, 'l', 'L');
    end
end

function [cosAvg, sinAvg] = ComputeRegCosAndSinAverages(ts, nu)
%- Finds the cosine and sine averages at the frequency indices in nu in a regular time series ts.
%- In:  ts[0..N-1]      Time series
%       nu[1..nFreqs]   Real-valued frequency indices, can be outside [0,N/2]
%- Out: cosAvg[1..nFreqs]
%       sinAvg[1..nFreqs]
n = length(ts);
oneON = 1/n;
twoPiON = 2*pi/n;
for i = 1:length(nu)
   c = 0;
   s = 0;
   twoPiNuOnN = twoPiON * nu(i);
   for tau = n-1:-1:0
      g = ts(tau+1);
      radians = twoPiNuOnN * tau;
      c = c + g * cos(radians);
      s = s + g * sin(radians);
   end
   cosAvg(i) = c * oneON;
   sinAvg(i) = s * oneON;
end
end
