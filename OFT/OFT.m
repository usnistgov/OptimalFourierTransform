function [freqs, OFT, fracErr] = OFT (ts, time, relativeAbsDev, bDoRecon, maxNFreqsPerStage, bWaitBar)
global kMaxNFreqsAtOnceOFT targetAbsDevW minChangeAbsDevW Figure
% Ported from VBA by Allen Goldstein, NIST from:
% http://jonova.s3.amazonaws.com/cfa/climate.xlsm
% written by: Dr David Evans
%             david.evans@sciencespeak.com
%****************************************************************************************************
%- Computes the Optimal Fourier Transform (OFT) of a real-valued time series with data points spaced
%  equally in time. Very slow.
%- Reconnaissance slows it, but it can find frequencies that are close enough not to produce separate
%  peaks in the DFT. Generally worth it, though sometimes it triples the exectution time for no gain.
%- First stage does most of the work. Following stages just wring out structure out of what is left
%  in the time series, until either the fit is near-perfect or the residue is just pointwise noise.
%- In:  ts  [0..N-1]            Time series.
%       time[0..N-1]            Times of data points (relevant iff "regular" is true).
%       relativeAbsDev          OFT stops finding new frequencies if the absolute deviation (after sinusoids found
%                                 so far are removed) drops to this fraction of the absolute deviation of ts.
%       bDoRecon                 Perform a reconasssiance stage first. Slower but often better results.
%       maxNFreqsPerStage       Max number of sinusoids found by the OFT in any one stage.
%       bwaitBar                Display a wait bar while the stages run
%- Out: Freqs[1..nFreqs]         Frequencies at which OFT detects sinusoids, in cycles per unit of time.
%       OFT[1..nFreqs]          Complex Optimal Fourier transform
%       fracErr                 Absolute deviation of synthesisTS from OFT of ts / absolute deviation of ts.
%===============================================================================================%
% The below code is not part of the function call.  It returns handles to
% all the local subfunctions for the purpose of unit testing if the value
% of the "nu" input (normally a double) is the string '-test'
if ischar(ts) && strcmp(ts,'-test')
    OFT = 0+0i;
    fracErr = 0;
    freqs = localfunctions;
    return
end  
%===============================================================================================%
% constants         'Shared quantities for the OFT (suffix "W" = "working")
kMinChangeTargetAbsDevFraction = 0.0001;
kMaxNFreqsAtOnceOFT = 10;               % OFT can only consider this many sinosiuds at one go.
kMaxNFreqsPerStage = 3 * kMaxNFreqsAtOnceOFT / 2;
kMaxNStages = 3;  % 3?, 4?                                       
%===============================================================================================%

%n = length(ts);
extent = time(end)-time(1);
% Initialize absolute deviations
absDevOfOriginalTS_W = AbsDevOfTS(ts,length(ts));
targetAbsDevW = relativeAbsDev * absDevOfOriginalTS_W;
minChangeAbsDevW = kMinChangeTargetAbsDevFraction * targetAbsDevW;


% ts for the first stage
tsStage = ts;

% don't do too many frequecies per stage
if maxNFreqsPerStage < 10
    maxNFreqsPerStage = 10;
end
freqStage = [];
nFreqs = 0;
j = maxNFreqsPerStage;
if j < kMaxNFreqsPerStage; j = kMaxNFreqsPerStage; end
nuStage = zeros(1,j);

% stage counter and progress bar
stageN = 0;
waitBarOFT(stageN,'Starting First Stage',bWaitBar);


while stageN < kMaxNStages
   stageN = stageN+1; 
   waitBarOFT(stageN/kMaxNStages,sprintf('Stage %d',stageN),bWaitBar);  
   
   [freqStage, MftStage] = OneStageOfOFT(tsStage, extent, bDoRecon, maxNFreqsPerStage, freqStage);  
   if ishandle(Figure)
      hold on
      plot (ts - tsStage, 'r');
      hold off
   end

    if stageN >=2
       ampBiggest = abs(MftStage(1));
       ampSmallest = abs(MftStage(length(freqStage)));
       ratio = ampSmallest / ampBiggest;
       % exit if spectrum of tsStage is "flat enough"
       if ratio >=0.4; break; end
    end
    
    % Append stage sinusoids to result
    j = nFreqs;
    nFreqs = nFreqs + length(freqStage);
    for i = 1:length(freqStage)
       j = j + 1;
       freqs(j) = freqStage(i);
       MFT_OFT(j) = MftStage(i);
       nuStage(i) = freqStage(i) * extent;
    end
    
    % Compute the TS for the next stage
    [tsStage, absDev] = SubtractMultipleSinusoidsFromTS (tsStage, real(MFT_OFT), imag(MFT_OFT), nuStage);
    if absDev < targetAbsDevW; break; end
           
end
waitBarOFT(100,'',bWaitBar);

if stageN == 1
    [fracErr] = FractionalError(ts, extent, freqs, real(MFT_OFT), imag(MFT_OFT));
else
    [freqs] = SortSinusoidsByAmplitude (freqs, MFT_OFT, removefraction);
    fracErr = absDev / absDevOfOriginalTS_W;
end

clearvars -global

end

%**************************************************************************%
function waitBarOFT (progress, msg, bWaitBar)
persistent wb
if bWaitBar
    if isempty(wb)
        wb = waitbar(progress,msg);
    else
        if progress > 1
            close(wb)
            wb = [];
        else
            waitbar(progress,wb,msg)
        end
    end
end
end
%-------------
function [freqStage, MftStage, fracErr] = OneStageOfOFT(tsStage, extent, bDoRecon, maxNFreqsPerStage, freqStage)
global kMaxNFreqsAtOnceOFT

% *** Reconnaisance part ***
if bDoRecon
    maxNFreqsAtOnce = kMaxNFreqsAtOnceOFT / 2;
    maxNFreqsTotal = 3 * kMaxNFreqsAtOnceOFT / 2;
    if maxNFreqsTotal > maxNFreqsPerStage
        maxNFreqsTotal = maxNFreqsPerStage;
        maxNFreqsAtOnce = ceil(maxNFreqsPerStage/3);
    end
    
    [freqStage, ~, tsStage] = FindFreqIxsInTS (tsStage, maxNFreqsAtOnce, maxNFreqsTotal, true, freqStage);
end

% *** Main Part ***
[freqStage, bracket, tsStage] = FindFreqIxsInTS (tsStage, maxNFreqsAtOnce, maxNFreqsPerStage, false, freqStage);

% *** Final Part ***
dT = extent/length(tsStage);
t = 0:dT:extent-dT;
[freqStage, MftStage, fracErr] = MFT (tsStage, t, freqStage, bracket);
end
%----------------
function [nu, bracket, tsStage] = FindFreqIxsInTS (tsStage, maxNFreqsAtOnce, maxNFreqsTotal, recon, nu)
%- Finds best frequency indices for sinusoids in ts.
%- Tries maxNFreqsAtOnce sinusoids at once, and finds up to maxNFreqsTotal.
%- In:  nu[1..nNu]     Guesses for nu's in the initial iteration, if present (ie if nNu > 0).
%- Out: nu[1..nNu]     Best estimates of frequency indices of sinusoids in ts.
global targetAbsDevW minChangeAbsDevW Figure
kNuConsolidate = 0.1;   % How close should frequencies be allowed to get before they are consolidated
tsW = tsStage;
NW = length(tsW);
nuMaxCW = floor(NW/2);
%cosPart_DFT=zeros(1,nuMaxCW);
%sinPart_DFT=zeros(1,nuMaxCW);
cosPart = zeros(1,maxNFreqsTotal);
sinPart = zeros(1,maxNFreqsTotal);
passNArr = zeros(1,maxNFreqsTotal);

absDev = AbsDevOfTS(tsW,length(tsW));
lastAbsDev = 0;

nNuInitial = length(nu);
nNu = 0;
passN = 0;
while nNu < maxNFreqsTotal && absDev > targetAbsDevW && abs(lastAbsDev - absDev) > minChangeAbsDevW
    passN = passN + 1;
    if ishandle(Figure)
        plot (tsStage); 
        hold on
        plot (tsStage - tsW, 'r');
        hold off
        drawnow;
    end
    if nNuInitial > 0
        nuGuessW=nu;
        nFreqsW = nNuInitial;
        nNuInitial=0;
    else
        [Xk]=dft(tsW);
        [nuGuessW] = FindFreqsOfAmplitudePeaks(Xk, maxNFreqsAtOnce, nuMaxCW);
        nFreqsW = length(nuGuessW);
    end
    if nFreqsW > maxNFreqsTotal - nNu
        nFreqsW = maxNFreqsTotal - nNu;
    end
    
    [nuGuessW] = MinimizeResidualByVaryingMultipleFreqs(tsW, nuGuessW, nuMaxCW);
    
    if ~recon        
        bRemoved = true;
        while bRemoved
            [nuGuessW,bRemoved] = ConsolidateFreqs(nuGuessW,kNuConsolidate);
            [nuGuessW] = MinimizeResidualByVaryingMultipleFreqs(tsW, nuGuessW, nuMaxCW);
            nFreqsW = length(nuGuessW);
        end  
    end
    
    [cosPartW, sinPartW] = EstimateContainedSinusoids(tsW, nuGuessW);
    if nNu + nFreqsW < maxNFreqsTotal
        lastAbsDev = absDev;
        [tsW, absDev] = SubtractMultipleSinusoidsFromTS (tsW, cosPartW, sinPartW, nuGuessW);
    end
    
    for i = 1:nFreqsW
       nNu = nNu +1;
       nu(nNu) = nuGuessW(i);
       cosPart(nNu) = cosPartW(i);
       sinPart(nNu) = sinPartW(i);
       passNArr(nNu) = passN;
    end   
end

if recon
    [nu] = FindBestFreqIxsForRecon(nu, cosPart, sinPart, kNuConsolidate);
    bracket = [];
else
    [nu, bracket] = BracketTheFreqIxs(nu, cosPart, sinPart, passNArr);
end

end
%----------------
function [Xk] = dft(xn)
N = length(xn);         % Length of the input
Xk=zeros(1,N);          % initialize the output
i=sqrt(-1);
for k=0:N-1
    for n=0:N-1
        Xk(k+1)=(Xk(k+1)+(xn(n+1)*exp((-i)*2*pi*k*n/N)));
    end
end
 Xk = Xk./N;

 % -fold the DFT:  
 % The first value is DC and is not duplicated; 
 % if the number of the rest of the values is even, then we can fold them,  
 % if odd, then we can fold all but the middle value (which will be the
 % highest frequency.
 
 N = length(Xk);
 if floor(N/2) ~= N/2
     N = N + 1;
 end
 count = floor((N)/2); 
 x = zeros(1,count);
 x(1) = Xk(1);
 for i = 2:count
    j = 2*(count)+2 - i;
    x(i) = Xk(i) + Xk(j)'; 
 end
 Xk = x;
 
end
%-------------------
function [nu] = FindFreqsOfAmplitudePeaks(Xk, maxNAtOnce, nuMaxC)
%- Returns frequency indicies of the distinct local peaks in the amplitude spectrum.
%- Keep parallel arrays nu and maxAmpSqW, ordered by maxAmps, of the maxNAtOnce biggest ampltudes so far.
%- In:  Xk[0..nuMaxC]    complex DFT
%- Out: nu[1..nPeaks]    distinct frequency indices of peak amp., nPeaks <= maxNAtOnce
kOfPeakThreshold = 0.005;
maxAmpSqW = -(ones(1,maxNAtOnce));
lastAmpSq = -1;
prevLastAmpSq = -2;
nu = zeros(1,maxNAtOnce);
for mu = nuMaxC:-1:1
    ampSq = abs(Xk(mu))^2;
    if ampSq < lastAmpSq && lastAmpSq > prevLastAmpSq
        [nu, maxAmpSqW]=UpdateMaxAmpFreqs (nu, mu, lastAmpSq, maxNAtOnce, maxAmpSqW);
    end
    prevLastAmpSq = lastAmpSq;
    lastAmpSq = ampSq;
end

if lastAmpSq > prevLastAmpSq
   [nu, maxAmpSqW]=UpdateMaxAmpFreqs (nu, 0, lastAmpSq, maxNAtOnce, maxAmpSqW); 
end

% count peaks
ampSqThreshold = kOfPeakThreshold * maxAmpSqW(1);
if ampSqThreshold <= 0
    nPeaks = 0;
else
    nPeaks = maxNAtOnce;
    for i = 2:maxNAtOnce
        if maxAmpSqW(i) < ampSqThreshold
            nPeaks = i - 1;
            break
        end
    end   
end

for i = 1:nPeaks
    if nu(i) < 0 || nu(i) > nuMaxC
       %msg =  sprintf ('OFT>FindFreqsOfAmplitudePeaks nu(%d) (%f) is too high or negative',i,nu(i));
       error('OFT>FindFreqsOfAmplitudePeaks nu(%d) (%f) is too high or negative',i,nu(i));        
    end
    if i < nPeaks
       if maxAmpSqW(i) < maxAmpSqW(i+1)
          error('OFT>FindFreqsOfAmplitudePeaks amplitudes not sorted')
       end
    end
    for j = 1:i-1
        if nu(j) == nu(i)
           error('OFT>FindFreqsOfAmplitudePeaks two nu values are the same')
        end
    end           
end

end
%---------------------
function [nu, maxAmpSqW] = UpdateMaxAmpFreqs (nu, newNu, newAmp, maxNAtOnce, maxAmpSqW)
%- maxAmpSqW(1) is the biggest, maxAmpSqW(maxNAtOnce) is the smallest.

% return if smaller than the smallest
if newAmp <= maxAmpSqW(maxNAtOnce)
    return
end

for i = maxNAtOnce:-1:1
   if newAmp <=  maxAmpSqW(i)
       maxAmpSqW(i+1) = newAmp;
      nu(i+1) = newNu;
      return
   else
       if i == 1
           maxAmpSqW(1) = newAmp;
           nu(1) = newNu;
           return
       end
       maxAmpSqW(i) = maxAmpSqW(i-1);
       nu(i) = nu(i-1);
   end    
end
end
%---------------------
function [nuGuessW] = MinimizeResidualByVaryingMultipleFreqs(tsStage, nuGuessW, nuMaxCW)
%- Minimizes ResidualForMultipleFreqs of nFreqsW variables.
%- Based on the "powell" routine from "Numerical Recipes in C", Press et al, 1988, pages 314-315,
%  with func = ResidualForMultipleFreqs.
%- In:   nuGuessW[1..nFreqsW]             Initial frequencies.
%- Out:  nuGuessW[1..nFreqsW]             Frequencies that minimize ResidualForMultipleFreqs
%- Note: dirnSetW[1..nFreqsW][1..nFreqsW] Latest set of directions. Starts with unit directions.
%        fRet                             ResidualForMultipleFreqs(latest nuGuessW), the best acheived so far.
%- Finish when an iteration fails to reduce ResidualForMultipleFreqs by ktolMinMulti, relatively.

% Constants
ktolMinMulti = 0.00000005;

nNu = length(nuGuessW);
[absDev] = ResidualForMultipleFreqs(tsStage, nuGuessW);
dirnSetW = zeros(nNu,nNu);
nuGuessAtStOfIterW = zeros(1,nNu);
for i = 1:nNu
   nuGuessAtStOfIterW(i) = nuGuessW(i);
   for j = 1:nNu
      dirnSetW(i,j)=0; 
   end
   dirnSetW(i,i)=1;
end

iter = 0;
while true
    iter = iter + 1;
    absDevAtStOfIter = absDev;
    biggestDecrease = 0;
    dirnIxOfBiggestDecrease = 0;
    dirnW = zeros(1,nNu);
    for i=1:nNu
       for j = 1:nNu
          dirnW(j) = dirnSetW(j,i); 
       end
       lastAbsDev = absDev;
       [dirnW, nuGuessW, ~] = MinimizeResidualAlongOneDirn (dirnW, nuGuessW,nuMaxCW, tsStage);
       decrease = abs(lastAbsDev - absDev);
       if decrease > biggestDecrease
           biggestDecrease = decrease;
           dirnIxOfBiggestDecrease = i;
       end
    end
    
    % finish?
    if abs(absDevAtStOfIter - absDev) <= ktolMinMulti * (abs(absDevAtStOfIter) + abs(absDev)) * 0.5
       return 
    end
    
    if iter >= 30
       %msg = sprintf ('OFT.MinimizeResidualByVaryingMultipleFreqs: %d iterations', iter); 
       warning('OFT.MinimizeResidualByVaryingMultipleFreqs: %d iterations', iter)
       return
    end
    
    extrapolatedNuGuessW = zeros(1,length(nuGuessW));
    for j=1:length(nuGuessW)
       dirnW(j) = nuGuessW(j) -  nuGuessAtStOfIterW(j); % direction of change
       extrapolatedNuGuessW(j) = nuGuessW(j) + dirnW(j);
       if extrapolatedNuGuessW(j) < 0; extrapolatedNuGuessW(j)= 0; end
       if extrapolatedNuGuessW(j) > nuMaxCW; extrapolatedNuGuessW(j)= nuMaxCW; end
       nuGuessAtStOfIterW(i) = nuGuessW(i);
    end
    
    [fExtrap] = ResidualForMultipleFreqs(tsStage, extrapolatedNuGuessW);
    if fExtrap < absDevAtStOfIter
        u = (absDevAtStOfIter - absDev) - biggestDecrease;
        v = absDevAtStOfIter - fExtrap;
        T = 2 * (absDevAtStOfIter - 2 * absDev + fExtrap) * u^2  - biggestDecrease * v^2;
        if T < 0    % reverse direction
            [dirnW, nuGuessW, ~] = MinimizeResidualAlongOneDirn (dirnW, nuGuessW,nuMaxCW, tsStage);
            for j = 1:length(nuGuessW)
                dirnW(j)=dirnSetW(j, dirnIxOfBiggestDecrease);
            end
        end
        
    end 
end

end
%---------------------
function [absDev] = ResidualForMultipleFreqs(ts, nu)
%- The m-function, i.e. the function to be minimized in searching for the nFreqsW sinusoids whose removal
%  most reduces the absolute deviation of the time series.
[cosPart, sinPart] = EstimateContainedSinusoids(ts, nu);
%fRet = AbsDevAfterSubtractMultipleSinusoidsFromTS(tsStage, nuGuessW, cosPart, sinPart);

n = length(ts);
twoPiOn = 2 * pi / n;
twoPiNuOn = nu * twoPiOn;

absDev = 0;
for tau = n-1:-1:0
    x = ts(tau+1);
    for i = 1:length(nu)
        radians = twoPiNuOn(i) * tau;
        x = x - (cosPart(i) * cos(radians) + sinPart(i) * sin(radians));        
    end
    x = abs(x);
    absDev = absDev + x;    
end
end
%--------------------
function [dirnW, nuGuessW, minSizeOfResidual] = MinimizeResidualAlongOneDirn (dirnW, nuGuessW, nuMaxCW, ts)
%- Minimizes absolute deviation along the line (dirnW) in frequency space that goes thru the starting nuGuessW.
%  Updates nuGuessW with the minimum found, and rescales dirnW using the minimum found.
%- Based on the "linmin" routine from "Numerical Recipes in C", Press et al, 1988, page 316.
%- Initial value of bx chosen so that on first run after using DFT peaks to estimate nuGuessW,
%  BracketGuess1D will initally guess cx = -1, or 1 distance in frequency space. Should be guaranteed by
%  DFT peak checking.
%- In:  nuGuessW[1..nFreqsW]  is a point in argument space
%       dirnW   [1..nFreqsW]  is a direction in argument space
%- Out: nuGuessW[1..nFreqsW]  minimizes ResidualForMultipleFreqs along line from nuGuessW in direction dirnW
%       dirnW   [1..nFreqsW]  = output nuGuessW - input nuGuessW
%       minSizeOfResidual     = ResidualForMultipleFreqs(output nuGuessW), the best acheived

[ax,bx,cx,foundMinimum,xMin,minSizeOfResidual] = BracketGuess1D(nuGuessW, dirnW, nuMaxCW, ts);
if ~foundMinimum
   [xMin, minSizeOfResidual] = MinimizeFn1D(ax, bx, cx, nuGuessW, dirnW, ts); 
end

for i = 1:length(nuGuessW)
    dirnW(i) = dirnW(i) * xMin;
    nuGuessW(i) = nuGuessW(i) + dirnW(i);
    if nuGuessW(i) < 0; nuGuessW(i) = 0; end
     if nuGuessW(i) > nuMaxCW; nuGuessW(i) = nuMacCW; end   
end
end
%--------------------
function [ax,bx,cx,foundMinimum,xMin,fxMin] = BracketGuess1D(nuGuessW,dirnW, nuMaxCW, ts)
%- Concerns argument x to ResidualAlongADirn, which then constructs a multivariable argument to
%  ResidualForMultipleFreqs, namely nuGuessW(i) + x * dirnW(i). (So x = 0 corresponds to nuGuessW.)
%- In:  -   (Implicitly: nuGuessW, ResidualAlongADirn, FindValidLimits1D.)
%  Out: ax, bx, cx      ax < bx < cx
%                       ResidualAlongADirn(ax) >= ResidualAlongADirn(bx) <= ResidualAlongADirn(cx).
%       foundMinimum    True iff found minimum instead of brackets for a minimum (due to limits on x)
%       xMin            Only set if foundMinimum. Value of x for the minimum.
%       fxMin           Only set if foundMinimum. The minimum value.
%- Based on the "mnbrak" routine from "Numerical Recipes in C", Press et al, 1988, page 297,
%  with func(x) = ResidualAlongADirn(x).
%- Searches downhill from x = 0, starting with small steps (relative to the range of allowable values),
%  taking larger and larger steps. If run into a limit, then search around a bit more near that limit.
%  If the limit still appears to be the minimum, declare that as the minimum (because could not find
%  a bracket).

% Constants
kTiny = 1e-20;
kGold = 1.618034;   % the golden ratio
kGoldLess1OP = 0.618034;
kScale = 1000;
kLimitTol = kScale * 0.0005;

xMin = 0;       % xMin, fxMin is only valid if foundMinimum is true.
fxMin = 0;

foundMinimum = false;
[lo, hi]=FindValidLimits1D(nuGuessW,dirnW,nuMaxCW);
if lo >= hi
   foundMinimum = true;
   xMin = lo;
   [fxMin] = ResidualAlongADirn(xMin, nuGuessW, dirnW, ts);
   return
end

scalor = (hi - lo)/kScale;   %Multiply by scalor to give a range of kScale

% initial guess for a and b
ax = 0;                     
bx = kGoldLess1OP * scalor;
% ensure a and b are valid directional moves
if ax < lo; ax = lo; end
if ax > hi; ax = hi; end
if bx < lo; bx = lo; end
if bx > hi; bx = hi; end
if ax == bx
    bx = ax + kGoldLess1OP * scalor;
    if bx > hi; bx = hi; end
    if ax == bx
        bx = ax - kGoldLess1OP * scalor;
        if bx < lo; bx = lo; end
    end
end

% ensure a -> b is downhill
fa = ResidualAlongADirn(ax, nuGuessW, dirnW, ts);
fb = ResidualAlongADirn(bx, nuGuessW, dirnW, ts);
if fb > fa      %a <--> b
    temp = ax;
    ax = bx;
    bx = temp;
    temp = fa;
    fa = fb;
    fb = temp;
end

%initial guess for c
cx = bx + kGold * (bx - ax);
if cx < lo; cx = lo; end
if cx > hi; cx = hi; end
fc = ResidualAlongADirn(cx, nuGuessW, dirnW, ts);

% have a new c (i.e. cx and fc)
while true
   if fc - fb > 1e-10 * abs(fc); return; end    % exit if b -> is uphill
   
   % exit is c is up against a limit
   if cx == lo || cx == hi
       while abs(cx - bx) > kLimitTol * scalor
          u = 0.5 * (bx + cx);
          fu = ResidualAlongADirn(u, nuGuessW, dirnW, ts);
          if fc - fu > 1e-9 * abs(fc)
             ax = bx;
             bx = u;
             return
          end
          if fu > fb
              cx = u;
              return
          end
          ax = bx;
          bx = u;
       end
       foundMinimum = true;
       xMin = cx;
       fxMin = fc;
       return
   end
   
   % compute a new cx
   r = (bx - ax) * (fb - fc);    % compute by parabolic extrapolation
   q = (bx - cx) * (fb - fa);
   diff = q - r;
   nzDiff = abs(diff);
   if nzDiff < kTiny; nzDiff = kTiny; end
   if diff < 0; nzDiff = -nzDiff; end
   u = bx - ((bx - cx) * q - (bx - ax) * r) / (2 * nzDiff);
   if u < lo; u = lo; end
   if u > hi; u = hi; end
   ulim = bx + 100 * (cx - bx);     % max step
   if ulim < lo; ulim = lo; end
   if ulim > hi; ulim = hi; end
   takeDefault = false;
   
   % u in (bx, cx)
   if (bx - u) * (u - cx) > 0
       fu =  ResidualAlongADirn(u, nuGuessW, dirnW, ts);
       if fu < fc  %b->a, u->b, return
           ax = bx;
           bx = u;
           return
       else
           if fu > fb %u->c, return
               cx = u;
               return
           end
       end
   else
       % u in (cx, ulim)
       if (cx - u) * (u - ulim) > 0
           fu = ResidualAlongADirn(u, nuGuessW, dirnW, ts);
           if fu < fc % c->b, u->c
               bx = cx;
               cx = u;
               fb = fc;
               fc = fu;
               takeDefault = true;
           end
       else
           % ulim in (cx, u)
           if (u - ulim) * (ulim - cx) >= 0
               u = ulim;
               fu = ResidualAlongADirn(u, nuGuessW, dirnW, ts);               
           else
               takeDefault = true;
           end
       end
   end
   
   if takeDefault
      u = cx + kGold * (cx - bx);
      if u < lo; u = lo; end
      if u > hi; u = hi; end
      fu = ResidualAlongADirn(u, nuGuessW, dirnW, ts);
   end
   
   ax = bx;
   bx = cx;
   cx = u;
   fa = fb;
   fb = fc;
   fc = fu;
end
end
%--------------------
function [lo, hi] = FindValidLimits1D(nuGuessW, dirnW, nuMaxCW)
lo = -9e99;
hi = 9e99;
for i = 1:length(nuGuessW)
   dirni = dirnW(i);
   if dirni ~= 0
      loi = -nuGuessW(i) / dirni;
      hii = (nuMaxCW - nuGuessW(i))/dirni;
      if loi > hii
         temp = loi;
         loi = hii;
         hii = temp;
      end
      if loi > lo
          lo = loi;
      end
      if hii < hi
          hi = hii;
      end  
   end
   if lo > hi
       error('OFT.FindValidLimits lo > hi');
   end
   if lo > 0
       error('OFT.FindValidLimits lo > 0')
   end 
end
end

%---------------------
function [fxMin]= ResidualAlongADirn (x, nuGuessW, dirnW, ts)
%- 1-D function minimized by MinimizeFn1D, along line in direction dirnW.
%- Based on the "f1Dim" routine from "Numerical Recipes in C", Press et al, 1988, page 317.
nuTryW = zeros(1,length(nuGuessW));
for i=1:length(nuGuessW)
   nuTryW(i) = nuGuessW(i) + x * dirnW(i); 
end
fxMin = ResidualForMultipleFreqs(ts, nuTryW);
end
%---------------------
function    [xMin, minSizeOfResidual] = MinimizeFn1D(ax, bx, cx, nuGuessW, dirnW, ts)
%- Brent's method for finding the minimum value of a function of one variable.
%- Sets xMin to the value of x between ax and cx that mimimizes ResidualAlongADirn.
%- Requires either ax < bx < cx or ax > bx > cx, and
%    ResidualAlongADirn(ax) >= ResidualAlongADirn(bx) <= ResidualAlongADirn(cx).
%- Based on the "brent" routine from "Numerical Recipes in C", Press et al, 1988, pages 299 - 302,
%  with f = MinFn1. Derviatives too expensive in our case.
%- Finds xMin to with in about +-ktolMin1D (bracketing interval of about 2 * ktolMin1D).
%  Set ktolMin1D to about the square root of machine precision (see p.294 of "Numerical Recipes in C").
%  Double is 64 bit with ~15 decimal points of precision, so setting ktolMin1D to less than 3 * 10^-8
%  is a waste of time.
%- x always remains in [ax, cx] by design.

% Constants
ktolMin1D = 0.0002;
kCGold = 0.381966;

% [a,b] bracket solution, i.e. minimum lies in [a,b]
if ax < cx
    a = ax;
    b = cx;
else
    a = cx;
    b = ax;
end

e = 0;
x = bx;
w = bx;
v = bx;
fx = ResidualAlongADirn(x, nuGuessW, dirnW, ts);
fv = fx;
fw = fx;
d = 0;

for iter = 1:100    %limit the number of iterations to 100
   xm = 0.5 * (a+b);
   tol1 = ktolMin1D * abs(x) + 1e-10;
   tol2 = 2 * tol1;
   if abs(x - xm) <= (tol2 - 0.5*(b-a))
      xMin = x;
      minSizeOfResidual = fx;
      return
   end
   
   if abs(e) > tol1
      r = (x - w) * (fx - fv);   % construct a trial parabilic fit 
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2 * (q - r);
      if q > r; p = -p; end
      q = abs(q);
      eTemp = e;
      e = d;
      if abs(p) >= abs(0.5 * q * eTemp) | p <= (q * (a - x)) | p >= (q * (b-x))
          % take the golden step into the larger of the segments
          if x > xm
             e = a - x;
          else
              e = b - x;
          end
          d = kCGold * e;
      else
          % parabolic step
          d = p / q;
          u = x + d;
          if u-a < tol2 || b-u < tol2
             d = abs(tol1);
             if x > xm; d = -d; end
          end
      end
   else
        % take the golden step into the larger of the segments
        if x > xm
           e = a - x;
        else
            e = b - x;
        end
        d = kCGold * e;
   end
   % make the move
   if abs(d) > tol1
      u = x + d;
   else
       if d > 0
           u = x + abs(tol1);
       else
           u = x - abs(tol1);
       end
   end
   [fu] = ResidualAlongADirn(u, nuGuessW, dirnW, ts);
   
   % decrease the size of the bracket
   if fu <= fx
       if u >= x
           a = x;
       else
           b = x;
       end
       v = w;
       w = x;
       x = u;
       fv = fw;
       fw = fx;
       fx = fu;
   else
       if u < x
          a = u;
       else
           b = u;
       end
       if fu<=fw | w==x
          v = w;
          w = u;
          fv = fw;
          fw = fu;
       else
           if fu <= fv | v == x
               v = u;
               fv = fu;
           end
       end
   end 
end

warning('OFT.MinimizeFn1D did not find minimum in 100 iterations')

xMin = x;
minSizeOfResidual = fx;
end
%-----------------
function [nu]=FindBestFreqIxsForRecon(nu,cosPart,sinPart, kNuConsolidate)
global kMaxNFreqsAtOnceOFT
kOfPeakThreshold = 0.005;
%- Puts frequency indices with highest amplitudes in nu[1..nNu]. Generally will lower nNu.
%- Need cos and sin parts as inputs, rather than amplitudes, in order to consolidate frequency indicies.
%- In:  nu     [1..nNu]  Frequency indicies.
%       cosPart[1..nNu]  Destroyed by this routine.
%       sinPart[1..nNu]  Destroyed by this routine.
%       kNuConsolidate   distance between two nu values at which they will be combined
%- Out: nu     [1..nNu]  nNu = min(nNu after frequency consolidations, kMaxNFreqsAtOnceOFT)
nNu = length(nu);
newNNu = nNu;
for i = 1:nNu
   if nu(i) >= 0
      for j = i+1:nNu
         if abs(nu(i)-nu(j)) < kNuConsolidate
             cosPart(i) = cosPart(i)+cosPart(j);
             sinPart(i) = sinPart(i)+sinPart(j);
             cosPart(j) = 0;
             sinPart(j) = 0;
             nu(j) = -1;
             newNNu = newNNu -1;
         end
      end
   end
end

if newNNu > kMaxNFreqsAtOnceOFT; newNNu=kMaxNFreqsAtOnceOFT; end

maxAmpSq = -1;
for i = 1:nNu
   c = cosPart(i);
   s = sinPart(i);
   a = c^2+s^2;
   if maxAmpSq < a; maxAmpSq = a; end
   cosPart(i) = a;
   sinPart(i)=0;
end

threshold = kOfPeakThreshold * maxAmpSq;

for i = 1:newNNu
    maxAmpSq = -1;
    for j = 1:nNu
        if maxAmpSq < cosPart(j)
           maxAmpSq = cosPart(j);
           jMax = j;
        end
    end
    if cosPart(jMax) > threshold; sinPart(jMax) = -1; end
    cosPart(jMax) = -1;
end

newNNu = 0;
for i = 1:nNu
   if sinPart(i) == -1
      newNNu = newNNu + 1;
      if newNNu < i; nu(newNNu) = nu(i); end
   end
end
% nnU = newNNu;     % not returned
end
%------------------
function [nu, bracket] = BracketTheFreqIxs(nu, cosPart, sinPart, passNArr)
%- Sorts nu and creates the bracket array for the nu's, based on the associated amplitudes.
%- Brackets are intended for the MFT, in the next stage. The MFT will estimate contained sinusoids
%  using the frequency indicies in each bracket, separately from those in other brackets.
%  Bracketing prevents sinusods with tiny amplitudes being
%  put on the same basis as sinusoids with larger amplitudes when estimating the sinusoids in the
%  time series, i.e. do the large-amplitude sinusoids first, separately.
%- Oh the other hand, if FindFreqIxsInTS is finding similar sized amplitudes each pass, the ones
%  found in later passes should not be in the same bracket as ones found in earlier passes.
%- In:  nu      [1..nNu]         Frequency indicies.
%       cosPart [1..nNu]
%       sinPart [1..nNu]
%       passNArr[1..nNu]         Pass number in FindFreqIxsInTS that the sinusoid was found in.
%- Out: nu      [1..nNu]         Order may have changed. nNu unchanged.
%       bracket [1..nBrackets]   nu[1..bracket(1)] are the freq. indicies in the first bracket, etc.

% Constants
kRatioHi = 1 / 10;
kRatioLo = 1 / 50;

nNu = length(nu);
ampSq = zeros(1,nNu);
for i = 1:nNu
   c = cosPart(i);
   s = sinPart(i);
   ampSq(i) = -(c^2+s^2);    % -ve to get sort order correct
end

ixArr = zeros(1,nNu);
for i=1:nNu
    ixArr(i) = i;
end

bracket = zeros(1,nNu);
nBrackets = 0;
passesEnIx = 0;
maxPassN = 0;
while passesEnIx < nNu
   maxPassN = maxPassN + 2; % Normally 2, maybe 3 or 4?. Need to set to 5 for ACRIM to get second low freq
                            % high amp sinusoid consolidated with first one.
   passesStIx = passesEnIx + 1;
   for i = passesStIx:nNu
      if passNArr(i) < maxPassN
         passesEnIx = i;
      else
          break  
      end
   end
end

% sort an index into ampsq
[ixArr] = SortArrayOfIndexesToDoubles(ixArr, ampSq, passesStIx, passesEnIx);

% use the sorted index to order nu amd ampsq
tempN = zeros(1,passesEnIx-passesStIx);
tempA = zeros(1,passesEnIx-passesStIx);
for i = passesStIx:passesEnIx
   tempN(i) = nu(i);
   tempA(i) = -ampSq(i); % return to +ve
end
for i = passesStIx:passesEnIx
   nu(i) = tempN(ixArr(i));
   ampSq(i) = tempN(ixArr(i));
end

% Loop by bracket (within a passN-group)
nNuStIx = passesStIx;
while nNuStIx <= passesEnIx
   thresholdHi = ampSq(nNuStIx) * kRatioHi;
   thresholdLo = ampSq(nNuStIx) * kRatioLo;
   
   for i = nNuStIx:passesEnIx
      if ampSq(i) < thresholdHi; break; end
      hiIx = i;
   end
   
   for i = hiIx:passesEnIx
      if ampSq(i) < thresholdLo; break; end
      loIx = i;
   end
   
   maxGap = -1;
   maxIx = passesEnIx;
   if loIx < nNu; enIx = loIx; else; enIx = loIx - 1; end
   for i = hiIx:enIx
      g = ampSq(i) / ampSq(i+1); %g >=1
      if maxGap < g
         maxGap = g;
         maxIx = i;
      end
   end
   nNuInBracket = maxIx - nNuStIx + 1;  % bracked ends in greatest gap
   nBrackets = nBrackets+1;
   bracket(nBrackets) = nNuInBracket;
   nNuStIx = nNuStIx + nNuInBracket;
end
   
end



