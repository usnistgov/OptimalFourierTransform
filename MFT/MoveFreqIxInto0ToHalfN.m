function [nuVal,mult] = MoveFreqIxInto0ToHalfN(n,nuVal)
%- Moves frequency index nu (any real number) into [0, N/2], such that the sinsoid is unchanged.
%  Sets mult to -1 or +1 if sine components need a change of sign or not (respectively).
%  Move nu into [0, N-1]
if nuVal >= 0
    mult = 1;
else
    nuVal=-nuVal;
    mult = -1;
end
if nuVal >= n
    nuVal = mod(nuVal,n);
end
if 2*nuVal > n
    nuVal = n-nuVal;
    mult = -mult;
end
end