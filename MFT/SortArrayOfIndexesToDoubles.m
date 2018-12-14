function [indexArr]=SortArrayOfIndexesToDoubles(indexArr, valueArr, minIx, maxIx)
% Takes in an array of indices to an array of double type values
% Sorts the indices so they are in order of increasing values in the
% valueArray
if minIx >= maxIx; return; end     % indexArr contains 0 or 1 item

i = floor((maxIx - minIx + 1) * rand + minIx);      % choose a random mid value
midValIx = indexArr(i);
midVal = valueArr(indexArr(i));

indexArr(i) = indexArr(minIx);
lo = minIx;
hi = maxIx;

while true
    while valueArr(indexArr(hi)) >= midVal
        hi = hi-1;
        if hi <= lo; break; end
    end
    if hi <= lo
       indexArr(lo) = midValIx;
       break
    end
    indexArr(lo) = indexArr(hi);
    lo = lo + 1;
    while valueArr(indexArr(lo)) < midVal
        lo = lo+1;
        if lo >= hi; break; end
    end
    if lo >= hi
       lo = hi;
       indexArr(hi) = midValIx;
       break
    end
    indexArr(hi) = indexArr(lo);
end
    [indexArr]=SortArrayOfIndexesToDoubles(indexArr, valueArr, minIx, lo-1);
    [indexArr]=SortArrayOfIndexesToDoubles(indexArr, valueArr, lo + 1, maxIx);
end
