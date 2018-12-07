function [bb] = SolveLinearEquations (aa, n, bb)
%- Solves the N linear equations aa * xx = bb for xx, returns xx in bb vector.
%- In: aa[1..N][1..N]   Array, N by N.
%       bb[1..N]         Vector, 1 by N.
%- Out: bb[1..N]         Vector, 1 by N, solution xx to aa * xx = bb.
%=========================================================================%
% The below code is not part of the function call.  It returns handles to
% all the local subfunctions for the purpose of unit testing if the value
% of the "aa" input (normally a double) is the string '-test'
if ischar(aa) && strcmp(aa,'-test')
    bb = localfunctions;
    return
end
%================================================================+========%

% This is the function call:
[aa, indx, ~] = ComputeLUDecomposition (aa, n);
[bb] = SolveLinearEquationsWithLUMatrix (aa, n, indx, bb);
end

%*************************************************************************%
% Local Functions
function [aa, indx, d] = ComputeLUDecomposition (aa, n)
global flags
%- Converts N by N matrix aa into LU form.
%- Based on the "ludcmp" routine from "Numerical Recipes in C", Press et al, 1988, page 43.
%- In:  A[1..N][1..N]   Array, N by N.
%- Out: A[1..N][1..N]   Array, N by N, LU decomposition of a row-wise permutation of input aa.
%       indx[1..N]      Records the row permutations effected by the partial pivoting.
%       d               +1 or -1 if number of row interchanges was even or odd, respectively.

d = 1;
vv = zeros(1,n);
for i = 1:n
    big=0;
    for j = 1:n
        temp = abs(aa(i,j));
        if big < temp
            big = temp;
        end
    end
    if big == 0
        error('ComputeLUDecomposition failed:  aa matrix has a row of all 0')
    end
    vv(i) = 1;
end

% Crouts method (loop over the columns)
indx = zeros(1,n);
for j = 1:n
    for i = 1:j-1
        sum = aa(i,j);
        for k = 1:i-1
            sum = sum - aa(i,k)*aa(k,j);
        end
        aa(i,j)=sum;
    end

    % Search for the largest pivot element
    big = 0;
    for i = j:n
        sum = aa(i,j);
        for k = 1:j-1
            sum = sum - aa(i,k) * aa(k, j);
        end
        aa(i,j) = sum;
        dum = vv(i) * abs(sum);
        if big <= dum
            big = dum;
            imax = i;
        end
    end
    
    % interchange rows
    if j ~= imax
        for k = 1:n
            dum = aa(imax, k);
            aa(imax,k) = aa(j,k);
            aa(j,k) = dum;
        end
        d = -d;
        vv(imax) = vv(j);
    end
    
    % divide by the pivot element
    indx(j) = imax;
    if aa(j,j) == 0
        aa(j,j) = 1e-20;
    end
    if j ~= n
        dum = 1/aa(j,j);
        %If the pivot element is 0 the matrix is singular (at least
        %to the precision of the algorithm). Sometimes better to
        %substitute a tiny value.
        for i = j+1:n
            aa(i,j) = dum * aa(i,j);
        end
    end
end

% if a 'flags' object exists, then flag for low matrix elements that could
% cause problems for the OFT or MFT
if isobject(flags)
    flags.FlagLowElementsInMatrix(aa, 'l' , 'L');
end
end

function [bb] = SolveLinearEquationsWithLUMatrix (aa, n, indx, bb)
%- Solves the n linear equations aa * xx = bb, where aa is in LU form as per ComputeLUDecomposition.
%  Returns xx in the bb vector.
%- Based on the "lubksb" routine from "Numerical Recipes in C", Press et al, 1988, page 44.
%- In:  aa[1..n][1..N]  Array, N by N.
%       indx[1..N]      The row permutations as per ComputeLUDecomposition.
%       bb[1..N]        Vector, 1 by N.
%- Out: bb[1..N]        Vector, 1 by N, solution xx to aa * xx = bb.

ii = 0;
for i = 1:n
    ip = indx(i);
    sum = bb(ip);
    bb(ip) = bb(i);
    if ii ~= 0
        for j = ii:i-1
            sum = sum - aa(i,j) * bb(j);
        end
    else
        if sum ~= 0
            ii = 1;
        end
    end
    bb(i) = sum;
end

for i = n:-1:1
    sum = bb(i);
    for j = i+1:n
        sum = sum - aa(i,j) * bb(j);
    end
    bb(i) = sum / aa(i,i);
end
end


        
        