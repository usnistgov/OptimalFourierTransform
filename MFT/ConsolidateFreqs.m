function [nu_MFT,bRemoved] = ConsolidateFreqs(nu_MFT,kNuConsolidate)
% Consolidates converging frequencies (those with indices very close together)
nu = nu_MFT;
nRemoved = 0;
bRemoved = false;
nFreqs = length(nu);

for i = 1:nFreqs 
    if nu(i) >= 0
        for j = i+1:nFreqs
            if abs(nu(i)-nu(j)) < kNuConsolidate
                nu(i) = (nu(i) + nu(j)) * 0.5;
                nu(j) = -1;
                nRemoved = nRemoved+1;
                bRemoved = true;
            end
        end
    end
end

if nRemoved > 0
    j = 0;
    for i = 1:nFreqs-nRemoved
        j = j+1;
        if nu(j) >= 0
            if i ~= j
                nu(i) = nu(j);
            end
        else
            while nu(j) < 0
                j = j+1;
                if j > nFreqs
                    break
                end
                nu(i) = nu(j);
            end
            if j > nFreqs; break; end
        end
    end
    nFreqs = nFreqs - nRemoved;
end
nu_MFT = nu(1:nFreqs);
end             