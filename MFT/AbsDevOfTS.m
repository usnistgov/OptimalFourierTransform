function [absDev] = AbsDevOfTS_Fn(ts, n)

absDev = 0;
for tau= n-1:-1:0
    x = abs(ts(tau+1));
    absDev = absDev + x;
end
end