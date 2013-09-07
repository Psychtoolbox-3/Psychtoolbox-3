function [peakWls,peakVals] = FindCmfPeaks(S_xxx,T_xxx)
%  [peakWls,peakVals] = FindCmfPeaks(S_xxx,T_xxx)
% 
% Return the peak wavelengths and peak values of a set of color matching functions (T_xxx).
%
% 8/12/13  dhb  Wrote it.

wls = SToWls(S_xxx);
for i = 1:size(T_xxx,1)
    [peakVals(i),index] = max(T_xxx(i,:));
    peakWls(i) = wls(index(1));
end

end