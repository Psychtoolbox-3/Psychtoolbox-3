function settings = GamutToSettingsExhaust(gammaInput, gammaTable, gamut)
% settings = GamutToSettingsExhaust(gammaInput, gammaTable, gamut)
%
% Find the best device settings to produce
% the passed linear device coordinates.
%
% This version works by exhaustive search of the fit gamma table,
% and returns the value in that table that produces output closest
% to desired.
%
% The passed coordinates should be in the range [0,1].
% The returned settings also run from [0,1], but after
% inversion of the device's gamma measurements.
%
% 5/26/12  dhb  Wrote it.

% Check dimensions and table sizes
[m,n] = size(gamut);
ng = size(gammaTable,2);
if (m > ng)
  error('Mismatch between device coordinate dimensions and gamma table');
end

% Use a search routine to find the best gamma function
settings = zeros(m,n);
for i = 1:m
    for j = 1:n
        [~,index] = min(abs(gammaTable(:,i)-gamut(i,j)));
        settings(i,j) = gammaInput(index);
    end
end
