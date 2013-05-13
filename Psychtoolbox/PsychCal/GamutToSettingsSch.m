function settings = GamutToSettingsSch(gammaInput, gammaTable, gamut)
% settings = GamutToSettingsSch(gammaInput, gammaTable, gamut)
%
% Find the best device settings to produce
% the passed linear device coordinates.
%
% This version works by linear interpolation on the fit values,
% with x as the output values and f(x) as the input values.
%  
% The passed coordinates should be in the range [0,1].
% The returned settings also run from [0,1], but after
% inversion of the device's gamma measurements.
%
% 9/26/93    dhb   Added calData argument.
% 10/19/93   dhb   Allow gamma table dimensions to exceed device settings.
% 11/11/93   dhb   Update for new calData routines.
% 8/4/96     dhb   Update for stuff bag routines.
%            dhb   Pulled out as a subroutine.
% 8/21/97	 dhb   Convert for structures.
% 11/16/06   dhb   Adjust for [0,1] world.
%            dhb   No more values return because we can't get at it in the [0,1] world
%            dhb   Pass input x values as well as y values.

% Check dimensions and table sizes
[m,n] = size(gamut);
[mg,ng] = size(gammaTable);
if (m > ng)
  error('Mismatch between device coordinate dimensions and gamma table');
end

% Use a search routine to find the best gamma function
settings = zeros(m,n);
% values = zeros(m,n);
for i = 1:m
  settings(i,:) = SearchGammaTable(gamut(i,:), gammaInput, gammaTable(:,i));
end
