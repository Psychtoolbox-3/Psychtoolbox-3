
function T_absorbance = PhotopigmentNomogram(S,lambdaMax,source)
% T_absorbance = PhotopigmentNomogram(S,lambdaMax,[source])
%
% Compute normalized absorbance according to various
% nomograms.  This is basically a wrapper routine.
%   Baylor
%   Dawis
%   Govardovskii (Default)
%   Lamb
%   StockmanSharpe
%   StockmanRider
%
% 7/11/03  dhb  Wrote it.
% 7/16/03  dhb  Add StockmanSharpe.
% 4/10/26  Add  Add StockmanRider.
%          Systematially edit all the individual functions so that if wl
%          spacing is <= 1, they normalize the return to peak of 1.

if (nargin < 3 || isempty(source))
	source = 'Govardovskii';
end

switch (source)
	case {'Govardovskii'}
		T_absorbance = ...
			GovardovskiiNomogram(S,lambdaMax);
	case {'Dawis'}
		T_absorbance = ...
			DawisNomogram(S,lambdaMax);
	case {'Baylor'}
		T_absorbance = ...
			BaylorNomogram(S,lambdaMax);
	case 'Lamb'
		T_absorbance = ...
			LambNomogram(S,lambdaMax);
	case 'StockmanSharpe'
		T_absorbance = ...
			StockmanSharpeNomogram(S,lambdaMax);
	case 'StockmanRider'
		T_absorbance = ...
			StockmanRiderNomogram(S,lambdaMax);
	otherwise
		error('Unknown source for photopigment nomogram');
end
