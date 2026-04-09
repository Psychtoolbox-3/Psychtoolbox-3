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
%
% 7/11/03  dhb  Wrote it.
% 7/16/03  dhb  Add StockmanSharpe.

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
		if (~exist('StockmanRiderNomogram.m','file'))
			error(['PhotopigmentNomogram: StockmanRider nomogram requested but ' ...
			       'StockmanRiderNomogram.m was not found on the MATLAB path. ' ...
			       'Add the BrainardLabToolbox StockmanRiderFormulae directory to your path.']);
		end
		T_absorbance = ...
			StockmanRiderNomogram(S,lambdaMax);
	otherwise
		error('Unknown source for photopigment nomogram');
end
