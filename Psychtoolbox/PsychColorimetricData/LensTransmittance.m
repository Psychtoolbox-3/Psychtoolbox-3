function [lensTransmit,lensDensity] = LensTransmittance(S,species,source,ageInYears,pupilDiameterMM)
% [lensTransmit,lensDensity] = LensTransmittance(S,[species],[source],[ageInYears],[pupilDiameterMM])
%
% Return an estimate of the transmittance of the lens.
%
% Allowable species:
%   Human (Default)
%
% Allowable sources:
%   StockmanSharpe (Default) - Stockman, Sharpe, & Fach (1999).
%   CIE                      - Formula from CIE 170-1:2006.
%   WyszeckiStiles           - W&S, Table 1(2.4.6), p. 109.  First data set in table.
%   None                     - Unity transmittance.
%
% The answer is returned in a row vector.  This function
% depends on data contained in directory
% PsychColorimetricData:PsychColorimetricMatFiles.
%
% The CIE source will take age in years and pupil size in mm.
% The default age is 32.  The acceptable range is 20-80.
% The default pupil size is 3 mm.  Sizes less than 3 are treated
% as 3.  Sizes greater than 7 are treated as 7.  The interpolation
% between 3 and 7 doesn't appear to be part of the standard but
% I coded it in anyway.
% 
%
% 7/8/03  dhb  Made this a separate function.
% 7/11/03 dhb  Species arg, change name.
% 7/23/03 dhb  Add Stockman estimate.
% 7/26/03 dhb  Extend functions, rather than zero truncate.
% 8/12/11 dhb  Start to write CIE version.  Return lensDensity too.
%         dhb  Finish. Add pupil size.
% 8/13/11 dhb  Linearly extrapolate read functions outside of range.
% 9/17/12 dhb  Return density for 'None' case as well.
% 8/9/13  dhb  More consistent returning of density for 'None' case.
% 8/11/13 dhb  Try to make dimensions of returned density match those of returned transmittance.

% Default
if (nargin < 2 || isempty(species))
	species = 'Human';
end
if (nargin < 3 || isempty(source))
	source = 'StockmanSharpe';
end
if (nargin < 4 || isempty(ageInYears))
	ageInYears = 32;
end
if (nargin < 5 || isempty(pupilDiameterMM))
	pupilDiameterMM = 3;
end

% Load correction for lens density
switch (species)
	case 'Human',
		switch (source)
			case 'None',
				lensTransmit = ones(S(3),1)';
                lensDensity = zeros(S(3),1)';
			case 'WyszeckiStiles',
				load den_lens_ws;
				lensDensity = SplineSrf(S_lens_ws,den_lens_ws,S,2)';
				lensTransmit = 10.^(-lensDensity);
			case 'StockmanSharpe',
				load den_lens_ssf;
				lensDensity = SplineSrf(S_lens_ssf,den_lens_ssf,S,2)';
				lensTransmit = 10.^(-lensDensity);
            case 'CIE'
                % Load CIE age dependent and age independent components
                load den_lens_cie_1
                load den_lens_cie_2
                lensDensity1 = SplineSrf(S_lens_cie_1,den_lens_cie_1,S,2)';
                lensDensity2 = SplineSrf(S_lens_cie_2,den_lens_cie_2,S,2)';
                
                % Combine them according to age using CIE formulae
                if (ageInYears < 20)
                    error('Specified age must be 20 or older');
                elseif (ageInYears <= 60)
                    lensDensity = lensDensity1*(1+0.02*(ageInYears-32))+lensDensity2;
                elseif (ageInYears <= 80)
                    lensDensity = lensDensity1*(1.56+0.0667*(ageInYears-60))+lensDensity2;
                else
                    error('Specified age must be 80 or younger');
                end
                
                % This is the answer for pupil size <= 3mm. But, can
                % correct for pupil diameter.
                %
                % The CIE report says for pupil size > 7 mm, multiply the density
                % values by 0.86207.  See note at Table 6.10.
                % So that's what we do for that case.  It is not clear that this
                % is actually part of the standard, and the effect is very small.
                %
                % It seems silly not to return anything for values between
                % 3 and 7 mm, so for this case I linearly interpolate the
                % factor.
                if (pupilDiameterMM > 3 && pupilDiameterMM < 7)
                    factor = (7-pupilDiameterMM)/4 + 0.86207*(pupilDiameterMM-3)/4;
                    lensDensity = factor*lensDensity;
                elseif (pupilDiameterMM >= 7)
                    lensDensity = 0.86207*lensDensity;
                end
                lensTransmit = 10.^(-lensDensity);
          
			otherwise,
				error('Unsupported lens density estimate specified');
		end

	otherwise,
		switch (source)
			case ('None'),
				lensTransmit = ones(S(3),1)';
                lensDensity = zeros(S(3),1)';
			otherwise,
				error('Unsupported species specified');
		end
end
