function [macTransmit,macDensity] = MacularTransmittance(S,species,source,fieldSizeDegrees)
% [macTransmit,macDensity] = MacularTransmittance(S,[species],[source],[fieldSizeDegrees])
%
% Return an estimate of the transmittance of the macular pigment transmittance
% as a function of wavelength.
%
% Allowable species:
%   Human (Default)
%
% Allowable sources:
%   CIE (Default)            - CIE 170-1:2006 values.
%   Bone                     - From Bone et al.  See CVRL database.
%   WyszeckiStiles           - From W&S, Table 2(2.4.6), p. 112.
%   Vos                      - From Vos.  See CVRL database.
%   None                     - Unity transmittance.
%
% For the CIE option, can pass fieldSizeDegrees [Default 2 degrees].
% This was buggy until the version of 5/8/12.
%
% The Bone values that we use a the basis for this calculation
% match those in  CIE 170-1:2006, Table 6.4 for a 2-degree observer.
% 
% The answer is returned in a row vector.  This function
% depends on data contained in directory
% PsychColorimetricData:PsychColorimetricMatFiles.
%
% 7/8/03  dhb  Made this a separate function.
% 7/11/03 dhb  Species arg, change name.
% 7/23/03 dhb  Change default.
% 7/26/03 dhb  Extend functions, rather than zero truncate.
% 8/12/11 dhb  Fixed default to match comments.
%         dhb  Add CIE option and made it default.
%         dhb  For CIE, can pass field size
%         dhb  Also return density
% 8/13/11 dhb  Linearly extrapolate read functions outside of range.
% 5/8/12  dhb  Fixed two bugs.  First, peak optical density correction is
%              multiplicative rather than additive.  Second, there was
%              an operator precedence grouping error in the computation
%              of the correction factor.
% 5/8/12  dhb  Removed comment that we can't reproduce CIE tabular 10 deg values.
% 9/17/12 dhb  Return density for 'None' case as well.
% 8/9/13  dhb  More consistent returning of density for 'None' case.
% 8/11/13 dhb  Try to make dimensions of returned density match those of returned transmittance.

% Default
if (nargin < 2 || isempty(species))
	species = 'Human';
end
if (nargin < 3 || isempty(source))
	source = 'CIE';
end
if (nargin < 4 || isempty(fieldSizeDegrees))
    fieldSizeDegrees = 2;
end

% Load correction for macular pigment density
switch (species)
	case 'Human',
		switch (source)
			case 'None',
				macTransmit = ones(S(3),1)';
                macDensity = zeros(S(3),1)';
			case 'WyszeckiStiles',
				load den_mac_ws;
				macDensity = SplineSrf(S_mac_ws,den_mac_ws,S,2)';
				macTransmit = 10.^(-macDensity);
			case 'Vos',
				load den_mac_vos;
				macDensity = SplineSrf(S_mac_vos,den_mac_vos,S,2)';
				macTransmit = 10.^(-macDensity);
			case 'Bone',
				load den_mac_bone;
				macDensity = SplineSrf(S_mac_bone,den_mac_bone,S,2)';
				macTransmit = 10.^(-macDensity);
            case 'CIE'
                load den_mac_bone;
				macDensity = SplineSrf(S_mac_bone,den_mac_bone,S,2)';
                
                % Adjust for field size by adjusting peak optical density.
                % This is a multiplicative adjustment, which is pretty
                % easy to derive from first principles.  See Rodieck, pp. 443-445.
                % Our Bone values have a peak of 0.35, but the CIE formula for peak
                % density produces coefficients for a peak of 1.  We handle this by
                % dividing the Bone values by their peak density and then applying the
                % formula.  And, the peak density of 0.35 is in fact obtained by
                % applying the CIE formula for a 2 degree field.  For tabulated
                % values to check this against, see CIE 170-1:2006, p.  29, Table 6.4.
                densityAdjustFieldSize = 0.485*exp(-fieldSizeDegrees/6.132) / (0.485*exp(-2/6.132));
                macDensity = macDensity*densityAdjustFieldSize;
                macDensity(macDensity < 0) = 0;
                
				macTransmit = 10.^(-macDensity);
			otherwise,
				error('Unsupported macular pigment density estimate specified');
		end

	otherwise,
		switch (source)
			case ('None'),
				macTransmit = ones(S(3),1)';
                macDensity = zeros(S(3),1)';
			otherwise,
				error('Unsupported species specified');
        end
end

