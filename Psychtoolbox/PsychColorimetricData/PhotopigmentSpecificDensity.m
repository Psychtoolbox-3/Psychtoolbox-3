function densities = PhotopigmentSpecificDensity(receptorTypes,species,source)
% densities = PhotopigmentSpecificDensity(receptorTypes,[species],[source])
%
% Return estimates of photopigment specific densities.
%
% Allowable receptor types depend on species and source, but the general
% list is:
% 	SCone, MCone, LCone, FovealSCone, FovealMCone, FovealLCone, Rod.
%
% The type argument may be a single string or a cell array of strings.  If it
% is an array, a column vector of values is returned.
%
% The foveal version of cone types is sensible only for primates.  Not all
% estimate sources support all receptor types.
%
% Note that the following three numbers are overdetermined: photopigment
% specific density (sd), photopigment axial density (ad), and outer segment
% length osl.  In particular, ad = sd*osl.  Depending on the measurement
% method, different sources provide different pairs of these numbers.
% We have attempted to enforce this consistency in the set of routines
% PhotopigmentSpecificDensity, PhotopigmentAxialDensity, and PhotoreceptorDimensions.
% That is to say, for the same source, species, and cone type, you should get
% a consistent triplet of numbers. 
%
% Supported species:
%		Human (Default), GuineaPig.
%
% Supported sources:
% 	Rodieck (Human) (Default).
%   Bowmaker (GuineaPig).
%   Generic
%   None (returns empty matrix as value)
%
% The Generic source returns a single number for all species and receptor types.
% This number is 0.015 /um.
%
% 7/11/03  dhb  Wrote it.
% 8/9/13   dhb  Comment clean up, allow 'None' to return empty as the value.

% Fill in defaults
if (nargin < 2 || isempty(species))
	species = 'Human';
end
if (nargin < 3 || isempty(source))
	source = 'Rodieck';
end

% Fill in specific density according to specified source
if (iscell(receptorTypes))
	densities = zeros(length(receptorTypes),1);
else
	densities = zeros(1,1);
end
for i = 1:length(densities)
	if (iscell(receptorTypes))
		type = receptorTypes{i};
	elseif (i == 1)
		type = receptorTypes;
	else
		error('Argument receptorTypes must be a string or a cell array of strings');
	end

	switch (source)
        case {'None'}
            densities = [];
		case {'Generic'}
			densities(i) = 0.015;

		case {'Bowmaker'}
			switch (species)
				case 'GuineaPig',
					% Specific density measured in Guinea Pig.
					% In Parry JWL, Bowmaker JK, 2002 paper
					% transverse measurement: OD in table Table 1
					% in the order of [M, S, Rod]
					%	
					% Note communication with Dr. J. Bowmaker:
					%	"The OD values measured in msp are always somewhat suspect
					%	snd this is more so for mammalian cones - the outer segments
					%	tend to collapse and because of their small size it is difficult
					% to ensure that the measuring beam is totally within the outer segment.
					%	
					%	"Because we use light polarized perpendicular to the long axis of
					%	the os, you can simply multiply the specific density by the length
					%	of the OS to give an axial OD.  The values will be in the right order of
					% magnitude, but could be out be a factor of 2 or 3."
					%
					% Also see http://cvrl.ucl.ac.uk, section on photopigment optical density.
					
					% Start with the transverse optical densities
					% and use to calculate an estimate of the specific densities.  The
					% diameter in the transverse measurement is roughly the pathlength,
					% so we divide by that to get the specific density.
					switch (type)
						case {'MCone'}
							transverseOD = 0.009;
							OSdiameter = PhotoreceptorDimensions('MCone','OSdiam','GuineaPig','SterlingLab');
							densities(i) = transverseOD/OSdiameter;
						case 'SCone',
							transverseOD = 0.008;
							OSdiameter = PhotoreceptorDimensions('SCone','OSdiam','GuineaPig','SterlingLab');
							densities(i) = transverseOD/OSdiameter;
						case 'Rod'
							transverseOD = 0.021;
							OSdiameter = PhotoreceptorDimensions('Rod','OSdiam','GuineaPig','SterlingLab');
							densities(i) = transverseOD/OSdiameter;
						otherwise,
							error(sprintf('Unsupported receptor type %s for %s estimates in %s',type,source,species));
					end
			otherwise,
				error(sprintf('%s estimates not available for species %s',source,species));
			end

	
		case {'Rodieck'}
			switch (species)
				case {'Human'},
					% Rodieck, The First Steps in Seeing, p. 472 gives a value for
					% axial specific density for rods and cones as about 0.015 /um.
					% Here we compute these numbers from his estimate of outersegment
					% length and axial optical density.  This enforces consistency
					% across different ways of getting the same number.
					switch (type)
						case {'FovealLCone','FovealMCone','FovealSCone','Rod'}
							OSlength = PhotoreceptorDimensions(type,'OSlength',species,source);
							axialOpticalDensity = PhotopigmentAxialDensity(type,species,source);
							densities(i) = axialOpticalDensity ./ OSlength;
						otherwise,
							error(sprintf('Unsupported receptor type %s for %s estimates in %s',type,source,species));
					end
				otherwise,
					error(sprintf('%s estimates not available for species %s',source,species));
			end	

		otherwise
			error(sprintf('Unknown source %s for specific density estimates',source));
	end
end
