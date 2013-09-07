function eyeLengthMM = EyeLength(species,source)
% eyeLengthMM = EyeLength(species,source)
%
% Return the length of the eye in mm.  Length is taken as distance
% between nodal point and fovea.  For foveal stimuli, this length
% may be used to convert between degrees of visual angle and mm/um
% of retina.  Since the nodal point isn't at the center of eye,
% the same conversion doesn't work for extrafoveal stimuli.
%
% Supported species:
%		Human (Default), Rhesus, Dog
%
% Supported sources:
%   LeGrand (Human, 16.6832 mm, Default)
%   Rodieck (Human, 16.1 mm)
%   Gulstrand (Human, 17 mm)
%   PerryCowey (Rhesus)
%   Packer (Rhesus)
%   PennDog (Dog)
%   None
%
% Passing a numeric value as source returns that value as the
% estimate, independent of species.  This is a hack that allows
% some debugging.
%
% Passing a string that is a number (that is, something that str2num
% will turn into a number) as the source will cause that number to
% be returned as the eye length.  This is a bit redundant with the
% numeric option above.
%
% Passing None is appropriate as an error check -- if a calculation
% uses the eye length when none is passed, NaN's will show up in
% the answer.
% 
% Finally, if you pass a decimal number as a string, this value
% will be returned.  Useful for passing arbitrary numbers through
% routines that rely on this function.
%
% 7/15/03  dhb  Wrote it.
% 2/27/13  dhb  Added 17 mm option
% 3/1/13   dhb  Changed '17' to 'Gulstrand'
%          dhb  Added option of passing a number as a string.
% 4/12/13  dhb  Fix bug introduced in previous update, which may have
%               been Matlab version dependent.
% 8/9/13   dhb  A few more comments.

% Fill in defaults
if (nargin < 1 || isempty(species))
	species = 'Human';
end
if (nargin < 2 || isempty(source))
	source = 'LeGrand';
end

% Check if a direct numerical value was passed as a string.  If so,
% ignore everything else and return it.
if (ischar(source))
    directVal = str2num(source);
    % If directVal is empty, don't do anything here and
    % go on to process the string.
    if (~isempty(directVal))
        eyeLengthMM = directVal;
        return;
    end
else
    % The source was passed as something not a string.  We assume
    % that it is the numerical values of the eye length in mm, and
    % return it.
    eyeLengthMM = source;
	return;
end

% Handle case where 'None' is passed.
if (streq(source,'None'))
	eyeLengthMM = NaN;
	return;
end

% Fill in length according to species, source
switch (source)
  % I took the LeGrand number from Wyszecki and Stiles' description of
  % his model eye.
	case {'LeGrand'}
		switch (species)
			case {'Human'}
				eyeLengthMM = 16.6832;
			otherwise,
				error(sprintf('%s estimates not available for species %s',source,species));
			end

	% Rodieck's standard observer, from Appendix B, The First Steps of Seeing.
	case {'Rodieck'}
		switch (species)
			case {'Human'}
				eyeLengthMM = 16.1;
			otherwise,
				error(sprintf('%s estimates not available for species %s',source,species));
        end
    
    % Gulstrand model eye, as described by Delori et. al., 2007, JOSA A, 24, pp. 1250-1265.
    case {'Gulstrand'}
        switch (species)
			case {'Human'}
				eyeLengthMM = 17;
			otherwise,
				error(sprintf('%s estimates not available for species %s',source,species));
        end

	% Orin Packer provided me with this information:
    % Monkey eye size varies a lot, so any particular number
    % is bound to be in error most of the time, sometimes substantially.
    % That said, I usually reference Perry & Cowey (1985) (Vis Res 25, 1795-1810).
    % For rhesus they report an eye diameter of 20 mm and a posterior nodal
    % distance of 12.8 mm.
	case {'PerryCowey'}
		switch (species)
			case {'Rhesus'}
				eyeLengthMM = 12.8;
			otherwise,
				error(sprintf('%s estimates not available for species %s',source,species));
			end

	% Orin Packer also said: For the fovea I usually use a conversion factor of
    % 210 um/degree which is the average of Perry & Cowey, Rolls & Cowey (1970)
    % (exp Br. Res., 10, 298) and deMonasterio et al (1985) (IOVS 26, 289-302).
    % This corresponds to 12.0324 mm.
	case {'Packer'}
		switch (species)
			case {'Rhesus'}
				mmPerDegree = 210*1e-3;
				eyeLengthMM = 0.5*mmPerDegree/atan((pi/180)*0.5);
			otherwise,
				error(sprintf('%s estimates not available for species %s',source,species));
        end

    % PennDog
    case {'PennDog'}
		switch (species)
			case {'Dog'}
				eyeLengthMM = 15;
			otherwise,
				error(sprintf('%s estimates not available for species %s',source,species));
        end
        
	otherwise
		error(sprintf('Unknown source %s for eye length estimates',source));
end
