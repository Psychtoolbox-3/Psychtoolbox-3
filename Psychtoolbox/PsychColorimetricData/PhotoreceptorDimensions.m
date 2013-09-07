function dimensions = PhotoreceptorDimensions(receptorTypes,whichDimension,species,source)
% dimensions = PhotoreceptorDimensions(receptorTypes,whichDimension,species,source)
%
% Return estimates of photoreceptor dimensions.
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
% That is to say, for the same source, species, and type, you should get
% a consistent triplet of numbers.
%
% Argument whichDimension may take on values:
% 	OSdiam, ISdiam, OSlength.
%
% Supported species:
%		Human (Default), GuineaPig, Dog
%
% Supported sources:
%   Rodeick (Default).
%  	CVRL (Human Cone OS length)
%   Webvision (Human Cone IS diameter)
%   Hendrickson (Human Rod OS length)
% 	SterlingLab (GuineaPig dimensions).
%   Generic
%   PennDog (Dog dimensions).
%   None (returns empty for the corresponding value)
%
% The Generic type returns a single number for all species/type.
%
% 7/11/03  dhb  Wrote it.
% 12/04/07 dhb  Added dog but with placeholder numbers.
% 8/9/13   dhb  Comment clean up, allow 'None' to return empty as the value.
% 8/10/13  dhb  Added Webvision source for IS diameter.

% Fill in defaults
if (nargin < 3 || isempty(species))
    species = 'Human';
end
if (nargin < 4 || isempty(source))
    source = 'Rodieck';
end

% Fill in dimensions according to specified source
if (iscell(receptorTypes))
    dimensions = zeros(length(receptorTypes),1);
else
    dimensions = zeros(1,1);
end
for i = 1:length(dimensions)
    if (iscell(receptorTypes))
        type = receptorTypes{i};
    elseif (i == 1)
        type = receptorTypes;
    else
        error('Argument receptorTypes must be a string or a cell array of strings');
    end

    switch (source)
        case {'None'}
            dimensions = [];
        case {'Generic'}
            % These are fairly generic dimensions
            switch (whichDimension)
                case 'OSlength'
                    dimensions(i) = 32;
                case 'OSdiam'
                    dimensions(i) = 2;
                case 'ISdiam'
                    dimensions(i) = 2;
                otherwise
                    error('Unsupported dimension requested');
            end 
        case {'Webvision'}
            % http://webvision.med.utah.edu
            switch (whichDimension)
                case 'ISdiam'
                    switch (type)
                        % http://webvision.med.utah.edu/book/part-ii-anatomy-and-physiology-of-the-retina/photoreceptors/
                        % Attributed to Helga Kolb
                        case {'FovealLCone', 'FovealMCone' 'FovealSCone'}
                            dimensions(i) = 1.5;
                        case {'LCone', 'MCone' 'SCone'}
                            dimensions(i) = 6;
                        case {'Rod'}
                            dimensions(i) = 2;
                        otherwise,
                            error(sprintf('Unsupported receptor type %s/%s for %s estimates in %s',...
                                type,whichDimension,source,species));
                    end
                otherwise
                    error('Unsupported dimension requested');
            end
        case ('PennDog')
            % Numbers we use for dog eyes at Penn.  Got these from
            % Gus Aguirre.  See emails sent about 12/5/07.
            switch (species)
                case {'Dog'}
                    switch (whichDimension)
                        case 'OSlength',
                            switch (type)
                                case {'LCone', 'SCone'}
                                    dimensions(i) = 13;
                                case {'Rod'}
                                    dimensions(i) = 13.5;
                                otherwise,
                                    error(sprintf('Unsupported receptor type %s/%s for %s estimates in %s',...
                                        type,whichDimension,source,species));
                            end
                        case 'ISdiam'
                            switch (type)
                                case {'LCone', 'SCone'}
                                    dimensions(i) = 2;
                                case {'Rod'}
                                    dimensions(i) = 2;
                                otherwise,
                                    error(sprintf('Unsupported receptor type %s/%s for %s estimates in %s',...
                                        type,whichDimension,source,species));
                            end
                        case 'OSdiam'
                            switch (type)
                                case {'LCone', 'SCone'}
                                    dimensions(i) = 1.25;
                                case {'Rod'}
                                    dimensions(i) = 1;
                                otherwise,
                                    error(sprintf('Unsupported receptor type %s/%s for %s estimates in %s',...
                                        type,whichDimension,source,species));
                            end
                        otherwise,
                            error(sprintf('Unsupported dimension %s requested',whichDimension));
                    end
                otherwise,
                    error(sprintf('%s estimates not available for species %s',source,species));
            end
     
        case ('Rodieck')
            % From Rodieck's "standard observer", Appendix B
            % in The First Steps of Seeing.
            switch (species)
                case {'Human'}
                    switch (whichDimension)
                        case 'OSlength',
                            switch (type)
                                case {'FovealLCone', 'FovealMCone', 'FovealSCone'}
                                    dimensions(i) = 33;
                                case {'Rod'}
                                    dimensions(i) = 31.2;
                                otherwise,
                                    error(sprintf('Unsupported receptor type %s/%s for %s estimates in %s',...
                                        type,whichDimension,source,species));
                            end
                        case 'ISdiam'
                            switch (type)
                                case {'FovealLCone', 'FovealMCone', 'FovealSCone'}
                                    dimensions(i) = 2.3;
                                case {'Rod'}
                                    dimensions(i) = 2.22;
                                otherwise,
                                    error(sprintf('Unsupported receptor type %s/%s for %s estimates in %s',...
                                        type,whichDimension,source,species));
                            end
                        otherwise,
                            error(sprintf('Unsupported dimension %s requested',whichDimension));
                    end
                otherwise,
                    error(sprintf('%s estimates not available for species %s',source,species));
            end

        case {'CVRL'}
            % These numbers are my encapsulations of CVRL's summary of a variety of data.
            % See CVRL summary text at:
            %   http://cvrl.ioo.ucl.ac.uk/database/text/intros/introlength.htm.
            switch (species)
                case {'Human'}
                    switch (whichDimension)
                        case 'OSlength',
                            switch (type)
                                case {'FovealLCone', 'FovealMCone'}
                                    dimensions(i) = 35.5;
                                case {'FovealSCone'}
                                    dimensions(i) = 35.5*0.95;
                                case {'LCone', 'MCone'}
                                    dimensions(i) = 18;
                                case {'SCone'}
                                    dimensions(i) = 18*0.82;
                                otherwise,
                                    error(sprintf('Unsupported receptor type %s/%s for %s estimates in %s',...
                                        type,whichDimension,source,species));
                            end
                        otherwise,
                            error(sprintf('Unsupported dimension %s requested',whichDimension));
                    end
                otherwise,
                    error(sprintf('%s estimates not available for species %s',source,species));
            end
       
        case {'Hendrickson'}
            % From Hendrickson and Drucker, numbers provided at CVRL database:
            % http://cvrl.ioo.ucl.ac.uk/database/text/outseg/length.htm.  40 um
            % is the number provided for mid-peripheral rods, 40-45 is cited for
            % parafoveal rods.  This routines returns 40.
            switch (species)
                case {'Human'}
                    switch (whichDimension)
                        case 'OSlength',
                            switch (type)
                                case {'Rod'}
                                    dimensions(i) = 40;
                                otherwise,
                                    error(sprintf('Unsupported receptor type %s/%s for %s estimates in %s',...
                                        type,whichDimension,source,species));
                            end
                        otherwise,
                            error(sprintf('Unsupported dimension %s requested',whichDimension));
                    end
                otherwise,
                    error(sprintf('%s estimates not available for species %s',source,species));
            end

        case {'SterlingLab'}
            % These are values that Lu Yin provided, based on unpublished
            % measurements used in the Sterling lab.
            switch (species)
                case {'GuineaPig'}
                    switch (whichDimension)
                        case 'OSdiam',
                            switch (type)
                                case {'LCone', 'MCone', 'SCone'}
                                    dimensions(i) = 2;
                                case 'Rod'
                                    dimensions(i) = 2;
                                otherwise,
                                    error(sprintf('Unsupported receptor type %s/%s for %s estimates in %s',...
                                        type,whichDimension,source,species));
                            end
                        case 'ISdiam',
                            switch (type)
                                case {'LCone', 'MCone', 'SCone'}
                                    dimensions(i) = 2.8;
                                case 'Rod'
                                    dimensions(i) = 2.4;
                                otherwise,
                                    error(sprintf('Unsupported receptor type %s/%s for %s estimates in %s',...
                                        type,whichDimension,source,species));
                            end
                        case 'OSlength',
                            switch (type)
                                case {'LCone', 'MCone', 'SCone'}
                                    dimensions(i) = 8;
                                case 'Rod'
                                    dimensions(i) = 16.2;
                                otherwise,
                                    error(sprintf('Unsupported receptor type %s/%s for %s estimates in %s',...
                                        type,whichDimension,source,species));
                            end
                        otherwise,
                            error(sprintf('Unsupported dimension %s requested',whichDimension));
                    end
                otherwise,
                    error(sprintf('%s estimates not available for species %s',source,species));
            end

            % Nope!
        otherwise
            error(sprintf('Unknown source %s for photoreceptor dimension estimates',source));
    end
end
