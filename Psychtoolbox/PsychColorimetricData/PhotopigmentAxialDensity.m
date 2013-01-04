function densities = PhotopigmentAxialDensity(receptorTypes,species,source,fieldSizeDegrees)
% densities = PhotopigmentAxialDensity(receptorTypes,[species],[source],[fieldSizeDegrees])
%
% Return estimates of photopigment axial density, sometimes called peak
% absorbance.
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
%		Human (Default).
%
% Supported sources:
% 	Rodieck (Human) (Default).
%   StockmanSharpe (Human).
%   CIE (Human).
%   Tsujimuar (Human, melanopsin gc's)
%
% The CIE method takes a field size argument.  This
% overrides the specified foveal or not part of the
% cone string.  If the field type is not passed and
% the method is CIE, it is set to 10-degrees for Scone, Mcone,
% and Lcone, and to 2-degrees for FovealScone, FovealMCone, and
% FovealLCone.
%
% The fieldSizeDegrees argument is ignored for sources other than
% CIE.
%
% 7/11/03  dhb  Wrote it.
% 8/12/11  dhb  Added CIE source, and allow passing of fieldSizeDegrees.
% 4/20/12  dhb  Add Tsujimura's estimate of melanopsin optical density in human.
% 12/16/12 dhb, ms Add Alpern's rod estimates from CVRL table.

% Fill in defaults
if (nargin < 2 || isempty(species))
	species = 'Human';
end
if (nargin < 3 || isempty(source))
	source = 'Rodieck';
end
if (nargin < 4 || isempty(fieldSizeDegrees))
    fieldSizeDegrees = [];
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
		case {'Rodieck'}
			switch (species)
				case {'Human'},
					% Rodieck, The First Steps in Seeing, Appendix B.
					switch (type)
						case {'FovealLCone','FovealMCone','FovealSCone'}
							densities(i) = 0.50;
						case 'Rod'
							densities(i) = 0.47;
						otherwise,
							error(sprintf('Unsupported receptor type %s for %s estimates in %s',type,source,species));
					end
				otherwise,
					error(sprintf('%s estimates not available for species %s',source,species));
            end	
            
        case {'Alpern'}
            % This value from the CVRL table of receptor optical density.  There are 3 papers from
            % Alpern's lab cited, with values of 0.342 (Alpern & Pugh, 1974), 0.342 (Zwas & Aloper, 1976), and
            % 0.318 (Zwas & Alpern, 1976). The value of 0.333 here is roughly the mean of these.
            switch (species)
                case {'Human'}
                    switch (type)
						case 'Rod'
							densities(i) = 0.333;
						otherwise,
							error(sprintf('Unsupported receptor type %s for %s estimates in %s',type,source,species));
                    end
                otherwise
                    error(sprintf('%s estimates not available for species %s',source,species));
            end

		case {'StockmanSharpe'}
			switch (species)
				case {'Human'},
					% Foveal values from Note c, Table 2, Stockman and Sharpe (2000), Vision Research.  These
					% are the values they used to produce a fit to their 2-degree fundamentals.  The peripheral
					% values were provided to me by Andrew Stockman, and were used to produce a fit to their
					% 10-degree fundamentals.
					switch (type)
						case {'FovealLCone','FovealMCone'}
							densities(i) = 0.50;
						case 'FovealSCone'
							densities(i) = 0.40;
						case {'LCone', 'MCone'}
							densities(i) = 0.38;
						case {'SCone'}
							densities(i) = 0.3;
						otherwise,
							error(sprintf('Unsupported receptor type %s for %s estimates in %s',type,source,species));
					end
				otherwise,
					error(sprintf('%s estimates not available for species %s',source,species));
            end	
            
        case {'CIE'}
			switch (species)
				case {'Human'},
					% These values computed according to formulae in CIE 170-1:2006.
					switch (type)
						case {'FovealLCone','FovealMCone'}
                            if (isempty(fieldSizeDegrees))
                                fieldSizeDegrees = 2;
                            end
                            densities(i) = 0.38+0.54*exp(-fieldSizeDegrees/1.333);
						case 'FovealSCone'
                             if (isempty(fieldSizeDegrees))
                                fieldSizeDegrees = 2;
                            end
							densities(i) = 0.30+0.45*exp(-fieldSizeDegrees/1.333);
						case {'LCone', 'MCone'}
                             if (isempty(fieldSizeDegrees))
                                fieldSizeDegrees = 10;
                            end
							densities(i) = 0.38+0.54*exp(-fieldSizeDegrees/1.333);
						case {'SCone'}
                             if (isempty(fieldSizeDegrees))
                                fieldSizeDegrees = 10;
                            end
							densities(i) = 0.30+0.45*exp(-fieldSizeDegrees/1.333);
						otherwise,
							error(sprintf('Unsupported receptor type %s for %s estimates in %s',type,source,species));
					end
				otherwise,
					error(sprintf('%s estimates not available for species %s',source,species));
            end	
            
        case {'Tsujimura'}
			switch (species)
				case {'Human'},
					switch (type)
						case {'Melanopsin'}
							densities(i) = 0.50;
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
