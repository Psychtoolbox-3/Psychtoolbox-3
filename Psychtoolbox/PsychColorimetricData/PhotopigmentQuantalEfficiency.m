function quantalEfficiencies = PhotopigmentQuantalEfficiency(receptorTypes,species,source)
% quantalEfficiencies = PhotopigmentQuantalEfficiency.(receptorTypes,[species],[source])
%
% Return estimate of the fraction of absorbed quanta that lead to an isomerization.
%
% Supported types:
%   Any
%
% Supported species:
%		Any.
%
% Supported sources:
% 	Generic (Default).
%   None (returns 1 for all photoreceptor types)
%
% The generic source returns the value 0.667 for any species or type.  This
% is currently the only interesting source implemented.  The value comes from
% Rodieck RW, The First Steps in Seeing, page 472.
%
% 7/25/03  dhb  Wrote it.
% 8/9/13   dhb  Add type of 1.

% Fill in defaults
if (nargin < 2 || isempty(species))
	species = 'Human';
end
if (nargin < 3 || isempty(source))
	source = 'Rodieck';
end

% Fill in specific density according to specified source
if (iscell(receptorTypes))
	quantalEfficiencies = zeros(length(receptorTypes),1);
else
	quantalEfficiencies = zeros(1,1);
end
for i = 1:length(quantalEfficiencies)
	if (iscell(receptorTypes))
		type = receptorTypes{i};
	elseif (i == 1)
		type = receptorTypes;
	else
		error('Argument receptorTypes must be a string or a cell array of strings');
	end

	switch (source)
        case {'None'};
            quantalEfficiencies(i) = 1;
            
		case {'Generic'}
			quantalEfficiencies(i) = 0.667;

		otherwise
			error(sprintf('Unknown source %s for quantal efficiency estimates',source));
	end
end
 
