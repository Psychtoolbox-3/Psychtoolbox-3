function photoreceptors = FillInPhotoreceptors(photoreceptors)
% photoreceptors = FillInPhotoreceptors(photoreceptors)
%
% Convert all source strings in a photoreceptors structures
% to numerical values, so that the result is ready to compute
% on.
%
% The typical usage of this routine would be:
%
%   clear photoreceptors
%   photoreceptors = DefaultPhotoreceptors('LivingHumanFovea');
%   ... statements here to override default values ...
%   photoreceptors = FillInPhotoreceptors;
%
% Computed fields that exist in the passed structure override
% passed values, with those computed last taking precedence
% over those computed earlier.
%
% NOTES: The override rules are pretty tricky, and this code
% should be improved to clarify.
%
% See also: DefaultPhotoreceptors, RetIrradianceToIsoRecSec
%   IsomerizationsInEyeDemo, IsomerizationsInDishDemo 
%
% 7/25/03  dhb  Wrote it.
% 8/14/11  dhb  Allow pass through of field size, pupil diameter, and age.
%               Try not to break old code in how this is handled.
% 4/26/12  dhb  Return density as well as transmittance for lens and macular pigment.

% Define common wavelength sampling for this function.
S = photoreceptors.nomogram.S;

% Some defaults to avoid breaking older code.
if (~isfield(photoreceptors,'fieldSizeDegrees'))
    photoreceptors.fieldSizeDegrees = [];
end
if (~isfield(photoreceptors,'ageInYears'))
    photoreceptors.ageInYears = [];
end
if (~isfield(photoreceptors,'pupilDiameter'))
    photoreceptors.pupilDiameter.source = '';
    photoreceptors.pupilDiameter.value = [];
elseif (~isfield(photoreceptors.pupilDiameter,'value'))
    photoreceptors.pupilDiameter.value = [];
end
    
% Consistency checks
if (length(photoreceptors.types) ~= length(photoreceptors.nomogram.lambdaMax))
	error('Mismatch between length of types and lambdaMax fields');
end

% Fill in photoreceptor dimensions
if (isfield(photoreceptors,'ISdiameter'))
	if (~isfield(photoreceptors.ISdiameter,'value'))
		photoreceptors.ISdiameter.value = ...
			PhotoreceptorDimensions(photoreceptors.types,'ISdiam', ...
				photoreceptors.species,photoreceptors.ISdiameter.source);
	end
end
if (isfield(photoreceptors,'OSdiameter'))
	if (~isfield(photoreceptors.OSdiameter,'value'))
		photoreceptors.OSdiameter.value = ...
			PhotoreceptorDimensions(photoreceptors.types,'OSdiam',...
				photoreceptors.species,photoreceptors.OSdiameter.source);
	end
end
if (isfield(photoreceptors,'OSlength'))
	if (~isfield(photoreceptors.OSlength,'value'))
		photoreceptors.OSlength.value = ...
			PhotoreceptorDimensions(photoreceptors.types,'OSlength',...
				photoreceptors.species,photoreceptors.OSlength.source);
	end
end

% Fill in specific density
if (isfield(photoreceptors,'specificDensity'))
	if (~isfield(photoreceptors.specificDensity,'value'))
		photoreceptors.specificDensity.value = ...
			PhotopigmentSpecificDensity(photoreceptors.types,...
				photoreceptors.species,photoreceptors.specificDensity.source);
	end
end

% Compute the axial optical density if it wasn't passed.  If it was passed,
% the source/value passed override the calculation.
if (~isfield(photoreceptors,'axialDensity'))
	[photoreceptors.axialDensity.value] = ComputeAxialDensity(photoreceptors.specificDensity.value,...
		photoreceptors.OSlength.value);
else
	if (~isfield(photoreceptors.axialDensity,'value'))
		photoreceptors.axialDensity.value = PhotopigmentAxialDensity(photoreceptors.types,...
				photoreceptors.species,photoreceptors.axialDensity.source,photoreceptors.fieldSizeDegrees);
	end
end

% Generate absorbance spectrum from specified nomogram
if (~isfield(photoreceptors,'absorbance'))
	photoreceptors.absorbance = ...
		PhotopigmentNomogram(photoreceptors.nomogram.S,photoreceptors.nomogram.lambdaMax, ...
		photoreceptors.nomogram.source);
end

% Convert the absorbance spectra of photoreceptors to absorbtance spectra
if (~isfield(photoreceptors,'absorbtance'))
	[photoreceptors.absorbtance] = AbsorbanceToAbsorbtance(...
		photoreceptors.absorbance,S,photoreceptors.axialDensity.value);
end

% Lens density.  Put in unity if there is none.
if (isfield(photoreceptors,'lensDensity'))
	if (isfield(photoreceptors.lensDensity,'source'))
		[photoreceptors.lensDensity.transmittance,photoreceptors.lensDensity.density] = ...
			LensTransmittance(S,photoreceptors.species,photoreceptors.lensDensity.source,...
            photoreceptors.ageInYears,photoreceptors.pupilDiameter.value);
	else
		if (~isfield(photoreceptors.lensDensity.transmittance))
			error('photoreceptors.lensDensity passed, but without source or transmittance');
		end
	end
else
	photoreceptors.lensDensity.transmittance = ones(S(3),1)';	
end

% Macular pigment density.  Put in unity if there is none.
if (isfield(photoreceptors,'macularPigmentDensity'))
	if (isfield(photoreceptors.macularPigmentDensity,'source'))
		[photoreceptors.macularPigmentDensity.transmittance,photoreceptors.macularPigmentDensity.density] = ...
			MacularTransmittance(S,photoreceptors.species,photoreceptors.macularPigmentDensity.source,photoreceptors.fieldSizeDegrees);
	else
		if (~isfield(photoreceptors.macularPigmentDensity.transmittance))
			error('photoreceptors.macularPigmentDensity passed, but without source or transmittance');
		end
	end
else
	photoreceptors.macularPigmentDensity.transmittance = ones(S(3),1)';	
end

% Compute overall pre-receptor transmittance
photoreceptors.preReceptoral.transmittance = photoreceptors.lensDensity.transmittance .* ...
	photoreceptors.macularPigmentDensity.transmittance;

% Compute effective absorbtance, which takes pre-receptor transmittance into account.
if (~isfield(photoreceptors,'effectiveAbsorbtance'))
	photoreceptors.effectiveAbsorbtance = photoreceptors.absorbtance .* ...
		(ones(length(photoreceptors.nomogram.lambdaMax),1)*photoreceptors.preReceptoral.transmittance);
end

% Get quantal efficiency of photopigment
if (isfield(photoreceptors,'quantalEfficiency'))
	if (~isfield(photoreceptors.quantalEfficiency,'value'))
		photoreceptors.quantalEfficiency.value = ...
			PhotopigmentQuantalEfficiency(photoreceptors.types,...
				photoreceptors.species,photoreceptors.quantalEfficiency.source);
	end
end

% Compute isomerizationAbsorbtance, which takes quantalEfficiency into account
if (~isfield(photoreceptors,'isomerizationAbsorbtance'))
	for i = 1:size(photoreceptors.effectiveAbsorbtance,1)
		photoreceptors.isomerizationAbsorbtance(i,:) = photoreceptors.quantalEfficiency.value(i) * ...
			photoreceptors.effectiveAbsorbtance(i,:);
	end
end

% Eye length
if (isfield(photoreceptors,'eyeLengthMM'))
	if (~isfield(photoreceptors.eyeLengthMM,'value'))
		photoreceptors.eyeLengthMM.value = EyeLength(photoreceptors.species,...
			photoreceptors.eyeLengthMM.source);
	end
end
