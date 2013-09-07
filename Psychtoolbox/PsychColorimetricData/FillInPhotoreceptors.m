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
%   ... statements here set value fields, to override default values
%       that are either filled in from source or that we directly set.
%   photoreceptors = FillInPhotoreceptors;
%
% Values supplied for fields that could be computed from other fields
% override the computed values, with those computed later in the code
% taking precedence over those computed earlier.  As of August, 2013
% a fair number of checks are in place to throw an error if you
% try to override a value that could have been computed from other passed
% parameters.
%
% This routine does not deal with filling in the pupil diameter field of the photoreceptors
% structure, which might be a mistake.  Pupil diameter is now used in the CIE calculations
% for lens density.  So perhaps this routine should fill in pupil diameter value of it
% doesn't exist. The problem is that computing pupil diameter from a source requires
% knowing the luminance of the stimulus, and this routine definitely does not want
% to know that.  So perhaps what should happen here is just a check that values needed
% exist.  See PupilDiameterFromLuminance if you are interested in knowing how pupil
% diameter might be determined given luminance.
%
% There are some useful comments towards the end of this routine about how to compute
% isomerizations from the absorbtances that this routine produces.
%
% See also: DefaultPhotoreceptors, RetIrradianceToIsoRecSec
%   IsomerizationsInEyeDemo, IsomerizationsInDishDemo, ComputeCIEConeFundamentals,
%   CIEConeFundamentalsTest, PrintPhotoreceptors.
%
% 7/25/03  dhb  Wrote it.
% 8/14/11  dhb  Allow pass through of field size, pupil diameter, and age.
%               Try not to break old code in how this is handled.
% 4/26/12  dhb  Return density as well as transmittance for lens and macular pigment.
% 8/9/13   dhb  Bulletproofing, by putting in a lot more consistency checks, and requiring
%               the calling program not to pass inconsistent information (e.g., you can't pass
%               a nomogram and an absorbance spectrum.)
% 8/11/13  dhb  More checking.  Add ability to adjust lens/macular density.  Return energy and quantal fundamentals (normalized to unity).
% 8/12/13  dhb  Fixed buglet resulting from forgetting to update after copy/paste.

%% Check that there is a nomogram field with an S subfield
%
% Long ago we tucked the wavelength sampling information into the nomogram field,
% never dreaming that there would come a day when we might pass absorbance or absorbtance
% directly.  So, we need to have the nomogram field because we need the wavelength
% sampling (S) field, and too much will break if we try to change to a more rational
% organization where the S field is hung directly onto the photoreceptors structure.
if (~isfield(photoreceptors,'nomogram'))
    error('You have to have a nomogram field, and it must contain an S subfield');
else
    if (~isfield(photoreceptors.nomogram,'S'))
        error('You have to have a nomogram field, and it must contain an S subfield');
    end
end
S = photoreceptors.nomogram.S;

%% Check that there is a types field.  These are just for people to read,
% but let's be draconian and make sure they are there.
if (~isfield(photoreceptors,'types'))
    error('The photoreceptors structure must have a types field')
end

%% Some defaults to avoid breaking older code.
% These parameters are passed to the lens and macular
% pigment density computing code, just in case the
% method specified is one that uses them.
if (~isfield(photoreceptors,'fieldSizeDegrees'))
    photoreceptors.fieldSizeDegrees = [];
end
if (~isfield(photoreceptors,'ageInYears'))
    photoreceptors.ageInYears = [];
end
if (~isfield(photoreceptors,'pupilDiameter'))
    photoreceptors.pupilDiameter.source = 'None';
elseif (~isfield(photoreceptors.pupilDiameter,'source') && ~isfield(photoreceptors.pupilDiameter,'value'))
    photoreceptors.pupilDiameter.source = 'None';
end

%% Consistency checks
if (isfield(photoreceptors,'nomogram') && isfield(photoreceptors.nomogram,'lambdaMax'))
    if (length(photoreceptors.types) ~= length(photoreceptors.nomogram.lambdaMax))
        error('Mismatch between length of types and lambdaMax fields');
    end
elseif (isfield(photoreceptors,'absorbance'))
    if (~isstr(photoreceptors.absorbance))
        if (length(photoreceptors.types) ~= size(photoreceptors.absorbance,1))
            error('Mismatch between length of types and absorbance fields');
        end
    end
elseif (isfield(photoreceptors,'absorbtance'))
    if (length(photoreceptors.types) ~= size(photoreceptors.absorbtance,1))
        error('Mismatch between length of types and absorbtance fields');
    end    
else
    error('Must specify either a photopigment nomogram or a photopigment absorbance spectrum');
end

%% Fill in photoreceptor dimensions
if (isfield(photoreceptors,'ISdiameter'))
    if (~isfield(photoreceptors.ISdiameter,'value'))
        photoreceptors.ISdiameter.value = ...
            PhotoreceptorDimensions(photoreceptors.types,'ISdiam', ...
            photoreceptors.species,photoreceptors.ISdiameter.source);
    else
        photoreceptors.ISdiameter.source = 'Value provided directly';
    end
end
if (isfield(photoreceptors,'OSdiameter'))
    if (~isfield(photoreceptors.OSdiameter,'value'))
        photoreceptors.OSdiameter.value = ...
            PhotoreceptorDimensions(photoreceptors.types,'OSdiam',...
            photoreceptors.species,photoreceptors.OSdiameter.source);
    else
        photoreceptors.OSdiameter.source = 'Value provided directly';
    end
end
if (isfield(photoreceptors,'OSlength'))
    if (~isfield(photoreceptors.OSlength,'value'))
        photoreceptors.OSlength.value = ...
            PhotoreceptorDimensions(photoreceptors.types,'OSlength',...
            photoreceptors.species,photoreceptors.OSlength.source);
    else
        photoreceptors.OSlength.source = 'Value provided directly';
    end
end

%% Pupil diameter
% The only place pupil diameter is used in this routine is to compute
% lens density when the source is CIE. In this case it must be directly
% provided as a value, since the various sources for pupil diameter
% compute the diameter from luminance, and this routine doesn't know
% about luminance.
%
% Pupil diameter is NOT used here to adjust the returend sensitivities
% for the manner in which the pupil affects the retinal illuminance.
if (isfield(photoreceptors,'pupilDiameter'))
    if (~isfield(photoreceptors.pupilDiameter,'value'))
        if (strcmp(photoreceptors.pupilDiameter.source,'None'))
            photoreceptors.pupilDiameter.value = [];
        end
    else
        photoreceptors.pupilDiameter.source = 'Value provided directly';
    end
end

%% Eye length
%
% This routine doesn't use eye length, but it will fetch it from the specified
% source so that it may be used elsewhere.  This was probably a design error
% but is maintained for backwards compatibility.
if (isfield(photoreceptors,'eyeLengthMM'))
    if (~isfield(photoreceptors.eyeLengthMM,'value'))
        photoreceptors.eyeLengthMM.value = EyeLength(photoreceptors.species,...
            photoreceptors.eyeLengthMM.source);
    else
        photoreceptors.eyeLengthMM.source = 'Value provided directly';
    end
end

%% Get quantal efficiency of photopigment 
% That is, the probability of an isomerization given an absorption.
% 
% If not passed, set it to 1.
if (isfield(photoreceptors,'quantalEfficiency'))
    if (~isfield(photoreceptors.quantalEfficiency,'value'))
        photoreceptors.quantalEfficiency.value = ...
            PhotopigmentQuantalEfficiency(photoreceptors.types,...
            photoreceptors.species,photoreceptors.quantalEfficiency.source);
    else
        photoreceptors.quantalEfficiency.source = 'Value provided directly';
    end
else
    photoreceptors.quantalEfficiency.source = 'None';
    for i = 1:size(photoreceptors.effectiveAbsorbtance,1)
        photoreceptors.quantalEfficiency.value(i) = 1;
    end
end

%% Fill in specific density
if (isfield(photoreceptors,'specificDensity'))
    if (~isfield(photoreceptors.specificDensity,'value'))
        photoreceptors.specificDensity.value = ...
            PhotopigmentSpecificDensity(photoreceptors.types,...
            photoreceptors.species,photoreceptors.specificDensity.source);
    else
        photoreceptors.specificDensity.source = 'Value provided directly';
    end
end

%% Determine the axial optical density.
%
% This can be computed from the specific density and OS length or obtained directly.
if (~isfield(photoreceptors,'axialDensity'))
    [photoreceptors.axialDensity.value] = ComputeAxialDensity(photoreceptors.specificDensity.value,...
        photoreceptors.OSlength.value);
    photoreceptors.axialDensity.source = 'Computed from specific density and OS length';
else
    if (~isfield(photoreceptors.axialDensity,'value'))
        photoreceptors.axialDensity.value = PhotopigmentAxialDensity(photoreceptors.types,...
            photoreceptors.species,photoreceptors.axialDensity.source,photoreceptors.fieldSizeDegrees);
    else
        photoreceptors.axialDensity.source = 'Value provided directly';
    end
    
    % Since we have ignored the specific density and OS length, check that 
    % this was intentional.  If there are value fields for either of the components,
    % throw an error.  The recommended fix is to remove these fields from the structure,
    % or set their source to 'None' or their value fields to empty, before the call to this routine.
    if (isfield(photoreceptors,'OSlength'))
        if (isfield(photoreceptors.OSlength,'value') && ~isempty(photoreceptors.OSlength.value))
            error('OS length provided but axial density overrides effect of that length');
        end
    end
    if (isfield(photoreceptors,'specificDensity'))
        if (isfield(photoreceptors.specificDensity,'value') && ~isempty(photoreceptors.specificDensity.value))
            error('Specific density provided by axial density overrides effect of that density');
        end
    end
end

%% Absorbance spectrum, either supplied or from nomogram
%
% The design of the fields here is a little inconsistent, because the source for the absorbance nomoogram is
% in the nomogram field.  More rational would have been to have the nomogram name as the source
% for the absorbance and get rid of the nomogram field, and then have a value field under
% absorbance to hold the absorbance.
%
% Similarly, given that we have a nomogram field, the lambdaMax field of the nomogram field would
% better be called value.
%
% But we've been living with this special case for so long that I am loathe to change it now for
% fear of all the work I'd have to do to unbreak calling routines, and these may exist far outside
% of the toolbox and my lab.
if (~isfield(photoreceptors,'absorbance') || isempty(photoreceptors.absorbance))
    if (~strcmp(photoreceptors.nomogram.source,'None'))
        photoreceptors.absorbance = ...
            PhotopigmentNomogram(photoreceptors.nomogram.S,photoreceptors.nomogram.lambdaMax, ...
            photoreceptors.nomogram.source);
    else
        error('No absorbance supplied, but the nomogram field source is ''None''.  That''s not going to work');
    end
else
    % Check for unhappiness.  Unhappiness occurs when there is both a nomogram specified and the absorbance has been
    % passed by the calling program.  The fix is to remove one or the other in the calling program.
    %
    % Note that there must abe a nomogram field, because we must have a .nomogram.S at call.  This is checked above.
    if (isfield(photoreceptors.nomogram,'source') && ~strcmp(photoreceptors.nomogram.source,'None'))
            error('There is a directly supplied absorbance and also a nomogram source field specified.  One or the other must go away');
    else
        % The PTB-style data file containing the absorbance may be specified and loaded.  Otherwise
        % it is assumed that the absorbance field has the right data at the right wavelenght sampling.
        if (isstr(photoreceptors.absorbance))
            theAbsorbanceStr = photoreceptors.absorbance;
            eval(['load(''T_' theAbsorbanceStr ''');']);
            eval(['photoreceptors.absorbance = 10.^SplineCmf(S_' theAbsorbanceStr ',T_' theAbsorbanceStr ',photoreceptors.nomogram.S,2);']);
            eval(['clear(''T_' theAbsorbanceStr ''',''S_' theAbsorbanceStr ''');']);
        end
    end
end

%% Absorbtance, either computed or supplied.
if (~isfield(photoreceptors,'absorbtance'))
    [photoreceptors.absorbtance] = AbsorbanceToAbsorbtance(...
        photoreceptors.absorbance,S,photoreceptors.axialDensity.value);
else
    if (isfield(photoreceptors,'absorbance'));
        error('There is a directly supplied absorbtance, but also either an absorbance or nomogram field specified.  Something must go away');
    end
end

%% Lens density.  Put in unity if there is no field yet.
if (isfield(photoreceptors,'lensDensity'))
    if (isfield(photoreceptors.lensDensity,'source'))
        if (isfield(photoreceptors.lensDensity,'transmittance') && ~isempty(photoreceptors.lensDensity.transmittance))
            error('Both source and transmittance passed in lens density field.  Choose one');
        end
        
        % You may wonder why age and pupil diameter are passed to a routine that computes
        % lens density.  This is because the CIE standard (which is one choice of source)
        % computes density in an age and pupil size dependent fashion.  These fields
        % are ignored for sources that don't use them, and can be empty in that case.  Indeed,
        % if you go back to the top of this routine you'll see that they are set to empty if not
        % passed, so that this code doesn't barf.
        [photoreceptors.lensDensity.transmittance,photoreceptors.lensDensity.density] = ...
            LensTransmittance(S,photoreceptors.species,photoreceptors.lensDensity.source,...
            photoreceptors.ageInYears,photoreceptors.pupilDiameter.value);
        
        % Adjust transmittance for a change in density from whatever standard was specified, if desired.
        if (isfield(photoreceptors.lensDensity,'adjustDen'))
            photoreceptors.lensDensity.transmittance = 10.^(-(photoreceptors.lensDensity.density + photoreceptors.lensDensity.adjustDen));
            photoreceptors.lensDensity.transmittance(photoreceptors.lensDensity.transmittance > 1) = 1;
        end
    else
        if (~isfield(photoreceptors.lensDensity.transmittance))
            error('photoreceptors.lensDensity field passed, but without source or transmittance');
        end
    end
else
    photoreceptors.lensDensity.source = 'None';
    photoreceptors.lensDensity.density = zeros(S(3),1)';
    photoreceptors.lensDensity.transmittance = ones(S(3),1)';
end

%% Macular pigment density.  Put in unity if there is none.
if (isfield(photoreceptors,'macularPigmentDensity'))
    if (isfield(photoreceptors.macularPigmentDensity,'source'))
        if (isfield(photoreceptors.macularPigmentDensity,'transmittance') && ~isempty(photoreceptors.macularPigmentDensity.transmittance))
            error('Both source and transmittance passed in macular density field.  Choose one');
        end
        
        % You may wonder why field size is passed to a routine that computes
        % macular density.  This is because the CIE standard (which is one choice of source)
        % computes density in a field size dependent fashion.  This field is
        % ignored for sources that don't use them, and can be empty in that case.  Indeed,
        % if you go back to the top of this routine you'll see that it is set to empty if not
        % passed, so that this code doesn't barf.
        [photoreceptors.macularPigmentDensity.transmittance,photoreceptors.macularPigmentDensity.density] = ...
            MacularTransmittance(S,photoreceptors.species,photoreceptors.macularPigmentDensity.source,photoreceptors.fieldSizeDegrees);
        
        % Adjust transmittance for a change in density from whatever standard was specified, if desired.
        if (isfield(photoreceptors.macularPigmentDensity,'adjustDen'))
            photoreceptors.macularPigmentDensity.transmittance = 10.^(-(photoreceptors.macularPigmentDensity.density + photoreceptors.macularPigmentDensity.adjustDen));
            photoreceptors.macularPigmentDensity.transmittance(photoreceptors.macularPigmentDensity.transmittance > 1) = 1;
        end
    else
        if (~isfield(photoreceptors.macularPigmentDensity.transmittance))
            error('photoreceptors.macularPigmentDensity field passed, but without source or transmittance');
        end
    end
else
    photoreceptors.macularPigmentDensity.source = 'None';
    photoreceptors.macularPigmentDensity.density = zeros(S(3),1)';
    photoreceptors.macularPigmentDensity.transmittance = ones(S(3),1)';
end

%% Compute overall pre-receptor transmittance
if (~isfield(photoreceptors,'preReceptoral'))
    photoreceptors.preReceptoral.source = 'Computed from lens and macular pigment transmittance';
    photoreceptors.preReceptoral.transmittance = photoreceptors.lensDensity.transmittance .* ...
        photoreceptors.macularPigmentDensity.transmittance;
elseif (~strcmp(photoreceptors.lensDensity.source,'None') | ~strcmp(photoreceptors.macularPigmentDensity.source,'None'))
    error('Pre-receptoral filtering specified directly, but non-unity lens or macular pigment transmittance also specified.');
end

%% Compute effective absorbtance, which takes pre-receptor transmittance into account.
% 
% This is the probability of a quantal absorption as a function of wavelength, referred to light
% entering the pupil.  That is, it accounts for pre-retinal absorbtions, but not pupil size or 
% eye length.  Typically, you'd compute retinal irradiance in quanta/[time-area], ignoring 
% pre-retinal absorptions, and then compute absorptions per receptor using this field and
% by multiplying by the desired integration time and collecting area.  We take the collecting area to be the 
% inner segment diameter in our computation routines, which is right for cones and I'm not
% sure about rods.
if (~isfield(photoreceptors,'effectiveAbsorbtance'))
    photoreceptors.effectiveAbsorbtance = photoreceptors.absorbtance .* ...
        (ones(size(photoreceptors.absorbtance,1),1)*photoreceptors.preReceptoral.transmittance);
else
    fprintf('If you''re passing the effectiveAbsorbtance to this routine, then you are ignoring essentially everything the routine does.\n')  
    fprintf('The only additional step it would take is to multiply by the quantal efficiency to get isomerizations from absorptions, and\n');
    fprintf('to set the eye length.  If you''ve gone so far as to do everything else by hand you should do that step by hand too.\n');
    error('We''re throwing an error to help you mend your ways.');
end

%% Compute isomerizationAbsorbtance, which takes quantalEfficiency into account
%
% When you want to compute isomerization rates from retinal irradiance (and who wouldn't), this
% is what you want.  See note just above about conventions with respect to conventions.
if (~isfield(photoreceptors,'isomerizationAbsorbtance'))
    for i = 1:size(photoreceptors.effectiveAbsorbtance,1)
        photoreceptors.isomerizationAbsorbtance(i,:) = photoreceptors.quantalEfficiency.value(i) * ...
            photoreceptors.effectiveAbsorbtance(i,:);
    end
else
    fprintf('If you''re passing the isomerizationAbsorbtance to this routine, then you are ignoring essentially everything the routine does.\n')  
    fprintf('If you''ve gone so far as to do everything by hand you should not be calling this routine.\n');
    error('We''re throwing an error to help you mend your ways.');
end

%% Compute normalized energy sensitivities (aka cone fundamentals in energy units)
%
% Note that we really do mean to call EnergyToQuanta here, because the PTB routines 
% implement the conversion for spectra.  The conversion for sensitivities is the
% inverse of the conversion for spectra.
photoreceptors.energyFundamentals = EnergyToQuanta(S,photoreceptors.isomerizationAbsorbtance')';
mx = max(photoreceptors.energyFundamentals,[],2);
photoreceptors.energyFundamentals = diag(1./mx)*photoreceptors.energyFundamentals;

%% Compute normalized quantal sensitivities (aka cone fundamentals in quantal units)
photoreceptors.quantalFundamentals = photoreceptors.isomerizationAbsorbtance;
mx = max(photoreceptors.quantalFundamentals,[],2);
photoreceptors.quantalFundamentals = diag(1./mx)*photoreceptors.quantalFundamentals;


