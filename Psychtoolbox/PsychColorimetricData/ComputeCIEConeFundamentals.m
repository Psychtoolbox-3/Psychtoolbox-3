function [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations,adjIndDiffParams,params,staticParams] = ...
    ComputeCIEConeFundamentals(S,fieldSizeDegrees,ageInYears,pupilDiameterMM,lambdaMax,whichNomogram,LserWeight, ...
    DORODS,rodAxialDensity,fractionPigmentBleached,indDiffParams)
% [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations,adjIndDiffParams,params,staticParams] = ...
%   ComputeCIEConeFundamentals(S,fieldSizeDegrees,ageInYears,pupilDiameterMM,[lambdaMax],[whichNomogram],[LserWeight], ...
%   [DORODS],[rodAxialDensity],[fractionPigmentBleached],[indDiffParams])
%
% Function to compute normalized cone quantal sensitivities
% from underlying pieces, as specified in CIE 170-1:2006.
%
% IMPORTANT: This routine returns quantal sensitivities.  You
% may want energy sensitivities.  In that case, use EnergyToQuanta to convert
%   T_energy = EnergyToQuanta(S,T_quantal')'
% and then renormalize.  (You call EnergyToQuanta because you're converting
% sensitivities, which go the opposite direction from spectra.)
% The routine also returns two quantal sensitivity functions. The first gives
% the probability that a photon will be absorbed.  The second is the probability
% that the photon will cause a photopigment isomerization.  It is the latter
% that is what you want to compute isomerization rates from retinal illuminance.
% See note at the end of function FillInPhotoreceptors for some information about
% convention.  In particular, this routine takes pre-retinal absorption into
% account in its computation of probability of absorptions and isomerizations,
% so that the relevant retinal illuminant is one computed without accounting for
% those factors.  This routine does not account for light attenuation due to
% the pupil, however.  The only use of pupil size here is becuase of its
% slight effect on lens density as accounted for in the CIE standard.
%
% This standard allows customizing the fundamentals for
% field size, observer age, and pupil size in mm.
%
% To get the Stockman-Sharpe/CIE 2-deg fundamentals, use
%   fieldSizeDegrees = 2;
%   ageInYears = 32;
%   pupilDiameterMM = 3;
% and don't pass the rest of the arguments.
%
% To get the Stockman-Sharpe/CIE 10-deg fundamentals, use
%   fieldSizeDegrees = 10;
%   ageInYears = 32;
%   pupilDiameterMM = 3;
% and don't pass the rest of the arguments.
%
% Although this routine will compute something over any wavelength
% range, I'd (DHB) recommend not going lower than 390 or above about 780 without
% thinking hard about how various pieces were extrapolated out of the range
% that they are specified in the standard.  Indeed, the lens optical
% density measurements only go down to 400 nm and these are extropolated
% to go below 400.
%
% This routine will compute from tabulated absorbance or absorbance based on a nomogram, where
% whichNomogram can be any source understood by the routine PhotopigmentNomogram.  To obtain
% the nomogram behavior, pass a lambdaMax vector. You can then also optionally pass a nomogram
% source (default: StockmanSharpe).  This option (using shifted nomograms) is not part of the
% CIE standard. See NOTE below for another way to handle individual differences 
%
% The nominal values of lambdaMax to fit the CIE 2-degree fundamentals with the
% Stockman-Sharpe nomogram are 558.9, 530.3, and 420.7 nm for the LMS cones respectively.
% These in fact do a reasonable job of reconstructing the CIE 2-degree fundamentals, although
% there are small deviations from what you get if you simply read in the tabulated cone
% absorbances.  Thus starting with these as nominal values and shifting is one way to
% produce fundamentals tailored to observers with different known photopigments.
%
% If you pass lambaMax and its length is 4, then first two values are treated as
% the peak wavelengths of the ser/ala variants of the L cone pigment, and these
% are then weighted according to LserWeight and (1-LserWeight).  The default
% for LserWeight is 0.56.  After travelling it for a distance to try to get better
% agreement between the nomogram based fundamentals and the tabulated fundamentals
% I (DHB) gave up and decided that using a single lambdaMax is as good as anything
% else I could come up with. If you are interested, see FitConeFundamentalsTest.
%
% NOTE 1: When we first implemented the CIE standard, adding this shifting feature
% seemed like a good idea to allow exploration of individual differences in photopigments.
% But, with 0 shift, none of the nomograms exactly reproduce the tabulated photopigment absorbance
% spectral sensitivities, and this is not so good.  We are phasing out our
% use of this feature in favor of simply shifting the tabulated
% photopigment absorbances, and indeed in favor of adopting the method
% published by Asano, Fairchild, & Blonde (2016), PLOS One, doi: 10.1371/journal.pone.0145671
% to tailor the CIE fundamentals to individual observers.  This is done by
% passing the argument indDiffParams, which is a structure as follows.
%   indDiffParams.dlens - deviation in % from CIE computed peak lens density
%   indDiffParams.dmac -  deviation in % from CIE peak macular pigment density
%   indDiffParams.dphotopigment - vector of deviations in % from CIE photopigment peak density.
%   indDiffParams.lambdaMaxShift - vector of values (in nm) to shift lambda max of each photopigment absorbance by.  
%   indDiffParams.shiftType - 'linear' (default) or 'log'. 'linear' gets the Asano et al. behavior
%
% You also can shift the absorbances along a wavenumber axis after you have
% obtained them.  To do this, pass argument lambdaMaxShift with the same
% number of entries as the number of absorbances that are used.
%
% The adjIndDiffParams output is a struct which is populated by ComputeRawConeFundamentals.
% It contains the actual parameter values for the parameters adjusted using the indDiffParams 
% input. It contains the following fields:
%    adjIndDiffParams.mac - the adjusted macular pigment transmittance as a function of wavelength
%                           as calculated in line 151 of ComputeRawConeFundamentals.
%    adjIndDiffParams.lens - the adjusted lens transmittance as a function of wavelength as calculated
%                            in line 41 of ComputeRawConeFundamentals.
%    adjIndDiffParams.dphotopigment - 3-vector of the adjusted photopigment axial density for
%                                     L, M and S cones (in that order), as calculated in lines
%                                     200-202 of ComputeRawConeFundamentals; or rods, as calculated
%                                     in line 216 of ComputeRawConeFundamentals if params.DORODS is true.
%    adjIndDiffParams.absorbance - Photopigment absorbance as given in line 188 of ComputeRawConeFundamentals
%    adjIndDiffParams.absorptance - Photopigment absorptance as given in line 230 of ComputeRawConeFundamentals
%    adjIndDiffParams.photopigmentQuantalEfficiency - Photopigment quantal efficiency as passed to ComptuteRawConeFundamentals.
%                                     We have this here in staticParams,
%                                     but useful to have in adjIndDiffParams as
%                                     well, sometimes.
%   adjIndDiffParams.ISdiameter - Inner segment diameter estimate (in microns), taken from
%                                    FillInPhotoreceptors return structure.
%
% For both adjIndDiffParams.mac and adjIndDiffParams.lens, the wavelength
% spacing is the same as in the S input variable of this function.
%
% The params and staticParams outputs are the argument strucutures that
% were passed to ComputeRawConeFundamentals by this routine to do the work.
% These can be useful if you'd like, say, to susequently use
% ComputeRawConeFundamentals to produce estimates for (e.g.) melanopsin or
% the rods, where you keep everything else as consistent as possible to
% what this routine does. Note that this is all a bit klugy for historical
% reasons, as there is redundancy between what you can/might do with
% adjIndDiffParams and with these two return outputs. In particular, these
% two return outputs would let you call ComputeRawConeFundamentals and get
% adjIndDiffParams directly from there.
%
% This function also has an option to compute rod spectral sensitivities, using
% the pre-retinal values that come from the CIE standard.  Set DORODS to true on
% call.  You then need to explicitly pass a single lambdaMax value.  You can
% also pass an optional rodAxialDensity value.  If you don't pass that, the
% routine uses the 'Alpern' estimate for 'Human'/'Rod' embodied in routine
% PhotopigmentAxialDensity.  The default nomogram for the rod spectral
% absorbance is 'StockmanSharpe', but you can override with any of the
% others available in routine PhotopigmentNomogram.  Use of this requires
% good choices for lambdaMax, rodAxialDensity, and the nomogram.  We are
% working on identifying those values more precisely.
%
% Finally, you can adjust the returned spectral sensitivities to account for
% the possibility that some of the pigment in the cones is bleached.  Pass
% a column vector with same length as number of spectral sensitivities beingt
% computed.  You need to estimate the fraction elsewhere.
%
% Relevant to individual differences, S & S (2000) estimate the wavelength difference
% between the ser/ala variants to be be 2.7 nm (ser longer).
%
% NOTE 2.  The CIE standard is specified for field sizes between 1 and 10
% degrees.  Our code will extrapolate using the given formulae to larger
% field sizes without complaining.  We think this is reasonable; see
% CIEConeFundamentalsFieldSizeTest and its header comments, but be aware
% that you have sailed into little charted territory if you do this.
%
% See also: ComputeRawConeFundamentals, CIEConeFundamentalsTest, CIEConeFundamentalsFieldSizeTest, 
% FitConeFundamentalsTest, FitConeFundamentalsWithNomogram, StockmanSharpeNomogram,
% ComputePhotopigmentBleaching.
%
% 8/13/11  dhb  Wrote it.
% 8/14/11  dhb  Clean up a little.
% 12/16/12 dhb, ms  Add rod option.
% 08/10/13 dhb  Test for consistency between what's returned by FillInPhotoreceptors and
%               what's returned by ComputeRawConeFundamentals.
% 05/24/14 dhb  Add fractionPigmentBleached optional arg.
% 05/26/14 dhb  Comment improvements.
% 02/08/16 dhb, ms  Add lambdaMaxShift argument.
%          ms   Don't do two way check when lambdaMax is shifted.
% 02/24/16 dhb, ms  Started to implement Asano et al. individual difference model
% 3/30/17  ms   Added output argument returning adjusted ind differences
% 8/1/17   dhb, ms  Added return of params and staticParams.

%% Are we doing rods rather than cones?
if (nargin < 8 || isempty(DORODS))
    DORODS = 0;
end

%% Check whether we'll adjust axial density for bleaching
if (nargin < 10 || isempty(fractionPigmentBleached))
    DOBLEACHING = 0;
else
    DOBLEACHING = 1;
end

%% Check for passed lambdaMaxShift
if (nargin < 11 || isempty(indDiffParams))
    params.indDiffParams = [];
else
    params.indDiffParams = indDiffParams;
end

%% Get some basic parameters.
%
%
% We start with default CIE parameters in 
% the photoreceptors structure, and then override
% as necessary.
% then override to match the CIE standard.
if (fieldSizeDegrees <= 4)
    whatCalc = 'CIE2Deg';
else
    whatCalc = 'CIE10Deg';
end
photoreceptors = DefaultPhotoreceptors(whatCalc);

%% Override default values so that FillInPhotoreceptors does
% our work for us.  The CIE standard uses field size, 
% age, and pupil diameter to computer other values.
% to compute other quantities.
photoreceptors.nomogram.S = S;
photoreceptors.fieldSizeDegrees = fieldSizeDegrees;
photoreceptors.pupilDiameter.value = pupilDiameterMM;
photoreceptors.ageInYears = ageInYears;

% Absorbance.  Use tabulated CIE values (which are in the
% default CIE photoreceptors structure) unless a nomogram and
% lambdaMax values are passed.
SET_ABSORBANCE = false;
if (nargin > 4 && ~isempty(lambdaMax))
    if (nargin < 6 || isempty(whichNomogram))
        whichNomogram = 'StockmanSharpe';
    end
    photoreceptors = rmfield(photoreceptors,'absorbance');
    photoreceptors.nomogram.source = whichNomogram;
    photoreceptors.nomogram.lambdaMax = lambdaMax;
    params.lambdaMax = lambdaMax;
    staticParams.whichNomogram = whichNomogram;
else
    % Absorbance is going to be specified directly.  We get
    % it after the call to FillInPhotoreceptors below,
    % which will convert a file containing the absorbance into
    % the needed data at the needed wavelength spacing.
    SET_ABSORBANCE = true;
end

%% Are we doing the rods?  In that case, a little more
% mucking is necessary.
if (DORODS)
    if (isempty(lambdaMax) || length(lambdaMax) ~= 1)
        error('When computing for rods, must specify exactly one lambda max');
    end
    photoreceptors.types = {'Rod'};
    photoreceptors.nomogram.lambdaMax = lambdaMax;
    photoreceptors.OSlength.source = 'None';
    photoreceptors.specificDensity.source = 'None';
    photoreceptors.axialDensity.source = 'Alpern';
    params.DORODS = true;
end

%% Pigment bleaching
%
% Hope for the best with respect to dimensionality of what is passed.
% FillInPhotoreceptors will throw an error if the dimension isn't
% matched to that of the axialDensity value field.
if (DOBLEACHING)
    photoreceptors.fractionPigmentBleached.value = fractionPigmentBleached;
end

%% Do the work.  Note that to modify this code, you'll want a good
% understanding of the order of precedence enforced by FillInPhotoreceptors.
% This is non-trivial, although the concept is that if a quantity that
% can be computed is specified directly in the passed structure is
% actually specified, the speciefied value overrides what could be computed.
photoreceptors = FillInPhotoreceptors(photoreceptors);
if (SET_ABSORBANCE)
    params.absorbance = photoreceptors.absorbance;
end

% Sometimes, we want to know things about how all the pieces came together.
% So we tuck awaay the photoreceptors structure, which has a lot of
% information in it.  We put this in a substructure just to make sure we
% don't collide in some structure name space somewhere else.
params.info.photoreceptors = photoreceptors;

%% Set up for call into the low level routine that computes the CIE fundamentals.
staticParams.S = photoreceptors.nomogram.S;
staticParams.fieldSizeDegrees = photoreceptors.fieldSizeDegrees;
staticParams.ageInYears = photoreceptors.ageInYears;
staticParams.pupilDiameterMM = photoreceptors.pupilDiameter.value;
staticParams.lensTransmittance = photoreceptors.lensDensity.transmittance;
staticParams.macularTransmittance = photoreceptors.macularPigmentDensity.transmittance;
staticParams.quantalEfficiency = photoreceptors.quantalEfficiency.value;
CHECK_FOR_AGREEMENT = true;
if (nargin < 7 || isempty(LserWeight))
    staticParams.LserWeight = 0.56;
else
    staticParams.LserWeight = LserWeight;
end
if (DORODS && nargin >= 9 && ~isempty(rodAxialDensity))
    params.axialDensity = rodAxialDensity;
    CHECK_FOR_AGREEMENT = false;
else
    params.axialDensity = photoreceptors.axialDensity.bleachedValue;
end

if (~isfield(params,'absorbance'))
    if (length(params.lambdaMax) ~= 3 && length(params.lambdaMax) ~= 1)
        CHECK_FOR_AGREEMENT = false;
    end
end

% Shift in lambda max bookkeeping.
if isfield(params,'indDiffParams')
    CHECK_FOR_AGREEMENT = false;
end

%% Drop into more general routine to compute
%
% See comment in ComputeRawConeFundamentals about the fact that
% we ought to unify this routine and what FillInPhotoreceptors does.
[T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations,adjIndDiffParams] = ComputeRawConeFundamentals(params,staticParams);
adjIndDiffParams.ISdiameter = photoreceptors.ISdiameter.value;

%% A little reality check.
%
% The call to FillInPhotoreceptors also computes what here is called
% T_quantal.  It is in the field effectiveAbsorptance.  For cases where
% we aren't playing games with the parameters after the call to 
% FillInPhotoreceptors, we can check for agreement.
if (CHECK_FOR_AGREEMENT)
    diffs = abs(T_quantalAbsorptions(:)-photoreceptors.effectiveAbsorptance(:));
    if (max(diffs(:)) > 1e-7)
        error('Two ways of computing absorption quantal efficiency referred to the cornea DO NOT AGREE');
    end
    diffs = abs(T_quantalIsomerizations(:)-photoreceptors.isomerizationAbsorptance(:));
    if (max(diffs(:)) > 1e-7)
        error('Two ways of computing isomerization quantal efficiency referred to the cornea DO NOT AGREE');
    end
end

end
