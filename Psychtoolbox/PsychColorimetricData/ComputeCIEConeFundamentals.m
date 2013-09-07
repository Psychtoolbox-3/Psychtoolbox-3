function [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations] = ComputeCIEConeFundamentals(S,fieldSizeDegrees,ageInYears,pupilDiameterMM,lambdaMax,whichNomogram,LserWeight,DORODS,rodAxialDensity)
% [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations] = ComputeCIEConeFundamentals(S,fieldSizeDegrees,ageInYears,pupilDiameterMM,[lambdaMax],[whichNomogram],[LserWeight],[DORODS],[rodAxialDensity])
%
% Function to compute normalized cone quantal sensitivities
% from underlying pieces, as specified in CIE 170-1:2006.
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
% Note that this routine returns quantal sensitivities.  You
% may want energy sensitivities.  In that case, use EnergyToQuanta to convert
%   T_energy = EnergyToQuanta(S,T_quantal')'
% and then renormalize.  (You call EnergyToQuanta because you're converting
% sensitivities, which go the opposite directoin from spectra.)
%
% The routine also returns two quantal sensitivity functions.  The first gives
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
% Although this routine will compute something over any wavelength
% range, I'd (DHB) recommend not going lower than 390 or above about 780 without
% thinking hard about how various pieces were extrapolated out of the range
% that they are specified in the standard.  Indeed, the lens optical
% density measurements only go down to 400 nm and these are extropolated
% to go below 400.
%
% This routine will compute from tabulated absorbance or absorbance based on a nomogram, where
% whichNomogram can be any source understood by the routine
% PhotopigmentNomogram.  To obtain the nomogram behavior, pass a lambdaMax vector.
% You can then also optionally pass a nomogram source (default: StockmanSharpe).
%
% The nominal values of lambdaMax to fit the CIE 2-degree fundamentals with the
% Stockman-Sharpe nomogram are 558.9, 530.3, and 420.7 nm for the LMS cones respectively.
% These in fact do a reasonable job of reconstructing the CIE 2-degree fundamentals, although
% there are small deviations from what you get if you simply read in the tabulated cone
% absorbances.  Thus starting with these as nominal values and shifting is a reasonable way to
% produce fundamentals tailored to observers with different known photopigments.
% 
% Relevant to that enterprise, S & S (2000) estimate the wavelength difference
% between the ser/ala variants to be be 2.7 nm (ser longer).
%
% If you pass lambaMax and its length is 4, then first two values are treated as
% the peak wavelengths of the ser/ala variants of the L cone pigment, and these
% are then weighted according to LserWeight and (1-LserWeight).  The default
% for LserWeight is 0.56.  After travelling it for a distance to try to get better
% agreement between the nomogram based fundamentals and the tabulated fundamentals
% I (DHB) gave up and decided that using a single lambdaMax is as good as anything
% else I could come up with. If you are interested, see FitConeFundametnalsTest.
%
% This function also has an option to compute rod spectral sensitivities, using
% the pre-retinal values that come from the CIE standard.  Set DORODS = true on
% call.  You then need to explicitly pass a single lambdaMax value.  You can
% also pass an optional rodAxialDensity value.  If you don't pass that, the
% routine uses the 'Alpern' estimate for 'Human'/'Rod' embodied in routine
% PhotopigmentAxialDensity.  The default nomogram for the rod spectral
% absorbance is 'StockmanSharpe', but you can override with any of the
% others available in routine PhotopigmentNomogram.  Use of this requires
% good choices for lambdaMax, rodAxialDensity, and the nomogram.  We are
% working on identifying those values more precisely.
%
% See also: ComputeRawConeFundamentals, CIEConeFundamentalsTest, 
% FitConeFundamentalsTest, FitConeFundamentalsWithNomogram, StockmanSharpeNomogram.
%
% 8/13/11  dhb  Wrote it.
% 8/14/11  dhb  Clean up a little.
% 12/16/12 dhb, ms  Add rod option.
% 08/10/13 dhb  Test for consistency between what's returned by FillInPhotoreceptors and
%               what's returned by ComputeRawConeFundamentals.

%% Are we doing rods rather than cones?
if (nargin < 8 || isempty(DORODS))
    DORODS = 0;
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

%% Do the work.  Note that to modify this code, you'll want a good
% understanding of the order of precedence enforced by FillInPhotoreceptors.
% This is non-trivial, although the concept is that if a quantity that
% can be computed is specified directly in the passed structure is
% actually specified, the speciefied value overrides what could be computed.
photoreceptors = FillInPhotoreceptors(photoreceptors);
if (SET_ABSORBANCE)
    params.absorbance = photoreceptors.absorbance;
end

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
    params.axialDensity = photoreceptors.axialDensity.value;
end

if (~isfield(params,'absorbance'))
    if (length(params.lambdaMax) ~= 3 & length(params.lambdaMax) ~= 1)
        CHECK_FOR_AGREEMENT = false;
    end
end

%% Drop into more general routine to compute
[T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations] = ComputeRawConeFundamentals(params,staticParams);

%% A little reality check.
%
% The call to FillInPhotoreceptors also computes what here is called
% T_quantal.  It is in the field effectiveAbsorbtance.  For cases where
% we aren't playing games with the parameters after the call to 
% FillInPhotoreceptors, we can check for agreement.
if (CHECK_FOR_AGREEMENT)
    diffs = abs(T_quantalAbsorptions(:)-photoreceptors.effectiveAbsorbtance(:));
    if (max(diffs(:)) > 1e-7)
        error('Two ways of computing absorption quantal efficiency referred to the cornea DO NOT AGREE\n');
    end
    diffs = abs(T_quantalIsomerizations(:)-photoreceptors.isomerizationAbsorbtance(:));
    if (max(diffs(:)) > 1e-7)
        error('Two ways of computing isomerization quantal efficiency referred to the cornea DO NOT AGREE\n');
    end
end

end

