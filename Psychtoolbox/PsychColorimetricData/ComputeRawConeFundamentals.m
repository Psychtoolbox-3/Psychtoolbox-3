function [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations] = ComputeRawConeFundamentals(params,staticParams)
% [T_quantalAbsorptionsNormalized,T_quantalAbsorptions,T_quantalIsomerizations] = ComputeRawConeFundamentals(params,staticParams)
%
% Function to compute normalized cone quantal sensitivities
% from underlying pieces and parameters.
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
% In the passed params structure, you can either pass
% the lambdaMax values for the photopigment, in which
% case the absorbance is computed from the specified
% nomogram, or you can pass the absorbance values
% directly in T_xxx format.  A typical choice in this
% case would be 10.^T_lgo10coneabsorbance_ss for the
% Stockman-Sharpe/CIE estimates.
%
% The typical use of this function is to be called by 
% ComputeCIEConeFundamentals, which sets up the
% passed structures acording to the CIE standard. 
% This routine, however, could in principle be used
% with a wide variety of choices of the component pieces.
%
% The other place this function is used is in attempts to
% fit a set of cone fundamentals by doing parameter search
% over the pieces.  It is this second use that led to the
% parameters being split over two separate structures, one
% that is held static during the fit and the other which
% contains the parameters that could be searched over by a calling
% routine.  For examples, see:
%   FitConeFundamentalsWithNomogram, FitConeFundamentalsTest.
% Looking around today (8/10/13), I (DHB) don't see any examples where
% this routine is called directly.  Rather, it is a subfunction
% called by ComputeCIEConeFundamentals.  The search routines above
% use ComputeCIEConeFundamentals, and only search over lambdaMax
% values.  I think I wrote this with the thought that I might one
% day search over more parameters, but lost motivation to carry it
% throught.
%
% The computations done here are very similar to those done in
% routine FillInPhotoreceptors.  I (DHB) think that I forgot about
% what FillInPhotoreceptors could do when I wrote this, which has
% led to some redundancy. FillInPhotoreceptors returns a field
% called effectiveAbsorbtance, which are the actual quantal efficiencies
% (not normalized) referred to the cornea.  FillInPhotoceptors also
% computes a field isomerizationAbsorbtance, which takes the quantal
% efficiency of isomerizations (probability of an isomerization given
% an absorption into acount.  This routine does not do that.
%
% It would probably be clever to unify the two sets of routines a
% little more, but we may never get to it.  The routine ComputeCIEConeFundamentals
% does contain a check that this routine and what is returned by FillInPhotoreceptors
% agree with each other, for cases where the parameters match.
%
% See also: ComputeCIEConeFundamentals, CIEConeFundamentalsTest, FitConeFundamentalsWithNomogram,
%           FitConeFundamentalsTest, DefaultPhotoreceptors, FillInPhotoreceptors.
%
% 8/12/11  dhb  Starting to make this actually work.
% 8/14/11  dhb  Change name, expand comments.
% 8/10/13  dhb  Expand comments.  Return unscaled quantal efficiencies too.

% Handle bad value
index = find(params.axialDensity <= 0.0001);
if (~isempty(index))
    params.axialDensity(index) = 0.0001;
end

% Handle optional values
if (~isfield(params,'extraLens'))
    params.extraLens = 0;
end
if (~isfield(params,'extraMac'))
    params.extraMac = 0;
end

% Prereceptor transmittance.  Sometimes adjustments of peak density
% recommended by various standards push the density less than
% zero at some wavelengths, so we need to make sure we never have
% transmittance greater than 1.
lens = 10.^-(-log10(staticParams.lensTransmittance)+params.extraLens);
lens(lens > 1) = 1;
mac = 10.^-(-log10(staticParams.macularTransmittance)+params.extraMac);
mac(mac > 1) = 1;

% Compute nomogram if absorbance wasn't passed directly.  We detect
% a direct pass by the existance of params.absorbance.
if (isfield(params,'absorbance'))
    absorbance = params.absorbance;
else
    absorbance = PhotopigmentNomogram(staticParams.S,params.lambdaMax,staticParams.whichNomogram);
end

% Compute absorbtance
%
% Handle special case where we deal with ser/ala polymorphism for L cone
if (size(absorbance,1) == 4)
    absorbtance = AbsorbanceToAbsorbtance(absorbance,staticParams.S,...
        [params.axialDensity(1) ; params.axialDensity(1) ; ...
        params.axialDensity(2) ; params.axialDensity(3)]);
elseif (size(absorbance,1) == 3)
    absorbtance = AbsorbanceToAbsorbtance(absorbance,staticParams.S,...
        [params.axialDensity(1) ; params.axialDensity(2) ; ...
        params.axialDensity(3)]);
elseif (size(absorbance,1) == 1 && params.DORODS)
    absorbtance = AbsorbanceToAbsorbtance(absorbance,staticParams.S,...
        params.axialDensity(1));
else
    error('Unexpected number of photopigment lambda max values passed');
end

%% Put together pre-receptor and receptor parts
for i = 1:size(absorbtance,1)
    absorbtance(i,:) = absorbtance(i,:) .* lens .* mac;
end

%% Put it into the right form
if (size(absorbtance,1) == 4)
    T_quantalAbsorptions = zeros(3,staticParams.S(3));
    T_quantalAbsorptions(1,:) = staticParams.LserWeight*absorbtance(1,:) + ...
        (1-staticParams.LserWeight)*absorbtance(2,:);
    T_quantalAbsorptions(2,:) = absorbtance(3,:);
    T_quantalAbsorptions(3,:) = absorbtance(4,:);
elseif (size(absorbtance,1) == 3)
    T_quantalAbsorptions = zeros(3,staticParams.S(3));
    T_quantalAbsorptions(1,:) = absorbtance(1,:);
    T_quantalAbsorptions(2,:) = absorbtance(2,:);
    T_quantalAbsorptions(3,:) = absorbtance(3,:);
elseif (size(absorbtance,1) == 1 && params.DORODS)
    T_quantalAbsorptions = zeros(1,staticParams.S(3));
    T_quantalAbsorptions(1,:) = absorbtance(1,:);
else
    error('Unexpected number of photopigment lambda max values passed');
end

%% Normalize to max of one for each receptor, and also compute isomerization quantal efficiency.
for i = 1:size(T_quantalAbsorptions,1)
    T_quantalIsomerizations = T_quantalAbsorptions*staticParams.quantalEfficiency(i);
    T_quantalAbsorptionsNormalized(i,:) = T_quantalAbsorptions(i,:)/max(T_quantalAbsorptions(i,:));
end


