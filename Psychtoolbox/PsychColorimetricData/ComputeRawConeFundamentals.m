function T_quantal = ComputeRawConeFundamentals(params,staticParams)
% T_quantal = ComputeRawConeFundamentals(params,staticParams)
%
% Function to compute normalized cone quantal sensitivities
% from underlying pieces and parameters.
%
% In the passed params structure, you can either pass
% the lambdaMax values for the photopigment, in which
% case the absorbance is computed from the specified
% nomogram, or you can pass the absorbance values
% directly in T_xxx format.  A typical choice in this
% case would be 10.^T_lgo10coneabsorbance_ss for the
% Stockman-Sharpe/CIE estimates.
%
% It's not clear that it is necessary to normalize; eventually
% this should probably be changed to provide fundamentals that
% give actual quantal efficiency.
%
% The typical use of this function is to be called by 
% ComputeCIEConeFundamentals, which sets up the
% passed structures acording to the CIE standard. 
% This routine, however, could in principle be used
% with a wide variety of choices of the component pieces.
%
% 8/12/11  dhb  Starting to make this actually work.
% 8/14/11  dhb  Change name, expand comment.

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
else
    error('Unexpected number of photopigment lambda max values passed');
end

%% Put together pre-receptor and receptor parts
for i = 1:size(absorbtance,1)
    absorbtance(i,:) = absorbtance(i,:) .* lens .* mac;
end

% Put it into the right form
T_quantal = zeros(3,staticParams.S(3));

if (size(absorbtance,1) == 4)
    T_quantal(1,:) = staticParams.LserWeight*absorbtance(1,:) + ...
        (1-staticParams.LserWeight)*absorbtance(2,:);
    T_quantal(2,:) = absorbtance(3,:);
    T_quantal(3,:) = absorbtance(4,:);
elseif (size(absorbtance,1) == 3)
    T_quantal(1,:) = absorbtance(1,:);
    T_quantal(2,:) = absorbtance(2,:);
    T_quantal(3,:) = absorbtance(3,:);
else
    error('Unexpected number of photopigment lambda max values passed');
end

% Normalize to max of one for each receptor
for i = 1:size(T_quantal,1)
    T_quantal(i,:) = T_quantal(i,:)/max(T_quantal(i,:));
end
