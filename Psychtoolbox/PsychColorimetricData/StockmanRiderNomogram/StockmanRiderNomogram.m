function T_absorbance = StockmanRiderNomogram(S, lambdaMax)
% T_absorbance = StockmanRiderNomogram(S, lambdaMax)
%
% Compute normalized absorbance using the Stockman-Rider Fourier polynomial
% nomogram, as described in:
%   Stockman, A. & Rider, A. T. (2023). Formulae for generating standard and
%   individual human cone spectral sensitivities. Color Research and Application,
%   48, 818-840.
%
% The appropriate L, M, or S cone Fourier template is selected for each
% lambdaMax value based on which cone type's range it falls in, using split
% points halfway between the nominal template peaks:
%   S template peak: 416.9 nm
%   M template peak: 529.8 nm
%   L template peak: 551.9 nm
%   S/M split:       473.35 nm  (midpoint of 416.9 and 529.8)
%   M/L split:       540.85 nm  (midpoint of 529.8 and 551.9)
%
% If lambdaMax is a 3-element column vector it is assumed to be [L; M; S]
% in that order. An error is thrown if the values do not fall in the
% expected ranges, since that likely indicates a calling error.
%
% The result is in quantal units (absorbance in the sense used by PTB
% nomograms). The output is normalized to a peak of 1 for each row,
% if wavelength spacing is <= 1 nm.
%
% To get sensitivity in energy units, apply EnergyToQuanta().
%
% Various utility routines from Stockman-Rider are include, with sr
% prepended to all of their names to avoid namespace collision with
% PTB routines that do similar things.
%
% Good lambda max values (L, M, S): 551.9, 529.8, 416.9 nm.
% These are the nominal template peaks built into the Stockman-Rider
% nomogram itself and reproduce the standard CIE photopigment absorbances
% well.
%
% See also: PhotopigmentNomogram, srLconelog, srMconelog, srSconelog,
%           srLMSconelog, srDemo

% History:
%   2026-04-09  dhb  Wrote it. The underlying code was translated by DHB 
%                    from the Stockman-Rider python version, with help from
%                    AI.

% Nominal template peaks and split points
Llmax_template = 551.9;
Mlmax_template = 529.8;
Slmax_template = 416.9;
splitSM = (Slmax_template + Mlmax_template) / 2;   % 473.35 nm
splitML = (Mlmax_template + Llmax_template) / 2;   % 540.85 nm

% Get wavelength vector from S spec
S = MakeItS(S);
wls = MakeItWls(S);
nWls = length(wls);
nT   = length(lambdaMax);

% Special check for 3-element input assumed to be [L; M; S]
if (nT == 3)
    if (lambdaMax(1) <= splitML)
        error(['StockmanRiderNomogram: lambdaMax is 3 elements, assumed [L; M; S] order.\n' ...
               'But lambdaMax(1) = %.1f nm is not in the L range (> %.2f nm).\n' ...
               'If you intended a different ordering or a non-LMS set, pass lambdaMax\n' ...
               'as a scalar or 2-element vector to bypass this check.'], ...
               lambdaMax(1), splitML);
    end
    if (lambdaMax(2) <= splitSM || lambdaMax(2) > splitML)
        error(['StockmanRiderNomogram: lambdaMax is 3 elements, assumed [L; M; S] order.\n' ...
               'But lambdaMax(2) = %.1f nm is not in the M range (%.2f to %.2f nm).\n' ...
               'If you intended a different ordering or a non-LMS set, pass lambdaMax\n' ...
               'as a scalar or 2-element vector to bypass this check.'], ...
               lambdaMax(2), splitSM, splitML);
    end
    if (lambdaMax(3) >= splitSM)
        error(['StockmanRiderNomogram: lambdaMax is 3 elements, assumed [L; M; S] order.\n' ...
               'But lambdaMax(3) = %.1f nm is not in the S range (< %.2f nm).\n' ...
               'If you intended a different ordering or a non-LMS set, pass lambdaMax\n' ...
               'as a scalar or 2-element vector to bypass this check.'], ...
               lambdaMax(3), splitSM);
    end
end

% Compute absorbance for each lambdaMax
T_absorbance = zeros(nT, nWls);
for i = 1:nT
    lmax = lambdaMax(i);
    if (lmax < splitSM)
        % S cone template
        shift = lmax - Slmax_template;
        logAbs = srSconelog(wls, shift);
    elseif (lmax <= splitML)
        % M cone template
        shift = lmax - Mlmax_template;
        logAbs = srMconelog(wls, shift);
    else
        % L cone template
        shift = lmax - Llmax_template;
        logAbs = srLconelog(wls, shift);
    end

    % Convert log10 to linear
    T_absorbance(i,:) = 10.^(logAbs(:)');
    
    % If wavelength spacing is sufficiently fine, normalize explicitly
    % to a max of 1.  Seems like a good idea since this should always
    % be true of absorbance functions.  Don't want to do it if
    % wavelength spacing is coarse, because the peak may be between the
    % samples.
    if (S(2) <= 1)
        T_absorbance(i,:) = T_absorbance(i,:)/max(T_absorbance(i,:));
    end
end

end
