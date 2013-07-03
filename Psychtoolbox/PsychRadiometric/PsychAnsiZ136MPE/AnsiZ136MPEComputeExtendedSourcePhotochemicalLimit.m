function [MPEPhotochemicalIntegratedRadiance_JoulesPerCm2Sr, ...
          MPEPhotochemicalRadiance_WattsPerCm2Sr, ...
          MPEPhotochemicalCornealIrradiance_WattsPerCm2, ...
          MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2] = ...
          AnsiZ136MPEComputeExtendedSourcePhotochemicalLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm,CONELIMITFLAG)
% [MPEPhotochemicalIntegratedRadiance_JoulesPerCm2Sr, ...
%  MPEPhotochemicalRadiance_WattsPerCm2S, ...
%  MPEPhotochemicalCornealIrradiance_WattsPerCm2, ...
%  MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2] = ...
%  AnsiZ136MPEComputeExtendedSourcePhotochemicalLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm,[CONELIMITFLAG])
%
% ****************************************************************************
% IMPORTANT: Before using the AnsiZ136 routines, please see the notes on usage
% and responsibility in PsychAnsiZ136MPE/Contents.m (type "help PsychAnsiZ136MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% Compute time photochemical MPE for an extended source
% ANSI Z136.1-2007, Table 5b, p. 75.
%
% Set CONELIMITFLAG to false (default true) to skip the asterisked 
% alternative computation desribed in Table 5 (see comments in code).
%
% 2/20/13  dhb  Wrote it.
% 3/2/13   dhb  Make limiting cone angle computation controllable.

%% Default arg 
if (nargin < 4 || isempty(CONELIMITFLAG))
    CONELIMITFLAG = true;
end

%% Check that we are in wavelength range
% Compute CB if so
%
% It's an error if we're below 400 nm, need to implement UV standard
if (stimulusWavelengthNm < 400)
    MPEPhotochemicalRadiance_WattsPerCm2Sr = NaN;
    MPEPhotochemicalIntegratedRadiance_JoulesPerCm2Sr = NaN;
    MPEPhotochemicalCornealIrradiance_WattsPerCm2 = NaN;
    MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2 = NaN;
    error('Photochemical MPE not yet implemented for wavelengths less than 400 nm');
end

% Above 600, other limits apply so we set photochemical limit to Inf
if (stimulusWavelengthNm >= 600)
    MPEPhotochemicalRadiance_WattsPerCm2Sr = Inf;
    MPEPhotochemicalIntegratedRadiance_JoulesPerCm2Sr = Inf;
    MPEPhotochemicalCornealIrradiance_WattsPerCm2 = Inf;
    MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2 = Inf;
    return
end
Cb = AnsiZ136MPEComputeCb(stimulusWavelengthNm);

%% Convert angle to mrad and get T2
stimulusSizeMrad = DegToMrad(stimulusSizeDeg);
stimulusAreaDeg2 = (pi/4)*stimulusSizeDeg^2;
T2Sec = AnsiZ136MPEComputeT2(stimulusSizeDeg);

%% Case of stimulus larger than 11 mradians
if (stimulusSizeMrad > 11)
    if (stimulusDurationSec < 0.7)
        % Below 0.7 seconds other limits apply and we set the photochemical limit to Inf
        MPEPhotochemicalRadiance_WattsPerCm2Sr = Inf;
        MPEPhotochemicalIntegratedRadiance_JoulesPerCm2Sr = Inf;
        MPEPhotochemicalCornealIrradiance_WattsPerCm2 = Inf;
        MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2 = Inf;
        return
    elseif (stimulusDurationSec >= 0.7 && stimulusDurationSec < 1e4)
        limitingConeAngleMrad = AnsiZ136MPEComputeLimitingConeAngle(stimulusDurationSec);
        
        % This next bit has to do with the asterisk for the case of angle greater than
        % 11 mrad in Table 5b.  This says that for this case, one can express the limit
        % in terms of the radiance measured through the limitingConeAngleMrad computed
        % just above.  When the stimulus is smaller than this limiting angle, this has
        % no effect.  But when the stimulus is larger, the effect is to move to a
        % measure that reflects the amount of light per unit area, I think.
        % 
        % The language in the table by the asterisk is not very clear -- it says "the limit may also be"
        % which I initially took to mean that this was another way of computing the same
        % thing.  But it doesn't produce the same answer, and the wording in the text itself
        % (top of p. 63) suggests that one should use the asterisked formula: "For photochemical
        % effects ... the MPE is provided in Table 5b as radiance or integrated radiance averaged
        % over a limiting cone angle.  Similarly the text in Section B.7.2 on p. 176.  
        %
        % The code as it is here and in parallel in the longer duration fork below has the
        % feature that it lets me reproduce Figures 10d and 10e of the standard.  Other 
        % variants of the same general idea, or omitting this part of the calculation,
        % cause a mismatch between what we compute for the cases of 10d and 10e and the
        % figures in the standard.  See AnsiZ136MPETest.
        if (limitingConeAngleMrad < stimulusSizeMrad || ~CONELIMITFLAG)
            MPEPhotochemicalIntegratedRadiance_JoulesPerCm2Sr = 100*Cb;
        else
            MPEPhotochemicalIntegratedRadiance_JoulesPerCm2Sr = 100*Cb*((limitingConeAngleMrad/stimulusSizeMrad)^2);
        end
        MPEPhotochemicalRadiance_WattsPerCm2Sr = MPEPhotochemicalIntegratedRadiance_JoulesPerCm2Sr/stimulusDurationSec;
        MPEPhotochemicalCornealIrradiance_WattsPerCm2 = ...
            RadianceAndDegrees2ToCornIrradiance(MPEPhotochemicalRadiance_WattsPerCm2Sr,stimulusAreaDeg2);
        MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2 = MPEPhotochemicalCornealIrradiance_WattsPerCm2*stimulusDurationSec;
    elseif (stimulusDurationSec >= 1e4 && stimulusDurationSec < 3e4)
        limitingConeAngleMrad = AnsiZ136MPEComputeLimitingConeAngle(stimulusDurationSec);
        if (limitingConeAngleMrad < stimulusSizeMrad || ~CONELIMITFLAG)
            MPEPhotochemicalRadiance_WattsPerCm2Sr = Cb*(1e-2);
        else
            MPEPhotochemicalRadiance_WattsPerCm2Sr = Cb*(1e-2)*((limitingConeAngleMrad/stimulusSizeMrad)^2);
        end
        MPEPhotochemicalIntegratedRadiance_JoulesPerCm2Sr = MPEPhotochemicalRadiance_WattsPerCm2Sr*stimulusDurationSec;
        MPEPhotochemicalCornealIrradiance_WattsPerCm2 = ...
            RadianceAndDegrees2ToCornIrradiance(MPEPhotochemicalRadiance_WattsPerCm2Sr,stimulusAreaDeg2);
        MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2 = MPEPhotochemicalCornealIrradiance_WattsPerCm2*stimulusDurationSec;
    else
        error('Limit not yet implemented for exposures greater than 3*10^4 seconds and size greater than 11 mrad');
    end
    
%% Case of stimulus smaller than or equal to 11 mradians
% The standard gives the limit as corneal irradiance for longer duration stimuli,
% and as radiant exposure for shorter durations, with the breakpoint being 100 seconds.
% The effect of this is to allow irradiance for shorter exposures, or equivalently
% higher radiant exposure for longer durations.  
%
% Since we have the stimulus duration, we can compute both quantities for both cases.
else
    if (stimulusDurationSec < 0.7)
        % Below 0.7 seconds other limits apply and we set the photochemical limit to Inf
        MPEPhotochemicalRadiance_WattsPerCm2Sr = Inf;
        MPEPhotochemicalIntegratedRadiance_JoulesPerCm2Sr = Inf;
        MPEPhotochemicalCornealIrradiance_WattsPerCm2 = Inf;
        MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2 = Inf;
        return
    elseif (stimulusDurationSec >= 0.7 && stimulusDurationSec < 100)
        MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2 = Cb * 10^-2;
        MPEPhotochemicalCornealIrradiance_WattsPerCm2 = MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2/stimulusDurationSec;
        MPEPhotochemicalRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPEPhotochemicalCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPEPhotochemicalIntegratedRadiance_JoulesPerCm2Sr = MPEPhotochemicalRadiance_WattsPerCm2Sr*stimulusDurationSec;
    elseif (stimulusDurationSec >= 100 && stimulusDurationSec < 3e4)
        MPEPhotochemicalCornealIrradiance_WattsPerCm2 = Cb * 10^-4;
        MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2 = MPEPhotochemicalCornealIrradiance_WattsPerCm2*stimulusDurationSec;
        MPEPhotochemicalRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPEPhotochemicalCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPEPhotochemicalIntegratedRadiance_JoulesPerCm2Sr = MPEPhotochemicalRadiance_WattsPerCm2Sr*stimulusDurationSec;
    else
        error('Limit not yet implemented for exposures greater than 3*10^4 seconds and size greater than 11 mrad');
    end
    
end

end


