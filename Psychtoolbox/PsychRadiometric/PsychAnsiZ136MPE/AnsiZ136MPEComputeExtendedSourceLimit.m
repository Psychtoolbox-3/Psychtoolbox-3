function [MPELimitIntegratedRadiance_JoulesPerCm2Sr, ...
    MPELimitRadiance_WattsPerCm2Sr, ...
    MPELimitCornealIrradiance_WattsPerCm2, ...
    MPELimitCornealRadiantExposure_JoulesPerCm2] = ...
    AnsiZ136MPEComputeExtendedSourceLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm,CONELIMITFLAG)
% [MPELimitIntegratedRadiance_JoulesPerCm2Sr, ...
%  MPELimitRadiance_WattsPerCm2S, ...
%  MPELimitCornealIrradiance_WattsPerCm2, ...
%  MPELimitCornealRadiantExposure_JoulesPerCm2] = ...
%  AnsiZ136MPEComputeExtendedSourceLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm,[CONELIMITFLAG])
%
% ****************************************************************************
% IMPORTANT: Before using the AnsiZ136 routines, please see the notes on usage
% and responsibility in PsychAnsiZ136MPE/Contents.m (type "help PsychAnsiZ136MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% Compute thermal MPE for an extended source ANSI Z136.1-2007, Table 5b, p. 75
%
% Set CONELIMITFLAG to false (default true) to skip the asterisked 
% alternative computation for the photochemical limit
% desribed in Table 5 (see comments in code).
%
% 2/20/13  dhb  Wrote it.
% 3/2/13   dhb  Make limiting cone angle computation controllable.

%% Default arg 
if (nargin < 4 || isempty(CONELIMITFLAG))
    CONELIMITFLAG = true;
end

%% Convert angle to mrad and get T2
stimulusSizeMrad = DegToMrad(stimulusSizeDeg);
stimulusAreaDeg2 = (pi/4)*stimulusSizeDeg^2;
T2Sec = AnsiZ136MPEComputeT2(stimulusSizeDeg);
Ce = AnsiZ136MPEComputeCe(stimulusSizeDeg);
Ca = AnsiZ136MPEComputeCa(stimulusWavelengthNm);

%% Outer conditional is on wavelength range, inner conditional
% by duration.
if (stimulusWavelengthNm < 400)
    error('MPE not yet implemented for wavelengths less than 400 nm');
    
elseif (stimulusWavelengthNm >= 400 && stimulusWavelengthNm < 700)
    % 400-700 nm.  Now split by time.
    if (stimulusDurationSec < 1e-13)
        error('Limit not yet implemented for exposures less than 1e-13 seconds');
    elseif (stimulusDurationSec >= 1e-13 && stimulusDurationSec < 1e-11)
        % At this duration, limit s directly specified
        MPELimitCornealRadiantExposure_JoulesPerCm2 = 1.5*Ce*1e-8;
        MPELimitCornealIrradiance_WattsPerCm2 = MPELimitCornealRadiantExposure_JoulesPerCm2/stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec;
    elseif (stimulusDurationSec >= 1e-11 && stimulusDurationSec < 1e-9)
        % At this duration, limit s directly specified
        MPELimitCornealRadiantExposure_JoulesPerCm2 = 2.7*Ce*(stimulusDurationSec^0.75);
        MPELimitCornealIrradiance_WattsPerCm2 = MPELimitCornealRadiantExposure_JoulesPerCm2/stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec;
    elseif (stimulusDurationSec >= 1e-9 && stimulusDurationSec < 18*1e-6)
        % At this duration, limit s directly specified
        MPELimitCornealRadiantExposure_JoulesPerCm2 = 5*Ce*1e-7;
        MPELimitCornealIrradiance_WattsPerCm2 = MPELimitCornealRadiantExposure_JoulesPerCm2/stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec;
    elseif (stimulusDurationSec >= 18*1e-6 && stimulusDurationSec < 0.7)
        % At this duration, limit s directly specified
        MPELimitCornealRadiantExposure_JoulesPerCm2 = 1.8*Ce*(stimulusDurationSec^0.75)*(1e-3);
        MPELimitCornealIrradiance_WattsPerCm2 = MPELimitCornealRadiantExposure_JoulesPerCm2/stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec; 
    elseif (stimulusDurationSec >= 0.7 && stimulusDurationSec < T2Sec)
        % In this range, we compute thermal limit, and compare it to the photochemical limit.  We
        % take whichever is smaller.
        MPEThermalCornealRadiantExposure_JoulesPerCm2 = 1.8*Ce*(stimulusDurationSec^0.75)*(1e-3);
        
        [~, ~, ~, MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2] = ...
          AnsiZ136MPEComputeExtendedSourcePhotochemicalLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm,CONELIMITFLAG);
        if (MPEThermalCornealRadiantExposure_JoulesPerCm2 <= MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2)
            MPELimitCornealRadiantExposure_JoulesPerCm2 = MPEThermalCornealRadiantExposure_JoulesPerCm2;
        else
            MPELimitCornealRadiantExposure_JoulesPerCm2 = MPEPhotochemicalCornealRadiantExposure_JoulesPerCm2;
        end
        MPELimitCornealIrradiance_WattsPerCm2 = MPELimitCornealRadiantExposure_JoulesPerCm2/stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec;
    elseif (stimulusDurationSec >= T2Sec && stimulusDurationSec < 3e4)
        % In this range, we compute thermal limit, and compare it to the photochemical limit.  We
        % take whichever is smaller.
        MPEThermalCornealIrradiance_WattsPerCm2 = 1.8*Ce*(T2Sec^(-0.25))*(1e-3);
        [~, ~, MPEPhotochemicalCornealIrradiance_WattsPerCm2, ~] = ...
          AnsiZ136MPEComputeExtendedSourcePhotochemicalLimit(stimulusDurationSec,stimulusSizeDeg,stimulusWavelengthNm,CONELIMITFLAG);
        if (MPEThermalCornealIrradiance_WattsPerCm2 <= MPEPhotochemicalCornealIrradiance_WattsPerCm2)
            MPELimitCornealIrradiance_WattsPerCm2 = MPEThermalCornealIrradiance_WattsPerCm2;
        else
            MPELimitCornealIrradiance_WattsPerCm2 = MPEPhotochemicalCornealIrradiance_WattsPerCm2;
        end
        MPELimitCornealRadiantExposure_JoulesPerCm2 = MPELimitCornealIrradiance_WattsPerCm2*stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec;
    else
        error('Limit not yet implemented for exposures greater than 3*10^4 seconds');
    end
    
elseif (stimulusWavelengthNm >= 700 && stimulusWavelengthNm < 1050)
    % 700-1050.  Now split by time.
    if (stimulusDurationSec < 1e-13)
        error('Limit not yet implemented for exposures less than 1e-13 seconds');
    elseif (stimulusDurationSec >= 1e-13 && stimulusDurationSec < 1e-11)
        % At this duration, limit s directly specified
        MPELimitCornealRadiantExposure_JoulesPerCm2 = 1.5*Ca*Ce*1e-8;
        MPELimitCornealIrradiance_WattsPerCm2 = MPELimitCornealRadiantExposure_JoulesPerCm2/stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec;
    elseif (stimulusDurationSec >= 1e-11 && stimulusDurationSec < 1e-9)
        % At this duration, limit s directly specified
        MPELimitCornealRadiantExposure_JoulesPerCm2 = 2.7*Ca*Ce*(stimulusDurationSec^0.75);
        MPELimitCornealIrradiance_WattsPerCm2 = MPELimitCornealRadiantExposure_JoulesPerCm2/stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec;
    elseif (stimulusDurationSec >= 1e-9 && stimulusDurationSec < 18*1e-6)
        % At this duration, limit s directly specified
        MPELimitCornealRadiantExposure_JoulesPerCm2 = 5*Ca*Ce*1e-7;
        MPELimitCornealIrradiance_WattsPerCm2 = MPELimitCornealRadiantExposure_JoulesPerCm2/stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec;
    elseif (stimulusDurationSec >= 18*1e-6 && stimulusDurationSec < T2Sec)
        MPELimitCornealRadiantExposure_JoulesPerCm2 = 1.8*Ca*Ce*(stimulusDurationSec^0.75)*(1e-3);
        MPELimitCornealIrradiance_WattsPerCm2 = MPELimitCornealRadiantExposure_JoulesPerCm2/stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec;
    elseif (stimulusDurationSec >= T2Sec && stimulusDurationSec < 3e4)
        MPELimitCornealIrradiance_WattsPerCm2 = 1.8*Ca*Ce*(T2Sec^(-0.25))*(1e-3);
        MPELimitCornealRadiantExposure_JoulesPerCm2 = MPELimitCornealIrradiance_WattsPerCm2*stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec;
    else
        error('Limit not yet implemented for exposures greater than 3*10^4 seconds');
    end
    
elseif (stimulusWavelengthNm >= 1050 && stimulusWavelengthNm < 1400)
    % 1050-1400.  Now split by time.
    Cc = AnsiZ136MPEComputeCc(stimulusWavelengthNm);
    if (stimulusDurationSec < 1e-13)
        error('Limit not yet implemented for exposures less than 1e-13 seconds');
    elseif (stimulusDurationSec >= 1e-13 && stimulusDurationSec < 1e-11)
        % At this duration, limit s directly specified
        MPELimitCornealRadiantExposure_JoulesPerCm2 = 1.5*Cc*Ce*1e-7;
        MPELimitCornealIrradiance_WattsPerCm2 = MPELimitCornealRadiantExposure_JoulesPerCm2/stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec;
    elseif (stimulusDurationSec >= 1e-11 && stimulusDurationSec < 1e-9)
        % At this duration, limit s directly specified
        MPELimitCornealRadiantExposure_JoulesPerCm2 = 27*Cc*Ce*(stimulusDurationSec^0.75);
        MPELimitCornealIrradiance_WattsPerCm2 = MPELimitCornealRadiantExposure_JoulesPerCm2/stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec;
    elseif (stimulusDurationSec >= 1e-9 && stimulusDurationSec < 50*1e-6)
        % At this duration, limit s directly specified
        MPELimitCornealRadiantExposure_JoulesPerCm2 = 5*Cc*Ce*1e-6;
        MPELimitCornealIrradiance_WattsPerCm2 = MPELimitCornealRadiantExposure_JoulesPerCm2/stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec;

    elseif (stimulusDurationSec >= 50*1e-6 && stimulusDurationSec < T2Sec)
        MPELimitCornealRadiantExposure_JoulesPerCm2 = 9*Cc*Ce*(stimulusDurationSec^0.75)*(1e-3);
        MPELimitCornealIrradiance_WattsPerCm2 = MPELimitCornealRadiantExposure_JoulesPerCm2/stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec;
    elseif (stimulusDurationSec >= T2Sec && stimulusDurationSec < 3e4)
        MPELimitCornealIrradiance_WattsPerCm2 = 9*Cc*Ce*(T2Sec^(-0.25))*(1e-3);
        MPELimitCornealRadiantExposure_JoulesPerCm2 = MPELimitCornealIrradiance_WattsPerCm2*stimulusDurationSec;
        MPELimitRadiance_WattsPerCm2Sr = CornIrradianceAndDegrees2ToRadiance(MPELimitCornealIrradiance_WattsPerCm2,stimulusAreaDeg2);
        MPELimitIntegratedRadiance_JoulesPerCm2Sr = MPELimitRadiance_WattsPerCm2Sr*stimulusDurationSec;
    else
        error('Limit not yet implemented for exposures greater than 3*10^4 seconds');
    end
    
else
    % Off the long wavelength end
    error('MPE not yet implemented for wavelengths greater than 1400 nm');
end


end


