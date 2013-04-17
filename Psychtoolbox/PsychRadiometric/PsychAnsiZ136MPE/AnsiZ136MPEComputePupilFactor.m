function P = AnsiZ136MPEComputePupilFactor(stimulusDurationSeconds,stimulusWavelengthNm)
% P = AnsiZ136MPEComputePupilFactor(stimulusDurationSeconds,stimulusWavelengthNm)
%
% ****************************************************************************
% IMPORTANT: Before using the AnsiZ136 routines, please see the notes on usage
% and responsibility in PsychAnsiZ136MPE/Contents.m (type "help PsychAnsiZ136MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% Delori et al. (2007) JOSA A, pp. 1250-1265 define a pupil
% factor (their Table 2) that they use to modify the ANSI
% standard as a function of stimulus duration and wavelength.
% This is applied when converting between the corneal irradiance
% of a a beam overfilling the pupil and the allowable power
% in the pupil.
%
% This routine computes that factor P.
%
% 3/2/13  dhb  Wrote it.

%% Implement formula.

% If wavelength < 400, it's not defined.
if (stimulusWavelengthNm < 400)
    error('Pupil factor not defined for wavelengths less than 400 nm');
end

% If wavelength > 700 nm, it's unity.
if (stimulusWavelengthNm > 700)
    P = 1;
    return;
end

% If duration less than or equal to 0.07 seconds, it's 1
if (stimulusDurationSeconds <= 0.07)
    P = 1;
    return;
end

% Compute for durations between 0.07 and 0.7 seconds
if (stimulusDurationSeconds > 0.07 && stimulusDurationSeconds <= 0.7)
    if (stimulusWavelengthNm >= 400 && stimulusWavelengthNm < 600)
        P = (stimulusDurationSeconds/0.07)^0.75;
    elseif (stimulusWavelengthNm >= 600 && stimulusWavelengthNm <= 700)
        P = ((stimulusDurationSeconds/0.07)^0.75)*10^(0.0074*(700-stimulusWavelengthNm));
        if (P < 1)
            P = 1;
        end        
    else
        error('Logic error in routine, should never get here.');
    end
    
% Compute for durations greater than 0.7 seconds
elseif (stimulusDurationSeconds > 0.7)
    if (stimulusWavelengthNm >= 400 && stimulusWavelengthNm < 600)
        P = 5.44;
    elseif (stimulusWavelengthNm >= 600 && stimulusWavelengthNm <= 700)
        P = 10^(0.0074*(700-stimulusWavelengthNm));
    else
        error('Logic error in routine, should never get here.');
    end
    
% This should not happen
else
    error('Logic error in routine, should never get here.');
end

end

