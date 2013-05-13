function Cc = AnsiZ136MPEComputeCc(stimulusWavelengthNm)
% Cc = AnsiZ136MPEComputeCc(stimulusSizeDeg)
%
% ****************************************************************************
% IMPORTANT: Before using the AnsiZ136 routines, please see the notes on usage
% and responsibility in PsychAnsiZ136MPE/Contents.m (type "help PsychAnsiZ136MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% Compute constant Cc, ANSI Z136.1-2007, Table 6, p. 76.
%
% This is only defined between 400 and 1400 nm.
%
% 2/20/13  dhb  Wrote it.

%% Implement formula.
if (stimulusWavelengthNm < 1050)
    error('Cc not defined for wavelengths less than 1050 nm');
elseif (stimulusWavelengthNm >=  1050 && stimulusWavelengthNm < 1150)
    Cc = 1;
elseif (stimulusWavelengthNm >= 1150 && stimulusWavelengthNm < 1200)
    Cc = 10^(18*(stimulusWavelengthNm/1000 - 1.150));
elseif (stimulusWavelengthNm >=  1200 && stimulusWavelengthNm < 1400)
    Cc = 8;
else
    error('Cc not defined for wavelengths greater than or equal to 1400 nm');
end

end


