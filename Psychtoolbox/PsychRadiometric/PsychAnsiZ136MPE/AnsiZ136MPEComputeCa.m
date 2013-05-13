function Ca = AnsiZ136MPEComputeCa(stimulusWavelengthNm)
% Ca = AnsiZ136MPEComputeCa(stimulusSizeDeg)
%
% ****************************************************************************
% IMPORTANT: Before using the AnsiZ136 routines, please see the notes on usage
% and responsibility in PsychAnsiZ136MPE/Contents.m (type "help PsychAnsiZ136MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% Compute constant Ca, ANSI Z136.1-2007, Table 6, p. 76.
%
% This is only defined between 400 and 1400 nm.
%
% 2/20/13  dhb  Wrote it.

%% Implement formula.
if (stimulusWavelengthNm < 400)
    error('Ca not defined for wavelengths less than 400 nm');
elseif (stimulusWavelengthNm >=  400 && stimulusWavelengthNm < 700)
    Ca = 1;
elseif (stimulusWavelengthNm >=  700 && stimulusWavelengthNm < 1050)
    Ca = 10^(2*(stimulusWavelengthNm/1000-0.700));
elseif (stimulusWavelengthNm >= 1050 && stimulusWavelengthNm < 1400)
    Ca = 5;
else
    error('Ca not defined for wavelengths greater than or equal to 1400 nm');
end

end


