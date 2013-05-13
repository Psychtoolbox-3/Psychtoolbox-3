function Cb = AnsiZ136MPEComputeCb(stimulusWavelengthNm)
% Cb = AnsiZ136MPEComputeCb(stimulusWavelengthNm)
%
% ****************************************************************************
% IMPORTANT: Before using the AnsiZ136 routines, please see the notes on usage
% and responsibility in PsychAnsiZ136MPE/Contents.m (type "help PsychAnsiZ136MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% Compute constant Cb, ANSI Z136.1-2007, Table 6, p. 76.
%
% This is only defined between 400 and 600 nm.
%
% 2/20/13  dhb  Wrote it.

%% Implement formula.  Factor only defined between 400 and 600 nm
if (stimulusWavelengthNm >= 400 && stimulusWavelengthNm < 450)
    Cb = 1;
elseif (stimulusWavelengthNm >= 450 && stimulusWavelengthNm < 600)
    Cb = 10^(20*(stimulusWavelengthNm/1000-0.450));
else
    Cb = NaN;
end

end


