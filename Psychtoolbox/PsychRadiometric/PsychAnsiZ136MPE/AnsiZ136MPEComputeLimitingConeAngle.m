function limitingConeAngleMrad = AnsiZ136MPEComputeLimitingConeAngle(stimulusDurationSec)
% limitingConeAngleMrad = AnsiZ136MPEComputeLimitingConeAngle(stimulusDurationSec)
%
% ****************************************************************************
% IMPORTANT: Before using the AnsiZ136 routines, please see the notes on usage
% and responsibility in PsychAnsiZ136MPE/Contents.m (type "help PsychAnsiZ136MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% Compute the limiting cone angle, ANSI Z136.1-2007, Table 5, p. 75.
%
% This is only defined for durations between 0.7 sec and 3*10^4 sec.
%
% 2/20/13  dhb  Wrote it.

%% Implement formula.  Factor only defined between 400 and 600 nm
if (stimulusDurationSec < 0.7)
    error('Limiting cone angle not defined for durations less than 0.7 sec');
elseif (stimulusDurationSec >= 0.7 && stimulusDurationSec < 100)
    limitingConeAngleMrad = 11;
elseif (stimulusDurationSec >= 100 && stimulusDurationSec < 1e4)
   limitingConeAngleMrad = 1.1 * (stimulusDurationSec^0.5);
elseif (stimulusDurationSec >= 1e4 && stimulusDurationSec < 3e4)
   limitingConeAngleMrad = 110;
else
    error('Limiting cone angle not defined for durations greater than 3e4 sec');
end

end


