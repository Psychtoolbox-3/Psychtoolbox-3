function T2Sec = AnsiZ136MPEComputeT2(stimulusSizeDeg)
% T2Sec = AnsiZ136MPEComputeT2(stimulusSizeDeg)
%
% ****************************************************************************
% IMPORTANT: Before using the AnsiZ136 routines, please see the notes on usage
% and responsibility in PsychAnsiZ136MPE/Contents.m (type "help PsychAnsiZ136MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% Compute time T2, ANSI Z136.1-2007, Table 6, p. 76
%
% 2/20/13  dhb  Wrote it.

%% Convert angle to mrad
stimulusSizeMrad = DegToMrad(stimulusSizeDeg);

%% Implement the table with its various special cases
if (stimulusSizeMrad < 1.5)
    T2Sec = 10;
elseif (stimulusSizeMrad > 100)
    T2Sec = 100;
else
    T2Sec = 10*10^((stimulusSizeMrad-1.5)/98.5);
end

end

