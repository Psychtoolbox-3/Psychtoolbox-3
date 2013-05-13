function AnsiZ136MPEPrintConditionalComparison(string,compareFormat,mainVal,compareVal,LOG10FLAG)
% AnsiZ136MPEPrintConditionalComparison(string,compareFormat,mainVal,compareVal,[LOG10FLAG])
%
% ****************************************************************************
% IMPORTANT: Before using the AnsiZ136 routines, please see the notes on usage
% and responsibility in PsychAnsiZ136MPE/Contents.m (type "help PsychAnsiZ136MPE"
% at the Matlab prompt.
% ****************************************************************************
%
% Printout utility for test scripts.
%
% 3/2/13  dhb  Wrote it.

% Optional arg
if (nargin < 5 || isempty(LOG10FLAG))
    LOG10FLAG = false;
end

% Handle log10 while preserving conditional info
if (LOG10FLAG)
    mainVal = log10(mainVal);
    if (compareVal > 0)
        compareVal = log10(compareVal);
    end
end

% Main print
fprintf(string,mainVal);

% Conditional print
if (compareVal > 0)
    fprintf([' (cf. ' compareFormat ')\n'],compareVal);
else
    fprintf('\n');
end
