function resultFlag = IsARM
% resultFlag = IsARM
%
% Returns true if the processor architecture is ARM.

% HISTORY
% 4.4.2013 mk   Wrote it.

persistent rc;
if isempty(rc)
     rc= ~isempty(strfind(computer, 'arm-'));
end

resultFlag = rc;
