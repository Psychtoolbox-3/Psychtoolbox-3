function PR670write(cmdStr, appendCR)
% PR670write - Write a string of characters to the PR-670.
%
% Syntax:
% PR670write(cmdStr)
% PR670write(cmdStr, appendCR)
%
% Description:
% Writes a string of characters to the PR-670.  By default, it appends a
% carriage return (CR) to the end of the string.  Note that if the CR is
% already there, it won't be appended. The CR can be disabled as
% some commands do not need it.  See the PR-670 documentation for which
% command need the CR.
%
% Input:
% cmdStr (1xN char) - String of characters to send.
% appendCR (scalar) - 1 to append a CR, 0 to not append.  Default: 1

global g_serialPort;

if nargin == 1
	appendCR = 1;
end

if appendCR
	% Check for the CR and add if necessary.
	if cmdStr(end) ~= char(13)
		cmdStr = [cmdStr char(13)];
	end
end

% Write sequence of chars to PR-670.
for i = 1:length(cmdStr)
    IOPort('write', g_serialPort, upper(cmdStr(i)));
    pause(0.05)
end
