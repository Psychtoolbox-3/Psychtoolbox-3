function varargout = PR705write(datastr, varargin)
% PR705write - Write a string of characters to the PR-705.
%
% Syntax:
% varargout = PR705write(cmdStr [, block=1])
%
% Description:
% This is a convenience wrapper to write data to the PR-705 device. Aside
% from appending a carriage return and not requiring an explicit port
% handle, this function is identical to IOPort's Write (see IOPort Write? 
% for more details).
%
% Inputs:
% datastr (1xN char) - data string to write to the PR-705.
% block (scalar) - enable (1, default) or disable (0) blocking writes
%
% Outputs:
% Refer to IOPort Write?
%
% 11/29/12    zlb   Wrote it.

global g_serialPort

if nargin < 1 || isempty(datastr) || nargin > 2 || nargout > 6
    error('Invalid input/output arguments. Please refer to ''help PR750write''.\n');
end

if datastr(end) ~= 13 && ~strcmp(datastr, 'PR705') % only command that doesn't require a CR
    datastr = [datastr 13];
end

% Write sequence of chars to PR-705.
[varargout{1:nargout}] = IOPort('Write', g_serialPort, upper(datastr), varargin{:});
