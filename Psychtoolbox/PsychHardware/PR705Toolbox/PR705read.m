function varargout = PR705read(varargin)
% PR705read - Read data from the PR-705.
%
% Syntax:
% serialData = PR705read([block=0] [, nbytes])
%
% Description:
% This is a convenience wrapper to read data from the PR-705 device. Aside
% from not requiring an explicit port handle, this function is identical to
% IOPort's Read (see IOPort Read? for more details).
%
% Output:
% serialData (1xN char).
%
% 11/29/12    zlb   Wrote it based on the PR670Toolbox.

global g_serialPort

if isempty(g_serialPort)
    error('Meter has not been initialized.');
end

if nargin > 2 || nargout > 3
    error('Invalid input/output arguments. Please refer to ''help PR750read''.\n');
end

[varargout{1:nargout}] = IOPort('Read', g_serialPort, varargin{:});

if ~isempty(varargout{1}) && (~nargin || varargin{1} == 0)
    % If data exists keep reading off the port until there's nothing left.
    while true
        WaitSecs(0.05);
        [tmpData varargout{2:nargout}] = IOPort('Read', g_serialPort);
        if isempty(tmpData), break; end
        varargout{1} = [varargout{1} tmpData];
    end
end

varargout{1} = char(varargout{1});
