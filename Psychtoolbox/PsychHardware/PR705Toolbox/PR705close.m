function PR705close()
% PR705close - Closes the PR-705 connection.
%
% Syntax:
% PR705close
%
% Description:
% Exits PR-705 remote mode, then closes the serial port connection.
%
% 11/29/12    zlb   Wrote it based on the PR670Toolbox.

global g_serialPort

% Close and clear serial port, if it is open.
% Also take the device out of remote mode.
if ~isempty(g_serialPort)
    IOPort('Purge', g_serialPort);
    PR705write('Q');
    IOPort('Close', g_serialPort);
    g_serialPort = [];
end
