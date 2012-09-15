function PR670close
% PR670close - Closes the PR-670 connection.
%
% Syntax:
% PR670close
%
% Description:
% Exits PR-670 remote mode, then closes the serial port connection.

global g_serialPort;

% Close and clear serial port, if it is open.
% Also take the device out of remote mode.
if ~isempty(g_serialPort)
   PR670write('Q', 0);
   IOPort('close', g_serialPort);
   g_serialPort = [];
end
