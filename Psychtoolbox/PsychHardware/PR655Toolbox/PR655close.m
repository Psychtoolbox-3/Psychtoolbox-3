function PR655close
% PR655close
%
% Exit PR655 Remote Mode, Close serial port, Reset serial global.
% 
% 01/16/09    tbc   Adapted from PR650Toolbox for use with PR655
%

global g_serialPort;

% Exit Remote Mode
PR655write('Q');

% Close and clear serial port
if ~isempty(g_serialPort)
   IOPort('close', g_serialPort);
   g_serialPort = [];
end

return
