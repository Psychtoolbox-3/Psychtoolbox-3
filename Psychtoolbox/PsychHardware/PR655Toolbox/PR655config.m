function configOUT = PR655config(configIN)
% Update PR-655 configuration according to configIN. If no input is provided, return the current
% configuration string.
%   [configIN] can be in form of the string returned from PR655config, or a cell containing values
%   ordered { Units, Exposure, N Samples, CIE Observer, Sync Mode, Frequency }. 
% 
% 02/09/09    tbc   Wrote it.
%

global g_serialPort;

% If no inputs, get current state of PR-655 and return string
if nargin == 0
    flush = PR655read;
    PR655write('D601')
    configOUT = PR655read;
    return
end

% otherwise, parse string input and configure PR
if ~iscell(configIN)
    configIN = textscan(configIN,'%s','delimiter',',');
    configIN = configIN{1}([end-7,end-4:end],:); 
end

% Setting suffix for PR setting commands
%  U = units, E = exposure, N = n samples, O = CIE observer, S = sync mode, K = freq.
sufx = 'UENOSK';

for j = 1:length(configIN)
    
    if ischar(configIN{j})
        cmdStr = ['S',sufx(j),configIN{j}];
    else
        cmdStr = ['S',sufx(j),num2str(configIN{j})];
    end
    
    PR655write(cmdStr)
    tmp = PR655read;
    if str2num(tmp) ~=0
        disp(['Error: could not write ',cmdStr,' to PR-655. Config aborted.'])
        configOUT = tmp;
        return
    end
end

PR655write('D601')
configOUT = PR655read;

return
