function configOUT = PR670config(configIN)
% PR670config - Update the PR-670 configuration.
%
% Syntax:
% configOUT = PR670config
% configOUT = PR670config(configIN)
%
% Description:
% Updates the PR-670 configuration according to 'configIN'. If no input is
% provided, returns the current configuration string. 'configIN' can be in
% the form of the string returned from PR670config, or a cell containing values
% ordered { Units, Exposure, N Samples, CIE Observer, Sync Mode, Frequency }.

% Note: I haven't checked if this works, I copied it from the PR-655
% toolbox.  CGB

% If no inputs, get current state of PR-655 and return string
if nargin == 0
    PR670read;
    PR670write('D601')
    configOUT = PR670read;
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
    
    PR670write(cmdStr)
    tmp = PR670read;
    if str2num(tmp) ~=0
        disp(['Error: could not write ',cmdStr,' to PR-670. Config aborted.'])
        configOUT = tmp;
        return
    end
end

PR670write('D601')
configOUT = PR670read;
