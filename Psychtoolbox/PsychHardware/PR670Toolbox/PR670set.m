function values = PR670setget(vargin)
% PR670setget - Set and get PR-670 attributes
%
% Syntax:
% configOUT = PR670setget(
% configOUT = PR670config(configIN)
%
% Description:
% Sets various attributes of the PR-670.
%
% Arguments are attribute, value pairs
%
% Supported attributes
%   'GetAperture' -- Get apperture.  Values are numeric in degrees
%   'SetAperture' -- Set apperture.  Values are numeric in degrees
%    Allowable values: 1 , 0.5, 0.25, 0.125

% ordered { Units, Exposure, N Samples, CIE Observer, Sync Mode, Frequency }.



%% Parse the input
p = inputParser;
addOptional(p,'Aperture',[],@isnumeric);
parse(p,varargin{:});

% Display parser results
disp(['Method: ', p.Results.method])
if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end
if ~isempty(p.UsingDefaults)
    disp('Using defaults: ')
    disp(p.UsingDefaults)
end
method = p.Results.method;

% If no inputs, get current state of PR-670 and return string
if nargin == 0
    PR670read;
    PR670write('D601')
    configOUT = PR670read;
    return
end

% PR-670 aperture codes
%      0 -> 1 deg
%      1 -> 1/2 deg
%      2 -> 1/4 deg
%      3 -> 1/8 deg

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
