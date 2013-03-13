function out = PR705config(varargin)
% PR705config - Update the PR-705 configuration.
%
% Syntax:
% config_str (string) = PR705config
% errcode (scalar) = PR705config('Option1', value1, 'Option2', value2, ...)
%
% Description:
% A wrapper for setting various options on the PR-705 using the 'S'
% command. With no input arguments, the function merely returns the device's
% current configuration string. In the multiple input arguments case, here
% are the case-insensitive options and the values they take:
% Lens         ID Number of PRIMARY accessory
% Add1         ID Number of 1st ADD ON accessory
% Add2         ID Number of 2nd ADD ON accessory
% Aperture     ID Number of Aperture
% Units        0 for English
%              1 for Metric
% ExposureTime 0 - Adaptive
%              25 ... 60000 ms
% CaptureMode  0 - Single Capture
%              1 - Continuous Capture
% Cycles       1 .. 99 - Number of Captures to average
% CalcMode     0 - Power
%              1 - Energy
% TriggerMode  0 - Manual
%              1 - External
% ViewShutter  0 - Open During Measurement
%              1 - Closed During Measurement
% CIEObserver  0 -  2 Degree
%              1 - 10 Degree
%
% Example:
% errcode = PR705config('Lens', 0, 'Units', 1, 'Cycles', 2, 'ViewShutter', 1);
%
% 12/06/12    zlb   Wrote it.

global g_serialPort

if ~nargin
    IOPort('Purge', g_serialPort);
    PR705write('D601');
    out = PR705read(1, 100);
    return
end

values = validate_sort_input(varargin);

setup_str = 'S%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d';
setup_str = sprintf(setup_str, values{:});

if length(setup_str) > 12 % i.e., 'values' had at least 1 non-empty element
    IOPort('Purge', g_serialPort);
    PR705write(setup_str);
    errcode = PR705read(1, 10);
    out = sscanf(errcode, '%d');
else % a noop is trivially successful
    out = 0;
end

function sorted_values = validate_sort_input(raw_input)
valid_sorted_options = {'Lens' 'Add1' 'Add2' 'Aperture' 'Units' 'ExposureTime' ...
    'CaptureMode' 'Cycles' 'CalcMode' 'TriggerMode' 'ViewShutter' 'CIEObserver'};

sorted_values = cell(size(valid_sorted_options));
if mod(numel(raw_input), 2)
    fprintf('Expecting an even number of arguments (option-value pairs)!\n');
    return
end

options = raw_input(1:2:end);
values = raw_input(2:2:end);
if ~all(cellfun(@ischar, options)) || ~all(cellfun(@isscalar, values))
    fprintf('Improper formatting of the option-value pairs!\n');
    return
end

[valid_members,location] = ismember(lower(options), lower(valid_sorted_options));
if ~all(valid_members)
    invalid_options = options(~valid_members);
    fprintf('The following options aren''t recognized and will be ignored:\n');
    fprintf('\t%s\n', invalid_options{:});
end

values = values(valid_members);
location = location(valid_members);

sorted_values(location) = values(:);
