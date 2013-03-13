function [rawspd,errcode] = PR705rawspd(timeout)
% PR705rawspd - Takes an spd measurement and returns the results.
%
% Syntax:
% rawspd = PR705rawspd(timeout)
%
% Input:
% timeout (scalar) - Timeout period in seconds.
%
% Output:
% rawspd (1xN char) - The raw character array resulting from a measurement.
% This will be an empty string if the timeout period was reached.
% errcode (scalar) - The error code reported by the meter (success = 0).
%
% 12/06/12   zlb   Wrote it.

global g_serialPort

if nargin < 1 || isempty(timeout)
    timeout = 300;
end

IOPort('Purge', g_serialPort);
PR705write('M5');

WaitSecs(0.5);
max_wait_time = GetSecs() + timeout;
while GetSecs() < max_wait_time && ~IOPort('BytesAvailable', g_serialPort)
    WaitSecs(0.5);
end

rawspd = PR705read();
errcode = sscanf(rawspd, '%d');
