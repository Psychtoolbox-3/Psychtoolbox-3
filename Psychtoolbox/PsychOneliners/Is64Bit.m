function rc = Is64Bit
% result = Is64Bit;
%
% Returns 1 (true) if the script is running inside a 64-Bit version of
% GNU/Octave or Matlab, 0 (false) otherwise.
%

% History:
% 3.09.2012  mk  Written.

rc = IsLinux(1) || IsOSX(1) || IsWin(1);

return;
