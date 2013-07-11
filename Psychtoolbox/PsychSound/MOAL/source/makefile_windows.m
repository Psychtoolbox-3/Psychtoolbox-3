function makefile_windows
% Makefile for building the moalcore MEX file for Matlab+OpenAL under
% MS-Windows. Call it while inside the .../MOAL/source folder. You'll
% have to install the Creative labs OpenAL SDK into ...
% C:\Program Files (x86)\OpenAL 1.1 SDK\
% ... which is the default install location on Windows-7.
%
% The SDK is currently available from:
% http://connect.creativelabs.com/openal
%

if ~IsWin
    error('This makefile is only for building moalcore on MS-Windows!');
end

if IsWin(1)
    % 64-Bit R2007a or later build:
    mex -v -outdir . -output moalcore -largeArrayDims -DWINDOWS -I"C:\Program Files (x86)\OpenAL 1.1 SDK\include" -L"C:\Program Files (x86)\OpenAL 1.1 SDK\libs\Win64" moalcore.c al_auto.c al_manual.c alm.c user32.lib -lOpenAL32
else
    % 32-Bit R2007a or later build:
    mex -v -outdir . -output moalcore -largeArrayDims -DWINDOWS -I"C:\Program Files (x86)\OpenAL 1.1 SDK\include" -L"C:\Program Files (x86)\OpenAL 1.1 SDK\libs\Win32" moalcore.c al_auto.c al_manual.c alm.c user32.lib -lOpenAL32
end

movefile(['moalcore.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);

return;
