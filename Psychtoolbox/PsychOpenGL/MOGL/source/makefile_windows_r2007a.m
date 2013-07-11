function makefile_windows_r2007a
% Build moglcore on MS-Windows for Matlab 32/64 Bit:
%
% Requires freeglut import libraries and header files which are bundled
% inside the PsychSourceGL\Cohorts\freeglut folder of the PTB source distro.
%
% Requires freeglut dll's included within Psychtoolbox distribution.
%

if ~IsWin
    error('This makefile is only for building moglcore on MS-Windows!');
end

if IsWin(1)
    % 64-Bit build:
    mex -v -outdir . -output moglcore -largeArrayDims -DWINR2007a -DWINDOWS -DGLEW_STATIC -I..\..\..\..\PsychSourceGL\Cohorts\freeglut\include -L..\..\..\..\PsychSourceGL\Cohorts\freeglut\lib\x64 -I. windowhacks.c gl_auto.c gl_manual.c mogl_rebinder.c moglcore.c glew.c ftglesGlue.c user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib -lfreeglut
else
    % 32-Bit build:
    mex -v -outdir . -output moglcore -largeArrayDims -DWINR2007a -DWINDOWS -DGLEW_STATIC -I..\..\..\..\PsychSourceGL\Cohorts\freeglut\include -L..\..\..\..\PsychSourceGL\Cohorts\freeglut\lib -I. windowhacks.c gl_auto.c gl_manual.c mogl_rebinder.c moglcore.c glew.c ftglesGlue.c user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib -lfreeglut
end

movefile(['moglcore.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);

return;
