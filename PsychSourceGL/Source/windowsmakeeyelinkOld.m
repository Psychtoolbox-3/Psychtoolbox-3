function windowsmakeeyelink(postR2007a)
% windowsmakeeyelink([postR2007a=1])
% windowsmakeeyelink -- Simple build-script for building
% Eyelink toolbox for M$-Windows:
% You will need to adapt the paths for your build system manually.
%
% postR2007a == 0 --> Build on Matlab R11. Unsupported.
% postR2007a == 1 --> Build on Matlab R2007a.

if nargin < 1
    postR2007a = [];
end

if isempty(postR2007a)
    postR2007a = 1;
end

if postR2007a
    % Build sequence for >= R2007a Matlabs on MPI HDR machine - build against
    % R2007a:

    % Copy our C++ PsychScriptingGlue.cc to a C style PsychScriptingGlue.c so the f%*§$!d M$-Compiler
    % can handle it:
    dos('copy T:\projects\OpenGLPsychtoolbox\trunk\PsychSourceGL\Source\Common\Base\PsychScriptingGlue.cc T:\projects\OpenGLPsychtoolbox\trunk\PsychSourceGL\Source\Common\Base\PsychScriptingGlue.c');

    % Build sequence for Eyelink.dll: Requires the freely downloadable (after registration) Eyelink-SDK for Windows.
    %mex -v -outdir T:\projects\OpenGLPsychtoolbox\trunk\PsychSourceGL\Projects\Windows\build\ -output Eyelink.dll  -ID:\install\QuickTimeSDK\CIncludes -I"C:\Program Files\SR Research\EyeLink\Includes\eyelink" -I"C:\Programme\Microsoft Visual Studio 8\VC\Include" -ICommon\Base -ICommon\Eyelink -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\Eyelink\*.c user32.lib gdi32.lib advapi32.lib winmm.lib "C:\Program Files\SR Research\EyeLink\libs\eyelink_core.lib" "C:\Program Files\SR Research\EyeLink\libs\eyelink_w32_comp.lib" "C:\Program Files\SR Research\EyeLink\libs\eyelink_exptkit20.lib"
    mex -v -outdir T:\projects\OpenGLPsychtoolbox\trunk\PsychSourceGL\Projects\Windows\build\ -output Eyelink -I"C:\Program Files\SR Research\EyeLink\Includes\eyelink" -I"C:\Programme\Microsoft Visual Studio 8\VC\Include" -ICommon\Base -ICommon\Eyelink -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\Eyelink\*.c user32.lib gdi32.lib advapi32.lib winmm.lib "C:\Program Files\SR Research\EyeLink\libs\eyelink_core.lib" "C:\Program Files\SR Research\EyeLink\libs\eyelink_w32_comp.lib" "C:\Program Files\SR Research\EyeLink\libs\eyelink_exptkit20.lib"

    % Move Eyelink.dll into its proper location:
    dos('copy T:\projects\OpenGLPsychtoolbox\trunk\PsychSourceGL\Projects\Windows\build\Eyelink.mexw32 T:\projects\OpenGLPsychtoolbox\trunk\Psychtoolbox\PsychBasic\MatlabWindowsFilesR2007a\');

    % Delete the temporary .c version of Scripting Glue:
    delete('T:\projects\OpenGLPsychtoolbox\trunk\PsychSourceGL\Source\Common\Base\PsychScriptingGlue.c');

else
    % Build sequence for pre R2007a Matlabs - Dead:
    error ('Building on pre R2007a is no longer supported.');
end

return;
