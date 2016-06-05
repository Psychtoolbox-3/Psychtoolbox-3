function windowsmakeit64_twisty(what, onoctave)
% Builds the 64-Bit Psychtoolbox on MS-Windows for Octave-4 and Matlab.
% As a bonus it also builds the 32-Bit Psychtoolbox for 32-Bit Octave-4.
% This script is customized for MK's build machine "twisty", building
% against a VirtualBox VM running MS-Windows 7 64-Bit.
%

if ~IsWin
    error('%s must be run on MS-Windows!', mfilename);
end

if ~Is64Bit && ~IsOctave
    error('%s must be run on MS-Windows within 32 or 64 Bit Octave or within 64-Bit Matlab!', mfilename);
end

if nargin < 1
    what = 0;
end

if nargin < 2
    onoctave = IsOctave;
end

% Rebuild all request?
if what == -1
    % Yes: Call ourselves recursively on all plugins/modes to rebuild
    % everything:
    tic;
    for what = 0:13
        windowsmakeit64_twisty(what);
    end
    elapsedsecs = toc;
    fprintf('Total rebuild time for all mex files was %f seconds. Bye.\n\n', elapsedsecs);
    return;
end

try

% Matlab or Octave build?
if onoctave == 0
    % Matlab build:
    if what == 0
        % Default: Build Screen with GStreamer-1 support: Needs the
        % www.gstreamer.net GStreamer-1.4 SDK for 64-Bit Windows.
        % Use this for verbose linker output: /VERBOSE:LIB
        mex -v -outdir ..\Projects\Windows\build -output Screen -DPTBMODULE_Screen -largeArrayDims -DPTB_USE_GSTREAMER -DGLEW_STATIC -LC:\gstreamer\1.0\x86_64\lib -IC:\gstreamer\1.0\x86_64\include -IC:\gstreamer\1.0\x86_64\include\gstreamer-1.0 -IC:\gstreamer\1.0\x86_64\include\glib-2.0 -IC:\gstreamer\1.0\x86_64\include\glib-2.0\include -IC:\gstreamer\1.0\x86_64\lib\glib-2.0\include -IC:\gstreamer\1.0\x86_64\lib\gstreamer-1.0\include -IC:\gstreamer\1.0\x86_64\include\libxml2 -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -ICommon\Base -ICommon\Screen -IWindows\Base -IWindows\Screen Windows\Screen\*.c Windows\Base\*.c Common\Base\*.c Common\Screen\*.c kernel32.lib user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib winmm.lib delayimp.lib -lgobject-2.0 -lglib-2.0 -lgstreamer-1.0 -lgstapp-1.0 -lgstvideo-1.0 -lgstpbutils-1.0 LINKFLAGS="$LINKFLAGS /DELAYLOAD:libgobject-2.0-0.dll /DELAYLOAD:libglib-2.0-0.dll /DELAYLOAD:libgstreamer-1.0-0.dll /DELAYLOAD:libgstapp-1.0-0.dll /DELAYLOAD:libgstvideo-1.0-0.dll /DELAYLOAD:libgstpbutils-1.0-0.dll"
        movefile(['..\Projects\Windows\build\Screen.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 1
        % Build WaitSecs
        mex -v -outdir ..\Projects\Windows\build -output WaitSecs -DPTBMODULE_WaitSecs -largeArrayDims  -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -ICommon\Base -ICommon\WaitSecs -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\WaitSecs\*.c kernel32.lib user32.lib winmm.lib
        movefile(['..\Projects\Windows\build\WaitSecs.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 2
        % Build PsychPortAudio
        mex -v -outdir ..\Projects\Windows\build -output PsychPortAudio -DPTBMODULE_PsychPortAudio -largeArrayDims -L..\Cohorts\PortAudio -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -ICommon\Base -ICommon\PsychPortAudio -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\PsychPortAudio\*.c kernel32.lib user32.lib winmm.lib delayimp.lib -lportaudio_x64 LINKFLAGS="$LINKFLAGS /DELAYLOAD:portaudio_x64.dll"
        movefile(['..\Projects\Windows\build\PsychPortAudio.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 3
        % Build GetSecs
        mex -v -outdir ..\Projects\Windows\build -output GetSecs -DPTBMODULE_GetSecs -largeArrayDims  -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -ICommon\Base -ICommon\GetSecs -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\GetSecs\*.c kernel32.lib user32.lib winmm.lib
        movefile(['..\Projects\Windows\build\GetSecs.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 4
        % Build IOPort
        mex -v -outdir ..\Projects\Windows\build -output IOPort -DPTBMODULE_IOPort -largeArrayDims  -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -ICommon\Base -ICommon\IOPort -IWindows\Base -IWindows\IOPort Windows\Base\*.c Common\Base\*.c Common\IOPort\*.c Windows\IOPort\*.c kernel32.lib user32.lib winmm.lib
        movefile(['..\Projects\Windows\build\IOPort.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 5
        % Build PsychCV
        % Disabled for now: As long as it contains 3rd party code, we can't
        % really distribute it in a precompiled version for Matlab under our
        % new license. Distribution of compiled mex files for octave would
        % be possible, but see below...
        % mex -v -outdir ..\Projects\Windows\build -output PsychCV -DPTBMODULE_PsychCV -largeArrayDims  -ID:\install\QuickTimeSDK\CIncludes -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -ICommon\Base -ICommon\PsychCV -IWindows\Base -I..\Cohorts\ARToolkit\include Windows\Base\*.c Common\Base\*.c Common\PsychCV\*.c kernel32.lib user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib winmm.lib delayimp.lib libARvideo.lib libARgsub.lib libARgsub_lite.lib libARgsubUtil.lib libARMulti.lib libAR.lib
        % movefile(['..\Projects\Windows\build\PsychCV.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 7
        % Build Eyelink:
        mex -v -outdir ..\Projects\Windows\build -output Eyelink -DPTBMODULE_Eyelink -largeArrayDims -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -I"C:\Program Files (x86)\SR Research\EyeLink\Includes\eyelink" -ICommon\Base -ICommon\Eyelink -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\Eyelink\*.c user32.lib gdi32.lib advapi32.lib winmm.lib "C:\Program Files (x86)\SR Research\EyeLink\libs\x64\eyelink_core64.lib" "C:\Program Files (x86)\SR Research\EyeLink\libs\x64\eyelink_w32_comp64.lib"
        movefile(['..\Projects\Windows\build\Eyelink.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 8
        % Build PsychKinectCore:
        mex -v -outdir ..\Projects\Windows\build -output PsychKinectCore -DPTBMODULE_PsychKinectCore -largeArrayDims -DWIN32 -L..\Cohorts\Kinect-v16-withsource\libusb\lib\msvc_x64 -I..\Cohorts\Kinect-v16-withsource -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -ICommon\Base -IWindows\Base -ICommon\PsychKinect Windows\Base\*.c Common\Base\*.c Common\PsychKinect\*.c ..\Cohorts\Kinect-v16-withsource\*.cpp kernel32.lib user32.lib winmm.lib -lusb
        movefile(['..\Projects\Windows\build\PsychKinectCore.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 9
        % Build PsychHID:
        mex -v -outdir ..\Projects\Windows\build -output PsychHID -DPTBMODULE_PsychHID -largeArrayDims -DWIN32 -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -L..\Cohorts\libusb1-win32\MS64\dll -I..\Cohorts\libusb1-win32\include\libusb-1.0 -ICommon\Base -IWindows\Base -ICommon\PsychHID Windows\PsychHID\*.cpp Windows\PsychHID\*.c Windows\Base\*.c Common\Base\*.c Common\PsychHID\*.c -ldinput8 kernel32.lib user32.lib winmm.lib -lusb-1.0 setupapi.lib
        movefile(['..\Projects\Windows\build\PsychHID.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end

    if what == 12
        % Build PsychOculusVRCore:
        % Needs the Oculus VR SDK v 0.5.01 installed in side-by-side to
        % the Psychtoolbox-3 folder, ie., in the same parent directory, and
        % renamed from OculusSDK to OculusSDKWin:
        mex -v -outdir ..\Projects\Windows\build -output PsychOculusVRCore -DPTBMODULE_PsychOculusVRCore -largeArrayDims -DWIN32 -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -I..\..\..\OculusSDKWin\LibOVR\Include -ICommon\Base -IWindows\Base -ICommon\PsychOculusVRCore Common\PsychOculusVRCore\*.c Windows\Base\*.c Common\Base\*.c kernel32.lib user32.lib winmm.lib ..\..\..\OculusSDKWin\LibOVR\Lib\Windows\x64\Release\VS2010\LibOVR.lib
        movefile(['..\Projects\Windows\build\PsychOculusVRCore.' mexext], [PsychtoolboxRoot 'PsychBasic\MatlabWindowsFilesR2007a\']);
    end
else
    % Octave-4 build:
    if Is64Bit
        target = [PsychtoolboxRoot 'PsychBasic\Octave4WindowsFiles64\'];
    else
        target = [PsychtoolboxRoot 'PsychBasic\Octave4WindowsFiles\'];
    end

    if what == 0
        % Default: Build Screen.mex
        if Is64Bit
            mexoctave -g -v --output ..\Projects\Windows\build\Screen.mex -DPTBMODULE_Screen -DPTB_USE_GSTREAMER -DPTBOCTAVE3MEX -DGLEW_STATIC -LC:\gstreamer\1.0\x86_64\lib -IC:\gstreamer\1.0\x86_64\include -IC:\gstreamer\1.0\x86_64\include\gstreamer-1.0 -IC:\gstreamer\1.0\x86_64\include\glib-2.0 -IC:\gstreamer\1.0\x86_64\include\glib-2.0\include -IC:\gstreamer\1.0\x86_64\lib\glib-2.0\include -IC:\gstreamer\1.0\x86_64\lib\gstreamer-1.0\include -IC:\gstreamer\1.0\x86_64\include\libxml2 -ICommon\Base -ICommon\Screen -IWindows\Base -IWindows\Screen Windows\Screen\*.c Windows\Base\*.c Common\Base\*.c Common\Screen\*.c kernel32.lib user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib winmm.lib delayimp.lib -lgobject-2.0 -lglib-2.0 -lgstreamer-1.0 -lgstapp-1.0 -lgstvideo-1.0 -lgstpbutils-1.0
        else
            mexoctave -g -v --output ..\Projects\Windows\build\Screen.mex -DPTBMODULE_Screen -DPTB_USE_GSTREAMER -DPTBOCTAVE3MEX -DGLEW_STATIC -LC:\gstreamer\1.0\x86\lib -IC:\gstreamer\1.0\x86\include -IC:\gstreamer\1.0\x86\include\gstreamer-1.0 -IC:\gstreamer\1.0\x86\include\glib-2.0 -IC:\gstreamer\1.0\x86\include\glib-2.0\include -IC:\gstreamer\1.0\x86\lib\glib-2.0\include -IC:\gstreamer\1.0\x86\lib\gstreamer-1.0\include -IC:\gstreamer\1.0\x86\include\libxml2 -ICommon\Base -ICommon\Screen -IWindows\Base -IWindows\Screen Windows\Screen\*.c Windows\Base\*.c Common\Base\*.c Common\Screen\*.c kernel32.lib user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib winmm.lib delayimp.lib -lgobject-2.0 -lglib-2.0 -lgstreamer-1.0 -lgstapp-1.0 -lgstvideo-1.0 -lgstpbutils-1.0
        end
        movefile(['..\Projects\Windows\build\Screen.' mexext], target);
    end

    if what == 1
        % Build WaitSecs.mex
        mexoctave -g -v --output ../Projects/Windows/build/WaitSecs.mex -DPTBMODULE_WaitSecs -DPTBOCTAVE3MEX -ICommon\Base -ICommon\WaitSecs -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\WaitSecs\*.c kernel32.lib user32.lib winmm.lib
        movefile(['..\Projects\Windows\build\WaitSecs.' mexext], target);
    end

    if what == 2
        % Build PsychPortAudio.mex
        % This needs a libportaudio_x86 compatible with mingw32 gcc. The way we do this
        % is by taking the libportaudio_x86.dll from the PsychSound folder and copying
        % it into the PsychSourceGL/Cohorts/PortAudio/MinGW32 folder, but with the file extension
        % changed from .dll to .lib, ie. libportaudio_x86.lib - Weird but true, this makes
        % mingws linker accept the file as input. Same procedure for libportaudio_x64.dll to
        % libportaudio_x64.lib into targetfolder MinGW64.
        if Is64Bit
            mexoctave -g -v --output ..\Projects\Windows\build\PsychPortAudio.mex -DPTBMODULE_PsychPortAudio -DPTBOCTAVE3MEX -L..\Cohorts\PortAudio\MinGW64 -ICommon\Base -ICommon\PsychPortAudio -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\PsychPortAudio\*.c kernel32.lib user32.lib winmm.lib delayimp.lib -lportaudio_x64
        else
            mexoctave -g -v --output ..\Projects\Windows\build\PsychPortAudio.mex -DPTBMODULE_PsychPortAudio -DPTBOCTAVE3MEX -L..\Cohorts\PortAudio\MinGW32 -ICommon\Base -ICommon\PsychPortAudio -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\PsychPortAudio\*.c kernel32.lib user32.lib winmm.lib delayimp.lib -lportaudio_x86
        end
        movefile(['..\Projects\Windows\build\PsychPortAudio.' mexext], target);
    end

    if what == 3
        % Build GetSecs.mex
        mexoctave -g -v --output ..\Projects\Windows\build\GetSecs.mex -DPTBMODULE_GetSecs -DPTBOCTAVE3MEX -ICommon\Base -ICommon\GetSecs -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\GetSecs\*.c kernel32.lib user32.lib winmm.lib
        movefile(['..\Projects\Windows\build\GetSecs.' mexext], target);
    end

    if what == 4
        % Build IOPort.mex
        mexoctave -g -v --output ..\Projects\Windows\build\IOPort.mex -DPTBMODULE_IOPort -DPTBOCTAVE3MEX -ICommon\Base -ICommon\IOPort -IWindows\Base -IWindows\IOPort Windows\Base\*.c Common\Base\*.c Common\IOPort\*.c Windows\IOPort\*.c kernel32.lib user32.lib winmm.lib
        movefile(['..\Projects\Windows\build\IOPort.' mexext], target);
    end

    if what == 5
        % Build PsychCV.mex
        % NOTE: Link is currently broken. Also we don't build and
        % distribute PsychCV.mex at the moment. Let's see if anybody
        % actually misses this mex file...
        % mexoctave -g -v --output ..\Projects\Windows\build\PsychCV.mex -DPTBMODULE_PsychCV -DPTBOCTAVE3MEX -ID:\install\QuickTimeSDK\CIncludes -ID:\MicrosoftDirectXSDK\Include -ICommon\Base -ICommon\PsychCV -IWindows\Base -I..\Cohorts\ARToolkit\include Windows\Base\*.c Common\Base\*.c Common\PsychCV\*.c kernel32.lib user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib winmm.lib delayimp.lib libARvideo.lib libARgsub.lib libARgsub_lite.lib libARgsubUtil.lib libARMulti.lib libAR.lib
        %movefile(['..\Projects\Windows\build\PsychCV.' mexext], target);
    end

    if what == 6
        % Build moglcore.mex: Please note that 32-Bit moglcore.mex for Octave
        % actually links against libfreeglut.a, *not* libfreeglut.lib in the
        % lib folder!
        curdir = pwd;
        cd('../../Psychtoolbox/PsychOpenGL/MOGL/source/')
        try
            if Is64Bit
                mexoctave -g -v --output moglcore.mex -DWINDOWS -DGLEW_STATIC -I..\..\..\..\PsychSourceGL\Cohorts\freeglut\include -L..\..\..\..\PsychSourceGL\Cohorts\freeglut\lib\x64 -I. windowhacks.c gl_manual.c mogl_rebinder.c moglcore.c glew.c ftglesGlue.c gl_auto.c user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib -lfreeglut
            else
                mexoctave -g -v --output moglcore.mex -DWINDOWS -DGLEW_STATIC -I..\..\..\..\PsychSourceGL\Cohorts\freeglut\include -L..\..\..\..\PsychSourceGL\Cohorts\freeglut\lib -I. windowhacks.c gl_manual.c mogl_rebinder.c moglcore.c glew.c ftglesGlue.c gl_auto.c user32.lib gdi32.lib advapi32.lib glu32.lib opengl32.lib -lfreeglut
            end
            movefile(['moglcore.' mexext], target);
        catch
            lasterr
        end
        % Remove stale object files:
        delete('*.o');
        cd(curdir);
    end

    if what == 7
        % Build Eyelink.mex
        try
            if Is64Bit
                copyfile('C:\Program Files (x86)\SR Research\EyeLink\libs\x64\*.dll', '..\..\..\');
                mexoctave -g -v -W'l,--enable-stdcall-fixup' --output ..\Projects\Windows\build\Eyelink.mex -DPTBMODULE_Eyelink -DPTBOCTAVE3MEX -I'C:\Program Files (x86)\SR Research\EyeLink\Includes\eyelink' -L..\..\.. -ICommon\Base -ICommon\Eyelink -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\Eyelink\*.c user32.lib gdi32.lib advapi32.lib winmm.lib -leyelink_core64 -leyelink_w32_comp64
            else
                copyfile('C:\Program Files (x86)\SR Research\EyeLink\libs\*.dll', '..\..\..\');
                mexoctave -g -v -W'l,--enable-stdcall-fixup' --output ..\Projects\Windows\build\Eyelink.mex -DPTBMODULE_Eyelink -DPTBOCTAVE3MEX -I'C:\Program Files (x86)\SR Research\EyeLink\Includes\eyelink' -L..\..\.. -ICommon\Base -ICommon\Eyelink -IWindows\Base Windows\Base\*.c Common\Base\*.c Common\Eyelink\*.c user32.lib gdi32.lib advapi32.lib winmm.lib -leyelink_core -leyelink_w32_comp
            end
            movefile(['..\Projects\Windows\build\Eyelink.' mexext], target);
        catch
        end
        delete('..\..\..\*.dll');
    end

    if what == 8
        % Build PsychKinectCore.mex:
        if Is64Bit
            mexoctave -g -v --output ..\Projects\Windows\build\PsychKinectCore.mex -DPTBMODULE_PsychKinectCore -DPTBOCTAVE3MEX -L..\Cohorts\Kinect-v16-withsource\libusb\lib\msvc_x64 -I..\Cohorts\Kinect-v16-withsource -ICommon\Base -IWindows\Base -ICommon\PsychKinect Windows\Base\*.c Common\Base\*.c Common\PsychKinect\*.c ..\Cohorts\Kinect-v16-withsource\*.cpp kernel32.lib user32.lib winmm.lib libusb.lib
        else
            mexoctave -g -v --output ..\Projects\Windows\build\PsychKinectCore.mex -DPTBMODULE_PsychKinectCore -DPTBOCTAVE3MEX -L..\Cohorts\Kinect-v16-withsource\libusb\lib\gcc -I..\Cohorts\Kinect-v16-withsource -ICommon\Base -IWindows\Base -ICommon\PsychKinect Windows\Base\*.c Common\Base\*.c Common\PsychKinect\*.c ..\Cohorts\Kinect-v16-withsource\*.cpp kernel32.lib user32.lib winmm.lib -lusb
        end
        movefile(['..\Projects\Windows\build\PsychKinectCore.' mexext], target);
    end;

    if what == 9
        % Build PsychHID.mex:
        if Is64Bit
            mexoctave -g -v --output ..\Projects\Windows\build\PsychHID.mex -DPTBMODULE_PsychHID -DPTBOCTAVE3MEX -L..\Cohorts\libusb1-win32\MinGW64\dll -I..\Cohorts\libusb1-win32\include\libusb-1.0 -ICommon\Base -IWindows\Base -ICommon\PsychHID Windows\PsychHID\*.cpp Windows\PsychHID\*.c Windows\Base\*.c Common\Base\*.c Common\PsychHID\*.c dinput8.lib kernel32.lib user32.lib winmm.lib libusb-1.0.lib setupapi.lib
        else
            mexoctave -g -v --output ..\Projects\Windows\build\PsychHID.mex -DPTBMODULE_PsychHID -DPTBOCTAVE3MEX -L..\Cohorts\libusb1-win32\MinGW32\dll -I..\Cohorts\libusb1-win32\include\libusb-1.0 -ICommon\Base -IWindows\Base -ICommon\PsychHID Windows\PsychHID\*.cpp Windows\PsychHID\*.c Windows\Base\*.c Common\Base\*.c Common\PsychHID\*.c dinput8.lib kernel32.lib user32.lib winmm.lib libusb-1.0.lib setupapi.lib
        end
        movefile(['..\Projects\Windows\build\PsychHID.' mexext], target);
    end

    if what == 10
        % Build moalcore 32-Bit:
        % For 32-Bit we copy the installed runtime DLL
        % C:\Windows\SysWOW64\OpenAL32dll into the build folder and abuse that
        % as import library for our mex file, as the import .lib of the SDK is
        % unuseable by the linker.
        % However, for 64-Bit strangely the .lib from the SDK works.
        curdir = pwd;
        cd('../../Psychtoolbox/PsychSound/MOAL/source/')
        if Is64Bit
            copyfile('C:\Program Files (x86)\OpenAL 1.1 SDK\libs\Win64\OpenAL32.lib', '.\');
        else
            copyfile('C:\Windows\SysWOW64\OpenAL32.dll', '.\');
        end

        try
            mexoctave -g -v --output moalcore.mex -DWINDOWS -I'C:\Program Files (x86)\OpenAL 1.1 SDK\include' -L. moalcore.c al_auto.c al_manual.c alm.c user32.lib -lOpenAL32
            movefile(['moalcore.' mexext], target);
        catch
            lasterr
        end
        % Remove stale object files:
        delete('*.o');
        delete('OpenAL32.lib');
        delete('OpenAL32.dll');
        cd(curdir);
    end

    if what == 11
        % Build pnet
        curdir = pwd;
        cd('../../Psychtoolbox/PsychHardware/iViewXToolbox/tcp_udp_ip/')
        try
            mexoctave -g -v --output pnet.mex -DWIN32 pnet.c ws2_32.lib winmm.lib
            movefile(['pnet.' mexext], target);
        catch
        end
        % Remove stale object files:
        delete('*.o');
        cd(curdir);
    end

    if what == 12
        % Build PsychOculusVRCore 32/64 Bit:
        % Needs the Oculus VR SDK v 0.5.01 installed in side-by-side to
        % the Psychtoolbox-3 folder, ie., in the same parent directory, and
        % renamed from OculusSDK to OculusSDKWin. Needs to Oculus VR runtime
        % for Windows installed as well. Then we copy the installed runtime DLL
        % C:\windows\system32\LibOVRRT64_0_5.dll or - for 32-Bit -
        % C:\Windows\SysWOW64\LibOVRRT32_0_5.dll into the SDK folder and abuse that
        % as import library for our mex file. Weird but true - this seems to work!
        if Is64Bit
            copyfile('C:\windows\system32\LibOVRRT64_0_5.dll', '..\..\..\OculusSDKWin\LibOVR\Lib\Windows\Win32\Release\VS2010\');
            mexoctave -g -v --output ..\Projects\Windows\build\PsychOculusVRCore.mex -DPTBMODULE_PsychOculusVRCore -DPTBOCTAVE3MEX -I..\..\..\OculusSDKWin\LibOVR\Include -ICommon\Base -IWindows\Base -ICommon\PsychOculusVRCore Common\PsychOculusVRCore\*.c Windows\Base\*.c Common\Base\*.c kernel32.lib user32.lib winmm.lib -L..\..\..\OculusSDKWin\LibOVR\Lib\Windows\Win32\Release\VS2010 -lOVRRT64_0_5
        else
            copyfile('C:\Windows\SysWOW64\LibOVRRT32_0_5.dll', '..\..\..\OculusSDKWin\LibOVR\Lib\Windows\Win32\Release\VS2010\');
            mexoctave -g -v --output ..\Projects\Windows\build\PsychOculusVRCore.mex -DPTBMODULE_PsychOculusVRCore -DPTBOCTAVE3MEX -I..\..\..\OculusSDKWin\LibOVR\Include -ICommon\Base -IWindows\Base -ICommon\PsychOculusVRCore Common\PsychOculusVRCore\*.c Windows\Base\*.c Common\Base\*.c kernel32.lib user32.lib winmm.lib -L..\..\..\OculusSDKWin\LibOVR\Lib\Windows\Win32\Release\VS2010 -lOVRRT32_0_5
        end
        movefile(['..\Projects\Windows\build\PsychOculusVRCore.' mexext], target);
    end

    % Remove stale object files:
    delete('*.o');
end

catch %#ok<*CTCH>
    % Empty. We just want to make sure the delete() call below is executed
    % in both success and failure case.
end

return;

% Special mex wrapper for Octave compile on Windows:
function mexoctave(varargin)
debugme = 0;
callmex = 1;

if (debugme), fprintf('nargin = %i\n', nargin); end;
myvararg = cell();
myvararg(end+1) = '-DPTBOCTAVE3MEX';
myvararg(end+1) = '-Wno-multichar';
myvararg(end+1) = '-Wno-unknown-pragmas';

outarg = '';
quoted = 0;
emitarg = 0;
for i=1:nargin
    curarg = char(varargin(i));

    if debugme, fprintf('Preparse Arg %i: %s\n', i, curarg); end
    if (length(strfind(curarg, '"')) > 0)
        if ~quoted
            % Start of quoted string:
            quoted = 1;
            outarg = [curarg];
            emitarg = 0;
        else
            % End of quoted string: Emit!
            quoted = 0;
            outarg = [ outarg curarg ];
            emitarg = 1;
        end
    else
        % Not start or end string of a quoted piece:
        if quoted
            % Within a quoted segment!
            outarg = [ outarg curarg ];
            emitarg = 0;
        else
            % Outside a quoted segment:
            % Expansion needed?
            ppos = strfind(curarg, '*.c');
            if length(ppos) > 0
                prefix  = curarg(1:ppos(1)-1);
                allfiles = dir(curarg);
                for j=1:length(allfiles)
                    expandedfiles = [ prefix allfiles(j).name ];
                    emitarg = 0;
                    myvararg(end+1) = expandedfiles;
                    if (debugme), fprintf('Emitted Arg %i : %s\n', length(myvararg), char(myvararg(end))); end;
                end
            else
                % Regular chunk: Emit it
                ppos = strfind(curarg, '.lib');
                if ~isempty(ppos)
                    curarg = [ '-l' curarg(1:ppos(1)-1) ];
                end

                outarg = curarg;
                emitarg = 1;
            end
        end
    end

    if emitarg
        emitarg = 0;
        myvararg(end+1) = outarg;
        if (debugme), fprintf('Emitted Arg %i : %s\n', length(myvararg), char(myvararg(end))); end;
    end
end

final = myvararg

if (debugme), outargtype = class(myvararg), end;

if (callmex)
    mex(myvararg{:});
end

return;
