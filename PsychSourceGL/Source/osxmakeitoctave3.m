function osxmakeitoctave3(mode)
% This is the MacOSX version of makeit: It is meant for building PTB for
% 64-Bit Octave-4 on OSX 10.10 Yosemite or later, using the 10.11 SDK.
%
% NOTE: This currently needs to be run on Octave 3.8, but builds for Octave-4.
% The reason is that in Octave-4 mkoctfile was converted from an easy to edit
% shell script into a C++ app, which needs tedious editing in C++ and a super
% expensive rebuilt of GNU/Octave from source code. Yay! We need to edit the
% mkoctfile because it does not handle compilation / linking of OSX -framework
% and OSX Objective-C files, but we obviously need these on OSX and for compiling
% our Objective-C based Cocoa glue in Screen.
%
% This build routine also encodes octave runtime library paths for octave
% version 4 into all mex files. The point of this is to make sure mex files
% don't fail to load due to libraries like liboctinterp.3.dylib not being
% found. This is a stop-gap measure: It would be much better to modify
% octave upstream to encode a -rpath = $OCTAVE_LIBS in the octave
% executable and libraries themselves, and to add a
% @rpath/liboctinterp.3.dylib etc. to the mkoctfile build function.
%
% This rpath encoding is done by the helper function osxsetoctaverpath().

if ~IsOSX(1) || ~IsOctave
    error('osxmakeitoctave3 only works with a 64-Bit version of Octave for OSX!');
end

if nargin < 1 || isempty(mode)
    mode = 0;
end

% Rebuild all request?
if mode == -1
    % Yes: Call ourselves recursively on all plugins/modes to rebuild
    % everything:
    tic;
    for mode = 0:14
        osxmakeitoctave3(mode);
    end
    elapsedsecs = toc;
    fprintf('Total rebuild time for all mex files was %f seconds. Bye.\n\n', elapsedsecs);
    return;
end

% Make our helper scripts available even if we change working directory:
addpath(pwd);

fprintf('Building plugin type %i ...\n\n', mode);

if mode==0
    % Build Screen:
    % Depends: GStreamer-1.4+, libdc1394-2, libusb-1.0

    % Need to build the PsychCocoaGlue.c separately as Objective-C code,
    % whereas the mex call below only compiles C-Code and links it with the
    % PsychCocoaGlue.o object code created here. This to avoid linker
    % trouble with unknown/illegal Objective-C language command line switches:
    mex -g  -v -c -DPTBMODULE_Screen -DPTBVIDEOCAPTURE_LIBDC -DPTBOCTAVE3MEX ...
        -I/usr/include -ICommon/Base -ICommon/Screen -IOSX/Base -IOSX/Screen -IOSX/Fonts -IOSX/EthernetAddress ...
        "-W, -x objective-c -mmacosx-version-min='10.10'" "OSX/Screen/PsychCocoaGlue.c"

    movefile('./OSX/Screen/PsychCocoaGlue.c', './OSX/Screen/PsychCocoaGlue.m');

    % Build with weak linking of GStreamer and libdc1394 via -weak_library
    % flag. This means that a missing GStreamer or libdc1394 installation
    % will not cause linker failure at Screen() load time -- Screen
    % continues to work as long as no GStreamer-dependent functions are
    % used. In that case it will abort with a helpful error message:
    % GStreamer 1.4.x+:
    mex -g -v --output ../Projects/MacOSX/build/Screen -DPTBMODULE_Screen -DPTB_USE_GSTREAMER -DPTBVIDEOCAPTURE_LIBDC -DPTBOCTAVE3MEX ...
    "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ \
    -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgstreamer-1.0.dylib \
    -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgstbase-1.0.dylib \
    -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgstapp-1.0.dylib \
    -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgstpbutils-1.0.dylib \
    -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgstvideo-1.0.dylib \
    -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgobject-2.0.dylib \
    -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgmodule-2.0.dylib \
    -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libxml2.dylib \
    -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgthread-2.0.dylib \
    -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libglib-2.0.dylib \
    -weak_library /usr/local/lib/libdc1394.dylib \
    -framework CoreServices -framework CoreFoundation -framework ApplicationServices -framework CoreAudio -framework OpenGL -framework CoreVideo \
    -framework IOKit -framework SystemConfiguration -framework Carbon -framework CoreText \
    -framework Cocoa,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" ...
    -I/usr/include -I/usr/local/include -I/Library/Frameworks/GStreamer.framework/Versions/Current/include/gstreamer-1.0 ...
    -I/Library/Frameworks/GStreamer.framework/Versions/Current/include/libxml2 -I/Library/Frameworks/GStreamer.framework/Versions/Current/include/glib-2.0 ...
    -I/Library/Frameworks/GStreamer.framework/Versions/Current/lib/glib-2.0/include ...
    -I/Library/Frameworks/GStreamer.framework/Versions/Current/lib/gstreamer-1.0/include/ ...
    -ICommon/Base -ICommon/Screen -IOSX/Base -IOSX/Screen -IOSX/Fonts -IOSX/EthernetAddress ...
    "OSX/Screen/*.c" "Common/Screen/*.c" "OSX/Base/*.c" "OSX/Fonts/*FontGlue*.c" "OSX/Fonts/FontInfo.c" "OSX/EthernetAddress/*.c" "Common/Base/*.c" ...
    -L/usr/lib -L/usr/local/lib ./PsychCocoaGlue.o

    movefile('./OSX/Screen/PsychCocoaGlue.m', './OSX/Screen/PsychCocoaGlue.c');

    osxsetoctaverpath('Screen');
    unix(['mv ../Projects/MacOSX/build/Screen.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);
end

if mode==1
    % Build GetSecs:
    mex -g -v --output ../Projects/MacOSX/build/GetSecs  -DPTBMODULE_GetSecs -DPTBOCTAVE3MEX "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/GetSecs "OSX/Base/*.c" "Common/Base/*.c" "Common/GetSecs/*.c"
    osxsetoctaverpath('GetSecs');
    unix(['mv ../Projects/MacOSX/build/GetSecs.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);
end

if mode==2
    % Build WaitSecs:
    mex -g -v --output ../Projects/MacOSX/build/WaitSecs  -DPTBMODULE_WaitSecs -DPTBOCTAVE3MEX "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/WaitSecs  "OSX/Base/*.c" "Common/Base/*.c" "Common/WaitSecs/*.c"
    osxsetoctaverpath('WaitSecs');
    unix(['mv ../Projects/MacOSX/build/WaitSecs.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);
end

if mode==3
    % Depends: PortAudio
    % Build PsychPortAudio:
    mex -g -v --output ../Projects/MacOSX/build/PsychPortAudio  -DPTBMODULE_PsychPortAudio -DPTBOCTAVE3MEX "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio -framework AudioToolbox -framework AudioUnit,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/PsychPortAudio  "OSX/Base/*.c" "Common/Base/*.c" "Common/PsychPortAudio/*.c" ../Cohorts/PortAudio/libportaudio_osx_64.a
    osxsetoctaverpath('PsychPortAudio');
    unix(['mv ../Projects/MacOSX/build/PsychPortAudio.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);
end

if mode==4
    % Depends: eyelink-SDK: eyelink_core.framework
    % Build Eyelink:
    mex -g -v --output ../Projects/MacOSX/build/Eyelink  -DPTBMODULE_Eyelink -DPTBOCTAVE3MEX "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio -framework eyelink_core,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -I/usr/local/include -ICommon/Base -IOSX/Base -ICommon/Eyelink  "OSX/Base/*.c" "Common/Base/*.c" "Common/Eyelink/*.c"
    osxsetoctaverpath('Eyelink');
    unix(['mv ../Projects/MacOSX/build/Eyelink.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);
end

if mode==5
    % Build IOPort:
    mex -g -v --output ../Projects/MacOSX/build/IOPort  -DPTBMODULE_IOPort -DPTBOCTAVE3MEX "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -I/usr/include  -ICommon/Base -IOSX/Base -ICommon/IOPort  "OSX/Base/*.c" "Common/Base/*.c" "Common/IOPort/*.c"
    osxsetoctaverpath('IOPort');
    unix(['mv ../Projects/MacOSX/build/IOPort.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);
end

if mode==6
    % Build moglcore:
    curdir = pwd;
    cd('../../Psychtoolbox/PsychOpenGL/MOGL/source/')
    try
        mex -g -v --output ./moglcore -DMACOSX -DGLEW_STATIC -DPTBOCTAVE3MEX "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework OpenGL -framework GLUT,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -I/usr/include moglcore.c gl_auto.c gl_manual.c glew.c mogl_rebinder.c ftglesGlue.c
    catch %#ok<*CTCH>
    end
    osxsetoctaverpath('moglcore', './');
    unix(['mv moglcore.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);
    cd(curdir);
end

if mode==7
    % Depends: libfreenect, libusb-1.0
    % Build PsychKinectCore:
    mex -g -v --output ../Projects/MacOSX/build/PsychKinectCore  -DPTBMODULE_PsychKinectCore -DPTBOCTAVE3MEX "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -I/usr/include -I/usr/local/include/libusb-1.0 -I/usr/local/include/libfreenect -ICommon/Base -IOSX/Base -ICommon/PsychKinect  "OSX/Base/*.c" "Common/Base/*.c" "Common/PsychKinect/*.c" -L/usr/local/lib -lfreenect -lusb-1.0
    osxsetoctaverpath('PsychKinectCore');
    unix(['mv ../Projects/MacOSX/build/PsychKinectCore.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);
end

if mode==8
    % Build PsychHID:
    % Dynamic link:   mex -g -v --output ../Projects/MacOSX/build/PsychHID  -DPTBMODULE_PsychHID -DPTBOCTAVE3MEX "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework ApplicationServices -framework Carbon -framework CoreAudio -framework IOKit,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/PsychHID -IOSX/PsychHID -I../Cohorts/HID_Utilities_64Bit/ -I../Cohorts/HID_Utilities_64Bit/IOHIDManager  "Common/PsychHID/*.c" "OSX/PsychHID/*.c" "OSX/Base/*.c" "Common/Base/*.c" -L../Cohorts/HID_Utilities_64Bit/build/Release -lHID_Utilities
    mex -g -v --output ../Projects/MacOSX/build/PsychHID  -DPTBMODULE_PsychHID -DPTBOCTAVE3MEX "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework ApplicationServices -framework CoreServices -framework CoreFoundation -framework Carbon -framework CoreAudio -framework IOKit,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/PsychHID -IOSX/PsychHID -I../Cohorts/HID_Utilities_64Bit/ -I../Cohorts/HID_Utilities_64Bit/IOHIDManager  "Common/PsychHID/*.c" "OSX/PsychHID/*.c" "OSX/Base/*.c" "Common/Base/*.c" ../Cohorts/HID_Utilities_64Bit/build/Release/libHID_Utilities64.a
    osxsetoctaverpath('PsychHID');
    unix(['mv ../Projects/MacOSX/build/PsychHID.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);
end

if mode==9
    % Build moalcore:
    curdir = pwd;
    cd('../../Psychtoolbox/PsychSound/MOAL/source/')
    try
        mex -g -v --output ./moalcore  -DMACOSX -DPTBOCTAVE3MEX "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework OpenAL -framework ApplicationServices -framework Carbon,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -I/usr/include moalcore.c al_auto.c al_manual.c alm.c
    catch
    end
    osxsetoctaverpath('moalcore', './');
    unix(['mv ./moalcore.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);

    cd(curdir);
end

if mode==10
    % Build FontInfo:
    mex -g -v --output ../Projects/MacOSX/build/FontInfo  -DPTBMODULE_FontInfo -DPTBOCTAVE3MEX "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework ApplicationServices -framework CoreAudio,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -I/usr/include -ICommon/Base -ICommon/Fonts -IOSX/Base -IOSX/Fonts  "Common/Base/*.c" "OSX/Base/*.c" "OSX/Fonts/*.c"
    osxsetoctaverpath('FontInfo');
    unix(['mv ../Projects/MacOSX/build/FontInfo.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);
end

if mode==11
    % Build MachAbsoluteTimeClockFrequency, MachGetPriorityMex,
    % MachSetPriorityMex, which are needed by the Priority() command on
    % OSX:

    % MachAbsoluteTimeClockFrequency
    mex -g -v --output ../Projects/MacOSX/build/MachAbsoluteTimeClockFrequency  -DPTBMODULE_MachAbsoluteTimeClockFrequency -DPTBOCTAVE3MEX "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/MachAbsoluteTimeClockFrequency  "OSX/Base/*.c" "Common/Base/*.c" "Common/MachAbsoluteTimeClockFrequency/*.c"
    osxsetoctaverpath('MachAbsoluteTimeClockFrequency');
    unix(['mv ../Projects/MacOSX/build/MachAbsoluteTimeClockFrequency.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);

    % MachGetPriorityMex
    mex -g -v --output ../Projects/MacOSX/build/MachGetPriorityMex -DPTBOCTAVE3MEX  "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -I/usr/include -ICommon/MachPriorityMex "Common/MachPriorityMex/MachPriorityCommonMex.c" "Common/MachPriorityMex/MachGetPriorityMex.c"
    osxsetoctaverpath('MachGetPriorityMex');
    unix(['mv ../Projects/MacOSX/build/MachGetPriorityMex.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);

    % MachSetPriorityMex
    mex -g -v --output ../Projects/MacOSX/build/MachSetPriorityMex -DPTBOCTAVE3MEX  "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -I/usr/include -ICommon/MachPriorityMex "Common/MachPriorityMex/MachPriorityCommonMex.c" "Common/MachPriorityMex/MachSetPriorityMex.c"
    osxsetoctaverpath('MachSetPriorityMex');
    unix(['mv ../Projects/MacOSX/build/MachSetPriorityMex.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);
end

if mode==12
    % Build Gestalt:
    mex -g -v --output ../Projects/MacOSX/build/Gestalt  -DPTBMODULE_Gestalt -DPTBOCTAVE3MEX "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -I/usr/include -ICommon/Base -IOSX/Base -IOSX/Gestalt -IOSX/OS9ToolboxFragments  "OSX/Base/*.c" "Common/Base/*.c" "OSX/Gestalt/*.c"
    osxsetoctaverpath('Gestalt');
    unix(['mv ../Projects/MacOSX/build/Gestalt.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);
end

if mode==13
    % Build pnet:
    curdir = pwd;
    cd('../../Psychtoolbox/PsychHardware/iViewXToolbox/tcp_udp_ip/')
    try
        mex -g -v --output ./pnet pnet.c "-Wl,-headerpad_max_install_names"
        osxsetoctaverpath('pnet', './');
        unix(['mv ./pnet.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);
    catch
    end
    cd(curdir);
end

if mode==14
    % Build PsychOculusVRCore.mex:
    % Depends on Oculus VR SDK v0.5
    try
        mex -g -v --output ../Projects/MacOSX/build/PsychOculusVRCore -DPTBMODULE_PsychOculusVRCore -DPTBOCTAVE3MEX "-Wl,-headerpad_max_install_names -F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio -framework LibOVR,-syslibroot,'/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk' -mmacosx-version-min='10.10'" -ICommon/Base -IOSX/Base -ICommon/PsychOculusVRCore -I/Library/Frameworks/LibOVR.framework/Versions/Current/Headers/ OSX/Base/*.c Common/Base/*.c Common/PsychOculusVRCore/*.c
    catch
        disp(psychlasterror);
    end

    osxsetoctaverpath('PsychOculusVRCore');
    unix(['mv ../Projects/MacOSX/build/PsychOculusVRCore.mex ' PsychtoolboxRoot 'PsychBasic/Octave4OSXFiles64/']);
end

delete('./*.o');

return;
