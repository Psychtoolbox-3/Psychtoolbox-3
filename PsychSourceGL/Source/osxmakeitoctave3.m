function osxmakeitoctave3(mode)
% This is the MacOS/X version of makeit: It is meant for building PTB for
% 64-Bit Octave on OSX IntelMac, 10.7 Lion or later:

if ~IsOSX(1) || ~IsOctave
    error('osxmakeit only works with a 64-Bit version of Octave for OSX!');
end

if nargin < 1
    mode = 0;
end

fprintf('Building plugin type %i ...\n\n', mode);

copyfile('./Common/Base/PsychScriptingGlue.cc', './Common/Base/PsychScriptingGlue.c');

if mode==0
    % Build Screen:
    % Depends: GStreamer-0.10, libdc1394-2, libusb-1.0

    % Need to build the PsychCocoaGlue.c separately as Objective-C code,
    % whereas the mex call below only compiles C-Code and links it with the
    % PsychCocoaGlue.o object code created here. This to avoid linker
    % trouble with unknown/illegal Objective-C language command line switches:
    mex -g  -v -c -DPTBMODULE_Screen -DPTB_USE_GSTREAMER -DPTBVIDEOCAPTURE_LIBDC -DPTBOCTAVE3MEX ...
        -I/usr/include -ICommon/Base -ICommon/Screen -IOSX/Base -IOSX/Screen -IOSX/Fonts -IOSX/EthernetAddress ...
        "-W, -x objective-c -mmacosx-version-min='10.5'" "OSX/Screen/PsychCocoaGlue.c"

    movefile('./OSX/Screen/PsychCocoaGlue.c', './OSX/Screen/PsychCocoaGlue.m');

    % Build with Weak linking of GStreamer and libdc1394 via -weak_library
    % flag. This means that a missing GStreamer or libdc1394 installation
    % will not cause linker failure at Screen() load time -- Screen
    % continues to work as long as no GStreamer-dependent functions are
    % used. In that case it will abort with a helpful error message:
    mex -g -v --output ../Projects/MacOSX/build/Screen -DPTBMODULE_Screen -DPTB_USE_GSTREAMER -DPTBVIDEOCAPTURE_LIBDC -DPTBOCTAVE3MEX ...
        "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgstreamer-0.10.0.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgstbase-0.10.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgstapp-0.10.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgstinterfaces-0.10.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgobject-2.0.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgmodule-2.0.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libxml2.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgthread-2.0.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libglib-2.0.dylib ...
        -weak_library /usr/local/lib/libdc1394.dylib ...
        -framework CoreServices -framework CoreFoundation -framework ApplicationServices -framework CoreAudio -framework OpenGL -framework CoreVideo ...
        -framework IOKit -framework AGL -framework SystemConfiguration -framework Carbon ...
        -framework Cocoa,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" ...
        -I/usr/include -I/Library/Frameworks/GStreamer.framework/Versions/Current/include/gstreamer-0.10 ...
        -I/Library/Frameworks/GStreamer.framework/Versions/Current/include/libxml2 -I/Library/Frameworks/GStreamer.framework/Versions/Current/include/glib-2.0 ...
        -I/Library/Frameworks/GStreamer.framework/Versions/Current/lib/glib-2.0/include ...
        -ICommon/Base -ICommon/Screen -IOSX/Base -IOSX/Screen -IOSX/Fonts -IOSX/EthernetAddress ...
        "OSX/Screen/*.c" "Common/Screen/*.c" "OSX/Base/*.c" "OSX/Fonts/*FontGlue*.c" "OSX/Fonts/FontInfo.c" "OSX/EthernetAddress/*.c" "Common/Base/*.c" ...
        "-W, -mmacosx-version-min='10.5'" ...
        -L/usr/lib -L/usr/local/lib ./PsychCocoaGlue.o
    
    movefile('./OSX/Screen/PsychCocoaGlue.m', './OSX/Screen/PsychCocoaGlue.c');
    unix(['mv ../Projects/MacOSX/build/Screen.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles64/']);
end

if mode==1
    % Build GetSecs:
    mex -g -v --output ../Projects/MacOSX/build/GetSecs  -DPTBMODULE_GetSecs -DPTBOCTAVE3MEX "-W, -mmacosx-version-min='10.5'" "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/GetSecs "OSX/Base/*.c" "Common/Base/*.c" "Common/GetSecs/*.c"
    unix(['mv ../Projects/MacOSX/build/GetSecs.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles64/']);
end

if mode==2
    % Build WaitSecs:
    mex -g -v --output ../Projects/MacOSX/build/WaitSecs  -DPTBMODULE_WaitSecs -DPTBOCTAVE3MEX "-W, -mmacosx-version-min='10.5'" "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/WaitSecs  "OSX/Base/*.c" "Common/Base/*.c" "Common/WaitSecs/*.c"
    unix(['mv ../Projects/MacOSX/build/WaitSecs.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles64/']);
end

if mode==3
    % Depends: PortAudio
    % Build PsychPortAudio:
    mex -g -v --output ../Projects/MacOSX/build/PsychPortAudio  -DPTBMODULE_PsychPortAudio -DPTBOCTAVE3MEX "-W, -mmacosx-version-min='10.5'" "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio -framework AudioToolbox -framework AudioUnit,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/PsychPortAudio  "OSX/Base/*.c" "Common/Base/*.c" "Common/PsychPortAudio/*.c" ../Cohorts/PortAudio/libportaudio_osx_64.a
    unix(['mv ../Projects/MacOSX/build/PsychPortAudio.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles64/']);
end

if mode==4
    % Depends: eyelink-SDK: eyelink_core.framework
    % Build Eyelink:
    mex -g -v --output ../Projects/MacOSX/build/Eyelink  -DPTBMODULE_Eyelink -DPTBOCTAVE3MEX "-W, -mmacosx-version-min='10.5'" "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio -framework eyelink_core,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/Eyelink  "OSX/Base/*.c" "Common/Base/*.c" "Common/Eyelink/*.c"
    unix(['mv ../Projects/MacOSX/build/Eyelink.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles64/']);
end

if mode==5
    % Build IOPort:
    mex -g -v --output ../Projects/MacOSX/build/IOPort  -DPTBMODULE_IOPort -DPTBOCTAVE3MEX "-W, -mmacosx-version-min='10.5'" "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" -I/usr/include  -ICommon/Base -IOSX/Base -ICommon/IOPort  "OSX/Base/*.c" "Common/Base/*.c" "Common/IOPort/*.c"
    unix(['mv ../Projects/MacOSX/build/IOPort.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles64/']);
end

if mode==6
    % Build moglcore:
    curdir = pwd;
    cd('../../Psychtoolbox/PsychOpenGL/MOGL/source/')
    try
        mex -g -v --output ./moglcore -DMACOSX -DGLEW_STATIC -DPTBOCTAVE3MEX "-W, -mmacosx-version-min='10.5'" "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework OpenGL -framework GLUT,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" -I/usr/include moglcore.c gl_auto.c gl_manual.c glew.c mogl_rebinder.c
    catch %#ok<*CTCH>
    end
    unix(['mv moglcore.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles64/']);
    cd(curdir);
end

if mode==7
    % Depends: libfreenect, libusb-1.0
    % Build PsychKinectCore:
    mex -g -v --output ../Projects/MacOSX/build/PsychKinectCore  -DPTBMODULE_PsychKinectCore -DPTBOCTAVE3MEX "-W, -mmacosx-version-min='10.5'" "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" -I/usr/include -I/usr/local/include/libusb-1.0 -I/usr/local/include/libfreenect -ICommon/Base -IOSX/Base -ICommon/PsychKinect  "OSX/Base/*.c" "Common/Base/*.c" "Common/PsychKinect/*.c" -L/usr/local/lib -lfreenect -lusb-1.0
    unix(['mv ../Projects/MacOSX/build/PsychKinectCore.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles64/']);
end

if mode==8
    % Build PsychHID:
    % Dynamic link:   mex -g -v --output ../Projects/MacOSX/build/PsychHID  -DPTBMODULE_PsychHID -DPTBOCTAVE3MEX "-W, -mmacosx-version-min='10.5'" "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework ApplicationServices -framework Carbon -framework CoreAudio -framework IOKit,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/PsychHID -IOSX/PsychHID -I../Cohorts/HID_Utilities_64Bit/ -I../Cohorts/HID_Utilities_64Bit/IOHIDManager  "Common/PsychHID/*.c" "OSX/PsychHID/*.c" "OSX/Base/*.c" "Common/Base/*.c" -L../Cohorts/HID_Utilities_64Bit/build/Release -lHID_Utilities
    mex -g -v --output ../Projects/MacOSX/build/PsychHID  -DPTBMODULE_PsychHID -DPTBOCTAVE3MEX "-W, -mmacosx-version-min='10.5'" "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework ApplicationServices -framework Carbon -framework CoreAudio -framework IOKit,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/PsychHID -IOSX/PsychHID -I../Cohorts/HID_Utilities_64Bit/ -I../Cohorts/HID_Utilities_64Bit/IOHIDManager  "Common/PsychHID/*.c" "OSX/PsychHID/*.c" "OSX/Base/*.c" "Common/Base/*.c" ../Cohorts/HID_Utilities_64Bit/build/Release/libHID_Utilities64.a
    unix(['mv ../Projects/MacOSX/build/PsychHID.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles64/']);
end

if mode==9
    % Build moalcore:
    curdir = pwd;
    cd('../../Psychtoolbox/PsychSound/MOAL/source/')
    try
        mex -g -v --output ./moalcore  -DMACOSX -DPTBOCTAVE3MEX "-W, -mmacosx-version-min='10.5'" "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework OpenAL -framework ApplicationServices -framework Carbon,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" -I/usr/include moalcore.c al_auto.c al_manual.c alm.c
    catch
    end
    unix(['mv ./moalcore.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles64/']);

    cd(curdir);
end

if mode==10
    % Build FontInfo:
    mex -g -v --output ../Projects/MacOSX/build/FontInfo  -DPTBMODULE_FontInfo -DPTBOCTAVE3MEX "-W, -mmacosx-version-min='10.5'" "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework ApplicationServices -framework CoreAudio,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" -I/usr/include -ICommon/Base -ICommon/Fonts -IOSX/Base -IOSX/Fonts  "Common/Base/*.c" "OSX/Base/*.c" "OSX/Fonts/*.c"
    unix(['mv ../Projects/MacOSX/build/FontInfo.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles64/']);
end

if mode==11
    % Build MachAbsoluteTimeClockFrequency, MachGetPriorityMex,
    % MachSetPriorityMex, which are needed by the Priority() command on
    % OSX:

    % MachAbsoluteTimeClockFrequency
    mex -g -v --output ../Projects/MacOSX/build/MachAbsoluteTimeClockFrequency  -DPTBMODULE_MachAbsoluteTimeClockFrequency -DPTBOCTAVE3MEX "-W, -mmacosx-version-min='10.5'" "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/MachAbsoluteTimeClockFrequency  "OSX/Base/*.c" "Common/Base/*.c" "Common/MachAbsoluteTimeClockFrequency/*.c"
    unix(['mv ../Projects/MacOSX/build/MachAbsoluteTimeClockFrequency.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles64/']);

    % MachGetPriorityMex
    mex -g -v --output ../Projects/MacOSX/build/MachGetPriorityMex -DPTBOCTAVE3MEX "-W, -mmacosx-version-min='10.5'"  "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" -I/usr/include -ICommon/MachPriorityMex "Common/MachPriorityMex/MachPriorityCommonMex.c" "Common/MachPriorityMex/MachGetPriorityMex.c"
    unix(['mv ../Projects/MacOSX/build/MachGetPriorityMex.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles64/']);

    % MachSetPriorityMex
    mex -g -v --output ../Projects/MacOSX/build/MachSetPriorityMex -DPTBOCTAVE3MEX "-W, -mmacosx-version-min='10.5'"  "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" -I/usr/include -ICommon/MachPriorityMex "Common/MachPriorityMex/MachPriorityCommonMex.c" "Common/MachPriorityMex/MachSetPriorityMex.c"
    unix(['mv ../Projects/MacOSX/build/MachSetPriorityMex.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles64/']);
end

if mode==12
    % Build Gestalt:
    mex -g -v --output ../Projects/MacOSX/build/Gestalt  -DPTBMODULE_Gestalt -DPTBOCTAVE3MEX "-W, -mmacosx-version-min='10.5'" "-Wl,-F/System/Library/Frameworks/ -F/Library/Frameworks/ -framework CoreServices -framework CoreFoundation -framework CoreAudio,-syslibroot,'/Developer/SDKs/MacOSX10.6.sdk' -mmacosx-version-min='10.5'" -I/usr/include -ICommon/Base -IOSX/Base -IOSX/Gestalt -IOSX/OS9ToolboxFragments  "OSX/Base/*.c" "Common/Base/*.c" "OSX/Gestalt/*.c"
    unix(['mv ../Projects/MacOSX/build/Gestalt.' mexext ' ' PsychtoolboxRoot 'PsychBasic/Octave3OSXFiles64/']);
end

delete('./Common/Base/PsychScriptingGlue.c');
delete('./*.o');

return;
