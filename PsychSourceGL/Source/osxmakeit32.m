function osxmakeit32(mode)
% This is the MacOS/X version of makeit for legacy 32-Bit Matlab:
%
% It is used for building against a 32-Bit Matlab R2007a on a legacy
% 10.4.x "Tiger" or 10.5.x "Leopard" OSX system.

if ~IsOSX || IsOSX(1) || IsOctave
    error('osxmakeit32 only works with a 32-Bit version of Matlab for OSX!');
end

if nargin < 1
    mode = 0;
end;

fprintf('Building plugin type %i ...\n\n', mode);

if mode==0
    % Build Screen:
    % Depends: GStreamer-0.10, libdc1394-2

    % Build with Weak linking of GStreamer and libdc1394 via -weak_library
    % flag. This means that a missing GStreamer or libdc1394 installation
    % will not cause linker failure at Screen() load time -- Screen
    % continues to work as long as no GStreamer-dependent functions are
    % used. In that case it will abort with a helpful error message:
    mex -v -outdir ../Projects/MacOSX/build -output Screen -largeArrayDims -DPTBMODULE_Screen -DPTB_USE_GSTREAMER -DPTBVIDEOCAPTURE_LIBDC ...
        LDFLAGS="\$LDFLAGS -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgstreamer-0.10.0.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgstbase-0.10.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgstapp-0.10.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgstinterfaces-0.10.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgstpbutils-0.10.dylib ...        
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgobject-2.0.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgmodule-2.0.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libxml2.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libgthread-2.0.dylib ...
        -weak_library /Library/Frameworks/GStreamer.framework/Versions/Current/lib/libglib-2.0.dylib ...
        -framework CoreServices -framework CoreFoundation -framework ApplicationServices -framework CoreAudio -framework OpenGL -framework CoreVideo ...
        -framework IOKit -framework AGL -framework SystemConfiguration -framework Carbon -framework Quicktime" ...
        -I/usr/include -I../Cohorts/libDC1394/libdc1394-2.0.1 ../Cohorts/libDC1394/libdc1394.a ....
        -I/Library/Frameworks/GStreamer.framework/Versions/Current/include/gstreamer-0.10 ...
        -I/Library/Frameworks/GStreamer.framework/Versions/Current/include/libxml2 -I/Library/Frameworks/GStreamer.framework/Versions/Current/include/glib-2.0 ...
        -I/Library/Frameworks/GStreamer.framework/Versions/Current/lib/glib-2.0/include ...
        -ICommon/Base -ICommon/Screen -IOSX/Base -IOSX/Screen -IOSX/Fonts -IOSX/EthernetAddress ...
        "Common/Base/*.cc" "OSX/Screen/*.c" "Common/Screen/*.c" "OSX/Base/*.c" "OSX/Fonts/*FontGlue*.c" "OSX/Fonts/FontInfo.c" "OSX/EthernetAddress/*.c" "Common/Base/*.c" ...
        -L/usr/local/lib

%     mex -v -outdir ../Projects/MacOSX/build -output Screen -largeArrayDims -DPTBMODULE_Screen -DPTBVIDEOCAPTURE_LIBDC ...
%         LDFLAGS="\$LDFLAGS ...
%         -framework CoreServices -framework CoreFoundation -framework ApplicationServices -framework CoreAudio -framework OpenGL -framework CoreVideo ...
%         -framework IOKit -framework AGL -framework SystemConfiguration -framework Carbon -framework Quicktime" ...
%         -I/usr/include -I../Cohorts/libDC1394/libdc1394-2.0.1 ../Cohorts/libDC1394/libdc1394.a ....
%         -ICommon/Base -ICommon/Screen -IOSX/Base -IOSX/Screen -IOSX/Fonts -IOSX/EthernetAddress ...
%         "Common/Base/*.cc" "OSX/Screen/*.c" "Common/Screen/*.c" "OSX/Base/*.c" "OSX/Fonts/*FontGlue*.c" "OSX/Fonts/FontInfo.c" "OSX/EthernetAddress/*.c" "Common/Base/*.c" ...
%         -L/usr/local/lib
    
    unix(['mv ../Projects/MacOSX/build/Screen.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end

if mode==1
    % Build GetSecs:
    mex -v -outdir ../Projects/MacOSX/build -output GetSecs -largeArrayDims -DPTBMODULE_GetSecs LDFLAGS="\$LDFLAGS -framework CoreServices -framework CoreFoundation -framework CoreAudio" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/GetSecs "Common/Base/*.cc" "OSX/Base/*.c" "Common/Base/*.c" "Common/GetSecs/*.c"
    unix(['mv ../Projects/MacOSX/build/GetSecs.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end

if mode==2
    % Build WaitSecs:
    mex -v -outdir ../Projects/MacOSX/build -output WaitSecs -largeArrayDims -DPTBMODULE_WaitSecs LDFLAGS="\$LDFLAGS -framework CoreServices -framework CoreFoundation -framework CoreAudio" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/WaitSecs "Common/Base/*.cc" "OSX/Base/*.c" "Common/Base/*.c" "Common/WaitSecs/*.c"
    unix(['mv ../Projects/MacOSX/build/WaitSecs.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end

if mode==3
    % Depends: PortAudio
    % Build PsychPortAudio:
    mex -v -outdir ../Projects/MacOSX/build -output PsychPortAudio -largeArrayDims -DPTBMODULE_PsychPortAudio LDFLAGS="\$LDFLAGS -framework CoreServices -framework CoreFoundation -framework CoreAudio -framework AudioToolbox -framework AudioUnit" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/PsychPortAudio "Common/Base/*.cc" "OSX/Base/*.c" "Common/Base/*.c" "Common/PsychPortAudio/*.c" ../Cohorts/PortAudio/libportaudio.a
    unix(['mv ../Projects/MacOSX/build/PsychPortAudio.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end

if mode==4
    % Depends: eyelink-SDK: eyelink_core.framework
    % Build Eyelink:
    mex -v -outdir ../Projects/MacOSX/build -output Eyelink -largeArrayDims -DPTBMODULE_Eyelink LDFLAGS="\$LDFLAGS -framework CoreServices -framework CoreFoundation -framework CoreAudio -framework eyelink_core" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/Eyelink "Common/Base/*.cc" "OSX/Base/*.c" "Common/Base/*.c" "Common/Eyelink/*.c"
    unix(['mv ../Projects/MacOSX/build/Eyelink.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end

if mode==5
    % Build IOPort:
    mex -v -outdir ../Projects/MacOSX/build -output IOPort -largeArrayDims -DPTBMODULE_IOPort LDFLAGS="\$LDFLAGS -framework CoreServices -framework CoreFoundation -framework CoreAudio" -I/usr/include  -ICommon/Base -IOSX/Base -ICommon/IOPort "Common/Base/*.cc" "OSX/Base/*.c" "Common/Base/*.c" "Common/IOPort/*.c"
    unix(['mv ../Projects/MacOSX/build/IOPort.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end

if mode==6
    % Build moglcore:
    curdir = pwd;
    cd('../../Psychtoolbox/PsychOpenGL/MOGL/source/')
    try
        mex -v -outdir ./ -output moglcore -DMACOSX -DGLEW_STATIC -largeArrayDims -f ../../../../PsychSourceGL/Source/mexopts.sh LDFLAGS="\$LDFLAGS -framework OpenGL -framework GLUT" -I/usr/include moglcore.c gl_auto.c gl_manual.c glew.c mogl_rebinder.c ftglesGlue.c
    catch %#ok<*CTCH>
    end
    unix(['mv moglcore.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
    cd(curdir);
end

if mode==7
    % Depends: libfreenect, libusb-1.0
    % Build PsychKinectCore:
    mex -v -outdir ../Projects/MacOSX/build -output PsychKinectCore -largeArrayDims -DPTBMODULE_PsychKinectCore LDFLAGS="\$LDFLAGS -framework CoreServices -framework CoreFoundation -framework CoreAudio" -I/usr/include -I/usr/local/include/libusb-1.0 -I/usr/local/include/libfreenect -ICommon/Base -IOSX/Base -ICommon/PsychKinect "Common/Base/*.cc" "OSX/Base/*.c" "Common/Base/*.c" "Common/PsychKinect/*.c" -L/usr/local/lib -lfreenect -lusb-1.0
    unix(['mv ../Projects/MacOSX/build/PsychKinectCore.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end

if mode==8
    % Build PsychHID:
    mex -v -outdir ../Projects/MacOSX/build -output PsychHID -largeArrayDims -DPTBMODULE_PsychHID LDFLAGS="\$LDFLAGS -framework ApplicationServices -framework Carbon -framework CoreAudio -framework IOKit" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/PsychHID -IOSX/PsychHID -I"../Cohorts/HID Utilities Source/" "Common/Base/*.cc" "Common/PsychHID/*.c" "OSX/PsychHID/*.c" "OSX/Base/*.c" "Common/Base/*.c" ../Projects/MacOSX/PsychToolbox/libHIDUtilities-Intel.a
    unix(['mv ../Projects/MacOSX/build/PsychHID.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end

if mode==9
    % Build moalcore:
    curdir = pwd;
    cd('../../Psychtoolbox/PsychSound/MOAL/source/')
    try
        mex -v -outdir ./ -output moalcore -largeArrayDims -DMACOSX -f ../../../../PsychSourceGL/Source/mexopts.sh LDFLAGS="\$LDFLAGS -framework OpenAL -framework ApplicationServices -framework Carbon" -I/usr/include moalcore.c al_auto.c al_manual.c alm.c
    catch
    end
    unix(['mv ./moalcore.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);

    cd(curdir);
end

if mode==10
    % Build FontInfo:
    mex -v -outdir ../Projects/MacOSX/build -output FontInfo -largeArrayDims -DPTBMODULE_FontInfo LDFLAGS="\$LDFLAGS -framework ApplicationServices -framework CoreAudio" -I/usr/include -ICommon/Base -ICommon/Fonts -IOSX/Base -IOSX/Fonts "Common/Base/*.cc" "Common/Base/*.c" "OSX/Base/*.c" "OSX/Fonts/*.c"
    unix(['mv ../Projects/MacOSX/build/FontInfo.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end

if mode==11
    % Build MachAbsoluteTimeClockFrequency, MachGetPriorityMex,
    % MachSetPriorityMex, which are needed by the Priority() command on
    % OSX:

    % MachAbsoluteTimeClockFrequency
    mex -v -outdir ../Projects/MacOSX/build -output MachAbsoluteTimeClockFrequency -largeArrayDims -DPTBMODULE_MachAbsoluteTimeClockFrequency LDFLAGS="\$LDFLAGS -framework CoreServices -framework CoreFoundation -framework CoreAudio" -I/usr/include -ICommon/Base -IOSX/Base -ICommon/MachAbsoluteTimeClockFrequency "Common/Base/*.cc" "OSX/Base/*.c" "Common/Base/*.c" "Common/MachAbsoluteTimeClockFrequency/*.c"
    unix(['mv ../Projects/MacOSX/build/MachAbsoluteTimeClockFrequency.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);

    % MachGetPriorityMex
    mex -v -outdir ../Projects/MacOSX/build -output MachGetPriorityMex -largeArrayDims LDFLAGS="\$LDFLAGS -framework CoreServices -framework CoreFoundation -framework CoreAudio" -I/usr/include -ICommon/MachPriorityMex "Common/MachPriorityMex/MachPriorityCommonMex.c" "Common/MachPriorityMex/MachGetPriorityMex.c"
    unix(['mv ../Projects/MacOSX/build/MachGetPriorityMex.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);

    % MachSetPriorityMex
    mex -v -outdir ../Projects/MacOSX/build -output MachSetPriorityMex -largeArrayDims LDFLAGS="\$LDFLAGS -framework CoreServices -framework CoreFoundation -framework CoreAudio" -I/usr/include -ICommon/MachPriorityMex "Common/MachPriorityMex/MachPriorityCommonMex.c" "Common/MachPriorityMex/MachSetPriorityMex.c"
    unix(['mv ../Projects/MacOSX/build/MachSetPriorityMex.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end

if mode==12
    % Build Gestalt:
    mex -v -outdir ../Projects/MacOSX/build -output Gestalt -largeArrayDims -DPTBMODULE_Gestalt LDFLAGS="\$LDFLAGS -framework CoreServices -framework CoreFoundation -framework CoreAudio" -I/usr/include -ICommon/Base -IOSX/Base -IOSX/Gestalt -IOSX/OS9ToolboxFragments "Common/Base/*.cc" "OSX/Base/*.c" "Common/Base/*.c" "OSX/Gestalt/*.c"
    unix(['mv ../Projects/MacOSX/build/Gestalt.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end

return;
