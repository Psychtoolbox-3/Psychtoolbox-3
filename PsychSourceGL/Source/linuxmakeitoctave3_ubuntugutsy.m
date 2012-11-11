function linuxmakeitoctave3_ubuntugutsy(mode)
% This is the GNU/Linux version of makeit to build the Linux
% version of PTB's binary plugins for GNU/Octave V3.2.0 and later.
% This version is adapted to build on the "funky" Laptop under Ubuntu Linux
% 10.10 Maverick or any Linux system with similar configuration.
%
% CAUTION: You *MUST* make sure that the MEX files are compiled *without*
% compiler optimization, ie., with a level of -O0 !!! Otherwise the files
% will be *miscompiled* !!! At first glance they seem to work, but the will
% malfunction and mysteriously crash in many ways, almost impossible to debug!!
%
% UPDATE: This seems to be no longer true for Octave 3.2.4 built on 10.10.
%         Our compiler settings are -O2 on this system and so far no bugs
%         or miscompiles were found during daily use and testing.
%
%         However, as -O2 never showed a significant speedup vs. -O0, it
%         pretty much doesn't matter if we do optimized builds or not.

if nargin < 1
    mode = 0;
end;

fprintf('Building plugin type %i ...\n\n', mode);

% Target folder depends if this is a 64 bit or 32 bit runtime:
if ~isempty(findstr(computer, '_64'))
    target = 'PsychBasic/Octave3LinuxFiles64/';
else
    target = 'PsychBasic/Octave3LinuxFiles/';
end

if mode==0
    % Build Screen.mex:
    mex -v -g "-W -std=gnu99" --output ../Projects/Linux/build/Screen.mex -DPTBMODULE_Screen -DPTB_USE_GSTREAMER -DPTBVIDEOCAPTURE_LIBDC -DPTBOCTAVE3MEX -D_GNU_SOURCE -I/usr/X11R6/include -I/usr/include/gstreamer-0.10 -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include -I/usr/lib/i386-linux-gnu/glib-2.0/include -I/usr/include/libxml2 -ICommon/Base -ICommon/Screen -ILinux/Base -ILinux/Screen -L/usr/X11R6/lib  Common/Base/*.cc Linux/Base/*.c Linux/Screen/*.c Common/Screen/*.c Common/Base/*.c -lc -ldl -lrt -lGL -lGLU -lX11 -lXext -lgstreamer-0.10 -lgstbase-0.10 -lgstapp-0.10 -lgstinterfaces-0.10 -lgobject-2.0 -lgmodule-2.0 -lxml2 -lgthread-2.0 -lglib-2.0 -lXxf86vm -ldc1394 -lusb-1.0 -lpciaccess -lXi -lXrandr -lXfixes
    unix(['mv ../Projects/Linux/build/Screen.mex ' PsychtoolboxRoot target]);
    striplibsfrommexfile([PsychtoolboxRoot target 'Screen.mex']);
end;

if mode==1
    % Build GetSecs.mex:
    mex -v -g --output ../Projects/Linux/build/GetSecs.mex -DPTBMODULE_GetSecs -DPTBOCTAVE3MEX -ICommon/Base -ILinux/Base -ICommon/GetSecs -ICommon/Screen Common/Base/*.cc Linux/Base/*.c Common/Base/*.c Common/GetSecs/*.c -lc -lrt 
    unix(['mv ../Projects/Linux/build/GetSecs.mex ' PsychtoolboxRoot target]);
    striplibsfrommexfile([PsychtoolboxRoot target 'GetSecs.mex']);
end;

if mode==2
    % Build WaitSecs.mex:
    mex -v -g --output ../Projects/Linux/build/WaitSecs.mex -DPTBMODULE_WaitSecs -DPTBOCTAVE3MEX -ICommon/Base -ILinux/Base -ICommon/WaitSecs -ICommon/Screen Common/Base/*.cc Linux/Base/*.c Common/Base/*.c Common/WaitSecs/*.c -lc -lrt 
    unix(['mv ../Projects/Linux/build/WaitSecs.mex ' PsychtoolboxRoot target]);
    striplibsfrommexfile([PsychtoolboxRoot target 'WaitSecs.mex']);
end;

if mode==3
    % Build PsychPortAudio.mex:
    mex -v -g --output ../Projects/Linux/build/PsychPortAudio.mex -DPTBMODULE_PsychPortAudio -DPTBOCTAVE3MEX -ICommon/Base -ILinux/Base -ICommon/PsychPortAudio -ICommon/Screen Common/Base/*.cc Linux/Base/*.c Common/Base/*.c Common/PsychPortAudio/*.c /usr/local/lib/libportaudio.a -lc -lrt -lasound
    unix(['mv ../Projects/Linux/build/PsychPortAudio.mex ' PsychtoolboxRoot target]);
    striplibsfrommexfile([PsychtoolboxRoot target 'PsychPortAudio.mex']);
end

if mode==4
    % Build Eyelink.mex:
    mex -v -g --output ../Projects/Linux/build/Eyelink.mex -DPTBMODULE_Eyelink -DPTBOCTAVE3MEX -ICommon/Base -ILinux/Base -ICommon/Eyelink -ICommon/Screen Common/Base/*.cc Linux/Base/*.c Common/Base/*.c Common/Eyelink/*.c -leyelink_core -lc -lrt
    unix(['mv ../Projects/Linux/build/Eyelink.mex ' PsychtoolboxRoot target]);
    striplibsfrommexfile([PsychtoolboxRoot target 'Eyelink.mex']);
end

if mode==5
    % Build IOPort.mex:
    mex -v -g --output ../Projects/Linux/build/IOPort.mex -DPTBMODULE_IOPort -DPTBOCTAVE3MEX -ICommon/Base -ILinux/Base -ICommon/IOPort -ICommon/Screen Common/Base/*.cc Linux/Base/*.c Common/Base/*.c Common/IOPort/*.c -lc -lrt
    unix(['mv ../Projects/Linux/build/IOPort.mex ' PsychtoolboxRoot target]);
    striplibsfrommexfile([PsychtoolboxRoot target 'IOPort.mex']);
end

if mode==6
    % Build moglcore.mex:
    curdir = pwd;
    cd('../../Psychtoolbox/PsychOpenGL/MOGL/source/')
    try
       mex -v -g --output moglcore.mex -DLINUX -DGLEW_STATIC -DPTBOCTAVE3MEX -I/usr/X11R6/include -L/usr/X11R6/lib -lc -lGL -lGLU -lglut moglcore.c gl_auto.c gl_manual.c glew.c mogl_rebinder.c
    catch
    end
    unix(['mv moglcore.mex ' PsychtoolboxRoot target]);
    cd(curdir);
    striplibsfrommexfile([PsychtoolboxRoot target 'moglcore.mex']);
end;

if mode==7
    % Build PsychKinectCore.mex:
    mex -v -g --output ../Projects/Linux/build/PsychKinectCore.mex -DPTBMODULE_PsychKinectCore -DPTBOCTAVE3MEX -I/usr/include/libusb-1.0 -I/usr/include/libfreenect -ICommon/Base -ILinux/Base -ICommon/PsychKinect -ICommon/Screen Common/Base/*.cc Linux/Base/*.c Common/Base/*.c Common/PsychKinect/*.c -lc -lrt -lfreenect -lusb-1.0
    unix(['mv ../Projects/Linux/build/PsychKinectCore.mex ' PsychtoolboxRoot target]);
    striplibsfrommexfile([PsychtoolboxRoot target 'PsychKinectCore.mex']);
end;

if mode==8
    % Build PsychHID.mex:
    mex -v -g --output ../Projects/Linux/build/PsychHID.mex -DPTBMODULE_PsychHID -DPTBOCTAVE3MEX -I/usr/include/libusb-1.0 -ICommon/Base -ILinux/Base -ICommon/PsychHID -ILinux/PsychHID -ICommon/Screen Common/Base/*.cc Linux/Base/*.c Common/Base/*.c Common/PsychHID/*.c Linux/PsychHID/*.c -lc -ldl -lrt -lusb-1.0 -lX11 -lXi -lutil
    unix(['mv ../Projects/Linux/build/PsychHID.mex ' PsychtoolboxRoot target]);
    striplibsfrommexfile([PsychtoolboxRoot target 'PsychHID.mex']);
end;

if mode==9
    % Build moalcore.mex:
    curdir = pwd;
    cd('../../Psychtoolbox/PsychSound/MOAL/source/')
    try
       mex -v -g --output moalcore.mex -DLINUX -DPTBOCTAVE3MEX -lc -lopenal moalcore.c al_auto.c al_manual.c alm.c 
    catch
    end
    unix(['mv moalcore.mex ' PsychtoolboxRoot target]);
    cd(curdir);
    striplibsfrommexfile([PsychtoolboxRoot target 'moalcore.mex']);
end;

return;
