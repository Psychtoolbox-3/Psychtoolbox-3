function linuxmakeit_ubuntugutsy(mode)
% This is the GNU/Linux version of makeit to build the Linux
% version of PTB's Screen - command.
% This version is adapted to build on the "funky" Laptop under Ubuntu Linux 10.10 Maverick against Matlab R2007a and later 
if nargin < 1
    mode = 0
end;

fprintf('Building plugin type %i ...\n\n', mode);

if mode==0
    % Build Screen.mexglx:
    mex CFLAGS='$CFLAGS -fPIC -std=gnu99 -fexceptions' -v -outdir ../Projects/Linux/build/ -output Screen -DPTBMODULE_Screen -DPTB_USE_GSTREAMER -DPTBVIDEOCAPTURE_LIBDC -D_GNU_SOURCE -I/usr/X11R6/include -I/usr/include/gstreamer-0.10 -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include -I/usr/lib/i386-linux-gnu/glib-2.0/include -I/usr/include/libxml2 -ICommon/Base -ICommon/Screen -ILinux/Base -ILinux/Screen -L/usr/X11R6/lib Common/Base/*.cc Linux/Base/*.c Linux/Screen/*.c Common/Screen/*.c Common/Base/*.c -lc -lrt -ldl -lGL -lGLU -lX11 -lXext -lgstreamer-0.10 -lgstbase-0.10 -lgstapp-0.10 -lgstinterfaces-0.10 -lgobject-2.0 -lgmodule-2.0 -lxml2 -lgthread-2.0 -lglib-2.0 -lXxf86vm -ldc1394 -lusb-1.0 -lpciaccess -lXi -lXrandr -lXfixes
    unix(['mv ../Projects/Linux/build/Screen.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end;

if mode==1
    % Build GetSecs.mexglx:
    mex CFLAGS='$CFLAGS -fPIC -std=gnu99 -fexceptions' -v -outdir ../Projects/Linux/build/ -output GetSecs -DPTBMODULE_GetSecs -ICommon/Base -ILinux/Base -ICommon/GetSecs -ICommon/Screen Common/Base/*.cc Linux/Base/*.c Common/Base/*.c Common/GetSecs/*.c -lc -lrt 
    unix(['mv ../Projects/Linux/build/GetSecs.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end;

if mode==2
    % Build WaitSecs.mexglx:
    mex CFLAGS='$CFLAGS -fPIC -std=gnu99 -fexceptions' -v -outdir ../Projects/Linux/build/ -output WaitSecs -DPTBMODULE_WaitSecs -ICommon/Base -ILinux/Base -ICommon/WaitSecs -ICommon/Screen Common/Base/*.cc Linux/Base/*.c Common/Base/*.c Common/WaitSecs/*.c -lc -lrt 
    unix(['mv ../Projects/Linux/build/WaitSecs.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end;

if mode==3
    % Build PsychPortAudio.mexglx:
    mex CFLAGS='$CFLAGS -fPIC -std=gnu99 -fexceptions' -v -outdir ../Projects/Linux/build/ -output PsychPortAudio -DPTBMODULE_PsychPortAudio -ICommon/Base -ILinux/Base -ICommon/PsychPortAudio -ICommon/Screen Common/Base/*.cc Linux/Base/*.c Common/Base/*.c Common/PsychPortAudio/*.c /usr/local/lib/libportaudio.a -lc -lrt -lasound
    unix(['mv ../Projects/Linux/build/PsychPortAudio.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end

if mode==4
    % Build Eyelink.mexglx:
    mex CFLAGS='$CFLAGS -fPIC -std=gnu99 -fexceptions' -v -outdir ../Projects/Linux/build/ -output Eyelink -DPTBMODULE_Eyelink -ICommon/Base -ILinux/Base -ICommon/Eyelink -ICommon/Screen Common/Base/*.cc Linux/Base/*.c Common/Base/*.c Common/Eyelink/*.c -leyelink_core -lc -lrt
    unix(['mv ../Projects/Linux/build/Eyelink.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end

if mode==5
    % Build IOPort.mexglx:
    mex CFLAGS='$CFLAGS -fPIC -std=gnu99 -fexceptions' -v -outdir ../Projects/Linux/build/ -output IOPort -DPTBMODULE_IOPort -ICommon/Base -ILinux/Base -ICommon/IOPort -ICommon/Screen Common/Base/*.cc Linux/Base/*.c Common/Base/*.c Common/IOPort/*.c -lc -lrt
    unix(['mv ../Projects/Linux/build/IOPort.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end

if mode==6
    % Build moglcore.mexglx:
    curdir = pwd;
    cd('../../Psychtoolbox/PsychOpenGL/MOGL/source/')
    try
    mex CFLAGS='$CFLAGS -fPIC -fexceptions' -v -outdir ./ -output moglcore -DLINUX -DGLEW_STATIC -I/usr/X11R6/include -L/usr/X11R6/lib -lc -lGL -lGLU -lglut moglcore.c gl_auto.c gl_manual.c glew.c mogl_rebinder.c
    catch
    end
    unix(['mv moglcore.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
    cd(curdir);
end

if mode==7
    % Build PsychKinectCore.mexglx:
    mex CFLAGS='$CFLAGS -fPIC -std=gnu99 -fexceptions' -v -outdir ../Projects/Linux/build/ -output PsychKinectCore -DPTBMODULE_PsychKinectCore -I/usr/include/libusb-1.0 -I/usr/include/libfreenect -ICommon/Base -ILinux/Base -ICommon/PsychKinect -ICommon/Screen Common/Base/*.cc Linux/Base/*.c Common/Base/*.c Common/PsychKinect/*.c -lc -lrt -lfreenect -lusb-1.0
    unix(['mv ../Projects/Linux/build/PsychKinectCore.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end;

if mode==8
    % Build PsychHID.mexglx:
    mex CFLAGS='$CFLAGS -fPIC -std=gnu99 -fexceptions' -v -outdir ../Projects/Linux/build/ -output PsychHID -DPTBMODULE_PsychHID -I/usr/include/libusb-1.0 -ICommon/Base -ILinux/Base -ICommon/PsychHID -ILinux/PsychHID -ICommon/Screen Common/Base/*.cc Linux/Base/*.c Common/Base/*.c Common/PsychHID/*.c Linux/PsychHID/*.c -lc -ldl -lrt -lusb-1.0 -lX11 -lXi -lutil
    unix(['mv ../Projects/Linux/build/PsychHID.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);
end;

if mode==9
    % Build moalcore.mexglx:
    curdir = pwd;
    cd('../../Psychtoolbox/PsychSound/MOAL/source/')
    try
       mex CFLAGS='$CFLAGS -fPIC -fexceptions' -v -outdir ../Projects/Linux/build/ -output moalcore -DLINUX -lc -lopenal moalcore.c al_auto.c al_manual.c alm.c 
    catch
    end
    unix(['mv ../Projects/Linux/build/moalcore.' mexext ' ' PsychtoolboxRoot 'PsychBasic/']);

    cd(curdir);
end;

return;
