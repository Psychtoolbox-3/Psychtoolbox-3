#!/bin/bash

g++ -g -fPIC -I. -I/usr/include/ -I/usr/include/freetype2/ -I/usr/include/GLES -L/usr/lib libptbdrawtext_ftgles.cpp -lGLESv1_CM -lfontconfig -lfreetype -pie -shared -Wl,-Bsymbolic -Wl,-Bsymbolic-functions -Wl,--version-script=linuxexportlist.txt -o ../../../Psychtoolbox/PsychBasic/PsychPlugins/libptbdrawtext_ftgles_arm.so.1 libftgles_arm32.a