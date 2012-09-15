% Psychtoolbox:PsychSound:MOAL Contents of MOAL Matlab-OpenAL toolbox
%
% Moal is a collection of M-File wrappers and a MEX file that allow to call
% all OpenAL commands from Matlab as one is used to from the C programming
% language.
%
% Directory structure is as follows: 
% 
%     moaldemo.m          -- demonstration of how to use the toolbox
% 
%     core/
% 
%         (first group:  main toolbox functions)
% 
%         moalcore.mexmac -- main MEX interface to OpenAL functions
%         oalconst.mat    -- constants used by OpenAL routines
% 
%     source/
% 
%         (first group:  files that generate interface code)
%         
%         al_auto_init.c  -- file used in generating al_auto.c;  contains
%                            top portion of file, i.e., #includes, etc.
%         oalconst.m      -- MATLAB script that searches through OpenAL header
%                            files for #defined constants, and writes them
%                            to oalconst.mat as variables
% 
%         (second group:  files that compile to produce moalcore.mexmac)
% 
%         al_auto.c       -- automatically generated interfaces to OpenAL functions
%         al_manual.c     -- manually generated interfaces to OpenAL functions
%         alm.c           -- ALM library of ALC like functions.
%         moalcore.c      -- main MEX interface function
%         moaltypes.h     -- useful data types
%         windowshacks.c  -- hacks needed for Windows compatibility.
%
%         (third group:   Makefiles and build scripts.)
%         makefile        -- makefile to compile C files into moalcore.mexmac on PPC.
%         makefile_intelmac -- makefile for IntelMac.
%         makefile_linux    -- makefile for GNU/Linux.
%         makefile_linuxoctave -- makefile for Linux + Octave.
%         makefile_windows.m -- makefile for M$-Windows.
% 
%     wrap/*          -- wrapper M-files that check arguments, etc., and
%                        then call to moalcore.mexmac to run OpenAL functions
% 
%
% The following three commands will completely regenerate moal.
% 
% >> autocode(1,[],1)     % generate al_auto.c and wrapper M-files
% >> !make                % compile C code to produce MEX files
% >> oalconst             % save constants from header files in a .mat file

% 06-Feb-2007 -- created (MK)
