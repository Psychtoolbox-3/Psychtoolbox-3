# setup.py -- Build-Script for building Psychtoolbox-3 "mex" files as Python extensions.
#
# (c) 2018 Mario Kleiner - Licensed under MIT license.
#

from distutils.core import setup, Extension # Build system.
import os, fnmatch, shutil                  # Directory traversal, file list building.
import platform                             # OS detection.
import sys                                  # cpu arch detection.
import numpy                                # To get include dir on macOS.

def get_sourcefiles(path):
    sources = [];
    pattern1 = '*.c';
    pattern2 = '*.cpp';
    for filename in os.listdir(path):
        if fnmatch.fnmatch(filename, pattern1) or fnmatch.fnmatch(filename, pattern2):
            sources += [os.path.join(path,filename)];

    # Fancy schmanzi, not needed atm. for recursive dir traversal:
    #    for root, dirs, files in os.walk(path):
    #        for filename in files:
    #            sources += [os.path.join(root,filename)];

    return(sources);

def get_basemacros(name, osname):
    return([('PTBMODULE_' + name, None), ('PTBMODULENAME', name)] + base_macros);

def get_baseincludedirs(name, osname):
    return(['PsychSourceGL/Source/Common/' + name] + baseincludes_common + ['PsychSourceGL/Source/' + osname + '/Base'] + ['PsychSourceGL/Source/' + osname + '/' + name]);

def get_basesources(name, osname):
    extrafiles = [];
    if os.access('./PsychSourceGL/Source/' + osname + '/' + name, os.F_OK):
         extrafiles = get_sourcefiles('./PsychSourceGL/Source/' + osname + '/' + name);

    return(basefiles_common + get_sourcefiles('./PsychSourceGL/Source/' + osname + '/Base/') + get_sourcefiles('./PsychSourceGL/Source/Common/' + name) + extrafiles);

# Treating some special cases like Octave seems to be the right thing to do,
# PSYCH_LANGUAGE setting is self-explanatory:
base_macros = [('PTBOCTAVE3MEX', None), ('PSYCH_LANGUAGE', 'PSYCH_PYTHON')];
# Disabled: This would enable Py_LIMITED_API, to allow to build one set of modules for all versions of
# Python >= 3.2. The downside is loss of important functionality [PsychRuntimeEvaluateString()) does not work!].
# Also, we have build failure on at least Ubuntu 18.04 LTS with Python 3.6, so it is a no-go on Linux for now!
#base_macros = [('PTBOCTAVE3MEX', None), ('PSYCH_LANGUAGE', 'PSYCH_PYTHON'), ('Py_LIMITED_API', None)];

# Common infrastructure and the scripting glue module for interfacing with the Python runtime:
basefiles_common = get_sourcefiles('./PsychSourceGL/Source/Common/Base') + ['PsychSourceGL/Source/Common/Base/PythonGlue/PsychScriptingGluePython.c'];
baseincludes_common = [numpy.get_include(), 'PsychSourceGL/Source/Common/Base', 'PsychSourceGL/Source/Common/Screen'];

is_64bits = sys.maxsize > 2**32;

# OS detection and file selection for the different OS specific backends:
print('Platform reported as: %s\n' % platform.system());
if platform.system() == 'Linux':
    # Linux specific backend code:
    print('Building for Linux...\n');
    osname = 'Linux';
    # All libraries to link to all modules:
    base_libs = ['c', 'rt', 'dl'];
    # No "no reproducible builds" warning:
    base_compile_args = ['-Wno-date-time'];
    # Extra OS specific libs for PsychPortAudio:
    audio_libdirs = [];
    audio_extralinkargs = [];
    audio_libs = ['portaudio', 'asound'];
    audio_objects = [];

    # Static linking for audio_objects aka libportaudio.a no longer used as of v3.0.15:
    #if is_64bits == True:
    #    audio_objects = ['PsychSourceGL/Cohorts/PortAudio/libportaudio64Linux.a'];
    #else:
    #    audio_objects = ['PsychSourceGL/Cohorts/PortAudio/libportaudio32Linux.a'];

    # libusb includes:
    usb_includes = ['/usr/include/libusb-1.0'];

    # Extra OS specific libs for PsychHID:
    psychhid_includes = usb_includes;
    psychhid_libdirs = [];
    psychhid_libs = ['dl', 'usb-1.0', 'X11', 'Xi', 'util'];
    psychhid_extra_objects = [];

    # Extra files needed, e.g., libraries:
    extra_files = {};

if platform.system() == 'Windows':
    print('Building for Windows...\n');
    osname = 'Windows';
    base_libs = ['kernel32', 'user32', 'advapi32', 'winmm'];
    base_compile_args = [];

    # Extra OS specific libs for PsychPortAudio:
    audio_libdirs = ['PsychSourceGL/Cohorts/PortAudio'];
    audio_extralinkargs = []; # No runtime delay loading atm. No benefit with current packaging method: ['/DELAYLOAD:portaudio_x64.dll'];
    audio_libs = ['delayimp', 'portaudio_x64'];
    audio_objects = [];

    # libusb includes:
    usb_includes = ['PsychSourceGL/Cohorts/libusb1-win32/include/libusb-1.0']

    # Extra OS specific libs for PsychHID:
    psychhid_includes = usb_includes;
    psychhid_libdirs = ['PsychSourceGL/Cohorts/libusb1-win32/MS64/dll'];
    psychhid_libs = ['dinput8', 'libusb-1.0', 'setupapi'];  
    psychhid_extra_objects = [];

    # Extra files needed, e.g., libraries.
    # The 64-Bit portaudio dll for PsychPortAudio and libusb1 dll for PsychHID:
    extra_files = {'psychtoolbox' : ['portaudio_x64.dll', 'libusb-1.0.dll']};

    # These extra files must be in the PsychPython folder that defines our
    # 'psychtoolbox' package. Copy them temporarily from their actual locations
    # in Psychtoolbox to that folder:
    shutil.copy('Psychtoolbox/PsychSound/portaudio_x64.dll', 'PsychPython/');
    shutil.copy('Psychtoolbox/PsychContributed/libusb-1.0.dll', 'PsychPython/');

if platform.system() == 'Darwin':
    print('Building for macOS...\n');
    osname = 'OSX';
    # These should go to extra_link_args in Extension() below, but apparently distutils
    # always appends the extra_link_args at the end of the linker command line, which is
    # wrong for -framework's, as they must be stated *before* the .o object files which
    # want to use functions from them. A solution to this problem doesn't exist in distutils
    # almost two decades after the debut of Mac OSX, because hey, take your time!
    #
    # The hack is to set the LDFLAGS environment variable to what we need, as LDFLAGS
    # apparently gets prepended to the linker invocation arguments, so -framework statements
    # precede the .o'bjects they should apply to and the linker is happy again.
    #
    # Downside is that now we have to pass the union of *all* -framework switches ever
    # used for *any* extension module, as os.environ can't be changed during the build
    # sequencing for a distribution package.
    #
    # Is this awful? Absolutely! And i thought the Octave/Matlab build process on macOS
    # sucked big time, but apparently somebody in Pythonland shouted "Hold my beer!" :(
    # Hopefully some new info about this issue will prove me wrong, and there is a sane
    # and elegant solution, but this is the best i could find after hours of Googling and
    # trying.
    #
    # The following would be the full list of frameworks, but apparently including -framework Carbon
    # is enough. Maybe a catch-all including all other frameworks?
    #
    # -framework Carbon -framework CoreServices -framework CoreFoundation -framework CoreAudio -framework AudioToolbox -framework AudioUnit
    # -framework ApplicationServices -framework OpenGL -framework CoreVideo -framework IOKit -framework SystemConfiguration
    # -framework CoreText -framework Cocoa
    #
    os.environ['LDFLAGS'] = '-framework Carbon -framework CoreAudio'
    base_libs = [];

    # No "no reproducible builds" warning:
    base_compile_args = ['-Wno-date-time', '-mmacosx-version-min=10.11'];

    # Extra OS specific libs for PsychPortAudio:
    audio_libdirs = [];
    audio_extralinkargs = [];
    audio_libs = [];
    # Include our statically linked on-steroids version of PortAudio:
    audio_objects = ['PsychSourceGL/Cohorts/PortAudio/libportaudio_osx_64.a'];

    # Include Apples open-source HID Utilities for all things USB-HID device handling:
    psychhid_includes = ['PsychSourceGL/Cohorts/HID_Utilities_64Bit/', 'PsychSourceGL/Cohorts/HID_Utilities_64Bit/IOHIDManager'];
    psychhid_libdirs = [];
    psychhid_libs = [];
    # Extra objects for PsychHID - statically linked HID utilities:
    psychhid_extra_objects = ['PsychSourceGL/Cohorts/HID_Utilities_64Bit/build/Release/libHID_Utilities64.a'];

    # Extra files needed, e.g., libraries:
    extra_files = {};

# GetSecs module: Clock queries.
name = 'GetSecs';
GetSecs = Extension(name,
                    extra_compile_args = base_compile_args,
                    define_macros = get_basemacros(name, osname),
                    include_dirs = get_baseincludedirs(name, osname),
                    sources = get_basesources(name, osname),
                    libraries = base_libs,
                   )

# WaitSecs module: Timed waits.
name = 'WaitSecs';
WaitSecs = Extension(name,
                     extra_compile_args = base_compile_args,
                     define_macros = get_basemacros(name, osname),
                     include_dirs = get_baseincludedirs(name, osname),
                     sources = get_basesources(name, osname),
                     libraries = base_libs
                    )

# PsychPortAudio module: High precision, high reliability, multi-channel, multi-card audio i/o.
name = 'PsychPortAudio';
PsychPortAudio = Extension(name,
                           extra_compile_args = base_compile_args,
                           define_macros = get_basemacros(name, osname),
                           include_dirs = get_baseincludedirs(name, osname),
                           sources = get_basesources(name, osname),
                           library_dirs = audio_libdirs,
                           libraries = base_libs + audio_libs,
                           extra_link_args = audio_extralinkargs,
                           extra_objects = audio_objects
                          )

# PsychHID module: Note the extra include_dirs and libraries:
name = 'PsychHID';
PsychHID = Extension(name,
                     extra_compile_args = base_compile_args,
                     define_macros = get_basemacros(name, osname),
                     include_dirs = get_baseincludedirs(name, osname) + psychhid_includes,
                     sources = get_basesources(name, osname),
                     library_dirs = psychhid_libdirs,
                     libraries = base_libs + psychhid_libs,
                     extra_objects = psychhid_extra_objects
                    )

# IOPort module:
name = 'IOPort';
IOPort = Extension(name,
                   extra_compile_args = base_compile_args,
                   define_macros = get_basemacros(name, osname),
                   include_dirs = get_baseincludedirs(name, osname),
                   sources = get_basesources(name, osname),
                   libraries = base_libs
                  )

setup (name = 'psychtoolbox',
       version = '3.0.15',
       description = 'Pieces of Psychtoolbox-3 ported to CPython.',
       author = 'Mario Kleiner',
       author_email = 'mario.kleiner.de@gmail.com',
       url = 'http://psychtoolbox.org',
       packages = ['psychtoolbox'],
       package_dir = {'psychtoolbox' : 'PsychPython'},
       package_data = extra_files,
       ext_package= 'psychtoolbox',
       ext_modules = [WaitSecs, GetSecs, IOPort, PsychHID, PsychPortAudio]
      )

if platform.system() == 'Windows':
    # Get rid of the now no longer needed copies of dll's inside PsychPython,
    # now that setup() has already copied them into the distribution.
    os.remove('PsychPython/portaudio_x64.dll');
    os.remove('PsychPython/libusb-1.0.dll');
