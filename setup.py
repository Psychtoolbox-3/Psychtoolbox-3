# setup.py -- Build-Script for building Psychtoolbox-3 "mex" files as Python extensions.
#
# (c) 2018-2023 Mario Kleiner
# (c) 2020-2022 Alex Forrence
# (c) 2019      Jonathan Peirce
# (c) 2019      Richard Höchenberger
# (c) 2019      Eric Larson
# (c) 2019      Bernhard M. Wiedemann
#
# Licensed under the MIT license.
#

# To build Python extensions/wheels, you'll want pip >= 10.0 for pyproject.toml support.
# This file will allow pip to automatically download the required build dependencies (i.e. numpy).
# To build this wheel locally into the dist/ directory, you could use the following
# (assuming the current working directory is next to the setup.py file):
#
# > pip wheel . -w dist
#
# (see https://pip.pypa.io/en/stable/cli/pip_wheel/ for all options)
#
# To build and install in the same step:
#
# > pip install .
#
# (see https://pip.pypa.io/en/stable/cli/pip_install/ for all options)
#
# If you wish to build wheels for distribution, you may want to compile with Python's Limited API
# enabled. The resulting wheel can be used across multiple Python releases (see
# https://docs.python.org/3/c-api/stable.html for details).
# This requires Python >= 3.7, and PTB_LIMITED_WHEEL must be set to 1. For example,
#
# > PTB_LIMITED_WHEEL=1 pip wheel . -w dist
#
# Will produce a wheel that will work for Python versions >= 3.7.
# It is recommended to compile these on the oldest supported Python version (i.e. 3.7).

# from distutils.core import setup, Extension # Build system.
from setuptools import setup, Extension, find_packages
import os, fnmatch, shutil                  # Directory traversal, file list building.
import platform                             # OS detection.
import sys                                  # cpu arch detection.
import numpy                                # To get include dir on macOS.

from wheel.bdist_wheel import bdist_wheel

# We borrowed the custom wheel configuration from https://github.com/joerick/python-abi3-package-sample
# It lives in a separate file for licensing simplicity
sys.path.append(os.path.join(os.path.dirname(__file__), 'PsychPython', 'psychtoolbox'))
from _abi3_wheel import abi3_wheel

is_64bits = sys.maxsize > 2**32

# unified version number, read from simple text file
def get_version():
    import re
    VERSIONFILE = "PsychPython/psychtoolbox/_version.py"
    with open(VERSIONFILE, "rt") as fid:
        verstrline = fid.read()
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    mo = re.search(VSRE, verstrline, re.M)
    if mo:
        verstr = mo.group(1)
    else:
        raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))
    return verstr
version = get_version()

def get_sourcefiles(path):
    sources = []
    pattern1 = '*.c'
    pattern2 = '*.cpp'
    for filename in sorted(os.listdir(path)):
        if fnmatch.fnmatch(filename, pattern1) or fnmatch.fnmatch(filename, pattern2):
            sources += [os.path.join(path,filename)]

    # Fancy schmanzi, not needed atm. for recursive dir traversal:
    #    for root, dirs, files in os.walk(path):
    #        dirs.sort()
    #        for filename in sorted(files):
    #            sources += [os.path.join(root,filename)]

    return(sources)

def get_basemacros(name, osname):
    return([('PTBMODULE_' + name, None), ('PTBMODULENAME', name)] + base_macros)

def get_baseincludedirs(name, osname):
    return(['PsychSourceGL/Source/Common/' + name] + baseincludes_common + ['PsychSourceGL/Source/' + osname + '/Base'] + ['PsychSourceGL/Source/' + osname + '/' + name])

def get_basesources(name, osname):
    extrafiles = []
    if os.access('./PsychSourceGL/Source/' + osname + '/' + name, os.F_OK):
         extrafiles = get_sourcefiles('./PsychSourceGL/Source/' + osname + '/' + name)

    return(basefiles_common + get_sourcefiles('./PsychSourceGL/Source/' + osname + '/Base/') + get_sourcefiles('./PsychSourceGL/Source/Common/' + name) + extrafiles)

# Treating some special cases like Octave seems to be the right thing to do,
# PSYCH_LANGUAGE setting is self-explanatory:
base_macros = [('PTBOCTAVE3MEX', None), ('PSYCH_LANGUAGE', 'PSYCH_PYTHON')]

# if desired by wheel-builder and building on py >= 3.7, build as a limited API wheel.
if sys.version_info >= (3, 7) and os.environ.get('PTB_LIMITED_WHEEL') == '1':
    base_macros.append(('Py_LIMITED_API', '0x03070000'))
    py_limited_api = True
    bdist = abi3_wheel
else:
    py_limited_api = False
    bdist = bdist_wheel

# Common infrastructure and the scripting glue module for interfacing with the Python runtime:
basefiles_common = get_sourcefiles('./PsychSourceGL/Source/Common/Base') + ['PsychSourceGL/Source/Common/Base/PythonGlue/PsychScriptingGluePython.c']
baseincludes_common = [numpy.get_include(), 'PsychSourceGL/Source/Common/Base', 'PsychSourceGL/Source/Common/Screen']

# OS detection and file selection for the different OS specific backends:
print('Platform reported as: %s\n' % platform.system())
if platform.system() == 'Linux':
    # Linux specific backend code:
    print('Building for Linux...\n')
    osname = 'Linux'
    # All libraries to link to all modules:
    base_libs = ['c', 'rt', 'dl']
    # No "no reproducible builds" warning:
    base_compile_args = ['-Wno-date-time']
    # Extra OS specific libs for PsychPortAudio:
    audio_libdirs = []
    audio_extralinkargs = []
    audio_libs = ['portaudio', 'asound']
    audio_objects = []

    # Extra OS specific libs for PsychHID:
    psychhid_includes = ['/usr/include/libusb-1.0']
    psychhid_libdirs = []
    psychhid_libs = ['dl', 'usb-1.0', 'X11', 'Xi', 'util']
    psychhid_extra_objects = []
    psychhid_extralinkargs = []

    # Extra files needed, e.g., libraries:
    extra_files = {}

if platform.system() == 'Windows':
    print('Building for Windows...\n')
    osname = 'Windows'
    base_libs = ['kernel32', 'user32', 'advapi32', 'winmm']
    base_compile_args = []

    # libusb includes:
    psychhid_includes = ['PsychSourceGL/Cohorts/libusb1-win32/include/libusb-1.0']
    psychhid_extralinkargs = []

    if is_64bits:
        # 64bit supports PsychPortAudio
        psychhid_libdirs = ['PsychSourceGL/Cohorts/libusb1-win32/MS64/dll']
        # and copy the files to the folder
        shutil.copy('PsychSourceGL/Cohorts/libusb1-win32/MS64/dll/libusb-1.0.dll',
                    'PsychPython/psychtoolbox/')
        shutil.copy('Psychtoolbox/PsychSound/portaudio_x64.dll',
                    'PsychPython/psychtoolbox/')
        # list them so they get packaged
        extra_files = {'psychtoolbox': ['portaudio_x64.dll', 'libusb-1.0.dll']}

        # Extra OS specific libs for PsychPortAudio:
        audio_libdirs = ['PsychSourceGL/Cohorts/PortAudio']
        audio_extralinkargs = [] # No runtime delay loading atm. No benefit with current packaging method: ['/DELAYLOAD:portaudio_x64.dll']
        audio_libs = ['delayimp', 'portaudio_x64']
        audio_objects = []
    else:
        # for win32 we use a different libusb dll and the portaudio dll is not supported
        psychhid_libdirs = ['PsychSourceGL/Cohorts/libusb1-win32/MS32/dll']
        shutil.copy('PsychSourceGL/Cohorts/libusb1-win32/MS32/dll/libusb-1.0.dll',
                    'PsychPython/psychtoolbox/')
        # list them so they get packaged
        extra_files = {'psychtoolbox': ['libusb-1.0.dll']}

    psychhid_libs = ['dinput8', 'libusb-1.0', 'setupapi']
    psychhid_extra_objects = []

if platform.system() == 'Darwin':
    machine = platform.machine()
    print(f'Building for macOS {machine}...\n')
    osname = 'OSX'
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
    base_libs = []

    # No "no reproducible builds" warning. macOS minimum version 10.9 selected by Jon Peirce.
    # May work for current modules, but is completely untested and unsupported by Psychtoolbox
    # upstream as of v3.0.15+, which only allows 10.11 as minimum version for Psychtoolbox mex
    # files and only tests with 10.13. Also note that already 10.11 is unsupported by Apple and
    # therefore a security risk.
    base_compile_args = ['-Wno-date-time', '-mmacosx-version-min=10.9']
    # Assuming the compiler is clang, this is how to turn warnings into errors
    base_link_args = ['-Wl,-fatal_warnings']

    # Extra OS specific libs for PsychPortAudio:
    audio_libdirs = []
    audio_extralinkargs = base_link_args
    audio_libs = []
    # Include our statically linked on-steroids version of PortAudio:
    # TODO: This needs to depend on intel vs arm64!
    audio_objects = ['PsychSourceGL/Cohorts/PortAudio/libportaudio_osx_64.a']

    # Include Apples open-source HID Utilities for all things USB-HID device handling, also libusb-1.0:
    psychhid_includes = ['PsychSourceGL/Cohorts/HID_Utilities_64Bit/', 'PsychSourceGL/Cohorts/HID_Utilities_64Bit/IOHIDManager',
                         'PsychSourceGL/Cohorts/libusb1-win32/include/libusb-1.0']
    psychhid_libdirs = []
    psychhid_libs = []
    # Weak-link libusb-1.0.dylib, so user does not need it on their local system as long as they
    # don't use PsychHID functions like USBControlTransfer/USBInterruptTransfer/USBBulkTransfer:
    psychhid_extralinkargs = ['-weak_library', 'PsychSourceGL/Cohorts/libusb1-win32/libusb-1.0.dylib'] + base_link_args

    # Extra objects for PsychHID - statically linked HID utilities:
    # TODO: This needs to depend on Intel vs arm64!
    psychhid_extra_objects = ['PsychSourceGL/Cohorts/HID_Utilities_64Bit/build/Release/libHID_Utilities64.a']

    # Extra files needed, e.g., libraries:
    extra_files = {}

ext_modules = []
# GetSecs module: Clock queries.
name = 'GetSecs'
GetSecs = Extension(name,
                    extra_compile_args = base_compile_args,
                    extra_link_args = base_link_args,
                    define_macros = get_basemacros(name, osname),
                    include_dirs = get_baseincludedirs(name, osname),
                    sources = get_basesources(name, osname),
                    libraries = base_libs,
                    py_limited_api = py_limited_api,
                   )
ext_modules.append(GetSecs)

# WaitSecs module: Timed waits.
name = 'WaitSecs'
WaitSecs = Extension(name,
                     extra_compile_args = base_compile_args,
                     define_macros = get_basemacros(name, osname),
                     include_dirs = get_baseincludedirs(name, osname),
                     sources = get_basesources(name, osname),
                     libraries = base_libs,
                     extra_link_args = base_link_args,
                     py_limited_api = py_limited_api,
                     )
ext_modules.append(WaitSecs)

# PsychPortAudio module: High precision, high reliability, multi-channel, multi-card audio i/o.
if is_64bits or platform.system() == 'Linux':
    # This won't compile on 32bit windows or macOS. Linux also has 32 Bit non-Intel variants, e.g., RaspberryPi
    name = 'PsychPortAudio'
    PsychPortAudio = Extension(name,
                               extra_compile_args = base_compile_args,
                               define_macros = get_basemacros(name, osname),
                               include_dirs = get_baseincludedirs(name, osname),
                               sources = get_basesources(name, osname),
                               library_dirs = audio_libdirs,
                               libraries = base_libs + audio_libs,
                               extra_link_args = audio_extralinkargs,
                               extra_objects = audio_objects,
                               py_limited_api = py_limited_api,
                               )
    ext_modules.append(PsychPortAudio)

# PsychHID module: Note the extra include_dirs and libraries:
name = 'PsychHID'
PsychHID = Extension(name,
                     extra_compile_args = base_compile_args,
                     define_macros = get_basemacros(name, osname),
                     include_dirs = get_baseincludedirs(name, osname) + psychhid_includes,
                     sources = get_basesources(name, osname),
                     library_dirs = psychhid_libdirs,
                     libraries = base_libs + psychhid_libs,
                     extra_link_args = psychhid_extralinkargs,
                     extra_objects = psychhid_extra_objects,
                     py_limited_api = py_limited_api,
                    )
ext_modules.append(PsychHID)

# IOPort module:
name = 'IOPort'
IOPort = Extension(name,
                   extra_compile_args = base_compile_args,
                   define_macros = get_basemacros(name, osname),
                   include_dirs = get_baseincludedirs(name, osname),
                   sources = get_basesources(name, osname),
                   libraries = base_libs,
                   extra_link_args = base_link_args,
                   py_limited_api = py_limited_api,
                  )
ext_modules.append(IOPort)

setup (name = 'psychtoolbox',
       version = version,
       description = 'Pieces of Psychtoolbox-3 ported to CPython.',
       author = 'Mario Kleiner',
       author_email = 'mario.kleiner.de@gmail.com',
       url = 'http://psychtoolbox.org',
       packages = ['psychtoolbox', 'psychtoolbox.demos'],
       package_dir = {'' : 'PsychPython',
                      'psychtoolbox' : 'PsychPython/psychtoolbox',
                      'psychtoolbox.demos' : 'PsychPython/demos'},
       package_data = extra_files,
       ext_package = 'psychtoolbox',
       ext_modules = ext_modules,
       include_package_data=True,  # Include files listed in MANIFEST.in
       install_requires = ['numpy>=1.13.3'], # Oldest supported numpy on Python 3.6
       cmdclass={"bdist_wheel": bdist},
      )

if platform.system() == 'Windows':
    # Get rid of the now no longer needed copies of dll's inside PsychPython,
    # now that setup() has already copied them into the distribution.
    if os.path.exists('PsychPython/psychtoolbox/portaudio_x64.dll'):
        os.remove('PsychPython/psychtoolbox/portaudio_x64.dll')
    os.remove('PsychPython/psychtoolbox/libusb-1.0.dll')
