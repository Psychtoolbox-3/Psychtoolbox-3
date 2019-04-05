# __init__.py
#
# Import Psychtoolbox-3 "mex files" for use in Python and make them available
# under the same names and syntax as people are used to from Octave and Matlab.
#
# Usage in Python script/module:
#
#           from psychtoolbox import *
#
# Copyright (c) 2018 Mario Kleiner. Licensed under MIT license.
#

# Import all PTB modules already ported to Python extension modules:
from . import WaitSecs, GetSecs, PsychPortAudio, PsychHID, IOPort

# Extract their "main" function, which is named just like the module,
# and assign its function handle to a variable named just like the module.
# This will cause the function handles to shadow the original module, so
# now the module can be called with (almost) identical syntax as one is
# used to from Octave or Matlab:
WaitSecs = getattr(WaitSecs, 'WaitSecs');
GetSecs = getattr(GetSecs, 'GetSecs');
PsychPortAudio = getattr(PsychPortAudio, 'PsychPortAudio');
PsychHID = getattr(PsychHID, 'PsychHID');
IOPort = getattr(IOPort, 'IOPort');
