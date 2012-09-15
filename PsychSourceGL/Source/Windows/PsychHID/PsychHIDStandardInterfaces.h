/*
  PsychToolbox3/Source/Windows/PsychHID/PsychHIDStandardInterfaces.h
  
  PROJECTS: PsychHID only.
  
  PLATFORMS:  Windows.
  
  AUTHORS:
  mario.kleiner@tuebingen.mpg.de    mk

  HISTORY:
  27.07.2011     mk     Created.
  
  TO DO:

*/

#ifndef PSYCH_IS_INCLUDED_PsychHIDStandardIntefaces
#define PSYCH_IS_INCLUDED_PsychHIDStandardIntefaces

// CAUTION: The exact order of include and define statements below matters!
// Define direct input version explicitely to avoid compiler warnings.
// This sequence of InitGuid, version define and dinput.h include defines the
// needed GUID's so we do no longer need to link against dinput.lib or dxguid.lib,
// which both have been removed from the current Windows-7 platform SDK, dinput.lib,
// because it is deprecated since 2007 -- it implemented DirectInput-7, and dxguid.lib,
// because the InitGuid.h method is the preferred method on Windows now.
#include <InitGuid.h>
#define DIRECTINPUT_VERSION 0x0800
#include <dinput.h>

#include "PsychHID.h"

#endif
