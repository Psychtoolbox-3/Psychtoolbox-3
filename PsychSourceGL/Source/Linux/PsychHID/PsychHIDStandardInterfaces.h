/*
  PsychToolbox3/Source/Linux/PsychHID/PsychHIDStandardInterfaces.h
  
  PROJECTS: PsychHID only.
  
  PLATFORMS:  Linux.
  
  AUTHORS:
  mario.kleiner@tuebingen.mpg.de    mk

  HISTORY:
  27.07.2011     mk     Created.
  
  TO DO:

*/

#ifndef PSYCH_IS_INCLUDED_PsychHIDStandardIntefaces
#define PSYCH_IS_INCLUDED_PsychHIDStandardIntefaces

// Needed for dlsym() RTLD_DEFAULT definition in PsychHIDStandardInterfaces.c
#define _GNU_SOURCE

#include "PsychHID.h"

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/extensions/XInput.h>
#include <X11/extensions/XInput2.h>

#endif
