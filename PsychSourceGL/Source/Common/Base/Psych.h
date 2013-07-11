/*
  PsychToolbox2/Source/Common/Psych.h		
  
  AUTHORS:
  Allen.Ingling@nyu.edu		awi 
  
  PLATFORMS: All OS
  
  PROJECTS:
  07/17/02	awi		Screen on MacOS9
   

  HISTORY:
  07/17/02  awi		Wrote it.  
  
  DESCRIPTION:
  
  T0 DO: 
	  

*/

//begin include once 
#ifndef PSYCH_IS_INCLUDED_Psych
#define PSYCH_IS_INCLUDED_Psych

#define __USE_UNIX98 1

//can include this alone
#include "PsychConstants.h"

#if defined(__cplusplus)
extern "C" {
#endif

//if using any of below include Psych.h
#include "MiniBox.h"
#include "ProjectTable.h"
#include "PsychError.h"
#include "PsychScriptingGlue.h"
#include "PsychStructGlue.h"
#include "PsychCellGlue.h"
#include "PsychRegisterProject.h"
#include "PsychAuthors.h"
#include "PsychVersioning.h"
#include "PsychHelp.h"
#include "PsychInit.h"
#include "PsychMemory.h"
#include "PsychTimeGlue.h"
#include "PsychInstrument.h"	

#ifndef PTBINSCRIPTINGGLUE
// This is provided by the project. We do not
// include it when building PsychScriptingGlue.cc
// in order to prevent build-failure due to some
// broken header dependencies...
#pragma once
#include "RegisterProject.h"
#include "TimeLists.h"
#endif

#if defined(__cplusplus)
}
#endif

//end include once
#endif

