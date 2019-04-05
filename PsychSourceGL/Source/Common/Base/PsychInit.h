/*
  PsychToolbox2/Source/Common/PsychInit.h		
  
  AUTHORS:
  Allen.Ingling@nyu.edu		awi 
  
  PLATFORMS: All
  
  PROJECTS:
  08/25/02	awi		Screen on MacOS9
   

  HISTORY:
  08/25/02  awi		wrote it.  
  
  DESCRIPTION:
  
	Holds the master initialization for the Psychtoolbox function.
	sections hold their own inits which are called by the master
	init function contained herein.  After the Psychtoolbox library
	has initialzed itself in then invokes the project init which 
	must be named PsychModuleInit(). 	 
   
  
	
*/
//begin include once 
#ifndef PSYCH_IS_INCLUDED_PsychInit
#define PSYCH_IS_INCLUDED_PsychInit

#include "Psych.h"

PsychError PsychInit(void);
PsychError PsychExit(void);

//end include once
#endif

