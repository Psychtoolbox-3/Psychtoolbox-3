/*
    PsychToolbox3/Source/Common/Screen/ScreenExit.cpp		

    AUTHORS:

        Allen.Ingling@nyu.edu           awi
        mario.kleiner@tuebingen.mpg.de  mk

    PLATFORMS:	

        All.

    HISTORY:

        12/20/01    awi     Created.
        1/25/04     awi     Added update provided by mk. It makes the ScreenCloseAllWindows call.  
        7/22/05     mk      Added call to CloseWindowBank() to free dynamic window bank array.

    DESCRIPTION:

        ScreenExitFunction is called before the Screen module is flushed.
*/

#include "Screen.h"

void PsychCleanupSCREENFillPoly(void);

PsychError ScreenExitFunction(void)
{
	//The timing array holds time values set by Screen internal diagnostics.  It allocates memory with 
	//malloc to hold the array of times.  This call frees the memory prior to unloading Screen
	ClearTimingArray();
  
	// Close all open onscreen windows and release their resources,
	// -> Perform exactly the same cleanup that Screen('CloseAll') would do.
	ScreenCloseAllWindows();
	CloseWindowBank();

    // Shutdown low-level display glue (Screens, displays, kernel-drivers et al.):
    PsychCleanupDisplayGlue();

	// Cleanup internal data structures of SCREEN('FillPoly');
	// This is defined in Common/Screen/SCREENFillPoly.c
	PsychCleanupSCREENFillPoly();

	// Release our internal locale object for character <-> unicode conversion:
	PsychSetUnicodeTextConversionLocale(NULL);

	return(PsychError_none);
}
