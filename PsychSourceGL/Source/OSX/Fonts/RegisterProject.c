/*
	PsychToolbox2/Source/OSX/FontInfo/RegisterProject.c		

	PROJECTS: 
	
		FontInfo only.  

	AUTHORS:
	
		Allen.Ingling@nyu.edu           awi
        mario.kleiner@tuebingen.mpg.de  mk

	PLATFORMS:	
	
		Mac OS X

	HISTORY:
	
		10/23/03  	awi		Created.
		3/7/06		awi		Changed references from "Font" to "FontInfo".  The function had been previously renamed, but not all references updated.

	DESCRIPTION:
		

	NOTES:

		FontInfo uses PsychFontGlue to enumerate fonts and store information in a native font list.  However, PsychFontGlue does not share the same
		list between modules.  So the font list maintained by Screen and that maintained by FontInfo will in fact be different lists, though 
		they should be identical lists and immutable.  This is wasteful but it should work until we consolidate the two by sharing
		memory between Psychtoolbox mex files.

  */ 

#include "RegisterProject.h"

PsychError PsychModuleInit(void)
{
	//register the project exit function
	PsychErrorExitMsg(PsychRegisterExit(&FontsExitFunction), NULL); 
	
	//register the project function which is called when the module
	//is invoked with no arguments:
	PsychErrorExitMsg(PsychRegister(NULL,  &PsychDisplayATSSynopsis), NULL);

	//register the module name
	PsychErrorExitMsg(PsychRegister("FontInfo", NULL),NULL);
	
	//register named subfunctions
	PsychErrorExit(PsychRegister("Fonts",  &FONTSFonts));
	PsychErrorExit(PsychRegister("NumFonts",  &FONTSNumFonts));
	PsychErrorExit(PsychRegister("Version",  &MODULEVersion));
        
	InitializeSynopsis();
	PsychSetModuleAuthorByInitials("awi");
    PsychSetModuleAuthorByInitials("mk");

	return(PsychError_none);
}
