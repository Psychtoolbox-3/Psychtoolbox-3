/*
	Common/Screen/PsychMovieWritingSupport.h
	
	PLATFORMS:	
	
		These are platform independent prototype definitions.
		
	AUTHORS:
	
		Mario Kleiner           mk              mario.kleiner@tuebingen.mpg.de

	HISTORY:
	
	DESCRIPTION:
	
		Psychtoolbox functions for dealing with movie editing. The implementation
		of the actual functions is elsewhere and platform dependent.
		
	AUTHORS:

        Mario Kleiner   mk  mario.kleiner at tuebingen.mpg.de
 
	HISTORY:
	04/18/10		mk		Wrote it. 
							
	DESCRIPTION:
	
*/

//include once
#ifndef PSYCH_IS_INCLUDED_PsychMovieWritingSupport
#define PSYCH_IS_INCLUDED_PsychMovieWritingSupport

// Maximum number of allowed simultaneously open capture devices. As structures are allocated
// statically, this should not be increased beyond reasonable limits. 100 is always a nice number:
#define PSYCH_MAX_MOVIEWRITERDEVICES 10

#include "Screen.h"

// These are the generic entry points, to be called by SCREENxxxx videocapture functions and
// other parts of screen. They dispatch into API specific versions, depending on users choice
// of capture system and support by OS:
void PsychMovieWritingInit(void);
void PsychExitMovieWriting(void);
void PsychDeleteAllMovieWriters(void);
int PsychCreateNewMovieFile(char* moviefile, int width, int height, double framerate, char* movieoptions);
int PsychFinalizeNewMovieFile(int movieHandle);
int PsychAddVideoFrameToMovie(int moviehandle, int frameDurationUnits, psych_bool isUpsideDown);
unsigned char*	PsychGetVideoFrameForMoviePtr(int moviehandle, unsigned int* twidth, unsigned int* theight);
psych_bool PsychAddAudioBufferToMovie(int moviehandle, unsigned int nrChannels, unsigned int nrSamples, double* buffer);

//end include once
#endif
