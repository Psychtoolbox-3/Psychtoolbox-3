/*
 
 Psychtoolbox3/Source/Common/SCREENGetMovieTimeIndex.c		
 
 AUTHORS:
 
 mario.kleiner at tuebingen.mpg.de   mk
 
 PLATFORMS:	
 
 This file should build on any platform. 
 
 HISTORY:
 
 10/23/05  mk		Created. 
 
 DESCRIPTION:
 
 TO DO:
 
 */

#include "Screen.h"

static char useString[] = "timeindex = Screen('GetMovieTimeIndex', moviePtr);";
static char synopsisString[] = "Return current time index for movie object 'moviePtr'.";
static char seeAlsoString[] = "CloseMovie PlayMovie GetMovieImage GetMovieTimeIndex SetMovieTimeIndex";

PsychError SCREENGetMovieTimeIndex(void) 
{
    int moviehandle = -1;
    
    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};
    
    PsychErrorExit(PsychCapNumInputArgs(1));            // Max. 1 input args.
    PsychErrorExit(PsychRequireNumInputArgs(1));        // Min. 1 input args required.
    PsychErrorExit(PsychCapNumOutputArgs(1));           // One output arg.
    
    // Get the movie handle:
    PsychCopyInIntegerArg(1, TRUE, &moviehandle);
    if (moviehandle==-1) {
        PsychErrorExitMsg(PsychError_user, "GetMovieTimeIndex called without valid handle to a movie object.");
    }
    
    // Retrieve and return current movie time index:
    PsychCopyOutDoubleArg(1, TRUE, PsychGetMovieTimeIndex(moviehandle));
    
    return(PsychError_none);
}
