/*

 Psychtoolbox3/Source/Common/SCREENPlayMovie.c		

 AUTHORS:

 mario.kleiner at tuebingen.mpg.de   mk

 PLATFORMS:	

 This file should build on any platform. 

 HISTORY:

 10/23/05  mk		Created. 

 DESCRIPTION:

 Start/Stop playback of given movie object or change its sound volume.

 TO DO:

 */

#include "Screen.h"

static char useString[] = "[droppedframes] = Screen('PlayMovie', moviePtr, rate, [loop], [soundvolume]);";
static char synopsisString[] = 
"Start playback of movie associated with movieobject 'moviePtr'. 'rate' defines the playback rate: 0 == Stop playback, "
"1 == Normal speed forward, -1 == Normal speed backward, ... . 'loop' Enable looped playback. 'soundvolume' Output volume "
"of associated soundtrack: 0 = Mute sound output, 0.01 - 1.0 Volume in percent. Returns number of frames "
"that needed to be dropped in order to keep video playback in sync with realtime and audio playback.\n";

static char seeAlsoString[] = "CloseMovie PlayMovie GetMovieImage GetMovieTimeIndex SetMovieTimeIndex";	 

PsychError SCREENPlayMovie(void) 

{
    int                                 moviehandle = -1;
    double                              rate = 0;
    int                                 loop = 0;
    double                              sndvolume = 1;
    double                              dropped = 0;
    
    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

    PsychErrorExit(PsychCapNumInputArgs(4));            // Max. 4 input args.
    PsychErrorExit(PsychRequireNumInputArgs(2));        // Min. 2 input args required.
    PsychErrorExit(PsychCapNumOutputArgs(1));           // Max. 1 output args.

    // Get the movie handle:
    PsychCopyInIntegerArg(1, TRUE, &moviehandle);
    if (moviehandle==-1) {
        PsychErrorExitMsg(PsychError_user, "PlayMovie called without valid handle to a movie object.");
    }

    // Get the requested playback rate.
    PsychCopyInDoubleArg(2, TRUE, &rate);

    // Get the 'loop' flag: If zero, we don't loop, otherwise we loop.
    PsychCopyInIntegerArg(3, FALSE, &loop);

    // Get the requested sound volume. Defaults to 1 == Full blast!
    PsychCopyInDoubleArg(4, FALSE, &sndvolume);
    if (sndvolume<0) sndvolume = 0;
    if (sndvolume>1) sndvolume = 1;

    // Start playback via low-level routines:
    dropped = (double) PsychPlaybackRate(moviehandle, rate, loop, sndvolume);

    // Return optional count of dropped frames:
    PsychCopyOutDoubleArg(1, FALSE, dropped);

    // Ready!
    return(PsychError_none);
}
