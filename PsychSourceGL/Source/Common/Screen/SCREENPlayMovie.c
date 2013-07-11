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
"Start playback of movie associated with movieobject 'moviePtr'. 'rate' defines the desired playback rate: 0 == Stop playback, "
"1 == Normal speed forward, -1 == Normal speed backward, ... . Not all movie files allow reverse playback or playback at other "
"rates than normal speed forward. 'loop' Enable looped playback if set to a value greater than zero. A value of 1 will enable "
"repetitive looped playback with the default strategy, or the special strategy selected via the optional 'specialFlags1' settings "
"provided during the call to Screen('OpenMovie', ...). Values greater than 1 select different repetition strategies when combined together. "
"Different strategies exist to handle different quirks with some movie file formats and encodings and with some versions of GStreamer: "
"A flag of 2 requests looped playback via gapless reloading of the movie instead of rewinding it to the start. A flag of 1+4 uses "
"so called segment seeks for rewinding, whereas adding 8 (e.g., 1+8 or 1+4+8) asks to flush the video pipeline during rewinding. "
"Your mileage with these looping strategies will differ, but usually the default settings are good enough for most purposes. "
"There are however movie files or network video streams that can't be automatically repeated at all, so called non-seekable streams.\n"
"'soundvolume' Select the output audio volume of an associated soundtrack: 0 = Mute sound output, 0.01 - 1.0 Volume in percent. "
"You can choose the sound volume with low overhead while playback is active by calling this function, as long as you provide the "
"same parameters for all other settings as the ones you used when starting playback, ie. only 'soundvolume' may differ.\n"
"If the function is called to stop playback, it will return the number of frames that needed to be dropped in order to keep "
"video playback in sync with realtime and audio playback. Otherwise it returns zero.\n";

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
