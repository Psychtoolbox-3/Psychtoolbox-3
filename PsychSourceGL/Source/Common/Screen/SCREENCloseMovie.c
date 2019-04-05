/*

    Psychtoolbox3/Source/Common/SCREENCloseMovie.c

    AUTHORS:

    mario.kleiner.de@gmail.com      mk

    PLATFORMS:

    All.

    HISTORY:

    10/23/05  mk        Created.

    DESCRIPTION:

    Close a previously opened movie file and release all associated ressources.

*/

#include "Screen.h"

static char useString[] = "Screen('CloseMovie' [, moviePtr=all]);";
static char synopsisString[] = "Close movie object specified by 'moviePtr' and release all associated ressources.\n"
                               "If 'moviePtr' is omitted then all open movies are closed.\n";
static char seeAlsoString[] = "PlayMovie GetMovieImage GetMovieTimeIndex SetMovieTimeIndex";

PsychError SCREENCloseMovie(void)
{
    int moviehandle;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if (PsychIsGiveHelp()) { PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(1));            // Max. 1 input args.
    PsychErrorExit(PsychRequireNumInputArgs(0));        // Min. 0 input args required.
    PsychErrorExit(PsychCapNumOutputArgs(0));           // No output args.

    // Get the movie handle:
    if (!PsychCopyInIntegerArg(1, kPsychArgOptional, &moviehandle))
        PsychDeleteAllMovies();
    else
        PsychDeleteMovie(moviehandle);

    return(PsychError_none);
}
