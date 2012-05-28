/*
	PsychToolbox3/Source/Common/Screen/PsychMovieSupportQuickTime.h
	
	PLATFORMS:	This is the version shared between OS/X and Windows.  
				
	AUTHORS:

        mario.kleiner@tuebingen.mpg.de   mk   Mario Kleiner

	HISTORY:

        10/31/05    mk      Wrote it.
	11/30/10    mk      Moved to Common/Screen, renamed functions to PsychQTxxx

	DESCRIPTION:
	
	Psychtoolbox functions for dealing with movies, using Apple's Quicktime
	SDK/Engine in the classic V6.x/7.x version.

*/

//include once
#ifndef PSYCH_IS_INCLUDED_PsychMovieSupportQuickTime
#define PSYCH_IS_INCLUDED_PsychMovieSupportQuickTime

#include "Screen.h"

#ifdef PSYCHQTAVAIL

#if PSYCH_SYSTEM == PSYCH_WINDOWS
#include <Movies.h>
#else
#include <Quicktime/Movies.h>
#endif

void PsychQTMovieInit(void);
int  PsychQTGetMovieCount(void);
void PsychQTCreateMovie(PsychWindowRecordType *win, const char* moviename, double preloadSecs, int* moviehandle, int specialFlags1);
void PsychQTGetMovieInfos(int moviehandle, int* width, int* height, int* framecount, double* durationsecs, double* framerate, int* nrdroppedframes);
void PsychQTDeleteMovie(int moviehandle);
void PsychQTDeleteAllMovies(void);
void PsychQTFreeMovieTexture(PsychWindowRecordType *win);
int PsychQTGetTextureFromMovie(PsychWindowRecordType *win, int moviehandle, int checkForImage, double timeindex, PsychWindowRecordType *out_texture, double *presentation_timestamp);
int PsychQTPlaybackRate(int moviehandle, double playbackrate, int loop, double soundvolume);
void PsychQTExitMovies(void);
double PsychQTGetMovieTimeIndex(int moviehandle);
double PsychQTSetMovieTimeIndex(int moviehandle, double timeindex, psych_bool indexIsFrames);

#endif

//end include once
#endif
