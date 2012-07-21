/*
	PsychSourceGL/Source/Common/Screen/PsychMovieSupportGStreamer.h
	
	PLATFORMS:	All with PTB_USE_GSTREAMER defined.
				
	AUTHORS:

	mario.kleiner@tuebingen.mpg.de		mk	Mario Kleiner

	HISTORY:

        28.11.2010    mk      Wrote it.

	DESCRIPTION:
	
	Psychtoolbox functions for dealing with movies.

*/

#ifdef PTB_USE_GSTREAMER

//include once
#ifndef PSYCH_IS_INCLUDED_PsychMovieSupportGStreamer
#define PSYCH_IS_INCLUDED_PsychMovieSupportGStreamer

#include "Screen.h"

void PsychGSMovieInit(void);
int  PsychGSGetMovieCount(void);
void PsychGSCreateMovie(PsychWindowRecordType *win, const char* moviename, double preloadSecs, int* moviehandle, int asyncFlag, int specialFlags1, int pixelFormat);
void PsychGSGetMovieInfos(int moviehandle, int* width, int* height, int* framecount, double* durationsecs, double* framerate, int* nrdroppedframes, double* aspectRatio);
void PsychGSDeleteMovie(int moviehandle);
void PsychGSDeleteAllMovies(void);
void PsychGSFreeMovieTexture(PsychWindowRecordType *win);
int PsychGSGetTextureFromMovie(PsychWindowRecordType *win, int moviehandle, int checkForImage, double timeindex, PsychWindowRecordType *out_texture, double *presentation_timestamp);
int PsychGSPlaybackRate(int moviehandle, double playbackrate, int loop, double soundvolume);
void PsychGSExitMovies(void);
double PsychGSGetMovieTimeIndex(int moviehandle);
double PsychGSSetMovieTimeIndex(int moviehandle, double timeindex, psych_bool indexIsFrames);

//end include once
#endif

#endif
