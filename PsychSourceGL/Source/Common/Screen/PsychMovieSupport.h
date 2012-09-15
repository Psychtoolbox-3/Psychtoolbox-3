/*
	PsychToolbox3/Source/Common/Screen/PsychMovieSupport.h
	
	PLATFORMS:	All.  
				
	AUTHORS:

	Mario Kleiner	mk	mario.kleiner@tuebingen.mpg.de

	HISTORY:

        10/31/05    mk      Wrote it.

	DESCRIPTION:
	
	Psychtoolbox functions for dealing with movies.

*/

//include once
#ifndef PSYCH_IS_INCLUDED_PsychMovieSupport
#define PSYCH_IS_INCLUDED_PsychMovieSupport

#include "Screen.h"

typedef struct PsychAsyncMovieInfo {
    unsigned char asyncstate;
    char* moviename;
    PsychWindowRecordType windowRecord;
    int moviehandle;
    double preloadSecs;
	int asyncFlag;
	int specialFlags1;
    int pixelFormat;
    int maxNumberThreads;
    psych_thread pid;
} PsychAsyncMovieInfo;

void PsychMovieInit(void);
int PsychGetMovieCount(void);
void* PsychAsyncCreateMovie(void* inmovieinfo);
void PsychCreateMovie(PsychWindowRecordType *win, const char* moviename, double preloadSecs, int* moviehandle, int asyncFlag, int specialFlags1, int pixelFormat, int maxNumberThreads);
void PsychGetMovieInfos(int moviehandle, int* width, int* height, int* framecount, double* durationsecs, double* framerate, int* nrdroppedframes, double* aspectRatio);
void PsychDeleteMovie(int moviehandle);
void PsychDeleteAllMovies(void);
void PsychFreeMovieTexture(PsychWindowRecordType *win);
int PsychGetTextureFromMovie(PsychWindowRecordType *win, int moviehandle, int checkForImage, double timeindex, PsychWindowRecordType *out_texture, double *presentation_timestamp);
int PsychPlaybackRate(int moviehandle, double playbackrate, int loop, double soundvolume);
void PsychExitMovies(void);
double PsychGetMovieTimeIndex(int moviehandle);
double PsychSetMovieTimeIndex(int moviehandle, double timeindex, psych_bool indexIsFrames);

//end include once
#endif

