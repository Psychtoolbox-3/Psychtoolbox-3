/*
	PsychToolbox3/Source/Common/Screen/PsychTextureSupport.h
	
	PLATFORMS:	This is the OS independent version.  
				
	AUTHORS:
	Allen Ingling		awi     Allen.Ingling@nyu.edu
	Mario Kleiner		mk      mario.kleiner@tuebingen.mpg.de

    HISTORY:

    3/9/04		awi		Wrote it 
							
	DESCRIPTION:
	
	Psychtoolbox functions for dealing with textures.

*/

//include once
#ifndef PSYCH_IS_INCLUDED_PsychTextureSupport
#define PSYCH_IS_INCLUDED_PsychTextureSupport

#include "Screen.h"

void PsychInitWindowRecordTextureFields(PsychWindowRecordType *winRec);
void PsychCreateTexture(PsychWindowRecordType *win);
void PsychFreeTextureForWindowRecord(PsychWindowRecordType *win);
void PsychBlitTextureToDisplay(PsychWindowRecordType *source, PsychWindowRecordType *target, double *sourceRect, double *targetRect,
                               double rotationAngle, int filterMode, double globalAlpha);
GLenum PsychGetTextureTarget(PsychWindowRecordType *win);
void PsychMapTexCoord(PsychWindowRecordType *tex, double* tx, double* ty);
void PsychDetectTextureTarget(PsychWindowRecordType *win);

//end include once
#endif
