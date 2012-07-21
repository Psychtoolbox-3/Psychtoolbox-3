/*
  PsychToolbox2/Source/Common/WindowHelpers.h		
  
  AUTHORS:
  Allen.Ingling@nyu.edu		awi 
  
  PLATFORMS: This file is included for all targets.
    
  PROJECTS:
  12/20/02	awi		Screen
   

  HISTORY:
  12/20/02  	awi		Wrote it. 
  
  DESCRIPTION:
  
  Convenience functions for extracting properties of windows from a window index or from a window record pointer.
  Sometimes doing that is a hassle because the properties of the window are actually properties of the associeate screen
  structure, hence we provide these conveience functions.
  
   
*/

//begin include once 
#ifndef PSYCH_IS_INCLUDED_WindowHelpers
#define PSYCH_IS_INCLUDED_WindowHelpers

int					PsychGetWindowDepthValueFromWindowRecord(PsychWindowRecordType *windowRecord);
void				PsychSetTextColorInWindowRecord(PsychColorType *textColor,  PsychWindowRecordType *winRec);
void				PsychSetTextBackgroundColorInWindowRecord(PsychColorType *textColor,  PsychWindowRecordType *winRec);
int					PsychGetDepthFromWindowRecord(PsychWindowRecordType *winRec);
void				PsychGetRectFromWindowRecord(double *rect, PsychWindowRecordType *windowRecord);
int					PsychGetNumBuffersFromWindowRecord(PsychWindowRecordType *windowRecord);
int					PsychGetNumPlanesFromWindowRecord(PsychWindowRecordType *windowRecord);
void				PsychSetWindowRecordValid(PsychWindowRecordType *winRec);
void				PsychCheckIfWindowRecordIsValid(PsychWindowRecordType *winRec);
psych_bool			PsychIsOnscreenWindow(PsychWindowRecordType *windowRecord);
psych_bool			PsychIsOffscreenWindow(PsychWindowRecordType *windowRecord);
psych_bool			PsychIsTexture(PsychWindowRecordType *windowRecord);
psych_bool			PsychIsMatlabOnscreenWindow(PsychWindowRecordType *windowRecord);

//end include once
#endif
