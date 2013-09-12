/*
 
 PsychToolbox3/Source/OSX/Screen/PsychCocoaGlue.h
 
 PLATFORMS:	
 
 OSX
 
 TARGETS:
 
 Screen
 
 AUTHORS:
 
 Mario Kleiner       mk      mario.kleiner@tuebingen.mpg.de
 
 DESCRIPTION:
 
 Glue code for window management, using Objective-C wrappers to use Cocoa.
 
 These functions are called by PsychWindowGlue.c.
 
*/

#ifndef PSYCH_COCOAGLUE_H
#define PSYCH_COCOAGLUE_H

// Create windowed onscreen window via Cocoa:
PsychError PsychCocoaCreateWindow(PsychWindowRecordType *windowRecord, int windowLevel, void** outWindow);
void PsychCocoaDisposeWindow(PsychWindowRecordType *windowRecord);
void PsychCocoaSetUserFocusWindow(void* window);
void* PsychCocoaGetUserFocusWindow(void);
void PsychCocoaGetWindowBounds(void* window, PsychRectType globalBounds);
void PsychCocoaShowWindow(void* window);
psych_bool PsychCocoaSetupAndAssignOpenGLContextsFromCGLContexts(void* window, PsychWindowRecordType *windowRecord);
void PsychCocoaSendBehind(void* window);
void PsychCocoaSetWindowLevel(void* window, int inLevel);
void PsychCocoaSetWindowAlpha(void* window, float inAlpha);
void PsychCocoaSetThemeCursor(int inCursor);

#endif
