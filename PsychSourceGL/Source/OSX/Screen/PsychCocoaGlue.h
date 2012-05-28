/*
 
 PsychToolbox3/Source/OSX/Screen/PsychCocoaGlue.h
 
 PLATFORMS:	
 
 This is the OS X 64-Bit Cocoa version only.
 
 TARGETS:
 
 Screen
 
 AUTHORS:
 
 Mario Kleiner       mk      mario.kleiner@tuebingen.mpg.de
 
 DESCRIPTION:
 
 Glue code for window management for the 64-Bit Psychtoolbox,
 using Objective-C wrappers to use Cocoa.
 
 These functions are called by PsychWindowGlue.c. Some are
 specific to Cocoa window management. Some are drop-in replacements
 for missing Carbon 32-Bit functions.
 
*/

#ifndef PSYCH_COCOAGLUE_H
#define PSYCH_COCOAGLUE_H

// Create windowed onscreen window via Cocoa:
PsychError PsychCocoaCreateWindow(PsychWindowRecordType *windowRecord,
                                  PsychRectType      screenRect,
                                  const Rect *       contentBounds,
                                  WindowClass        wclass,
                                  WindowAttributes   addAttribs,
                                  WindowRef *        outWindow,
                                  psych_bool         enableTransparentGL);

// Close Cocoa onscreen window, release OpenGL resources:
void PsychCocoaDisposeWindow(PsychWindowRecordType *windowRecord);

// Imposer for Carbon function: Returns windowref of current window which has
// keyboard input focus (for GetMouse()):
WindowRef GetUserFocusWindow(void);

// Show the window onscreen:
void PsychCocoaShowWindow(WindowRef window);

// Legacy path for < 10.6: Setup minimally viable OpenGL contexts for Cocoa window:
psych_bool PsychCocoaSetupAndAssignLegacyOpenGLContext(WindowRef window, PsychWindowRecordType *windowRecord);

// Modern >= 10.6 path: Setup all OpenGL contexts and CGL contexts:
psych_bool PsychCocoaSetupAndAssignOpenGLContextsFromCGLContexts(WindowRef window, PsychWindowRecordType *windowRecord);

#endif
