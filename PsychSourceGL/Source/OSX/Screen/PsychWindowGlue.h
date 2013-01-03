/*
	PsychToolbox3/Source/OSX/Screen/PsychWindowGlue.h
	
	PLATFORMS:	
	
		This is the OS-X Coregraphics version only.  
				
	AUTHORS:
	
		Allen Ingling		awi		Allen.Ingling@nyu.edu

	HISTORY:
	
		12/20/02		awi		Wrote it mostly by modifying SDL-specific refugees (from an experimental SDL-based Psychtoolbox).
		11/16/04		awi		Added description.  
		04/03/05                mk              Support for stereo display output and enhanced Screen('Flip') behaviour.
                05/09/05                mk              New function PsychGetMonitorRefreshInterval -- queries (and measures) monitor refresh.
	
        DESCRIPTION:
	
		Functions in this file comprise an abstraction layer for probing and controlling window state, except for window content.  
		
		Each C function which implements a particular Screen subcommand should be platform neutral.  For example, the source to SCREENPixelSizes() 
		should be platform-neutral, despite that the calls in OS X and Windows to detect available pixel sizes are different.  Platform 
		specificity is abstracted out in C files which end it "Glue", for example PsychScreenGlue, PsychWindowGlue, PsychWindowTextClue.

	NOTES:
	
	TO DO: 
	
		¥ The "glue" files should should be suffixed with a platform name.  The original (bad) plan was to distingish platform-specific files with the same 
		name by their placement in a directory tree.
		
		¥ PsychFlushGL() and PsychGLFlush() need to be merged.

        

*/

//include once
#ifndef PSYCH_IS_INCLUDED_PsychWindowGlue
#define PSYCH_IS_INCLUDED_PsychWindowGlue

#include "Screen.h"

psych_bool PsychRealtimePriority(psych_bool enable_realtime);
psych_bool PsychOSOpenOnscreenWindow(PsychScreenSettingsType *screenSettings, PsychWindowRecordType *windowRecord, int numBuffers, int stereomode, int conserveVRAM);
psych_bool PsychOSOpenOffscreenWindow(double *rect, int depth, PsychWindowRecordType **windowRecord);
void	PsychOSCloseOnscreenWindow(PsychWindowRecordType *windowRecord);
void	PsychOSCloseWindow(PsychWindowRecordType *windowRecord);
void	PsychOSCloseOffscreenWindow(PsychWindowRecordType *windowRecord);
void    PsychOSFlipWindowBuffers(PsychWindowRecordType *windowRecord);
void    PsychOSSetVBLSyncLevel(PsychWindowRecordType *windowRecord, int swapInterval);
void	PsychOSSetGLContext(PsychWindowRecordType *windowRecord);
void	PsychOSUnsetGLContext(PsychWindowRecordType *windowRecord);
double  PsychOSGetVBLTimeAndCount(PsychWindowRecordType *windowRecord, psych_uint64* vblCount);
void	PsychOSSetUserGLContext(PsychWindowRecordType *windowRecord, psych_bool copyfromPTBContext);
psych_bool PsychOSSetupFrameLock(PsychWindowRecordType *masterWindow, PsychWindowRecordType *slaveWindow);
psych_int64 PsychOSScheduleFlipWindowBuffers(PsychWindowRecordType *windowRecord, double tWhen, psych_int64 targetMSC, psych_int64 divisor, psych_int64 remainder, unsigned int specialFlags);
psych_int64 PsychOSGetSwapCompletionTimestamp(PsychWindowRecordType *windowRecord, psych_int64 targetSBC, double* tSwap);
void	PsychOSProcessEvents(PsychWindowRecordType *windowRecord, int flags);
//end include once
#endif

