/*
	SCREENFlip.c		
  
	AUTHORS:

		Allen.Ingling@nyu.edu               awi
		mario.kleiner at tuebingen.mpg.de   mk
  
	PLATFORMS:	
	
		All.

	HISTORY:

		09/10/02	awi		Created.
		10/12/04	awi		Changed "SCREEN" to "Screen" in useString.
		04/03/05    mk      Add optional sync/nosync to VBL, don't clear fb on flip, flip after deadline, and return timestamps.
		05/16/05	mk      Add optional flag "dontsync" and some more timestamps.
		06/09/05	mk      Add optional flag "multiflip" for experimental multiflip support.
 
	DESCRIPTION:
  
		Interchanges front and back buffers of a double-buffered display.
		It can be selected if buffer should be cleared to background color after the flip or if
		it is left intact. One can also select to flip after a specific deadline -
		useful for easier timing of stimulus presentation. And Flip returns multiple timestamps
		that help to check proper presentation timing.

	TO DO:
  
*/

#include "Screen.h"

// SCREENFlip internally implements both Screen('Flip') and Screen('FlipAsyncBegin') and Screen('FlipAsyncEnd'):
PsychError SCREENFlip(void) 
{
	// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
	static char useString0[] = "[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('Flip', windowPtr [, when] [, dontclear] [, dontsync] [, multiflip]);";
	static char synopsisString0[] = 
	"Flip front and back display surfaces in sync with vertical retrace and return completion timestamps.\n"
	"\"windowPtr\" is the id of the onscreen window whose content should be shown at flip time. "
	"\"when\" specifies when to flip: If set to zero (default), it will flip on the next possible video retrace. "
	"If set to a value when > 0, it will flip at the first video retrace after system time 'when' has been reached. "
	"\"dontclear\" If set to 1, flip will not clear the framebuffer after Flip - this allows incremental drawing "
	"of stimuli. The default is zero, which will clear the framebuffer to background color after each flip. "
	"A value of 2 will prevent Flip from doing anything to the framebuffer after flip. This leaves the job "
	"of setting up the buffer to you - the framebuffer is in an undefined state after flip. "
	"\"dontsync\" If set to zero (default), Flip will sync to the vertical retrace and will pause execution of your script "
	"until the Flip has happened. If set to 1, Flip will still synchronize stimulus onset to the vertical retrace, but "
	"will *not* wait for the flip to happen: Flip returns immediately and all returned timestamps are invalid. A value of 2 "
	"will cause Flip to show the stimulus *immediately* without waiting/syncing to the vertical retrace. "
	"\"multiflip\" defaults to zero: If set to a value greater than zero, Flip will flip *all* onscreen windows instead of "
	"just the specified one. This allows to synchronize stimulus onset on multiple displays, e.g., for multidisplay stereo "
	"setups or haploscopes. You need to (somehow) synchronize all attached displays for this to operate tear-free. "
	"Flip (optionally) returns a high-precision estimate of the system time (in seconds) when the actual flip has happened "
	"in the return argument 'VBLTimestamp'. An estimate of Stimulus-onset time is returned in 'StimulusOnsetTime'. "
	"Beampos is the position of the monitor scanning beam when the time measurement was taken (useful for correctness tests). " 
	"FlipTimestamp is a timestamp taken at the end of Flip's execution. Use the difference between FlipTimestamp and VBLTimestamp "
	"to get an estimate of how long Flips execution takes. This is useful to get a feeling for the timing error if you try to sync "
	"script execution to the retrace, e.g., for triggering acquisition devices like EEG, fMRI, or for starting playback of a sound. "
	"\"Missed\" indicates if the requested presentation deadline for your stimulus has been missed. A negative value means that dead- "
	"lines have been satisfied. Positive values indicate a deadline-miss. The automatic detection of deadline-miss is not fool-proof "
	"- it can report false positives and also false negatives, although it should work fairly well with most experimental setups. "
	"If you are picky about timing, please use the provided timestamps or additional methods to exercise your own tests. ";

	static char useString1[] = "[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('AsyncFlipBegin', windowPtr [, when] [, dontclear] [, dontsync] [, multiflip]);";
	static char synopsisString1[] = 
	"Schedule an asynchronous flip of front and back display surfaces for given onscreen window. "
	"\"windowPtr\" is the id of the onscreen window whose content should be shown at flip time. "
	"\"when\" is the requested stimulus onset time, a value of zero or no argument asks for flip "
	"at next possible vertical retrace.\n"
	"For the meaning and explanation of all other parameters, see the help for 'Screen Flip?'.\n\n"
    "If this command is called while a previously scheduled asynchronous flip is still in progress, it will "
    "wait for that pending async flip to finish and return its results (timestamps etc.). If no such "
    "operation is in progress, it will return results of the most recently finished async or sync flip. "
    "Waiting for previous flips to complete and returning their results is a convenience function. "
    "In most cases, in order to have more control over "
    "the execution of your script and your flip timing, you will rather want to use one of the finalizer "
	"commands Screen('AsyncFlipCheckEnd') or Screen('AsyncFlipEnd') mentioned below "
	"to collect information about the final result and timing of the asynchronous flip operation.\n\n"
	"The difference between Screen('AsyncFlipBegin',...); and the more commonly used Screen('Flip', ...); "
	"is that Screen('Flip') operates synchronously: Execution of your code is paused until the flip operation "
	"has finished, ie. at least until the requested onset deadline 'when' has passed.\n\n"
	"Screen('AsyncFlipBegin') will prepare everything for a flip at the requested time 'when' - or "
	"at next vertical retrace if 'when' is omitted - but then immediately return control to your "
	"code. Your code can continue to execute and do other things, e.g., schedule flips for other onscreen "
	"windows, perform keyboard or mouse queries, etc.\n\n"
	"You can check the progress state of asynchronous flips or wait for them to finish and collect timing "	
	"information for the finalized flip by use of the commands Screen('AsyncFlipCheckEnd') and Screen('AsyncFlipEnd');\n\n"
	"In general you should avoid using asynchronous flips and instead use conventional 'Flip' unless you have a "
	"good reason to do otherwise, because async flips come with a couple of strings attached:\n\n"
	"- You are restricted in what you can do with Screen() or OpenGL while async\n"
	"  flips are in progress: You can not do anything with textures or offscreen\n"
	"  windows while their parent-onscreen window is in async flip state. You can\n"
	"  only access onscreen windows which are not participating in an async flip\n"
	"  operation.\n\n"
    "- If you enable the Psychtoolbox image processing pipeline, most restrictions\n"
    "  on drawing during async flips are relaxed: You can draw into any windows while\n"
    "  async flips are pending, even the window for which the flip is pending. Only\n"
    "  use of the Screen('GetImage') command is forbidden on async flipping onscreen\n"
    "  windows and potentially problematic on offscreen windows while an async flip\n"
    "  is in progress. However, this is somewhat theoretical. In practice many\n"
    "  operating systems, graphics drivers and graphics cards can't really handle the\n"
    "  load of parallel drawing and async flipping, due to system bugs or design\n"
    "  constraints. On such systems you may observe inconsistent timing, degraded\n"
    "  performance, and on some systems even visual stimulus corruption, malfunctions\n"
    "  or hard system crashes [e.g., Apple MacOS/X 10.4.11 with ATI Radeon X1600]!\n\n"
    "- Even the restricted set of allowed reliably working Screen/OpenGL commands\n"
	"  should be avoided, because some graphics hardware and drivers may not be\n"
	"  able to handle such concurrent graphics operations without degraded stimulus\n"
	"  onset timing accuracy, ie. you may experience more missed stimulus deadlines\n"
	"  and timing glitches -- or inconsistent behaviour accross different computers\n"
	"  and graphics cards or operating system releases. In the end it allows you to\n"
	"  do non-Screen related things like sound, I/O, keyboard checks...\n\n"
	"- Parallel processing of flips puts additional burden onto your CPU,\n"
	"  GPU and operating system, so it incurs additional overhead and may degrade\n"
	"  absolute drawing performance and cause more timing issues and glitches if\n"
	"  your system is not reliably able to handle the concurrent load.\n\n"
	"- Code with async flips - as any piece of parallely executing code - is harder\n"
	"  to implement correctly and more challenging to debug for you.\n\n"
	"- Using a non-zero \"multiflip\" argument is not allowed.\n\n"
	"- Asynchronous updates of gamma tables will likely not work reliably.\n\n"
	"- Stereo stimulus display in stereomode 10 (two separate onscreen windows) will\n"
    "  likely not work with reliable timing or have possible tearing artifacts.\n\n"
	"- Use of the 'UserspaceBufferDrawingPrepare' hook-chain of the imaging\n"
	"  pipeline is not allowed.\n"
	"\n\n"
	"Our general stance is that most code can be written efficiently without need for async flips, so this feature is "
	"provided for the few demanding special cases where this is not the case and the benefits outweight the costs.";

	static char useString2[] = "[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('AsyncFlipEnd', windowPtr);";
	static char synopsisString2[] = 
	"Wait for completion of a previously scheduled asynchronous flip operation (see Screen AsyncFlipBegin? for help).\n"
	"This command will wait for completion on onscreen window \"windowPtr\", then return the result of the operation, "
	"ie. all the stimulus onset timestamps and other diagnostic information. See help for Screen Flip? for explanation "
	"of the returned info. If you call this function without having called 'AsyncFlipBegin' before, it will return "
    "immediately with the results of the last completed flip, regardless if it was a synchronous Screen('Flip') or an "
    "asynchronous flip long finished ago.\n"
	"Screen('AsyncFlipCheckEnd') provides a non-blocking, polling version of this command -- one that doesn't pause if "
	"the referenced operation hasn't completed yet.";

	static char useString3[] = "[VBLTimestamp StimulusOnsetTime FlipTimestamp Missed Beampos] = Screen('AsyncFlipCheckEnd', windowPtr);";
	static char synopsisString3[] = 
	"Check if a previously scheduled asynchronous flip operation has completed (see Screen AsyncFlipBegin? for help).\n"
	"This command will check for completion on onscreen window \"windowPtr\", then return the result of the operation, "
	"ie. all the stimulus onset timestamps and other diagnostic information. See help for Screen Flip? for explanation "
	"of the returned info. If the operation hasn't completed yet, it will return a 'VBLTimestamp' of zero "
	"and you'll have to retry later. If you call this function without having called 'AsyncFlipBegin' before, "
    "it will return the results of the last completed flip, regardless if it was a synchronous Screen('Flip') or an "
    "asynchronous flip long finished ago.\n"
	"Screen('AsyncFlipEnd') provides a blocking version of this command -- one that pauses until the operation completes if "
	"the referenced operation hasn't completed yet.";

	static char seeAlsoString[] = "DrawingFinished WaitUntilAsyncFlipCertain AsyncFlipBegin AsyncFlipCheckEnd AsyncFlipEnd Flip";

	PsychFlipInfoStruct*  flipRequest;
	PsychWindowRecordType *windowRecord;
	int dont_clear;
	int vbl_synclevel;
	int multiflip;
	double flipwhen, tNow;
	double vbl_timestamp;
	int beamposatflip;
	double miss_estimate;
	double time_at_flipend;
	double time_at_onset;
	unsigned int opmode;
	psych_bool flipstate;

	// Change our "personality" depending on the name with which we were called:
	if (PsychMatch(PsychGetFunctionName(), "AsyncFlipBegin")) {
		// Async flip invocation:
		opmode = 1;

		// Push usage string and/or give online help:
		PsychPushHelp(useString1, synopsisString1, seeAlsoString);
	}
	else if (PsychMatch(PsychGetFunctionName(), "AsyncFlipEnd")) {
		// Finalize async flip operation:
		opmode = 2;

		// Push usage string and/or give online help:
		PsychPushHelp(useString2, synopsisString2, seeAlsoString);
	}
	else if (PsychMatch(PsychGetFunctionName(), "AsyncFlipCheckEnd")) {
		// Poll for finalization of async flip operation:
		opmode = 3;

		// Push usage string and/or give online help:
		PsychPushHelp(useString3, synopsisString3, seeAlsoString);
	}
	else {
		// Standard synchronous flip:
		opmode = 0;

		// Push usage string and/or give online help:
		PsychPushHelp(useString0, synopsisString0, seeAlsoString);
	};
	
	// Give online help, if requested:
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

	PsychErrorExit(PsychCapNumInputArgs((opmode < 2)  ? 5 : 1));		// The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(1));						// The required number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(5));							// The maximum number of outputs
	
	// Get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);
	
	if (!PsychIsOnscreenWindow(windowRecord)) PsychErrorExitMsg(PsychError_user, "Flip called on something else than an onscreen window. You can only flip onscreen windows.");
	if(windowRecord->windowType!=kPsychDoubleBufferOnscreen) PsychErrorExitMsg(PsychError_user, "Flip called on window without backbuffers. Specify numberOfBuffers=2 in Screen('OpenWindow') if you want to use Flip.");
	
	// Only retrieve additional arguments if this isn't a finish on an async flip:
	if ((opmode != 2) && (opmode != 3)) {
		// Query optional dont_clear argument: 0 (default) clear backbuffer to background color after flip.
		// 1 == Restore backbuffer to state before flip - this allows incremental drawing/updating of stims.
		dont_clear=0;
		PsychCopyInIntegerArg(3,FALSE,&dont_clear);
		if(dont_clear < 0 || dont_clear > 2) {
			PsychErrorExitMsg(PsychError_user, "Only 'dontclear' values 0 (== clear after flip), 1 (== don't clear) and 2 (== don't do anything) are supported");
		}
		
		// Query optional vbl_synclevel argument: 0 (default) Synchronize flip to VBL and wait (block) for VBL onset.
		// 1 == Synchronize flip to VBL, but don' wait for flip to happen -- timestamps are invalid.
		// 2 == Don't sync to the retrace *at all* Flip as soon as stimulus is ready.
		// EXPERIMENTAL modes:
		// 3 == Like 0, but use spin-waiting to busy wait for VBL onset -> Non-blocking -> Low delay -> But dangerous!
		vbl_synclevel=0;
		PsychCopyInIntegerArg(4,FALSE,&vbl_synclevel);
		if(vbl_synclevel < 0 || vbl_synclevel > 3) {
			PsychErrorExitMsg(PsychError_user, "Only 'dontsync' values 0 (== fully synchronize with VBL), 1 (== don't wait for VBL) and 2 (== Ignore VBL) are supported");
		}
		
		// Query optional multiflip argument: 0 (default) Only flip the requested onscreen window.
		// 1 = Flip *all* onscreen windows simultaneously, but individually synced to their respective VBL onset.
		// -> This is tear- and flicker-free and ressource efficient due to its blocking nature, but it requires
		// displays that retrace in sync, otherwise stimulus onset on each display can happen phase-shifted to other
		// displays.
		// 2 = Like 1, but we wait for VBL onset on requested onscreen window and then immediately flip all other
		// windows without waiting for their respective retrace -> All stimuli show up at the same time/with minimal
		// time delay, but if monitors are not in-phase, some of them will tear/flicker.
		multiflip=0;
		PsychCopyInIntegerArg(5,FALSE,&multiflip);
		if(multiflip < 0 || multiflip > 2)
			PsychErrorExitMsg(PsychError_user, "Only 'multiflip' values 0 (== don't multiflip), 1 (== multiflip method 1) and 2 (== multiflip method 2) are supported");
		if (multiflip==2 && vbl_synclevel!=0 && vbl_synclevel!=3) {
			PsychErrorExitMsg(PsychError_user, "Only 'dontsync' values 0 and 3 are allowed when multiflip is set to 2.");
		}
		
		// Current multiflip > 0 implementation is not thread-safe, so we don't support this in async flip mode:
		if ((multiflip != 0) && (opmode != 0))  PsychErrorExitMsg(PsychError_user, "Using a non-zero 'multiflip' flag while starting an asynchronous flip! Sorry, this is currently not possible.\n");

		// Query optional flipwhen argument:
		// 0 (default value) == Flip at next vertical retrace and sync to VBL. This is the old PTB behaviour as of PTB 1.0.42.
		// flipwhen>0 == Sleep/Wait until system time "flipwhen" and then flip at the next VBL after time "flipwhen".
		flipwhen=0;
		PsychCopyInDoubleArg(2,FALSE,&flipwhen);
		if(flipwhen < 0) {
			PsychErrorExitMsg(PsychError_user, "Only 'when' values greater or equal to 0 are supported");
		}

		PsychGetAdjustedPrecisionTimerSeconds(&tNow);
		if (flipwhen - tNow > 1000) {
			PsychErrorExitMsg(PsychError_user, "\nYou specified a 'when' value to Flip that's over 1000 seconds in the future?!? Aborting, assuming that's an error.\n\n");
		}

		// Pack all parameters of the fliprequest into the flipinfo struct:
		// At least our part of the struct. Initial setup of threds and locks etc. is done by the actual
		// PsychFlipWindowBuffersIndirect() routine:
		flipRequest = windowRecord->flipInfo;

		if (flipRequest->asyncstate != 0) {
			// Started, executing or finalized async flip in progress. We can't trigger a new flip request
			// before the current one has finished. Perform a blocking wait for flip completion, basically
			// a Screen('AsyncFlipEnd') op, collect its results for return to usercode, then continue with
			// scheduling the new flip request:
			flipRequest->opmode = 2;
			flipstate = PsychFlipWindowBuffersIndirect(windowRecord);

			// Reset state to zero, ie. ready for new adventures ;-)
			if (flipstate) flipRequest->asyncstate = 0;

			// Ok, async flip completed. Its completion data is stored in flipRequest.
		}
		else {
			// No async flip -- Fake a "success", so cached results from
			// previous flips can be returned:
			flipstate = TRUE;
		}

		// If this is a Screen('AsyncFlipBegin') aka opmode 1, we can return the
		// completion data of previous async flips to usercode.
		// In opmode 0 aka Screen('Flip') returning data here would not make
		// sense as userspace expects the results from the synchronous flip we will
		// schedule next.
		if (opmode == 1) {
			vbl_timestamp	= (flipstate) ? flipRequest->vbl_timestamp : 0.0;
			time_at_onset	= flipRequest->time_at_onset;
			time_at_flipend = flipRequest->time_at_flipend;
			miss_estimate	= flipRequest->miss_estimate;
			beamposatflip	= flipRequest->beamPosAtFlip;
            
			// Return timestamp at start of VBL time:
			PsychCopyOutDoubleArg(1, FALSE, vbl_timestamp);
			// Return timestamp at stimulus onset time:
			PsychCopyOutDoubleArg(2, FALSE, time_at_onset);
			// Return time when Flip ended:
			PsychCopyOutDoubleArg(3, FALSE, time_at_flipend);
			// Return current estimate of deadline miss, if any:
			PsychCopyOutDoubleArg(4, FALSE, miss_estimate);
			// Return beam position at VBL time:
			PsychCopyOutDoubleArg(5, FALSE, (double) beamposatflip);
		}

		// This info needs to be provided for flip mechanism:
		flipRequest->opmode		= opmode;
		flipRequest->dont_clear		= dont_clear;
		flipRequest->flipwhen		= flipwhen;
		flipRequest->multiflip		= multiflip;
		flipRequest->vbl_synclevel	= vbl_synclevel;		

		// This field is set to -1. It will be set by the flipperThread at flip
		// completion. We do this, so the 'WaitUntilAsyncFlipCertain' subroutine
		// knows if a swap is still pending, or really finished. If this value is
		// not -1, then the routine can simply skip its wait op and take values from
		// the initialized timestamps of the flipperThread, as we know that that
		// thread already has detected swap-completion and done all the timestamping
		// work:
		flipRequest->vbl_timestamp = -1;
		
		// Ok, the struct is filled with spec for a synchronous or asynchronous flip...
		
		// Store current preflip GPU graphics surface addresses, if supported:
		PsychStoreGPUSurfaceAddresses(windowRecord);
	}
	else {
		// opmode == 2 or 3 - 'AsyncFlipEnd' or 'AsyncFlipCheckEnd':
		flipRequest = windowRecord->flipInfo;
		if (flipRequest->asyncstate == 0) {
			// No started, executing or finalized async flip in progress! No async flip operation triggered
			// which we could finalize. This is fine. We basically no-op and return the cached last known
			// values from previous async flips or sync flips, or all zeros if no flip was ever executed:
			vbl_timestamp	= flipRequest->vbl_timestamp;
			time_at_onset	= flipRequest->time_at_onset;
			time_at_flipend = flipRequest->time_at_flipend;
			miss_estimate	= flipRequest->miss_estimate;
			beamposatflip	= flipRequest->beamPosAtFlip;
            
			// Return timestamp at start of VBL time:
			PsychCopyOutDoubleArg(1, FALSE, vbl_timestamp);
			// Return timestamp at stimulus onset time:
			PsychCopyOutDoubleArg(2, FALSE, time_at_onset);
			// Return time when Flip ended:
			PsychCopyOutDoubleArg(3, FALSE, time_at_flipend);
			// Return current estimate of deadline miss, if any:
			PsychCopyOutDoubleArg(4, FALSE, miss_estimate);
			// Return beam position at VBL time:
			PsychCopyOutDoubleArg(5, FALSE, (double) beamposatflip);
            
			// We're done:
			return(PsychError_none);
		}

		// Ok, there's an async flip going on and we have a handle to it...
		flipRequest->opmode	= opmode;
	}

	// Execute the flip request with the indirect routine: In opmode == 0 (Synchronous flip)
	// this will block until flip is done, then return with the results. In async mode it
	// will dispatch the flip to the helper thread, then return immediately. In async end
	// or async poll mode, it will either block until async op finished, or return after
	// polling with a FALSE result, telling that not yet finished:
	flipstate = PsychFlipWindowBuffersIndirect(windowRecord);
	
	// Only have return args in synchronous mode or in return path from end/successfull poll of async flip:
	if (opmode != 1) {
		// Async flip is either zero in synchronous mode, or it's 2 if an async flip
		// successfully finished:
		if ((flipRequest->asyncstate!=0) && (flipRequest->asyncstate!=2) && (flipstate)) {
			printf("PTB-ERROR: flipRequest->asyncState has impossible value %i at end of flipop! This is a PTB DESIGN BUG!", flipRequest->asyncstate);
			PsychErrorExitMsg(PsychError_internal, "flipRequest->asyncState has impossible value at end of flipop! This is a PTB DESIGN BUG!");
		}
		
		// Reset it to zero, ie. ready for new adventures ;-)
		if (flipstate) flipRequest->asyncstate=0;
		
		// Return return arguments from flip: We return a zero vbl_timestamp in case a poll for flip completion failed.
		// That indicates that flip not yet finished and all other return values are invalid:
		vbl_timestamp	= (flipstate) ? flipRequest->vbl_timestamp : 0.0;
		time_at_onset	= flipRequest->time_at_onset;
		time_at_flipend = flipRequest->time_at_flipend;
		miss_estimate	= flipRequest->miss_estimate;
		beamposatflip	= flipRequest->beamPosAtFlip;
		
		// Return timestamp at start of VBL time:
		PsychCopyOutDoubleArg(1, FALSE, vbl_timestamp);
		// Return timestamp at stimulus onset time:
		PsychCopyOutDoubleArg(2, FALSE, time_at_onset);
		// Return time when Flip ended:
		PsychCopyOutDoubleArg(3, FALSE, time_at_flipend);
		// Return current estimate of deadline miss, if any:
		PsychCopyOutDoubleArg(4, FALSE, miss_estimate);
		// Return beam position at VBL time:
		PsychCopyOutDoubleArg(5, FALSE, (double) beamposatflip);
        
		// Execute hook chain for preparation of user space drawing ops:
		PsychPipelineExecuteHook(windowRecord, kPsychUserspaceBufferDrawingPrepare, NULL, NULL, FALSE, FALSE, NULL, NULL, NULL, NULL);        
	}
	
	return(PsychError_none);
}

PsychError SCREENWaitUntilAsyncFlipCertain(void) 
{
	// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
	static char useString[] = "[VBLTimestamp StimulusOnsetTime swapCertainTime] = Screen('WaitUntilAsyncFlipCertain', windowPtr);";
	static char synopsisString[] =
	"Wait until it is certain that a previously initiated Screen('AsyncFlipBegin',...); "
	"operation for onscreen window 'windowPtr' will happen at next vertical retrace or has happened "
	"already, then return timestamps of when flip operation has happened or when the flip operation will happen.\n"
	"This function only works on Linux and MacOS/X and only with recent AMD/ATI hardware of the "
	"Radeon X1000/HD2000/... series or later (or the equivalent FireGL cards). On MacOS/X, the "
	"PsychtoolboxKernelDriver needs to be loaded for this to work (see help PsychtoolboxKernelDriver).\n"
	"The method tries to \"read the mind\" of your graphics card to find out if the card will swap the "
	"display buffers at the next vertical retrace to interrogate if the swap is certain, ie., not possibly "
	"delayed by graphics card overload, other overload conditions or system timing jitter/delays.\n"
	"Ideally it will allow your script to know when stimulus onset will happen a few milliseconds before "
	"stimulus onset. This allows, e.g., to initiate some other operations that need to be tightly and "
	"reliably synchronized to stimulus onset and that have some small startup latency, so they need some "
	"headstart. Examples would be sound output or output of some slow trigger signal which takes multiple "
	"milliseconds to take action.\n"
	"This feature only works in fullscreen mode and is considered experimental - It may or may not work "
	"reliably on your setup. Good luck!\n"
	"Returns a high-precision estimate of the system time (in seconds) when the actual flip will- or has happened "
	"in the return argument 'VBLTimestamp'. An estimate of Stimulus-onset time is returned in 'StimulusOnsetTime'. "
	"'swapCertainTime' is the system time when certainty of bufferswap was detected. ";
	static char seeAlsoString[] = "'AsyncFlipBegin' 'AsyncFlipCheckEnd' 'AsyncFlipEnd'";

#if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_LINUX	
	PsychWindowRecordType *windowRecord;
	double timestamp, stimulusOnset;
	int	   beamposition;
	double vbl_startline, vbl_endline, vbl_lines_elapsed, onset_lines_togo;
	double vbl_time_elapsed, onset_time_togo, currentrefreshestimate;
	psych_bool   swappending;
	long   scw, sch;
#endif

	// Push usage string and/or give online help:
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	
	// Give online help, if requested:
	if(PsychIsGiveHelp()) {PsychGiveHelp();return(PsychError_none);};

	PsychErrorExit(PsychCapNumInputArgs(1));		// The maximum number of inputs
	PsychErrorExit(PsychRequireNumInputArgs(1));	// The required number of inputs
	PsychErrorExit(PsychCapNumOutputArgs(3));		// The maximum number of outputs

#if PSYCH_SYSTEM == PSYCH_OSX || PSYCH_SYSTEM == PSYCH_LINUX
	
	// Get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);

	// Child protection:
	if (!PsychIsOnscreenWindow(windowRecord)) PsychErrorExitMsg(PsychError_user, "WaitUntilAsyncFlipCertain called on something else than an onscreen window.");
	if(windowRecord->windowType!=kPsychDoubleBufferOnscreen) PsychErrorExitMsg(PsychError_user, "WaitUntilAsyncFlipCertain called on window without backbuffers.");
	
	// Just if we are called on a window for which an async flip operation is active.
	// The routine can only be used for active async flips, so bail out on anything else.
	if ((windowRecord->flipInfo == NULL) || (windowRecord->flipInfo->asyncstate == 0)) {
		// No async flip operation active: Either no flip triggered at all, or at least not an async one:
		PsychErrorExitMsg(PsychError_user, "WaitUntilAsyncFlipCertain only works for async flips: May only be called between Screen('AsyncFlipBegin') and Screen('AsyncFlipEnd') or Screen('AsyncFlipCheckEnd')!");		
	}
	
	// Just need to check if GPU low-level access is supported:
	if (!PsychOSIsKernelDriverAvailable(windowRecord->screenNumber)) {
		PsychErrorExitMsg(PsychError_user, "WaitUntilAsyncFlipCertain needs the PsychtoolboxKernelDriver to be loaded on supported hardware!");
	}

	// Swap still pending? A non -1 vbl_timestamp would mean that the async flip helper
	// thread has already detected swap completion and calculated all relevant timestamps,
	// so we can skip all this and just take the timestamps provided by the flipper thread:
	if (windowRecord->flipInfo->vbl_timestamp == -1)  {
		// Swap still pending - Do your work:
		
		// Retrieve estimate of monitor refresh interval:
		if (windowRecord->VideoRefreshInterval > 0) {
			currentrefreshestimate = windowRecord->VideoRefreshInterval;
		}
		else {
			currentrefreshestimate=0;
			// We abort - This is too unsafe...
			PsychErrorExitMsg(PsychError_user,"Estimate of monitor refresh interval is INVALID -> Aborting!");
		}
		
		if (windowRecord->VBL_Endline == -1) {
			PsychErrorExitMsg(PsychError_user,"Estimate of VBL endline is INVALID -> Aborting!");
		}
		
		// Retrieve scheduled swap time: We use it as baseline for a timeout for the mechanism:
		timestamp = windowRecord->flipInfo->flipwhen;
		// If scheduled time is "on next retrace", we set baseline to "now":
		if (timestamp == 0) PsychGetAdjustedPrecisionTimerSeconds(&timestamp);

		// Should flip within 1 refresh of timestamp, ie., in reality at least within 2 or a few
		// refreshes. Certainly within 1 second. We generously set a deadline of 5 seconds for swap detection:
		timestamp+=5.0;

		// Wait for bufferswap completion or bufferswap certain (or until timeout time 'timestamp' elapsed):
		swappending = PsychWaitForBufferswapPendingOrFinished(windowRecord, &timestamp, &beamposition);
		if (timestamp == -1) PsychErrorExitMsg(PsychError_user, "Malfunctioned: Aborted due to timeout exceeded. Seems your graphics card doesn't support 'WaitUntilAsyncFlipCertain' properly. Sorry.");
		
		// Copy out optional timestamp of bufferswap detection:
		PsychCopyOutDoubleArg(3, kPsychArgOptional, timestamp);
		
		// Compute timestamps:
		PsychGetScreenSize(windowRecord->screenNumber, &scw, &sch);
		vbl_startline = (int) sch;
		vbl_endline = windowRecord->VBL_Endline;
		
		if (swappending) {
			// Swap is pending:
			onset_lines_togo = vbl_endline - beamposition + 1;
			vbl_lines_elapsed = vbl_startline - beamposition + 1; 
			
			// From the elapsed number we calculate the remaining time till swap et al.:
			vbl_time_elapsed = vbl_lines_elapsed / vbl_endline * currentrefreshestimate; 
			onset_time_togo = onset_lines_togo / vbl_endline * currentrefreshestimate;
			
			// Add to basetime:
			stimulusOnset = timestamp + onset_time_togo;
			timestamp = timestamp + vbl_time_elapsed;
			
			// printf("PENDING\n");
		}
		else {
			// Swap happened already: Do same calculations as in PsychFlipWindowBuffers():
			
			if (beamposition >= vbl_startline) {
				vbl_lines_elapsed = beamposition - vbl_startline;
				onset_lines_togo = vbl_endline - beamposition + 1;
			}
			else {
				vbl_lines_elapsed = vbl_endline - vbl_startline + 1 + beamposition;
				onset_lines_togo = -1.0 * beamposition;
			}
			
			// From the elapsed number we calculate the elapsed time since VBL start:
			vbl_time_elapsed = vbl_lines_elapsed / vbl_endline * currentrefreshestimate; 
			onset_time_togo = onset_lines_togo / vbl_endline * currentrefreshestimate;
			
			// Compute of stimulus-onset, aka time when retrace is finished:
			stimulusOnset = timestamp + onset_time_togo;
			// Now we correct our time_at_vbl by this correction value:
			timestamp = timestamp - vbl_time_elapsed;
		}
		
		// Copy out optional vbl timestamp of bufferswap:
		PsychCopyOutDoubleArg(1, kPsychArgOptional, timestamp);
		
		// Copy out optional timestamp of stimulus onset:
		PsychCopyOutDoubleArg(2, kPsychArgOptional, stimulusOnset);
	}
	else {
		// Swap already finished and timestamped by flipper thread.
		// Just retrieve and return its values:
		// Copy out optional vbl timestamp of bufferswap:
		PsychCopyOutDoubleArg(1, kPsychArgOptional, windowRecord->flipInfo->vbl_timestamp);
		
		// Copy out optional timestamp of stimulus onset:
		PsychCopyOutDoubleArg(2, kPsychArgOptional, windowRecord->flipInfo->time_at_onset);
		
		// Copy out optional timestamp of bufferswap detection:
		PsychGetAdjustedPrecisionTimerSeconds(&timestamp);
		PsychCopyOutDoubleArg(3, kPsychArgOptional, timestamp);
	}
#else
	PsychErrorExitMsg(PsychError_user, "Screen('WaitUntilAsyncFlipCertain') is not supported on Microsoft Windows.");
#endif

	// Done.
	return(PsychError_none);	
}
