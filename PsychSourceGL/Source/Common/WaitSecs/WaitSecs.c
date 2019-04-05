/*
	PsychToolbox3/Source/Common/WaitSecs/WaitSecs.c
	
	PLATFORMS:	
	
		All.

	AUTHORS:

		Allen Ingling		awi		Allen.Ingling@nyu.edu
		Mario Kleiner		mk		mario.kleiner@tuebingen.mpg.de

	HISTORY:

		1/20/02			awi		wrote it.
		4/6/05			awi		Use mach_wait_until() instead of looping.  Mario's suggestion.  
		4/7/05			awi		Relocate mach_wait_until() call within PsychWaitIntervalSeconds().
		1/2/08			mk		Add subfunction for waiting until absolute time, and return of wakeup time. 
		

	NOTES: 
	
		We add a small error here by doing some stuff before calling PsychWaitIntervalSeconds.  The error should be small.  
		If we need to do better, mark time from the entry point into WAITSECSWaitSecs().
		
		We add a compile flag to the project MachTimebase settings: "-Wno-long-double" turns off a warning caused by using type long double.  
		The warning is:
			use of `long double' type; its size may change in a future release (Long double usage is reported only once for each file.
			To disable this warning, use -Wno-long-double.)
			
		Turning off the warning is probably ok because we do not depend on long doubles being any particular size or precision; if does change,
		that's not going to break anything.  
			
		The -Wno-long-double flag is appended to the project setting "Other Warning Flags".
  
*/


#include "WaitSecs.h"

#define MAX_SYNOPSIS_STRINGS 500
static const char *synopsisSYNOPSIS[MAX_SYNOPSIS_STRINGS];

const char** InitializeSynopsis(void)
{
    int i = 0;
    const char **synopsis = synopsisSYNOPSIS;  // abbreviate the long name

    synopsis[i++] = "";
    synopsis[i++] = "WaitSecs - Timed waits:";
    synopsis[i++] = "-----------------------";
    synopsis[i++] = "";
    synopsis[i++] = "[realWakeupTimeSecs] = WaitSecs(waitPeriodSecs);              -- Wait for at least 'waitPeriodSecs' seconds. Try to be precise.";
    synopsis[i++] = "[realWakeupTimeSecs] = WaitSecs('UntilTime', whenSecs);       -- Wait until at least time 'whenSecs'.";
    synopsis[i++] = "[realWakeupTimeSecs] = WaitSecs('YieldSecs', waitPeriodSecs); -- Wait for at least 'waitPeriodSecs' seconds. Be more sloppy.";
    synopsis[i++] = "\nThe optional 'realWakeupTimeSecs' is the real system time when WaitSecs finished waiting,";
    synopsis[i++] = "just as if you'd call realWakeupTimeSecs = GetSecs; after calling WaitSecs. This for your";
    synopsis[i++] = "convenience and to reduce call overhead and drift a bit for this common combo of commands.";
    synopsis[i++] = NULL;

    return(synopsisSYNOPSIS);
}

PsychError WAITSECSSynopsis(void)
{
    int i;

    for (i = 0; synopsisSYNOPSIS[i] != NULL; i++)
        printf("%s\n",synopsisSYNOPSIS[i]);

    return(PsychError_none);
}

PsychError WAITSECSWaitSecs(void)
{
    double	waitPeriodSecs;
    double	now;

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(1));

    if (!PsychCopyInDoubleArg(1, FALSE, &waitPeriodSecs)) {
        // Called without arguments. Output synopsis:
        WAITSECSSynopsis();
        return(PsychError_none);
    }

    // Wait for requested interval:
    PsychWaitIntervalSeconds(waitPeriodSecs);

    // Return current system time at end of sleep:
    PsychGetAdjustedPrecisionTimerSeconds(&now);
    PsychCopyOutDoubleArg(1, FALSE, now);

    return(PsychError_none);
}

PsychError WAITSECSWaitUntilSecs(void)
{
    static char useString[] = "[realWakeupTimeSecs] = WaitSecs('UntilTime', whenSecs);";
    //                                                                      1 
    static char synopsisString[] = 
    "Wait until at least system time \"whenSecs\" has been reached. "
    "Optionally, return the real wakeup time \"realWakeupTimeSecs\".\n"
    "This allows conveniently waiting until an absolute point in time "
    "has been reached, or to allow drift-free waiting for a well defined "
    "interval, more accurate than the standard WaitSecs() call.\n"
    "Example:\n"
    "Wait until 0.6 secs after last stimulus onset, if vbl=Screen('Flip', window); "
    "was the onset timestamp vbl from a previous flip:\n"
    "realwakeup = WaitSecs('UntilTime', vbl + 0.6);\n\n"
    "In a perfect world, realwakeup == vbl + 0.6, in reality it will be\n"
    "realwakeup == vbl + 0.6 + randomjitter; with randomjitter being the "
    "hopefully small scheduling delay of your operating system. If the "
    "delay is high or varies a lot between trials then your system has "
    "noisy timing or real timing problems.\n";
	
    static char seeAlsoString[] = "";	

    double	waitUntilSecs;
    double	now;

    //all sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(1));
    
    PsychCopyInDoubleArg(1,TRUE,&waitUntilSecs);
    PsychWaitUntilSeconds(waitUntilSecs);

    // Return current system time at end of sleep:
    PsychGetAdjustedPrecisionTimerSeconds(&now);
    PsychCopyOutDoubleArg(1, FALSE, now);

    return(PsychError_none);	
}

PsychError WAITSECSYieldSecs(void)
{
    static char useString[] = "[realWakeupTimeSecs] = WaitSecs('YieldSecs', waitPeriodSecs);";
    //                                                                      1 
    static char synopsisString[] = 
    "Wait for at least \"waitPeriodSecs\", don't care if it takes a few milliseconds longer. "
    "Optionally, return the real wakeup time \"realWakeupTimeSecs\".\n"
	"This call is useful if you want your code to release the cpu for a few milliseconds, "
	"e.g., to avoid overloading the cpu in a spinning loop, and you don't care if the "
	"wait takes a few msecs longer than specified. If you do care, use one of the other "
	"WaitSecs() variants! The other variants emphasize accuracy of timed waits, even if "
	"this causes a high load on the processor.\n";
	
    static char seeAlsoString[] = "";	

    double	waitPeriodSecs;
    double	now;

    //all sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(1));
    
    PsychCopyInDoubleArg(1,TRUE,&waitPeriodSecs);
    PsychYieldIntervalSeconds(waitPeriodSecs);

    // Return current system time at end of sleep:
    PsychGetAdjustedPrecisionTimerSeconds(&now);
    PsychCopyOutDoubleArg(1, FALSE, now);

    return(PsychError_none);	
}
