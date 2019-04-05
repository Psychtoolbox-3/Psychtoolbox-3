/*
 *    PsychSourceGL/Source/Common/GetSecs/GetSecs.c
 *
 *    PROJECTS:
 *
 *      GetSecs only
 *
 *    AUTHORS:
 *
 *      Allen.Ingling@nyu.edu           awi
 *      mario.kleiner.de@gmail.com      mk
 *
 *    PLATFORMS:
 *
 *      All.
 *
 *    HISTORY:
 *
 *      1/20/02         awi     Derived the GetSecs project from Screen .
 *      4/6/05          awi     Updated header comments.
 *
 *   DESCRIPTION:
 *
 *      Return the time in seconds with high precision. On MS-Windows, it allows to
 *      trigger some debug operations and returns some diagnostic values that allow
 *      to diagnose timer problems. Not so on Linux or OS/X: They don't have such problems afaik ;-)
 */

#include "GetSecs.h"

#define MAX_SYNOPSIS_STRINGS 500
static const char *synopsisSYNOPSIS[MAX_SYNOPSIS_STRINGS];

const char** InitializeSynopsis(void)
{
    int i = 0;
    const char **synopsis = synopsisSYNOPSIS;
    synopsis[i++] = "[GetSecsTime, WallTime, syncErrorSecs] = GetSecs('AllClocks' [, maxError=0.000020]);";
    synopsis[i++] = NULL;

    return(synopsisSYNOPSIS);
}

#if PSYCH_SYSTEM != PSYCH_WINDOWS
#include <sys/time.h>
#endif

PsychError GETSECSGetSecs(void)
{
    double returnValue;
    double referenceValue, realValue;
    int healthy, opmode;

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(5));
    PsychErrorExit(PsychCapNumInputArgs(1));

    PsychGetAdjustedPrecisionTimerSeconds(&returnValue);
    PsychCopyOutDoubleArg(1, FALSE, returnValue);

    // Special code for diagnosing problems with TimeGlue on systems that
    // are broken by design(TM) aka MS-Windows:
    #if PSYCH_SYSTEM == PSYCH_WINDOWS
        if (PsychCopyInIntegerArg(1, FALSE, &opmode) && (opmode!=0)) {
            referenceValue = PsychGetTimeGetTimeValueAtLastTimeQuery(&realValue);
            healthy = (int) PsychGetTimeBaseHealthiness();
            PsychCopyOutDoubleArg(2, FALSE, referenceValue);
            PsychCopyOutDoubleArg(3, FALSE, realValue);
            PsychCopyOutDoubleArg(4, FALSE, (double) healthy);

            // A positive opmode value allows to change the thread affinity mask of the PTB Thread.
            // The old affinity mask is returned in that case. This to check for TSC sync across cpu cores:
            if (opmode > 0)
                PsychCopyOutDoubleArg(5, FALSE, (double) SetThreadAffinityMask(GetCurrentThread(), (DWORD) opmode));

            // An opmode setting of smaller than -1 will try to cancel our timeBeginPeriod(1) requests, as
            // automatically done by the PsychTimeGlue: We try to reset the low-res timer to something
            // like its normal 10 or 15 msecs duty cycle. Could help to spot timers that are actually
            // broken, but do work with the new PTB due to the increased IRQ load and therefore the
            // reduced power management. Post-hoc test if timing was reliable in earlier PTB releases.
            if (opmode < -1)
                timeEndPeriod(1);
        }
    #else
        // For MacOS/X and Linux, we return gettimeofday() as reference value and
        // some default "no error" healthy flag. Please note that gettimeofday() is just a
        // drift/NTP corrected version of the system timebase mach_absolute_time(), ie the
        // same value as returned as primary time, just with some offset. Therefore its not
        // possible to check for timebase malfunction, only for clock drift against some
        // external timebase, e.g., UTC as determined by NTP protocol.
        if (PsychCopyInIntegerArg(1, FALSE, &opmode) && (opmode!=0)) {
            struct timeval tv;
            gettimeofday(&tv, NULL);
            referenceValue = ((double) tv.tv_sec) + (((double) tv.tv_usec) / 1000000.0);
            healthy = 0;
            realValue = returnValue;
            PsychCopyOutDoubleArg(2, FALSE, referenceValue);
            PsychCopyOutDoubleArg(3, FALSE, realValue);
            PsychCopyOutDoubleArg(4, FALSE, (double) healthy);

            // Copy out a fake thread affinity mask of the PTB Thread if requested:
            if (opmode > 0)
                PsychCopyOutDoubleArg(5, FALSE, 0);
        }
    #endif

    return(PsychError_none);
}

PsychError GETSECSAllClocks(void)
{
    static char useString[] = "[GetSecsTime, WallTime, syncErrorSecs] = GetSecs('AllClocks' [, maxError=0.000020]);";
    //                          1            2         3                                       1
    static char synopsisString[] =
    "Return current time in seconds according to all supported clocks.\n\n"
    "'GetSecsTime' is the usual GetSecs() clock, as returned by GetSecs(), in the timebase "
    "used by all other Psychtoolbox functions, e.g., Screen('Flip') or PsychPortAudio timestamps.\n\n"
    "'WallTime' is real-world system time in an operating system specific timebase. This clock is "
    "expected to be subject to time adjustments by the system administrator or by automated mechanisms "
    "like NTP time adjustments. Useful for synchronizing clocks across multiple machines on a local "
    "network, as it is possible for this clock to get automatically corrected for drift.\n"
    "On Linux and OSX, the timebase is gettimeofday(), seconds since 1. January 1970 00:00:00 UTC, "
    "with about 1 microsecond precision.\n"
    "On MS-Windows the timebase is also UTC time, with about 1 millisecond granularity, measuring "
    "elapsed seconds since 1. January 1601 00:00:00 UTC via the GetSystemTimeAsFileTime() function.\n\n"
    "'syncErrorSecs' How tightly together did the returned clock times get queried? A measure of "
    "confidence as to how much all returned times actually denote the same point in physical time.\n\n"
    "The input argument 'maxError' allows to set an allowable upper bound to 'syncErrorSecs'. The "
    "default value is 20 microseconds. The function will try up to 10 times to get a result no worse "
    "than 'maxError', and output a warning if it doesn't manage, e.g., due to some severely "
    "overloaded or deficient system.";
    static char seeAlsoString[] = "";

    static psych_bool firstTimeWarning = TRUE;
    int maxRetries = 10;
    double getSecsClock, getSecsClock2;
    double maxError = 20 * 1e-6;                // Default to 20 usecs max clock mapping error.
    double wallClock;

    // All sub functions should have these two lines:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if (PsychIsGiveHelp()) { PsychGiveHelp(); return(PsychError_none); };

    // Check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(3));
    PsychErrorExit(PsychCapNumInputArgs(1));

    PsychCopyInDoubleArg(1, FALSE, &maxError);
    if (maxError < 0.000001)
        PsychErrorExitMsg(PsychError_user, "Invalid 'maxError' argument supplied. Lower than minimum allowed value of 1 microsecond.\n");

    // Repeat clock queries up to maxRetries times if they can't be completed within
    // maxError seconds and thereby "clock sync" between returned values is not good
    // enough:
    do {
        PsychGetAdjustedPrecisionTimerSeconds(&getSecsClock);
        wallClock = PsychGetWallClockSeconds();
        PsychGetAdjustedPrecisionTimerSeconds(&getSecsClock2);
    } while ((maxRetries-- > 0) && (getSecsClock2 - getSecsClock > maxError));

    // Warn on excessive error, at least once:
    if (firstTimeWarning && (getSecsClock2 - getSecsClock > 2 * maxError)) {
        firstTimeWarning = FALSE;
        printf("PTB-WARNING: GetSecs('AllClocks') sync margin %f secs > 2 times maxError %f even after multiple retries! System timing problems?1?\n",
               getSecsClock2 - getSecsClock, 2 * maxError);
        printf("PTB-WARNING: This one-time warning will not repeat, even on successive failure to reach good precision. Check your system.\n");
    }

    // Return arg 1 is as usual GetSecs time:
    PsychCopyOutDoubleArg(1, FALSE, getSecsClock);

    // Return arg 2 is Wall clock, system specific but usually related to real-world
    // time, e.g., UTC, NTP time, or Unix time, something that translates in normal
    // clocks like humans use it. On Unix (Linux and OSX) this is gettimeofday() Unix
    // time, system time, measured in seconds since 1. January 1970. On Windows this is
    // UTC time, measured in seconds since whatever. This clock is subject to changes
    // of system time, manually by the system administrator, or automatically by NTP
    // and similar clock sync services.
    PsychCopyOutDoubleArg(2, FALSE, wallClock);

    // Return arg3 = Confidence interval for sync between clock queries:
    PsychCopyOutDoubleArg(3, FALSE, getSecsClock2 - getSecsClock);

    return(PsychError_none);
}
