/*
    PsychToolbox3/Source/OSX/Base/PsychTimeGlue.c

    AUTHORS:

    Allen.Ingling@nyu.edu           awi
    mario.kleiner.de@gmail.com      mk

    PLATFORMS: OSX Only

    PROJECTS:

    All

    HISTORY:

    1/20/03     awi     Wrote it.

    5/27/05     mk      Add while-loops around mach_wait_until, so no problems with interruptions...
                        New routine PsychWaitUntilSeconds() for waiting until a specific time.
    1/03/09     mk      Add generic Mutex locking support as service to ptb modules. Add PsychYieldIntervalSeconds().

    DESCRIPTION:

    TO DO:

*/

#include "Psych.h"
#include <mach/mach_time.h>
#include <sys/sysctl.h>
#include <sched.h>

#include <mach/mach_init.h>
#include <mach/task_policy.h>
#include <mach/thread_policy.h>

static double       precisionTimerAdjustmentFactor = 1;
static psych_bool   isKernelTimebaseFrequencyHzInitialized = FALSE;
static long double  kernelTimebaseFrequencyHz;

void PsychWaitUntilSeconds(double whenSecs)
{
    kern_return_t   waitResult;
    uint64_t        deadlineAbsTics;

    // Compute deadline for wakeup in mach absolute time units:
    deadlineAbsTics= (uint64_t) (kernelTimebaseFrequencyHz * ((long double) whenSecs));

    if (!(deadlineAbsTics > 0 && whenSecs > 0)) return;

    // Call mach_wait_unit in an endless loop, because it can fail with retcode>0.
    // In that case we just restart...
    while(mach_wait_until(deadlineAbsTics));
}

void PsychWaitIntervalSeconds(double delaySecs)
{
    long double     waitPeriodTicks;
    kern_return_t   waitResult;
    uint64_t        startTimeAbsTics, deadlineAbsTics;

    if (delaySecs <= 0) return;

    startTimeAbsTics = mach_absolute_time();
    waitPeriodTicks= kernelTimebaseFrequencyHz * delaySecs;
    deadlineAbsTics= startTimeAbsTics + (uint64_t) waitPeriodTicks;
    while(mach_wait_until(deadlineAbsTics));
}

/* PsychYieldIntervalSeconds() - Yield the cpu for given 'delaySecs'
 *
 * PsychYieldIntervalSeconds() differs from PsychWaitIntervalSeconds() in that
 * it is supposed to release the cpu to other threads or processes for *at least*
 * the given amount of time 'delaySecs', instead of *exactly* 'delaySecs'.
 *
 * If one wants to wait an exact amount of time, one uses PsychWaitIntervalSeconds().
 * If one just "has nothing to do" for some minimum amount of time, and wants to
 * play nice to other threads/processes and exact timing is not crucial, then
 * this is the routine of choice. Typical use is within polling loops, where one
 * wants to pause between polling cycles and it doesn't matter if the pause takes
 * a bit longer.
 *
 * A 'delaySecs' of <= zero will just release the cpu for the remainder of
 * the current scheduling timeslice. If you don't know what to do, choose a
 * zero setting.
 *
 */
void PsychYieldIntervalSeconds(double delaySecs)
{
    if (delaySecs <= 0) {
        // Yield cpu for remainder of this timeslice:
        sched_yield();
    }
    else {
        // On OS/X we use standard wait ops - they're good for us:
        PsychWaitIntervalSeconds(delaySecs);
    }
}

double	PsychGetKernelTimebaseFrequencyHz(void)
{
    long double                 clockPeriodNSecs;
    mach_timebase_info_data_t   tbinfo;

    if(!isKernelTimebaseFrequencyHzInitialized){
        // Retrieve the mach absolute time timebase.  The kernel expresses the period in two integers, the ratio of which is the clock period.
        mach_timebase_info(&tbinfo);

        // Calculate the mach timebase period from values reported from the mach kernel.
        clockPeriodNSecs = ((long double) tbinfo.numer) / ((long double) tbinfo.denom);

        // Convert the mach timebase period from awkward units into frequency in Hz.
        // Frequency in Hz is a convenient form because it makes converting from a period in seconds into a period in mach timebase units easy:
        //  time_interval_in_mach_units= time_interval_in_seconds * clockFrequencyHz;
        kernelTimebaseFrequencyHz = 1000000000.0 / clockPeriodNSecs;
        isKernelTimebaseFrequencyHzInitialized=TRUE;
    }

    return((double)kernelTimebaseFrequencyHz);
}

/* Called at Module init time: */
void PsychInitTimeGlue(void)
{
    // Force initialization of the kernel timebase frequency:
    PsychGetKernelTimebaseFrequencyHz();
    return;
}

/* Called at module shutdown/jettison time: */
void PsychExitTimeGlue(void)
{
    // Nothing to do on OS/X yet:
    return;
}

void PsychGetPrecisionTimerTicks(psych_uint64 *ticks)
{
    *ticks=AudioGetCurrentHostTime();
}

void PsychGetPrecisionTimerTicksPerSecond(double *frequency)
{
    *frequency = (double)AudioGetHostClockFrequency();
}

void PsychGetPrecisionTimerTicksMinimumDelta(psych_uint32 *delta)
{
    *delta=AudioGetHostClockMinimumTimeDelta();
}

void PsychGetPrecisionTimerSeconds(double *secs)
{
    *secs= mach_absolute_time() / kernelTimebaseFrequencyHz;
}

void PsychGetAdjustedPrecisionTimerSeconds(double *secs)
{
    double  rawSecs;

    PsychGetPrecisionTimerSeconds(&rawSecs);
    *secs=rawSecs * precisionTimerAdjustmentFactor;
}

void PsychGetPrecisionTimerAdjustmentFactor(double *factor)
{
    *factor=precisionTimerAdjustmentFactor;
}

void PsychSetPrecisionTimerAdjustmentFactor(double *factor)
{
    precisionTimerAdjustmentFactor=*factor;
}

/*
    PsychEstimateGetSecsValueAtTickCountZero()
 */
void PsychEstimateGetSecsValueAtTickCountZero(void)
{
    // Dead as of PTB 3.0.12, as Apple deprecated TickCount(), and it seems
    // GetSecs() zero == TickCount zero, so no point here anymore.
    return;
}

double PsychGetEstimatedSecsValueAtTickCountZero(void)
{
    return(0.0);
}

/* Init a Mutex: */
int	PsychInitMutex(psych_mutex* mutex)
{
    return(pthread_mutex_init(mutex, NULL));
}

/* Deinit and destroy a Mutex: */
int	PsychDestroyMutex(psych_mutex* mutex)
{
    return(pthread_mutex_destroy(mutex));
}

/* Lock a Mutex, blocking until mutex is available if it isn't available: */
int PsychLockMutex(psych_mutex* mutex)
{
    return(pthread_mutex_lock(mutex));
}

/* Try to lock a Mutex, returning immediately, with a return code that tells if mutex could be locked or not: */
int PsychTryLockMutex(psych_mutex* mutex)
{
    return(pthread_mutex_trylock(mutex));
}

/* Unlock a Mutex: */
int PsychUnlockMutex(psych_mutex* mutex)
{
    return(pthread_mutex_unlock(mutex));
}

/* Create a parallel thread of execution, invoke its main routine: */
int PsychCreateThread(psych_thread* threadhandle, void* threadparams, void *(*start_routine)(void *), void *arg)
{
    // threadparams not yet used, this line just to make compiler happy:
    (void) threadparams;

    // Return result code of pthread_create - We're a really thin wrapper around this Posix call:
    return( pthread_create(threadhandle, NULL, start_routine, arg) );
}

/* Join a parallel thread - Wait for its termination, then return its result code: */
int PsychDeleteThread(psych_thread* threadhandle)
{
    // Join on the thread, wait for termination:
    int rc = pthread_join(*threadhandle, NULL);
    // Null out now invalid thread handle of dead thread:
    *threadhandle = NULL;
    // Return return code of joined thread:
    return(rc);
}

/* Send abort request to thread: */
int PsychAbortThread(psych_thread* threadhandle)
{
    return( pthread_cancel(*threadhandle) );
}

/* Check for abort request to thread: Exit thread gracefully if abort requested: */
void PsychTestCancelThread(psych_thread* threadhandle)
{
    // threadhandle unused on POSIX: This line just to make compiler happy:
    (void) threadhandle;

    // Test for cancellation, cancel if so:
    pthread_testcancel();
}

/* Return handle of calling thread: */
psych_threadid PsychGetThreadId(void)
{
    return( pthread_self() );
}

/* Check if two given thread handles do refer to the same thread: */
int PsychIsThreadEqual(psych_thread threadOne, psych_thread threadTwo)
{
    return( pthread_equal(threadOne, threadTwo) );
}

/* Check if current (invoking) thread has an id equal to given threadid: */
int PsychIsCurrentThreadEqualToId(psych_threadid threadId)
{
    return( pthread_equal(PsychGetThreadId(), threadId) );
}

/* Check if current (invoking) thread is equal to given threadhandle: */
int PsychIsCurrentThreadEqualToPsychThread(psych_thread threadhandle)
{
    return( pthread_equal(PsychGetThreadId(), threadhandle) );
}

// Helper for PsychSetThreadPriority(): Setup of Mach RT scheduling.
// We want / "promise to not use more" than "computation" cycles out of every "period" cycles.
// Once we started execution, we want to finish our "computation" cycles within at most "constraint" cycles.
// We allow / or don't allow to be "isPreemptible" preempted - to finish our "computation" cycles split up into
// multiple pieces, but finishing within at most "constraint" cycles.
int set_realtime(task_t threadID, int period, int computation, int constraint, psych_bool isPreemptible) {
    struct thread_time_constraint_policy ttcpolicy;
    int ret;

    // Set realtime scheduling with following parameters:
    ttcpolicy.period = period;
    ttcpolicy.computation = computation;
    ttcpolicy.constraint = (constraint >= computation) ? constraint : computation;
    ttcpolicy.preemptible = (isPreemptible) ? 1 : 0;

    ret = thread_policy_set(threadID, THREAD_TIME_CONSTRAINT_POLICY, (thread_policy_t) &ttcpolicy, THREAD_TIME_CONSTRAINT_POLICY_COUNT);
    return(ret);
}

/* Change priority for thread 'threadhandle', or for the calling thread if 'threadhandle' == NULL.
 * threadhandle == 0x1 means "Main Psychtoolbox thread" and may incur special treatment.
 * 'basePriority' can be 0 for normal scheduling, 1 for higher priority and 2 for highest priority.
 * 'tweakPriority' modulates more fine-grained within the category given by 'basepriority'. It
 * can be anywhere between 0 and some big value where bigger means more priority.
 *
 * Returns zero on success, non-zero on failure to set new priority.
 */
int PsychSetThreadPriority(psych_thread* threadhandle, int basePriority, int tweakPriority)
{
    int                         rc = 0;
    pthread_t                   thread;
    int                         kernError;
    task_t                      threadID;
    thread_policy_t             threadPolicy;
    mach_msg_type_number_t      policyCount, policyCountFilled;
    boolean_t                   isDefault;

    if ((NULL != threadhandle) && ((psych_thread*) 0x1 != threadhandle)) {
        // Retrieve thread handle of thread to change:
        thread = *threadhandle;
    }
    else {
        // Retrieve handle of calling thread:
        thread = pthread_self();
    }

    // Map Posix thread handle to Mach thread handle:
    threadID = pthread_mach_thread_np(thread);

    // Get timebase:
    double ticksPerSec = PsychGetKernelTimebaseFrequencyHz();
    double baseQuantum = 0.010;

    // tweakPriority <= 0 -> 10% cpu. Can go up to 90% at level >=8 in 10% increments.
    if (tweakPriority < 0) tweakPriority = 0;
    if (tweakPriority > 8) tweakPriority = 8;

    switch(basePriority) {
        case 0:	// Normal priority: Drop to standard scheduling.
            threadPolicy = (thread_policy_t) malloc(sizeof(thread_standard_policy_data_t));
            policyCount  = THREAD_STANDARD_POLICY_COUNT;
            policyCountFilled = policyCount;
            isDefault = TRUE;
            kernError = thread_policy_get(threadID, THREAD_STANDARD_POLICY, threadPolicy, &policyCountFilled, &isDefault);
            if (kernError == 0) kernError = thread_policy_set(threadID, THREAD_STANDARD_POLICY, threadPolicy, policyCountFilled);
            free(threadPolicy);
            rc = (int) kernError;
        break;

        case 1: // High priority: Up to 90% cpu utilization, but preemptible for urgent tasks, with an allowable total time to completion of baseQuantum.
            // This basically says: "I am more important than bog-standard threads, and i want to have x msecs of 10 msecs very 10 msecs, but i don't care
            // about startup delay (reaction times) or interruptions, as long as i don't lose more than 10 msecs. Good for high priority compute tasks
            // with lots of wiggle room wrt. when stuff happens, e.g., some data producer thread or i/o thread which needs to deliver/handle a certain
            // amount of data processing/shuffling/io within a certain time quantum, because it is feeding some other realtime thread or hw process,
            // due to things like intermediate fifo buffering, itself can tolerate a certain lag.
            // This may become useful in the future for i/o functions in IOPort/PsychHID, movie playback/recording helper threads etc...
            rc = set_realtime(threadID, baseQuantum * ticksPerSec, (((double) tweakPriority + 1) / 10) * baseQuantum * ticksPerSec, baseQuantum * ticksPerSec, TRUE);
        break;

        case 2: // Realtime priority: Up to (tweakPriority + 1) msecs out of 10 msecs of *uninterrupted* computation (non-preemptible after start).
            // However, after the thread becomes runnable, its actual start of uninterrupted execution can be delayed by up to 1 msec, e.g., if more
            // important (basePriority 10) threads are executing, or a high priority or lower priority thread needs some computation time.
            // This is our most common use-case: Most of our realtime threads are completely triggered (= reactive to) by external hardware input events.
            // They wait on the arrival of some external event, e.g., a user key press or mouse click, some trigger signal from some I/O device like
            // response box, serial port or parallel port, DAQ board etc., or for some timer going off at a certain time. Most often they have to respond
            // to some trigger event by either executing some action, or by simply timestamping the event, like a button press of a subject, or some TTL
            // trigger from some equipment. Executing the actual action, or timestamping, or storing the received data in some queue, is usually fast,
            // the computation finishes quickly. As timestamping or external hardware control can be involved, we don't want to get preempted once running,
            // to avoid impairing precision of timestamps or clock-sync algorithms or hw control actions. However for typical neuro-science experiments,
            // we can tolerate a random time delay (or imprecision in acquired timestamps) of 1 msec.
            // Typical consumers of this setup: IOPort, PsychHID, Movie playback or video capture high-level control.
            rc = set_realtime(threadID, baseQuantum * ticksPerSec, (((double) tweakPriority + 1) / 10) * baseQuantum * ticksPerSec, (((double) tweakPriority + 1 + 1) / 10) * baseQuantum * ticksPerSec, FALSE);
        break;

        case 10:  // Critical priority: Up to (tweakPriority + 1) msecs out of 10 msecs of *uninterrupted* computation (non-preemptible after start),
            // must run as soon as possible and then complete without distraction. This is good for timestamping operations that must not be interrupted
            // in the wrong moment, because that would impair timestamps significantly, and for time-based triggering of execution of operations with
            // the highest possible timing precision.
            // Out main client of this is currently the OpenGL flipperThread used by Screen for async flip scheduling and timestamping, and for
            // frame-sequential stereo fallback. For those apps, uninterrupted low latency is crucial. flipperThread uses tweakPriority == 2, so could
            // run for up to 3 msecs uninterrupted, something it usually won't do (closer to << 1 msec is expected), but can do in a worst case scenario,
            // where various workarounds for broken GPU drivers are active and screen resolution/refresh rate settings are especially suboptimal.
            // The other client is video refresh rate calibration during Screen('GetFlipInterval') active calibration or during Screen('Openwindow')
            // default calibration.
            //
            // Future clients may be found in the IOPort async-task framework for highly timing sensitive i/o operations.
            rc = set_realtime(threadID, baseQuantum * ticksPerSec, (((double) tweakPriority + 1) / 10) * baseQuantum * ticksPerSec, 0, FALSE);
        break;

        default:
            printf("PTB-CRITICAL: In call to PsychSetThreadPriority(): Invalid/Unknown basePriority %i provided!\n", basePriority);
            rc = 2;
    }

    // Try to apply new priority and scheduling method:
    if (rc != 0) {
        printf("PTB-WARNING: In call to PsychSetThreadPriority(): Failed to set new basePriority %i, tweakPriority %i, effective %i [%s] for thread %p provided!\n",
                basePriority, tweakPriority, tweakPriority, (basePriority > 0) ? "REALTIME" : "NORMAL", (void*) threadhandle);
        printf("PTB-WARNING: This can lead to timing glitches and odd performance behaviour.\n");
    }

    // rc is either zero for success, or 2 for invalid arg, or some other non-zero failure code:
    return(rc);
}

/* Assign a name to a thread, for debugging: */
void PsychSetThreadName(const char *name)
{
    // OSX interface only allows assigning name to current thread, different to Linux:
    pthread_setname_np(name);
}

/* Initialize condition variable:
 * CAUTION: Use of condition_attribute is non-portable! Code using it will not work properly
 * on MS-Windows as this attribute is unsupported there! Pass NULL for this argument for
 * portable operation!
 */
int PsychInitCondition(psych_condition* condition, const pthread_condattr_t* condition_attribute)
{
    return(pthread_cond_init(condition, condition_attribute));
}

/* Destroy condition variable: */
int PsychDestroyCondition(psych_condition* condition)
{
    return(pthread_cond_destroy(condition));
}

/* Signal/wakeup exactly one thread waiting on the given condition variable: */
int PsychSignalCondition(psych_condition* condition)
{
    return(pthread_cond_signal(condition));
}

/* Signal/Wakeup all threads waiting on the given condition variable:
 * CAUTION: Use of this function is non-portable to MS-Windows for now! Code
 * using it will malfunction if used on MS-Windows!
 */
int PsychBroadcastCondition(psych_condition* condition)
{
    return(pthread_cond_broadcast(condition));
}

/* Atomically release the 'mutex' lock and go to sleep, waiting for the 'condition' variable
 * being signalled, then waking up and trying to re-lock the 'mutex'. Will return with
 * mutex locked.
 */
int PsychWaitCondition(psych_condition* condition, psych_mutex* mutex)
{
    return(pthread_cond_wait(condition, mutex));
}

/* Atomically release the 'mutex' lock and go to sleep, waiting for the 'condition' variable
 * being signalled, then waking up and trying to re-lock the 'mutex'. Will return with
 * mutex locked.
 *
 * Like PsychWaitCondition, but function will timeout if it fails being signalled before
 * timeout interval 'maxwaittimesecs' expires. In any case, it will only return after
 * reacquiring the mutex. It will retun zero on successfull wait, non-zero (ETIMEDOUT) if
 * timeout was triggered without the condition being signalled.
 */
int PsychTimedWaitCondition(psych_condition* condition, psych_mutex* mutex, double maxwaittimesecs)
{
    struct timespec abstime;
    struct timeval gtod_time;
    double tnow;

    // Convert relative wait time to absolute system time. As pthread_cond_timedwait()
    // uses gettimeofday() time as reference, we can't query or regular GetSecs clock,
    // but need to use gettimeofday():
    gettimeofday(&gtod_time, NULL);

    // Convert gtod_time and maxwaittimesecs into timespec format, add it...
    abstime.tv_sec  = (time_t) maxwaittimesecs + gtod_time.tv_sec;
    abstime.tv_nsec = (long) (((double) maxwaittimesecs - (double) ((time_t) maxwaittimesecs)) * (double) (1e9));
    abstime.tv_nsec+= (long) (gtod_time.tv_usec * 1000);

    // ... (Re-)split into seconds and nanoseconds:
    while (abstime.tv_nsec >= 1e9) {
        abstime.tv_nsec-= 1e9;
        abstime.tv_sec+= 1;
    }

    // Perform wait for signalled condition with a timeout at absolute system time abstime:
    return(pthread_cond_timedwait(condition, mutex, &abstime));
}

/* Set thread affinity mask of calling thread to the modules global cpuMask:
 * 
 * 'curCpuMask' is an in/out pointer. If NULL, it is completely ignored. If non-NULL,
 * the target variable of the pointer will contain the new cpu mask after a change
 * of mask. If the target variable already contains a valid (non-zero) current cpu mask
 * and it matches the new target cpu mask, then the actual mask change is omitted, assuming it
 * is redundant, thereby saving some system call overhead.
 *
 * Threads can avoid redundant switches if they keep track of their current cpu mask
 * by caching it in the *curCpuMask pointer target. They can pass in a zero value if
 * unsure, or a NULL pointer if they are neither interested in caching, nor in the old
 * or new value.
 *
 * Returns the old pre-switch affinity mask as a 64-Bit bitfield.
 * Comparing the return value (previous mask) with the *curCpuMask value (new mask)
 * allows the caller to check if the affinity mask was actually changed, resulting
 * in a possible thread migration.
 *
 * If this function is called without the time lock held, ie., from outside
 * of other timeglue functions, a small race condition exists which may cause
 * deferred updated to the real new affinity mask due to 
 *
 */
psych_uint64 PsychAutoLockThreadToCores(psych_uint64* curCpuMask)
{
    // No-op on OSX:
    return(INT64_MAX);
}

/* Query / derive / return OSX minor version from Darwin kernel major version.
 * This is a makeshift replacement for Gestalt(), which was sadly deprecated by
 * the iPhone company. It only gives us the x in OSX 10.x.y, but that's usually
 * all we need.
 */
int PsychGetOSXMinorVersion(void)
{
    int mib[2] = { CTL_KERN, KERN_OSRELEASE };
    int minorVersion;
    char tempStr[256];
    size_t tempStrSize = sizeof(tempStr);

    // Query kernel version string:
    if (sysctl(mib, 2, tempStr, &tempStrSize, NULL, 0)) {
        printf("PTB-WARNING: Could not query Darwin kernel version! This will end badly...\n");
    }

    // Parse out major version: That - 4 == OSX minor version:
    if (1 != sscanf(tempStr, "%i", &minorVersion)) {
        printf("PTB-WARNING: Could not parse Darwin kernel major version! This will end badly...\n");
    }

    minorVersion = minorVersion - 4;

    // Return minorVersion of the OSX version number: 10.minorVersion
    return(minorVersion);
}

/* Report official support status for this operating system release.
 * The string "Supported" means supported.
 * Other strings describe lack of support.
 */
const char* PsychSupportStatus(void)
{
    // Operating system minor version:
    int osMinor;

    // Init flag to -1 aka unknown:
    static int  isSupported = -1;
    static char statusString[256];

    if (isSupported == -1) {
        // First call: Do the query!

        // Query OS/X version:
        osMinor = PsychGetOSXMinorVersion();

        // Only OSX 10.11 is officially supported:
        isSupported = (osMinor == 11) ? 1 : 0;

        if (isSupported) {
            sprintf(statusString, "OSX 10.%i Supported and tested.", osMinor);
        }
        else {
            sprintf(statusString, "OSX version 10.%i is not supported or tested anymore.", osMinor);
        }
    }

    return(statusString);
}
