/*

  	PsychToolbox3/Source/Linux/Base/PsychTimeGlue.h

	AUTHORS:

  	mario.kleiner@tuebingen.mpg.de		mk

	PLATFORMS: Linux Only

 
  	PROJECTS:

  	HISTORY:

  	2/20/06	mk		Wrote it.  
	1/03/09	mk		Add generic Mutex locking support as service to ptb modules.
	
  	DESCRIPTION:

	TO DO:

*/

//begin include once 
#ifndef PSYCH_IS_INCLUDED_PsychTimeGlue
#define PSYCH_IS_INCLUDED_PsychTimeGlue

#include "Psych.h"

void PsychWaitUntilSeconds(double whenSecs);
void PsychWaitIntervalSeconds(double seconds);
void PsychYieldIntervalSeconds(double seconds);
double	PsychGetKernelTimebaseFrequencyHz(void);
void PsychGetPrecisionTimerTicks(psych_uint64 *ticks);
void PsychGetPrecisionTimerTicksPerSecond(double *frequency);
void PsychGetPrecisionTimerTicksMinimumDelta(psych_uint32 *delta);
void PsychGetPrecisionTimerSeconds(double *secs);
void PsychGetAdjustedPrecisionTimerSeconds(double *secs);
void PsychGetPrecisionTimerAdjustmentFactor(double *factor);
void PsychSetPrecisionTimerAdjustmentFactor(double *factor);
void PsychInitTimeGlue(void);
void PsychExitTimeGlue(void);
void PsychEstimateGetSecsValueAtTickCountZero(void);
double PsychGetEstimatedSecsValueAtTickCountZero(void);
int	PsychInitMutex(psych_mutex* mutex);
int	PsychDestroyMutex(psych_mutex* mutex);
int PsychLockMutex(psych_mutex* mutex);
int PsychTryLockMutex(psych_mutex* mutex);
int PsychUnlockMutex(psych_mutex* mutex);
int PsychCreateThread(psych_thread* threadhandle, void* threadparams, void *(*start_routine)(void *), void *arg);
int PsychDeleteThread(psych_thread* threadhandle);
int PsychAbortThread(psych_thread* threadhandle);
void PsychTestCancelThread(psych_thread* threadhandle);
psych_threadid PsychGetThreadId(void);
int PsychIsThreadEqual(psych_thread threadOne, psych_thread threadTwo);
int PsychIsCurrentThreadEqualToId(psych_threadid threadId);
int PsychSetThreadPriority(psych_thread* threadhandle, int basePriority, int tweakPriority);
int PsychInitCondition(psych_condition* condition, const pthread_condattr_t* condition_attribute);
int PsychDestroyCondition(psych_condition* condition);
int PsychSignalCondition(psych_condition* condition);
int PsychBroadcastCondition(psych_condition* condition);
int PsychWaitCondition(psych_condition* condition, psych_mutex* mutex);
int PsychTimedWaitCondition(psych_condition* condition, psych_mutex* mutex, double maxwaittimesecs);
#define PsychIsMSVista() 0
int PsychIsCurrentThreadEqualToPsychThread(psych_thread threadhandle);
psych_uint64 PsychAutoLockThreadToCores(psych_uint64* curCpuMask);

// Linux specific: CLOCK_MONOTONIC time in seconds -- Usually the system uptime:
double PsychOSGetLinuxMonotonicTime(void);
double PsychOSMonotonicToRefTime(double monotonicTime);
//end include once
#endif
