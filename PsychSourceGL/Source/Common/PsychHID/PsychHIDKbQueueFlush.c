/*
	PsychtoolboxGL/Source/Common/PsychHID/PsychHIDKbQueueFlush.c		
  
	PROJECTS: 
	
		PsychHID only.
  
	PLATFORMS:  
	
		All.  
  
	AUTHORS:
	
		rwoods@ucla.edu		rpw 
        mario.kleiner@tuebingen.mpg.de      mk
      
	HISTORY:
		8/23/07  rpw        Created.
  
	NOTES:
	
		The routines PsychHIDKbQueueCreate, PsychHIDKbQueueStart, PsychHIDKbQueueCheck, PsychHIDKbQueueStop
		and PsychHIDKbQueueRelease comprise a replacement for PsychHIDKbCheck, providing the following
		advantages:
		
			1) Brief key presses that would be missed by PsychHIDKbCheck are reliably detected
			2) The times of key presses are recorded more accurately
			3) Key releases are also recorded
		
		Requires Mac OS X 10.3 or later. The Matlab wrapper functions (KbQueueCreate, KbQueueStart,
		KbQueueCheck, KbQueueStop and KbQueueRelease screen away Mac OS X 10.2 and earlier, and the C code 
		does nothing to verify the Mac OS X version.
		
		Only a single device can be monitored at any given time. The deviceNumber can be specified only
		in the call to PsychHIDKbQueueCreate. The other routines then relate to that specified device. If
		deviceNumber is not specified, the first device is the default (like PyschHIDKbCheck). If
		PsychHIDKbQueueCreate has not been called first, the other routines will generate an error 
		message. Likewise, if PsychHIDKbQueueRelease has been called more recently than PsychHIDKbQueueCreate,
		the other routines will generate error messages.
		
		It is acceptable to cal PsychHIDKbQueueCreate at any time (e.g., to switch to a new device) without
		calling PsychKbQueueRelease.
		
		PsychHIDKbQueueCreate:
			Creates the queue for the specified (or default) device number
			No events are delivered to the queue until PsychHIDKbQueueStart is called
			Can be called again at any time
			
		PsychHIDKbQueueStart:
			Starts delivering keyboard events from the specified device to the queue
			
		PsychHIDKbQueueStop:
			Stops delivery of new keyboard events from the specified device to the queue.
			Data regarding events already queued is not cleared and can be recovered by PsychHIDKbQueueCheck
			
		PsychHIDKbQueueCheck:
			Obtains data about keypresses on the specified device since the most recent call to
			this routine or to PsychHIDKbQueueStart
			
			Clears all currently scored events (unscored events may still be in the queue)
			
		PsychHIDKbQueueFlush:
			Flushes unscored events from the queue and zeros all previously scored events
			
		PsychHIDKbQueueRelease:
			Releases queue-associated resources; once called, PsychHIDKbQueueCreate must be invoked
			before using any of the other routines
			
			This routine is called automatically at clean-up and can be omitted at the potential expense of
			keeping memory allocated unnecesarily
				

		---

*/

#include "PsychHID.h"

static char useString[]= "[navail] = PsychHID('KbQueueFlush' [, deviceIndex][, flushType=1])";
static char synopsisString[] = 
		"Flushes all scored and unscored keyboard events from a queue.\n"
		"Returns number of events 'navail' in keyboard event buffer before the flush takes place.\n"
		"If 'flushType' is 0, only the number of currently queued events will be returned.\n"
		"If 'flushType' is 1, only events returned by KbQueueCheck will be flushed. This is the default.\n"
		"If 'flushType' is 2, only events returned by KbQueueGetEvent will be flushed.\n"
		"If 'flushType' is 3, events returned by both KbQueueCheck and KbQueueGetEvent will be flushed.\n"
        "If 'flushType' & 4, only the number of key-press events with valid, mapped ASCII CookedKey field will be returned.\n"
        "PsychHID('KbQueueCreate') must be called before this routine.\n"
        "On Linux, the optional 'deviceIndex' is the index of the device whose queue should be flushed. "
        "If omitted, the default devices queue will be flushed. On other systems, the last queue will be flushed.\n";
        
static char seeAlsoString[] = "KbQueueCreate, KbQueueStart, KbQueueStop, KbQueueCheck, KbQueueRelease";
 
PsychError PSYCHHIDKbQueueFlush(void) 
{
    int deviceIndex, flushType;
    
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(2));

    deviceIndex = -1;
    PsychCopyInIntegerArg(1, kPsychArgOptional, &deviceIndex);

	flushType = 1;
    PsychCopyInIntegerArg(2, kPsychArgOptional, &flushType);
	
	// Return current count of contained events pre-flush:
	PsychCopyOutDoubleArg(1, FALSE, (double) PsychHIDAvailEventBuffer((PSYCH_SYSTEM != PSYCH_OSX) ? deviceIndex : 0, (flushType & 4) ? 1 : 0));

    if (flushType & 1) PsychHIDOSKbQueueFlush(deviceIndex);
    if (flushType & 2) PsychHIDFlushEventBuffer((PSYCH_SYSTEM != PSYCH_OSX) ? deviceIndex : 0);

    return(PsychError_none);	
}

#if PSYCH_SYSTEM == PSYCH_OSX
#include "PsychHIDKbQueue.h"

extern UInt32 modifierKeyState;
extern AbsoluteTime *psychHIDKbQueueFirstPress;
extern AbsoluteTime *psychHIDKbQueueFirstRelease;
extern AbsoluteTime *psychHIDKbQueueLastPress;
extern AbsoluteTime *psychHIDKbQueueLastRelease;
extern HIDDataRef hidDataRef;

extern pthread_mutex_t psychHIDKbQueueMutex;
    
void PsychHIDOSKbQueueFlush(int deviceIndex)
{
	if(!hidDataRef || !psychHIDKbQueueFirstPress || !psychHIDKbQueueFirstRelease || !psychHIDKbQueueLastPress || !psychHIDKbQueueLastRelease){
		PsychErrorExitMsg(PsychError_user, "Queue has not been created.");
	}
	
	// Drain the queue of any unprocessed events
	{
		IOHIDEventStruct event;
		AbsoluteTime zeroTime= {0,0};
		
		while(  (*hidDataRef->hidQueueInterface)->getNextEvent(hidDataRef->hidQueueInterface, &event, zeroTime, 0) == kIOReturnSuccess){
			if ((event.longValueSize != 0) && (event.longValue != NULL)) free(event.longValue);
		}
	}
	
	pthread_mutex_lock(&psychHIDKbQueueMutex);

	// Zero out the scored values
	{
		int i;
		for(i=0; i<256; i++){
			
			psychHIDKbQueueFirstPress[i].hi=0;
			psychHIDKbQueueFirstPress[i].lo=0;

			psychHIDKbQueueFirstRelease[i].hi=0;
			psychHIDKbQueueFirstRelease[i].lo=0;

			psychHIDKbQueueLastPress[i].hi=0;
			psychHIDKbQueueLastPress[i].lo=0;

			psychHIDKbQueueLastRelease[i].hi=0;
			psychHIDKbQueueLastRelease[i].lo=0;
		}
        
        modifierKeyState = 0;        
	}
	pthread_mutex_unlock(&psychHIDKbQueueMutex);
}

#endif
