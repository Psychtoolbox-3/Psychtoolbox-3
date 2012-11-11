/*
	PsychtoolboxGL/Source/Common/PsychHID/PsychHIDKbQueueRelease.c		
  
	PROJECTS: 
	
		PsychHID only.
  
	PLATFORMS:  
	
		All.  
  
	AUTHORS:
	
		rwoods@ucla.edu		rpw 
        mario.kleiner@tuebingen.mpg.de      mk
      
	HISTORY:
		8/19/07  rpw		Created.
		8/23/07  rpw        Added PsychHIDQueueFlush to documentation; removed call to PsychHIDVerifyInit()
  
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

static char useString[]= "PsychHID('KbQueueRelease' [, deviceIndex])";
static char synopsisString[] = 
		"Releases queue and other resources allocated by PsychHID('KbQueueCreate').\n"
        "On Linux and Windows, the optional 'deviceIndex' is the index of the device whose queue should be released. "
        "If omitted, the default devices queue will be released. On other systems, the last queue will be released.\n";
        
static char seeAlsoString[] = "KbQueueCreate, KbQueueStart, KbQueueStop, KbQueueCheck, KbQueueFlush";

PsychError PSYCHHIDKbQueueRelease(void) 
{
    int deviceIndex;
    
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    deviceIndex = -1;
    PsychCopyInIntegerArg(1, kPsychArgOptional, &deviceIndex);
    
    PsychHIDOSKbQueueRelease(deviceIndex);

    return(PsychError_none);
}

#if PSYCH_SYSTEM == PSYCH_OSX
#include "PsychHIDKbQueue.h"

extern AbsoluteTime *psychHIDKbQueueFirstPress;
extern AbsoluteTime *psychHIDKbQueueFirstRelease;
extern AbsoluteTime *psychHIDKbQueueLastPress;
extern AbsoluteTime *psychHIDKbQueueLastRelease;
extern HIDDataRef hidDataRef;
extern pthread_mutex_t psychHIDKbQueueMutex;
extern CFRunLoopRef psychHIDKbQueueCFRunLoopRef;
extern pthread_t psychHIDKbQueueThread;

void PsychHIDOSKbQueueRelease(int deviceIndex)
{
    pthread_mutex_lock(&psychHIDKbQueueMutex);    
    
	// Remove the source from the CFRunLoop so queue transitions from empty
	// to non-empty cannot trigger a new callout on the CFRunLoop thread
	if(psychHIDKbQueueCFRunLoopRef){
		if(CFRunLoopContainsSource(psychHIDKbQueueCFRunLoopRef, hidDataRef->eventSource, kCFRunLoopDefaultMode))
			CFRunLoopRemoveSource(psychHIDKbQueueCFRunLoopRef, hidDataRef->eventSource, kCFRunLoopDefaultMode);
	}

	// Stop, drain, dispose, release the queue
	if(hidDataRef){
		IOHIDQueueInterface **queue=(hidDataRef->hidQueueInterface);
		if(queue){
			HRESULT result=(*queue)->stop(queue);
			{
				IOHIDEventStruct event;
				AbsoluteTime zeroTime= {0,0};
				while(  (*hidDataRef->hidQueueInterface)->getNextEvent(hidDataRef->hidQueueInterface, &event, zeroTime, 0) == kIOReturnSuccess){
					if ((event.longValueSize != 0) && (event.longValue != NULL)) free(event.longValue);
				}
			}
			result = (*queue)->dispose(queue);
			(*queue)->Release(queue);
		}
		if(hidDataRef->hidElementDictionary) CFRelease(hidDataRef->hidElementDictionary);
		free(hidDataRef);
		hidDataRef=NULL;
	}

    pthread_mutex_unlock(&psychHIDKbQueueMutex);    
    
	// Stop the CFRunLoop, which will allow its thread to exit
	// The mutex will be automatically unlocked and destroyed by the CFRunLoop thread
	// so it isn't even declared in this routine
	if(psychHIDKbQueueCFRunLoopRef){
		CFRunLoopStop(psychHIDKbQueueCFRunLoopRef);

        pthread_join(psychHIDKbQueueThread, NULL);
        psychHIDKbQueueThread = NULL;

		CFRelease(psychHIDKbQueueCFRunLoopRef);
		psychHIDKbQueueCFRunLoopRef=NULL;

		// Destroy the mutex
		pthread_mutex_destroy(&psychHIDKbQueueMutex);
	}
	
	// Free and null the memory where presses and releases are stored
	if(psychHIDKbQueueFirstPress){
		free(psychHIDKbQueueFirstPress);
		psychHIDKbQueueFirstPress=NULL;
	}
	if(psychHIDKbQueueFirstRelease){
		free(psychHIDKbQueueFirstRelease);
		psychHIDKbQueueFirstRelease=NULL;
	}
	if(psychHIDKbQueueLastPress){
		free(psychHIDKbQueueLastPress);
		psychHIDKbQueueLastPress=NULL;
	}
	if(psychHIDKbQueueLastRelease){
		free(psychHIDKbQueueLastRelease);
		psychHIDKbQueueLastRelease=NULL;
	}

	// Release kbqueue event buffer:
	PsychHIDDeleteEventBuffer(0);
}

#endif
