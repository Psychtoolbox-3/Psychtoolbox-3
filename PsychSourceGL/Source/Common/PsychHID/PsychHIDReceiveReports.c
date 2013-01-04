/*
	PsychToolbox/Source/Common/PsychHID/PsychHIDReceiveReports.c		

	PROJECTS: PsychHID

	PLATFORMS:  All

	AUTHORS:
	denis.pelli@nyu.edu                 dgp
    mario.kleiner@tuebingen.mpg.de      mk

	HISTORY:
	4/7/05  dgp	Wrote it, based on PsychHIDGetReport.c

	READ:
	bugs in mac os x retrieval of reports.
	http://lists.apple.com/archives/usb/2004/Jul/msg00003.html
	$900 USB bus analyzer
	http://lists.apple.com/archives/usb/2003/Dec/msg00011.html
	guido kolber reports problem with getting port 
	http://lists.apple.com/archives/usb/2004/Nov/msg00062.html
	usb mailing list:
	http://lists.apple.com/mailman/listinfo/usb
	header file with helpful comments:
	http://darwinsource.opendarwin.org/10.3.6/IOHIDFamily-86.21/IOHIDLib/IOHIDDeviceClass.h
	my query on the Apple usb list:
	http://lists.apple.com/archives/usb/2005/Apr/msg00050.html
	HID documentation:
	http://developer.apple.com/documentation/DeviceDrivers/Conceptual/HID/
	How to set up notification callback when a device is attached
	or removed:
	file:///Developer/Examples/IOKit/usb/Another%20USB%20Notification%20Example/USBNotificationExample.c
	http://developer.apple.com/samplecode/USBPrivateDataSample/USBPrivateDataSample.html
 
	NOTE ON GETREPORT: 
	4 April 2005. Implementing SetReport, to send a report to the
	device, was straightforward, just a matter of calling Apple's HID
	Utilities function HIDSetReport. Implementing "GetReport", to
	receive a report from the device, was arduous. I quickly
	implemented a call to Apple's HIDGetReport, and it worked with my
	keyboard and mouse, but it never worked with my PMD-1208FS. The
	GetReport request always resulted in a "stall", indicating that the
	PMD-1208FS was rejecting the request. This baffled me. Last friday,
	the firmware engineer, Nick, who programmed the PMD-1208FS,
	explained to me that I shouldn't use the GetReport command because
	that goes out on the USB bus, making a request to the device, which
	is too slow a way of communicating, so the device doesn't support
	it. This surprised me. I had supposed until then that GetReport
	merely grabbed whatever report had already been received by the
	host. Following a suggestion of Allen Ingling's I looked again in
	Apple's HID Manager documentation and discovered that the command
	one should use is in fact setInterruptReportHandlerCallback.
	However, the text explaining what this does and how to call it is
	very skimpy, not enough for me to figure it out. Fortunately
	googling setInterruptReportHandlerCallback turned up several helpful
	discussions on apple's usb developer forum. Remarkably, last
	november, the CEO of Code Mercenaries, who make IOWarrior, presented
	an abbreviated listing of exactly the code I needed, about five
	lines, asking why it broke in Mac OS 10.3.6. With apple's answer I
	was able to finish writing my code. I then had to read about a Mach
	feature in Mac OS X called a CFRunLoop, since
	setInterruptReportHandlerCallback implements the callback by
	installing a "source" in that loop. After about a day's fiddling it
	now works.
 */

/*
file:///Developer/ADC%20Reference%20Library/documentation/Darwin/Reference/IOKit/IOHIDLib/Classes/IOHIDDeviceInterface122/CompositePage.html#//apple_ref/doc/compositePage/c/func/setInterruptReportHandlerCallback
 
IOHIDDeviceInterface::
 
 IOHIDReportCallbackFunction

typedef void (*IOHIDReportCallbackFunction) (
											 void * target, 
											 IOReturn result, 
											 void * refcon, 
											 void * sender, 
											 UInt32 bufferSize); 


Type and arguments of callout C function that is used when a completion routine is called, see IOHIDLib.h:setReport().


Parameter Descriptions

target
void * pointer to your data, often a pointer to an object.

result
Completion result of desired operation.

refcon
void * pointer to more data.

sender
Interface instance sending the completion routine.

bufferSize
Size of the buffer received upon completion.
*/

#include "PsychHID.h"

// MK TODO FIXME: Question to former self: Why the hell does this exist?
void mexMakeMemoryPersistent(void *ptr);
#ifndef PTBOCTAVE
void mxFree(void *ptr);
void *mxCalloc(size_t n, size_t size);
#endif

// internal prototypes
PsychError ReceiveReportsStop(int deviceIndex);  // function is below.
void CountReports(char *string);

typedef struct ReportStruct{
	int deviceIndex;
	long int error;
	psych_uint32 bytes;
	double time;
	//int type; // 1=input, 2=output, 3=feature
	struct ReportStruct *next;
	psych_uint8 *report;
} ReportStruct;

static psych_bool firstTimeInit = TRUE;         // TRUE at PsychHID load & init time and after shutdown. FALSE during operation.

// These are out here for easy access by several routines in this file.
static psych_bool oneShotRealloc = FALSE;
static psych_bool ready[MAXDEVICEINDEXS];
static psych_bool optionsPrintReportSummary=0;	// options.print: Enable diagnostic print of report by ReportCallback.
static psych_bool optionsPrintCrashers=0;		// options.printCrashers
static psych_bool optionsConsistencyChecks=0;	// options.consistencyChecks
static int optionsMaxReports=10000;			// options.maxReports
static int optionsMaxReportSize=65;			// options.maxReportSize
static double optionsSecs=0.010;			// options.secs

// These are out here for easy access by my report callback function: ReportCallback.
static ReportStruct *freeReportsPtr[MAXDEVICEINDEXS];   // Per device linked list of free hid input reports.
static ReportStruct *deviceReportsPtr[MAXDEVICEINDEXS]; // Per device linked list of filled hid input reports.

static ReportStruct *allocatedReports[MAXDEVICEINDEXS]; // Per device linked list storage - list is tightly packed in memory.
static psych_bool reportsHaveBeenAllocated[MAXDEVICEINDEXS]; // Allocated flag.
static int MaxDeviceReports[MAXDEVICEINDEXS];           // Per device number of total reports.
static int MaxDeviceReportSize[MAXDEVICEINDEXS];        // Per device max size of each report.
psych_uint8 * reportData[MAXDEVICEINDEXS];              // Per device buffer for all reports databuffers, tightly packed.

// Set by PsychHIDSetReport, read by ReportCallback solely for the optionsPrintReportSummary.
double AInScanStart=0;

#if PSYCH_SYSTEM == PSYCH_OSX

#include <IOKit/HID/IOHIDLib.h>

void CheckRunLoopSource(int deviceIndex,char *caller,int line);
static CFRunLoopSourceRef source[MAXDEVICEINDEXS]; 
static CFStringRef        myRunLoopMode=NULL;		// CFSTR("myMode");

void ReportCallback(void *target,IOReturn result,void *refcon,void *sender,psych_uint32 bufferSize); // function is below.

void ReportCallback(void *target,IOReturn result,void *refcon,void *sender,psych_uint32 bufferSize)
{
	int deviceIndex,i,n,m;
	unsigned char *ptr;
	ReportStruct *r;
	
	CountReports("ReportCallback beginning.");
	
	deviceIndex = (long int) refcon;
	if(deviceIndex < 0 || deviceIndex >= MAXDEVICEINDEXS) {
		printf("ReportCallback received out-of-range deviceIndex %d. Aborting.\n", deviceIndex);
		return;
	}
	
	// take report from free list.
	if(freeReportsPtr[deviceIndex] == NULL){
		// Darn. We're full. It might be elegant to discard oldest report, but for now, we'll just ignore the new one.
		printf("ReportCallback warning. No more free reports. Discarding new report.\n");
		return;
	}
    
	r = freeReportsPtr[deviceIndex];
	freeReportsPtr[deviceIndex] = r->next;
	r->next=NULL;
	
	// install report into the device's list.
	r->next = deviceReportsPtr[deviceIndex];
	deviceReportsPtr[deviceIndex] = r;
	
	// fill in the rest of the report struct
	r->error=result;
	r->bytes=bufferSize;
	r->deviceIndex=deviceIndex;
	ptr=target;
    
    // Clamp amount of returned data to global and per-device limit:
    if (bufferSize > MAXREPORTSIZE) bufferSize = MAXREPORTSIZE;
    if (bufferSize > MaxDeviceReportSize[deviceIndex]) bufferSize = MaxDeviceReportSize[deviceIndex];
    
    // Copy data:
	for(i = 0; i < bufferSize; i++) r->report[i] = *(ptr+i);
    
	PsychGetPrecisionTimerSeconds(&r->time);
	if(optionsPrintReportSummary){
		// print diagnostic summary of the report
		int serial;
		
		serial=r->report[62]+256*r->report[63]; // 32-bit serial number at end of AInScan report from PMD-1208FS
		printf("Got input report %4d: %2ld bytes, dev. %d, %4.0f ms. ",serial,(long)r->bytes,deviceIndex,1000*(r->time-AInScanStart));
		if(r->bytes>0){
			printf(" report ");
			n=r->bytes;
			if(n>6)n=6;
			for(i=0;i<n;i++)printf("%3d ",(int)r->report[i]);
			m=r->bytes-2;
			if(m>i){
				printf("... ");
				i=m;
			}
			for(;i<r->bytes;i++)printf("%3d ",(int)r->report[i]);
		}
		printf("\n");
	}
	CountReports("ReportCallback end.");
	return;
}

PsychError ReceiveReports(int deviceIndex)
{
	long error=0;
	pRecDevice device;
	IOHIDDeviceInterface122** interface=NULL;
	int reason; // kCFRunLoopRunFinished, kCFRunLoopRunStopped, kCFRunLoopRunTimedOut, kCFRunLoopRunHandledSource

	PsychHIDVerifyInit();

	if(deviceIndex < 0 || deviceIndex >= MAXDEVICEINDEXS) PrintfExit("Sorry. Can't cope with deviceNumber %d (more than %d). Please tell denis.pelli@nyu.edu",deviceIndex,(int)MAXDEVICEINDEXS-1);

    // Allocate report buffers if needed:
    PsychHIDAllocateReports(deviceIndex);

	CountReports("ReceiveReports beginning.");
	if (freeReportsPtr[deviceIndex] == NULL) PrintfExit("No free reports.");

	device=PsychHIDGetDeviceRecordPtrFromIndex(deviceIndex);
	if(!HIDIsValidDevice(device))PrintfExit("PsychHID: Invalid device.\n");

	interface = PsychHIDGetDeviceInterfacePtrFromIndex(deviceIndex);
	if(interface==NULL)PrintfExit("PsychHID: No interface for device.\n");
	CheckRunLoopSource(deviceIndex,"ReceiveReports",__LINE__);
	if(!ready[deviceIndex]){
		// setInterruptReportHandlerCallback
		static unsigned char buffer[MAXREPORTSIZE];
		psych_uint32 bufferSize=MAXREPORTSIZE;
		psych_bool createSource;

		createSource=(source[deviceIndex]==NULL);
		if(createSource){
			if(optionsPrintCrashers && createSource)printf("%d: createAsyncEventSource\n",deviceIndex);
			error=(*interface)->createAsyncEventSource(interface,&(source[deviceIndex]));
			if(error)PrintfExit("ReceiveReports - createAsyncEventSource error 0x%lx.",error);
			if(0 && optionsPrintCrashers && createSource)
				printf("%d: source %4.4lx validity %d, CFRunLoopContainsSource is %d.\n",deviceIndex,(unsigned long)source[deviceIndex]
					   ,CFRunLoopSourceIsValid(source[deviceIndex])
					   ,CFRunLoopContainsSource(CFRunLoopGetCurrent(),source[deviceIndex],myRunLoopMode));
		}
		if(optionsPrintCrashers && createSource)printf("%d: getAsyncEventSource\n",deviceIndex);
		CheckRunLoopSource(deviceIndex,"ReceiveReports",__LINE__);
		if(optionsPrintCrashers && createSource)printf("%d: CFRunLoopAddSource\n",deviceIndex);
		CFRunLoopAddSource(CFRunLoopGetCurrent(),source[deviceIndex],myRunLoopMode);
		if(0 && optionsPrintCrashers && createSource)printf("%d: source %4.4lx validity %d, CFRunLoopContainsSource is %d.\n",deviceIndex,(unsigned long)source[deviceIndex]
			   ,CFRunLoopSourceIsValid(source[deviceIndex])
			   ,CFRunLoopContainsSource(CFRunLoopGetCurrent(),source[deviceIndex],myRunLoopMode));
		ready[deviceIndex]=1;
		CheckRunLoopSource(deviceIndex,"ReceiveReports",__LINE__);
		if(optionsPrintCrashers && createSource)printf("%d: setInterruptReportHandlerCallback\n",deviceIndex);
		error=(*interface)->setInterruptReportHandlerCallback(interface,buffer,bufferSize,ReportCallback,buffer,(void *)(long int) deviceIndex);
		if(error)PrintfExit("ReceiveReports - setInterruptReportHandlerCallback error 0x%lx.",error);
		if(optionsPrintCrashers && createSource)printf("%d: CFRunLoopRunInMode.\n",deviceIndex);
	}
	//printf("%d: CFRunLoopRunInMode\n",deviceIndex);
	reason=CFRunLoopRunInMode(myRunLoopMode,optionsSecs,false);
	if(reason!=kCFRunLoopRunTimedOut && reason!=kCFRunLoopRunHandledSource){
		char *s;
		switch(reason){
			case kCFRunLoopRunFinished: s="kCFRunLoopRunFinished"; break;
			case kCFRunLoopRunStopped: s="kCFRunLoopRunStopped"; break;
			case kCFRunLoopRunTimedOut: s="kCFRunLoopRunTimedOut"; break;
			case kCFRunLoopRunHandledSource: s="kCFRunLoopRunHandledSource"; break;
			default: s="of unknown reason.";
		}
		printf("RunLoop ended at %.3f s because %s.\n",CFAbsoluteTimeGetCurrent()-AInScanStart,s);			
	}
	CountReports("ReceiveReports end.");
	return error;
}

PsychError ReceiveReportsStop(int deviceIndex)
{
	CheckRunLoopSource(deviceIndex,"ReceiveReportsStop",__LINE__);
	if(ready[deviceIndex]){
		// Rob Yepez, at Apple, suggested that it might be better to call CFRunLoopRemoveSource than CFRunLoopSourceInvalidate.
		// He's right. There's no problem in re-enabling the callback after using CFRunLoopRemoveSource.
		// The source remains valid and can be added again to the run loop. Don't create it again.
		if(myRunLoopMode==NULL)myRunLoopMode=CFSTR("myMode");; // kCFRunLoopDefaultMode
		if(0 && optionsPrintCrashers)printf("%d: CFRunLoopRemoveSource\n",deviceIndex);
		CFRunLoopRemoveSource(CFRunLoopGetCurrent(),source[deviceIndex],myRunLoopMode);// kCFRunLoopDefaultMode
		if(0 && optionsPrintCrashers)printf("%d: source %4.4lx validity %d, CFRunLoopContainsSource is %d.\n",deviceIndex,(unsigned long)source[deviceIndex]
			   ,CFRunLoopSourceIsValid(source[deviceIndex])
			   ,CFRunLoopContainsSource(CFRunLoopGetCurrent(),source[deviceIndex],myRunLoopMode));
		ready[deviceIndex]=0;
	}
	CheckRunLoopSource(deviceIndex,"ReceiveReportsStop",__LINE__);
	return 0;
}

/* PsychHIDReceiveReportsCleanup(void) -- Called at PsychHID shutdown time:

1. one must call CFRunLoopSourceInvalidate to kill the callback established by setInterruptReportHandlerCallback 
 before calling HIDReleaseDeviceList.

2. after calling setInterruptReportHandlerCallback to enable the callback and CFRunLoopSourceInvalidate to disable 
it, it is then impossible to re-enable the callback with that source. To later re-enable, simply remove the source, instead of
invalidating it. Once i've called CFRunLoopSourceInvalidate it appears that my only option is to release the interface by calling 
 HIDReleaseDeviceList and start over.

*/
PsychError PsychHIDReceiveReportsCleanup(void) 
{
	int deviceIndex;
	
	//printf("Clean up before PsychHID is flushed.\n");
	for(deviceIndex=0;deviceIndex<MAXDEVICEINDEXS;deviceIndex++) if(source[deviceIndex]!=NULL) {
        CheckRunLoopSource(deviceIndex,"PsychHIDReceiveReportsCleanup",__LINE__);
		CFRunLoopRemoveSource(CFRunLoopGetCurrent(),source[deviceIndex],myRunLoopMode);// kCFRunLoopDefaultMode
        if (optionsPrintCrashers) printf("%d: CFRunLoopSourceInvalidate\n",deviceIndex);
        CFRunLoopSourceInvalidate(source[deviceIndex]);
        if(optionsPrintCrashers) printf("%d: source %4.4lx validity %d, CFRunLoopContainsSource is %d.\n",deviceIndex,(unsigned long)source[deviceIndex]
                                           ,CFRunLoopSourceIsValid(source[deviceIndex])
                                           ,CFRunLoopContainsSource(CFRunLoopGetCurrent(),source[deviceIndex],myRunLoopMode));
        ready[deviceIndex]=0;
        CheckRunLoopSource(deviceIndex,"PsychHIDReceiveReportsCleanup",__LINE__);
        source[deviceIndex]=NULL;
	}
    
    // Release all report linked lists, memory buffers etc.:
    PsychHIDReleaseAllReportMemory();

	return 0;
}

void CheckRunLoopSource(int deviceIndex,char *caller,int line){
	CFRunLoopSourceRef currentSource;
	pRecDevice device;
	IOHIDDeviceInterface122** interface = NULL;
	
    // Skip this function if hardcore debugging is not enabled:
    if (!optionsPrintCrashers) return;

	device=PsychHIDGetDeviceRecordPtrFromIndex(deviceIndex);
	if(!HIDIsValidDevice(device))PrintfExit("PsychHID: Invalid device.\n");

	interface = PsychHIDGetDeviceInterfacePtrFromIndex(deviceIndex);
	if(interface==NULL)PrintfExit("PsychHID: No interface for device.\n");
	currentSource=(*interface)->getAsyncEventSource(interface);
	if(source[deviceIndex] != currentSource)printf("%s (%d): source[%d] %4.4lx != current source %4.4lx.\n"
			,caller,line,(int)deviceIndex,(unsigned long)source[deviceIndex],(unsigned long)currentSource);

	if(ready[deviceIndex] && (source[deviceIndex]!=NULL) && !CFRunLoopContainsSource(CFRunLoopGetCurrent(),source[deviceIndex],myRunLoopMode))
		printf("%d: %s(%d): \"ready\" but source not in CFRunLoop.\n",(int)deviceIndex,caller,line);
	if(!ready[deviceIndex] && (source[deviceIndex]!=NULL) && CFRunLoopContainsSource(CFRunLoopGetCurrent(),source[deviceIndex],myRunLoopMode))
		printf("%d: %s(%d): \"!ready\" yet source is in CFRunLoop.\n",(int)deviceIndex,caller,line);
}

#else
// NON OSX CODE (Linux, Windows):
// ==============================

extern hid_device* source[MAXDEVICEINDEXS]; 
extern hid_device* last_hid_device;

/* Do all the report processing for all devices: Iterates in a fetch loop
 * until error condition, or a maximum allowable processing time of
 * optionSecs seconds has been exceeded.
 *
 * Calls hidlib function hid_read() to get reports, one at a time, enqueues
 * it in our own reports lists for later retrieval by 'GiveMeReports' or
 * 'GiveMeReport'.
 *
 */
PsychError ReceiveReports(int deviceIndex)
{
    int rateLimit[MAXDEVICEINDEXS] = { 0 };
    double deadline, now;
    pRecDevice device;
    int n, m;
    unsigned int i;
    ReportStruct *r;
    long error = 0;

    PsychHIDVerifyInit();

    if(deviceIndex < 0 || deviceIndex >= MAXDEVICEINDEXS) PrintfExit("Sorry. Can't cope with deviceNumber %d (more than %d). Please tell denis.pelli@nyu.edu",deviceIndex, (int) MAXDEVICEINDEXS-1);

    // Allocate report buffers if needed:
    PsychHIDAllocateReports(deviceIndex);

    CountReports("ReceiveReports beginning.");
    if (freeReportsPtr[deviceIndex] == NULL) PrintfExit("No free reports.");

    // Enable this device for hid report reception:
    ready[deviceIndex] = TRUE;

    PsychGetAdjustedPrecisionTimerSeconds(&now);
    deadline = now + optionsSecs;
    
    // Iterate until deadline reached or no more pending reports to process:
    while ((error == 0) && (now <= deadline)) {
        // Iterate over all active devices:
        for (deviceIndex = 0; deviceIndex < MAXDEVICEINDEXS; deviceIndex++) {            
            // Test for timeout:
            PsychGetAdjustedPrecisionTimerSeconds(&now);
            if (now > deadline) break;

            // Skip this device if it isn't enabled to receive hid reports:
            if (!ready[deviceIndex]) continue;
            
            // Free target report buffers?
            if(freeReportsPtr[deviceIndex] == NULL) {
                // Darn. We're full. It might be elegant to discard oldest report, but for now, we'll just ignore the new one.
                if (!rateLimit[deviceIndex]) printf("PsychHID: WARNING! ReportCallback warning. No more free reports for deviceIndex %i. Discarding new report.\n", deviceIndex);
                rateLimit[deviceIndex] = 1;
                continue;
            }
            
            // Handle one report for this device:
            CountReports("ReportCallback beginning.");

            device = PsychHIDGetDeviceRecordPtrFromIndex(deviceIndex);
            last_hid_device = (hid_device*) device->interface;

            // Get a report struct to fill in:
            r = freeReportsPtr[deviceIndex];

            // Fetch the actual data: Bytes fetched, or zero for no reports available, or
            // -1 for error condition.
            r->error = hid_read((hid_device*) device->interface, &(r->report[0]), MaxDeviceReportSize[deviceIndex]);

            // Skip remainder if no data received:
            if (r->error == 0) continue;

            // Ok, we got something, even if it is only an error code. Need
            // to move the (r)eport from the free list to the received list:
            freeReportsPtr[deviceIndex] = r->next;
            r->next=NULL;
            
            // install report into the device's list.
            r->next = deviceReportsPtr[deviceIndex];
            deviceReportsPtr[deviceIndex] = r;
            
            // fill in the rest of the report struct
            r->deviceIndex = deviceIndex;
            
            // Timestamp processing:
            PsychGetPrecisionTimerSeconds(&r->time);

            // Success or error?
            if (r->error > 0) {
                // Success: Reset error, assign size of retrieved report:
                r->bytes = r->error;
                r->error = 0;
            }
            else {
                // Error: No data assigned.
                r->bytes = 0;
                
                // Signal error return code -1:
                error = -1;
                
                // Abort fetch loop:
                break;
            }

            if (optionsPrintReportSummary) {
                // print diagnostic summary of the report
                int serial;
                
                serial = r->report[62] + 256 * r->report[63]; // 32-bit serial number at end of AInScan report from PMD-1208FS
                printf("Got input report %4d: %2ld bytes, dev. %d, %4.0f ms. ", serial, (long) r->bytes, deviceIndex, 1000 * (r->time - AInScanStart));
                if(r->bytes>0) {
                    printf(" report ");
                    n = r->bytes;
                    if (n > 6) n=6;
                    for(i=0; i < (unsigned int) n; i++) printf("%3d ", (int) r->report[i]);
                    m = r->bytes - 2;
                    if (m > (int) i) {
                        printf("... ");
                        i = m;
                    }
                    for(; i < r->bytes; i++) printf("%3d ", (int) r->report[i]);
                }
                printf("\n");
            }
            CountReports("ReportCallback end.");
        }
    }
    
    CountReports("ReceiveReports end.");
    return error;
}

PsychError ReceiveReportsStop(int deviceIndex)
{
	pRecDevice device;

	PsychHIDVerifyInit();
    
    // Disable HID report reception:
    ready[deviceIndex] = FALSE;

	device = PsychHIDGetDeviceRecordPtrFromIndex(deviceIndex);
	last_hid_device = (hid_device*) device->interface;

	if (device->interface) hid_close((hid_device*) device->interface);
	device->interface = NULL;
    
	return 0;
}

PsychError PsychHIDReceiveReportsCleanup(void) 
{	
    // Release all report linked lists, memory buffers etc.:
    PsychHIDReleaseAllReportMemory();

	return 0;
}

#endif

// OS INDEPENDENT CODE:
// ====================

void PsychHIDReleaseAllReportMemory(void)
{
    int deviceIndex;    
    for(deviceIndex = 0; deviceIndex < MAXDEVICEINDEXS; deviceIndex++) {        
        if (!firstTimeInit && reportsHaveBeenAllocated[deviceIndex]) {
            free(allocatedReports[deviceIndex]);
            free(reportData[deviceIndex]);
        }
        
        // Reset all stuff that needs to be reset at PsychHID init and shutdown:
        freeReportsPtr[deviceIndex] = NULL;
        deviceReportsPtr[deviceIndex] = NULL;
        allocatedReports[deviceIndex] = NULL;
        reportData[deviceIndex] = NULL;        
        MaxDeviceReports[deviceIndex] = 0;
        MaxDeviceReportSize[deviceIndex] = 0;
        reportsHaveBeenAllocated[deviceIndex] = FALSE;
        source[deviceIndex] = NULL;
        ready[deviceIndex] = FALSE;
    }

    // Reset defaults:
    optionsMaxReports = 10000; // options.maxReports
    optionsMaxReportSize = 65; // options.maxReportSize
    optionsSecs = 0.010;       // options.secs    
    oneShotRealloc = FALSE;
    
    // Toggle firstTimeInit - If this was a shutdown, next call will be init,
    // if it was an init, next call will be a shutdown:
    firstTimeInit = !firstTimeInit;
}

void PsychHIDAllocateReports(int deviceIndex)
{
    int i;
	ReportStruct *r;

    // Reallocation of report buffers requested by caller?
    if (oneShotRealloc) {
        // Yes. Reset flag, as this is a one-shot query:
        oneShotRealloc = FALSE;
        
        // Anything allocated that needs to be reallocated?
        if(reportsHaveBeenAllocated[deviceIndex]) {
            // Yes. Device stopped? Otherwise this is a no-go:
            if (ready[deviceIndex]) {
                // No-No:
                printf("PTB-WARNING:PsychHID:ReceiveReports: Tried to set new option.maxReportSize or option.maxReports on deviceIndex %i while report\n", deviceIndex);
                printf("PTB-WARNING:PsychHID:ReceiveReports: processing is active. Call PsychHID('ReceiveReportsStop', %i); first to release old reports!\n", deviceIndex);
            } else {
                // Release all databuffers, so they get reallocated below:
                free(allocatedReports[deviceIndex]);
                free(reportData[deviceIndex]);
                freeReportsPtr[deviceIndex] = NULL;
                deviceReportsPtr[deviceIndex] = NULL;
                allocatedReports[deviceIndex] = NULL;
                reportData[deviceIndex] = NULL;        
                MaxDeviceReports[deviceIndex] = 0;
                MaxDeviceReportSize[deviceIndex] = 0;

                // Done. Code below will realloc with current settings:
                reportsHaveBeenAllocated[deviceIndex] = FALSE;
            }
        }
    }
    
    // Reports for device already allocated? Do so if not:
	if (!reportsHaveBeenAllocated[deviceIndex]) {
		// Initial set up. Allocate free reports.

        // Allocate common buffer to store linked list of all
        // ReportStruct's, tightly packed:
        allocatedReports[deviceIndex] = (ReportStruct*) calloc(optionsMaxReports, sizeof(ReportStruct));
        if (NULL == allocatedReports[deviceIndex]) PsychErrorExitMsg(PsychError_outofMemory, "Out of memory while trying to allocate hid reports!");

        // Allocate common buffer to store actual report data
        // referenced by ReportStruct's, tightly packed:
        reportData[deviceIndex] = (psych_uint8*) calloc(optionsMaxReports, optionsMaxReportSize);
        if (NULL == reportData[deviceIndex]) {
            // Failed. Free previous allocations:
            free(allocatedReports[deviceIndex]);
            allocatedReports[deviceIndex] = NULL;
            PsychErrorExitMsg(PsychError_outofMemory, "Out of memory while trying to allocate hid report data buffers!");
        }

        // Ok, we have memory allocated. Set it up:

        // Store max number of reports and max size of reports
        // for this device, as defined and allocated from global setting:
        MaxDeviceReports[deviceIndex] = optionsMaxReports;
        MaxDeviceReportSize[deviceIndex] = optionsMaxReportSize;
        
        // Setup pointer mappings to create the linked-list in the memory buffers:
		freeReportsPtr[deviceIndex] = allocatedReports[deviceIndex];
		for(i = 0; i < optionsMaxReports; i++) {
            // Setup linked-list pointers for linked list of ReportStruct's
            // insie the allocatedReports[deviceIndex] buffer:
			r=&(allocatedReports[deviceIndex][i]);
			r->next=&(allocatedReports[deviceIndex][i+1]);
            
            // Setup pointer to associated actual HID report data buffer
            // inside the reportData[deviceIndex] buffer:
            r->report = &(reportData[deviceIndex][i * optionsMaxReportSize]);
		}
		r->next=NULL;        

		reportsHaveBeenAllocated[deviceIndex] = TRUE;
	}
}

void CountReports(char *string)
{
	int i, l1, l2;
	ReportStruct *r;

	// First time init at first invocation after PsycHID load time:
    #if PSYCH_SYSTEM == PSYCH_OSX
	if (myRunLoopMode==NULL) myRunLoopMode=CFSTR("myMode"); // kCFRunLoopDefaultMode
    #endif

	// Optional consistency check, disabled by default. Do the numbers of
	// reports enqueued in the different device lists and the free list
	// sum up to the total number of allocated reports? Print warning and
	// current numbers if this is not the case:
	if (optionsConsistencyChecks > 0) {
		for(i = 0; i < MAXDEVICEINDEXS; i++) {
			r = deviceReportsPtr[i];
            l1 = 0;
			while(r != NULL) {
				r = r->next;
				l1++;
			}
			
			r = freeReportsPtr[i];
            l2 = 0;
			while(r != NULL) {
				r = r->next;
				l2++;
			}
		
            if((l1 + l2) != MaxDeviceReports[i]) {
                printf("%s", string);
                printf(" device:reports. free:%3d, %2d:%3d",l2, i,l1);
                printf("\n");
            }
        }
	}
}

// GiveMeReports is called solely by PsychHIDGiveMeReports, but the code resides here
// in PsychHIDReceiveReports because it uses the typedefs and static variables that
// are defined solely in this file. The linked lists of reports are unknown outside of this file.
PsychError GiveMeReports(int deviceIndex,int reportBytes)
{
	mwSize dims[]={1,1};
	mxArray **outReports;
	ReportStruct *r,*rTail;
	const char *fieldNames[]={"report", "device", "time"};
	mxArray *fieldValue;
	unsigned char *reportBuffer;
	int i,n;
    unsigned int j;
	long error=0;
	double now;
	
	CountReports("GiveMeReports beginning.");

	outReports=PsychGetOutArgMxPtr(1); 
	r=deviceReportsPtr[deviceIndex];
	n=0;
	while(r!=NULL){
		n++;
		rTail=r;
		r=r->next;
	}
	*outReports=mxCreateStructMatrix(1,n,3,fieldNames);
	r=deviceReportsPtr[deviceIndex];
	PsychGetPrecisionTimerSeconds(&now);
	for(i=n-1;i>=0;i--){
		// assert(r!=NULL);
		if(r->error)error=r->error;
		dims[0]=1;
		//printf("%2d: r->bytes %2d, reportBytes %4d, -%4.1f s\n",i,(int)r->bytes,(int)reportBytes, now-r->time);
		if(r->bytes> (unsigned int) reportBytes)r->bytes=reportBytes;
		dims[1]=r->bytes;
		fieldValue=mxCreateNumericArray(2,(void *)dims,mxUINT8_CLASS,mxREAL);
		reportBuffer=(void *)mxGetData(fieldValue);
		for(j=0;j<r->bytes;j++)reportBuffer[j]=r->report[j];
		if(fieldValue==NULL)PrintfExit("Couldn't allocate report array.");
		mxSetField(*outReports,i,"report",fieldValue);
		fieldValue=mxCreateDoubleMatrix(1,1,mxREAL);
		*mxGetPr(fieldValue)=(double)r->deviceIndex;
		mxSetField(*outReports,i,"device",fieldValue);
		fieldValue=mxCreateDoubleMatrix(1,1,mxREAL);
		*mxGetPr(fieldValue)=r->time;
		mxSetField(*outReports,i,"time",fieldValue);
		r=r->next;
	}
	if(deviceReportsPtr[deviceIndex]!=NULL){
		// transfer all these now-obsolete reports to the free list
		rTail->next=freeReportsPtr[deviceIndex];
		freeReportsPtr[deviceIndex]=deviceReportsPtr[deviceIndex];
		deviceReportsPtr[deviceIndex]=NULL;
	}
	CountReports("GiveMeReports end.");
	return error;
}

// Called solely by PsychHIDGetReport, but resides here in order to access the linked list of reports.
PsychError GiveMeReport(int deviceIndex,psych_bool *reportAvailablePtr,unsigned char *reportBuffer,psych_uint32 *reportBytesPtr,double *reportTimePtr)
{
	ReportStruct *r,*rOld;
	long error;
	unsigned int i;
	
	CountReports("GiveMeReport beginning.");

	r=deviceReportsPtr[deviceIndex];
	if(r!=NULL){ // report available?
				 // grab the oldest report for this device
		*reportAvailablePtr=1;
		if(r->next==NULL){
			deviceReportsPtr[deviceIndex]=NULL;
		}else{
			while(r->next->next!=NULL)r=r->next;
			rOld=r;
			r=r->next;
			rOld->next=NULL;
		}
		if(*reportBytesPtr > r->bytes)*reportBytesPtr=r->bytes;
		for(i=0;i<*reportBytesPtr;i++)reportBuffer[i]=r->report[i];
		*reportTimePtr=r->time;
		error=r->error;
		
		// add it to the free list
		r->next=freeReportsPtr[deviceIndex];
		freeReportsPtr[deviceIndex]=r;
	}else{
		*reportAvailablePtr=0;
		*reportBytesPtr=0;
		*reportTimePtr=0.0;
		error=0;
	}
	CountReports("GiveMeReport end.");
	return error;
}				  

/* OS independent entry point from PsychHID(): */

static char useString[]= "err=PsychHID('ReceiveReports',deviceNumber[,options])";
static char synopsisString[]= 
"Receive and save, internally, all reports from the specified USB HID device, now and forever more (unless stopped).\n"
"Some parameters are persistent. When you don't explicitly supply a value they retain whatever value they had last time. "
"They do have an initial default (built-in) that is present when PsychHID is first called after being CLEARed. "
"Non-persistent parameters have a fixed default that applies every time you call PsychHID.\n"
"The returned value \"err.n\" is zero upon success and a nonzero error code upon failure, "
"as spelled out by \"err.name\" and \"err.description\". "
"\"deviceNumber\" specifies which device.\n"
"\"options.print\" =1 (initial default 0) enables diagnostic printing of a summary of each report when our callback routine receives it.\n"	
"\"options.printCrashers\" =1 (initial default 0) enables diagnostic printing of the creation of the callback source and its addition to the CFRunLoop.\n"
"\"options.consistencyChecks\" =1 (initial default 0) enables diagnostic printing of the consistency of all report structs. Very time consuming!\n"
"\"options.maxReports\" (initial default 10000) allocate space for at least this many reports for the given device.\n"
"\"options.maxReportSize\" (initial default 65) allocate this many bytes per report. Most HID devices don't send HID reports of more than "
"64 Bytes, so allowing for one extra byte for the reportID, a default of 65 Bytes is usually sufficient. If you need more, you can increase "
"this value up to 8192 Bytes. If you need even more, contact us, because likely you are doing something wrong. Smaller values than 65 may "
"make sense if you are very tight on memory.\n"
"\"options.secs\" (initial default 0.010 s) is how long to allow the function to process reports received from all active HID devices. "
"The operating system receives reports all the time after the first call to 'ReceiveReports' or 'GetReport'. "
"It has a small buffer capacity, discarding the oldest received reports if its small buffer is full. When requested by PsychHID, the OS "
"tranfers reports to PsychHID (for all devices for which ReceiveReports is still active). "
"Thus reports are received from the OS only during your call to ReceiveReports or GetReport (GetReport implies an automatic call to ReceiveReports). "
"You should call ReceiveReports frequently to avoid losing reports. "
"Reports can be received from multiple devices during a single call to ReceiveReports. "
"Calling ReceiveReports enables callbacks (forever) for the incoming reports from that device; "
"call ReceiveReportsStop to halt acquisition of further reports for a device; "
" you can resume acquisition for a device by calling ReceiveReports again. "
"Call GiveMeReports to get all the received reports and empty PsychHID's internal store for a device. "
"PsychHID can hold up to options.maxReports reports, and discards new incoming reports when it has no room to hold them. "
"For prolonged data acquisition you may need to call GiveMeReports periodically, emptying PsychHID's store before it becomes full.\n"
"PsychHID was enhanced by adding HID commands to send and receive HID reports to support the PMD-1208FS. "
"PsychHID is likely to work with other HID-compliant USB devices as well.\n"
"The device-specific programming for the PMD-1208FS is entirely in the MATLAB M files of the Daq Toolbox; "
"what is specific to the different operating systems is (hopefully) entirely in PsychHID. "
"PsychHID is entirely generic, following the HID standard; none of the internal code is specific to any device, "
"but we tested it primarily with the PMD-1208FS, "
"as well as keyboards, mice, and gamepads, so working with new devices may be an adventure. ";

static char seeAlsoString[]="SetReport, ReceiveReportsStop, GiveMeReports";

PsychError PSYCHHIDReceiveReports(void)
{
	long error=0;
	int deviceIndex;
	mxArray **mxErrPtr;
	const mxArray *mxOptions,*mx;

    PsychPushHelp(useString,synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(2));
	PsychCopyInIntegerArg(1,TRUE,&deviceIndex);

    if(deviceIndex < 0 || deviceIndex >= MAXDEVICEINDEXS) PrintfExit("Sorry. Can't cope with deviceNumber %d (more than %d). Please tell denis.pelli@nyu.edu",deviceIndex, (int) MAXDEVICEINDEXS-1);
    
	/*
	 "\"options.print\" =1 (default 0) enables diagnostic printing of a summary of each report when our callback routine receives it. "	
	 "\"options.printCrashers\" =1 (default 0) enables diagnostic printing of the creation of the callback source and its addition to the CFRunLoop. "
	 "\"options.maxReports\" (default 10000) allocate space for at least this many reports, shared among all devices. "
	 "\"options.maxReportSize\" (default 65) allocate this many bytes per report. "
	 */
	//optionsPrintReportSummary=0;	// options.print
	//optionsPrintCrashers=0;		// options.printCrashers
	//optionsMaxReports=10000;		// options.maxReports
	//optionsMaxReportSize=65;		// options.maxReportSize
	//optionsSecs=0.010;			// options.secs
    
	mxOptions=PsychGetInArgMxPtr(2);
	if(mxOptions!=NULL){
		mx=mxGetField(mxOptions,0,"print");
		if(mx!=NULL)optionsPrintReportSummary=(psych_bool)mxGetScalar(mx);
		mx=mxGetField(mxOptions,0,"printCrashers");
		if(mx!=NULL)optionsPrintCrashers=(psych_bool)mxGetScalar(mx);
		mx=mxGetField(mxOptions,0,"secs");
		if(mx!=NULL)optionsSecs=mxGetScalar(mx);
		mx=mxGetField(mxOptions,0,"consistencyChecks");
		if(mx!=NULL)optionsConsistencyChecks=(psych_bool)mxGetScalar(mx);

        // Changing maxReports or maxReportSize triggers a reallocation of
        // buffer memory:
		mx=mxGetField(mxOptions,0,"maxReports");
		if(mx!=NULL) {
            oneShotRealloc = TRUE;
            optionsMaxReports = (int) mxGetScalar(mx);
        }

		mx=mxGetField(mxOptions,0,"maxReportSize");
		if(mx!=NULL) {
            oneShotRealloc = TRUE;
            optionsMaxReportSize = (int) mxGetScalar(mx);
        }
	}

    // Sanity check:
	if(optionsMaxReports < 1) PsychErrorExitMsg(PsychError_user, "PsychHID ReceiveReports: Sorry, requested maxReports count must be at least 1!");
	if(optionsMaxReportSize < 1) PsychErrorExitMsg(PsychError_user, "PsychHID ReceiveReports: Sorry, requested maxReportSize must be at least 1 byte!");
	if(optionsMaxReportSize > MAXREPORTSIZE) {
        printf("PsychHID ReceiveReports: Sorry, requested maximum report size %d bytes exceeds built-in maximum of %d bytes.\n", optionsMaxReportSize, (int) MAXREPORTSIZE);
        PsychErrorExitMsg(PsychError_user, "Invalid option.maxReportSize provided!");
    }
    
    // Start reception of reports: This will also allocate memory for the reports
    // on first invocation for this deviceIndex:
	error = ReceiveReports(deviceIndex);

	mxErrPtr=PsychGetOutArgMxPtr(1);
	if(mxErrPtr!=NULL){
		const char *fieldNames[]={"n", "name", "description"};
		char *name="",*description="";
		mxArray *fieldValue;

		PsychHIDErrors(NULL, error,&name,&description); // Get error name and description, if available.
		*mxErrPtr=mxCreateStructMatrix(1,1,3,fieldNames);
		fieldValue=mxCreateString(name);
		if(fieldValue==NULL)PrintfExit("PSYCHHIDReceiveReports: Couldn't allocate \"err\".");
		mxSetField(*mxErrPtr,0,"name",fieldValue);
		fieldValue=mxCreateString(description);
		if(fieldValue==NULL)PrintfExit("PSYCHHIDReceiveReports: Couldn't allocate \"err\".");
		mxSetField(*mxErrPtr,0,"description",fieldValue);
		fieldValue=mxCreateDoubleMatrix(1,1,mxREAL);
		if(fieldValue==NULL)PrintfExit("PSYCHHIDReceiveReports: Couldn't allocate \"err\".");
		*mxGetPr(fieldValue)=(double)error;
		mxSetField(*mxErrPtr,0,"n",fieldValue);
	}
    return(PsychError_none);	
}
