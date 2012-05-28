/*
	PsychToolbox/Source/Common/PsychHID/PsychHIDGetReport.c		

	PROJECTS: PsychHID

	PLATFORMS:  All.

	AUTHORS:
	denis.pelli@nyu.edu                 dgp
    mario.kleiner@tuebingen.mpg.de      mk

	HISTORY:
	3/15/05  dgp	Wrote it.
	4/3/05   dgp	Got it working with the PMD-1208FS for the first time.
	4/7/05   dgp	Removed the guts of it, to create
					PsychHIDReceiveReports.c. We now call a routine there to do the real work.

	READ:
	bugs in mac os x retrieval of reports.
	http://lists.apple.com/archives/usb/2004/Jul/msg00003.html
	$900 USB bus analyzer
	http://lists.apple.com/archives/usb/2003/Dec/msg00011.html
	guido kolber reports problem with getting port 
	http://lists.apple.com/archives/usb/2004/Nov/msg00062.html

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

#include "PsychHID.h"

#if PSYCH_SYSTEM != PSYCH_OSX
extern hid_device* last_hid_device;
#endif

static char useString[]= "[report,err]=PsychHID('GetReport',deviceNumber,reportType,reportID,reportBytes)";
static char synopsisString[]=
	"Get a report from the connected USB HID device.\n"
	"\"deviceNumber\" specifies which device. \"reportType\" is 1=input, 3=feature (0 to just echo arguments).\n"
	"\"reportID\" is either zero or an integer (1 to 255) specifying the topic, e.g. read analog, read digital, write analog, etc.\n"
	"\"reportBytes\" is maximum size of report to retrieve in bytes.\n"
	"\"report\" is a uint8 vector with the data of the retrieved report. "
	"If your device uses reportID for this report then the first byte of the report is the reportID, "
    "otherwise the first byte will be already the first byte of the retrieved report.\n"
	"The returned value \"err.n\" is zero upon success and a nonzero error code upon failure, "
	"as spelled out by \"err.name\" and \"err.description\". "
	"When nothing has been received, \"report\" will be an empty matrix. "
	"However, some devices may send a zero-length report, which will also result in \"report\" being an empty matrix. "
	"Resolving this ambiguity, \"err.reportLength\" is -1 when no report was received, and the report length in bytes when a report was received.\n "
	"Each time that you use a new deviceNumber, PsychHID:GetReport enables callbacks for the incoming reports from that device. "
	"If you use many devices, all their reports will cause callbacks, as documented in the diagnostic printout. "
	"You can disable a callback that you no longer want by calling ReceiveReportsStop. This saves computation time. ";
static char seeAlsoString[]="SetReport, ReceiveReports, ReceiveReportsStop, GiveMeReports.";

PsychError PSYCHHIDGetReport(void) 
{
	long error=0;
	int deviceIndex;
	pRecDevice device;    
	int reportType; // 1=input, 2=output, 3=feature
	unsigned char *reportBuffer;
	psych_uint32 reportBytes=0;
	unsigned int reportBufferSize=0;
	int reportID=0;
	mwSize dims[]= {1,1};
	mxArray **outReport,**outErr;
	char *name="",*description="",string[256];
	psych_bool reportAvailable;
	double reportTime;

    PsychPushHelp(useString,synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
    PsychErrorExit(PsychCapNumOutputArgs(2));
    PsychErrorExit(PsychCapNumInputArgs(4));
	PsychCopyInIntegerArg(1,TRUE,&deviceIndex);
	PsychCopyInIntegerArg(2,TRUE,&reportType);
	PsychCopyInIntegerArg(3,TRUE,&reportID);
	PsychCopyInIntegerArg(4,TRUE,(int*) &reportBufferSize);
    if (reportBufferSize < 1) PsychErrorExitMsg(PsychError_user, "Size of receive buffer 'reportBytes' must be at least 1!");

	outReport=PsychGetOutArgMxPtr(1); 
	outErr=PsychGetOutArgMxPtr(2); // outErr==NULL if optional argument is absent.
	dims[0]=1;
	dims[1]=reportBufferSize;
	*outReport=mxCreateNumericArray(2,(void *)dims,mxUINT8_CLASS,mxREAL);
	if(*outReport==NULL)PrintfExit("Couldn't allocate report array.");
	reportBuffer=(void *)mxGetData(*outReport);
	if(reportBuffer==NULL)PrintfExit("Couldn't allocate report buffer.");
	reportBytes=reportBufferSize;

	PsychHIDVerifyInit();

    // No invalid numbers, no output reports for 'GetReport':
    if (reportType < 0 || reportType > 3 || reportType == 2) PsychErrorExitMsg(PsychError_user, "Invalid 'reportType' for this function provided!");
    
	if (reportType==0) {
        // Echo report:
		printf("GetReport(reportType %d, reportID %d, reportBytes %d)\n",reportType,reportID,(int)reportBytes);
	} else {
        #if PSYCH_SYSTEM == PSYCH_OSX
            if (reportType == 3) {
                // Feature report: Use special method for this:
                device = PsychHIDGetDeviceRecordPtrFromIndex(deviceIndex);
                if (!HIDIsValidDevice(device)) PrintfExit("PsychHID:GetReport: Invalid device.\n");

                IOHIDDeviceInterface122** interface = NULL;
                interface = PsychHIDGetDeviceInterfacePtrFromIndex(deviceIndex);
                if (interface==NULL) PrintfExit("PsychHID:GetReport: No interface for device.\n");                

                // Apple defines constants for the reportType with values (0,1,2) that are one less that those specified by USB (1,2,3).
                error = (*interface)->getReport(interface, reportType-1, reportID, reportBuffer, &reportBytes, -1, nil, nil, nil);
                if (error == 0) {
                    reportAvailable = 1;
                } else {
                    reportAvailable = 0;
                }
                
                PsychGetAdjustedPrecisionTimerSeconds(&reportTime);
            }
            else { 
                // Input report: Using setInterruptReportHandlerCallback and CFRunLoopRunInMode
                error = ReceiveReports(deviceIndex);
                if (!error) error = GiveMeReport(deviceIndex,&reportAvailable,reportBuffer,&reportBytes,&reportTime);
            }
        #else
            if (reportType == 3) {
                // Feature report: Use special HIDLIB method for this:
                device = PsychHIDGetDeviceRecordPtrFromIndex(deviceIndex);

                // 1st byte always contains reportID:
                reportBuffer[0] = (unsigned char) reportID;

                // Get it: error == -1 would mean error, otherwise it is number of bytes retrieved.
                last_hid_device = (hid_device*) device->interface;
                error = hid_get_feature_report((hid_device*) device->interface, reportBuffer, (size_t) reportBytes);
                if (error >= 0) {
                    reportBytes = error;
                    error = 0;
                    reportAvailable = 1;
                } else {
                    reportBytes = 0;
                    reportAvailable = 0;
                }
                PsychGetAdjustedPrecisionTimerSeconds(&reportTime);
            }
            else {
                // Regular Input reports:
                error=ReceiveReports(deviceIndex);
                if (!error) error = GiveMeReport(deviceIndex,&reportAvailable,reportBuffer,&reportBytes,&reportTime);
            }
        #endif
	}
    
	if(outReport==NULL) PrintfExit("Output argument is required."); // I think MATLAB always provides this.

	if(error==0 && reportBytes>reportBufferSize){
		error=2;
		name="Warning";
		description=string;
		sprintf(description,"GetReport overflow. Expected %ld but received %ld bytes.",(long)reportBufferSize,(long)reportBytes); 
	}

	if(error==0 && reportBytes<reportBufferSize){
		// Reduce the declared size of the array to that of the received report.
		if(reportBytes>0)error=3;
		name="Warning";
		description=string;
		sprintf(description,"GetReport expected %ld but received %ld bytes.",(long)reportBufferSize,(long)reportBytes);
		mxSetN(*outReport,reportBytes);
	}

	if(outErr!=NULL){
		const char *fieldNames[]={"n", "name", "description", "reportLength", "reportTime"};
		mxArray *fieldValue;

		if ((error != 2) && (error != 3)) PsychHIDErrors(NULL, error,&name,&description); // Get error name and description, if available.
		*outErr=mxCreateStructMatrix(1,1,5,fieldNames);
		fieldValue=mxCreateString(name);
		if(fieldValue==NULL)PrintfExit("Couldn't allocate \"err\".");
		mxSetField(*outErr,0,"name",fieldValue);
		fieldValue=mxCreateString(description);
		if(fieldValue==NULL)PrintfExit("Couldn't allocate \"err\".");
		mxSetField(*outErr,0,"description",fieldValue);
		fieldValue=mxCreateDoubleMatrix(1,1,mxREAL);
		if(fieldValue==NULL)PrintfExit("Couldn't allocate \"err\".");
		*mxGetPr(fieldValue)=(double)error;
		mxSetField(*outErr,0,"n",fieldValue);
		fieldValue=mxCreateDoubleMatrix(1,1,mxREAL);
		if(fieldValue==NULL)PrintfExit("Couldn't allocate \"err\".");
		if(reportAvailable)*mxGetPr(fieldValue)=(double)reportBytes;
		else *mxGetPr(fieldValue)=-1.0;
		mxSetField(*outErr,0,"reportLength",fieldValue);
		fieldValue=mxCreateDoubleMatrix(1,1,mxREAL);
		if(fieldValue==NULL)PrintfExit("Couldn't allocate \"err\".");
		*mxGetPr(fieldValue)=reportTime;
		mxSetField(*outErr,0,"reportTime",fieldValue);
	}
    return(PsychError_none);	
}
