/*
 Source/Common/PsychHID/PsychHIDSetReport.c		
 
 PROJECTS: PsychHID
 
 PLATFORMS:  All. 
 
 AUTHORS:
    denis.pelli@nyu.edu                 dgp
    mario.kleiner@tuebingen.mpg.de      mk
 
 HISTORY:
 3/15/05  dgp		Wrote it.
 4/23/05  dgp Now use low-level setReport instead of high-level HIDSetReport so that we can specify a time out (50 ms) instead of hanging forever when there is trouble.
 */

#include "PsychHID.h"

static char useString[]= "err=PsychHID('SetReport',deviceNumber,reportType,reportID,report)";
static char synopsisString[] = "Send a report to the connected USB HID device.\n"
"\"deviceNumber\" specifies which device.\n"
"\"reportType\" is 2=output, 3=feature (0 to just echo arguments).\n"
"\"reportID\" is either zero or an integer (1 to 255) specifying the topic, e.g. read analog, read digital, write analog, etc. "
"If you provide a non-zero reportID, the first byte of your \"report\" will be overwritten with this reportID. You have to "
"take this into account, ie., leave a leading byte of space for the reportID to avoid corrupting your actual report data. "
"If reportID is zero, then your \"report\" will be sent as-is, without any special treatment of the first byte.\n"
"\"report\" must be an array of char or integer (8-, 16-, or 32-bit, signed or unsigned) holding "
"the correct total number of bytes.\n"
"The returned value \"err.n\" is zero upon success and a nonzero error code upon failure, "
"as spelled out by \"err.name\" and \"err.description\".\n";

static char seeAlsoString[] = "GetReport";

/*
 The HID firmware notes I have to program the Measurement Computing PMD1208FS indicate that the reportID is one byte. It specifies
 what the report is for (0x03 read digital, 0x04 write digital, 0x10 read analog, 0x14 write analog, etc.).
 */

extern double AInScanStart;

PsychError PSYCHHIDSetReport(void) 
{
    static unsigned char scratchBuffer[MAXREPORTSIZE+1];
    
    const char *fieldNames[]={"n", "name", "description"};
    char *name="",*description="";
    mxArray *fieldValue;
    
	long error;
	pRecDevice device;
	int deviceIndex;
	int reportType; // kIOHIDReportTypeInput=0, kIOHIDReportTypeOutput=1, or kIOHIDReportTypeFeature=2
	int reportID;
	unsigned char *reportBuffer;
	int reportSize;
	const mxArray *report;
	mxArray **outErr;
	int i;
	
    PsychPushHelp(useString,synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(4));
	outErr=PsychGetOutArgMxPtr(1);
	PsychCopyInIntegerArg(1,TRUE,&deviceIndex);
	PsychCopyInIntegerArg(2,TRUE,&reportType);
	PsychCopyInIntegerArg(3,TRUE,&reportID);
	report=PsychGetInArgMxPtr(4); 
	reportBuffer=(void *)mxGetData(report);

	switch(mxGetClassID(report)){
		case mxCHAR_CLASS:    
		case mxINT8_CLASS: 
		case mxUINT8_CLASS:   
		case mxINT16_CLASS:  
		case mxUINT16_CLASS:  
		case mxINT32_CLASS:  
		case mxUINT32_CLASS: 
			break;
		default:
			PrintfExit("\"report\" array must be char or integer (8-, 16-, or 32-bit).");
			break;
	}

	reportSize=mxGetElementSize(report)*mxGetNumberOfElements(report);
    if (reportSize > MAXREPORTSIZE) PsychErrorExitMsg(PsychError_user, "Tried to send a HID report which exceeds the maximum allowable size! Aborted.");
    if (reportSize < 1) PsychErrorExitMsg(PsychError_user, "Tried to send an empty HID report! Aborted.");

    // No invalid numbers, no input reports for 'SetReport':
    if (reportType < 0 || reportType > 3 || reportType == 1) PsychErrorExitMsg(PsychError_user, "Invalid 'reportType' for this function provided!");

	PsychHIDVerifyInit();
    
    device=PsychHIDGetDeviceRecordPtrFromIndex(deviceIndex);

    // For non-zero reportID, the first byte of buffer is always overwritten with reportID on all os platforms:
	if (reportID!=0) *reportBuffer= 0xff & reportID; // copy low byte of reportID to first byte of report.

	if (reportType==0) {
        // Echo report:
		printf("SetReport(reportType %d, reportID %d, report ",reportType,reportID);
		for(i=0;i<reportSize;i++)printf("%d ",(int)reportBuffer[i]);
		printf(")\n");
		error=0;        
	} else {
#if PSYCH_SYSTEM == PSYCH_OSX
        IOHIDDeviceInterface122 **interface = NULL;
        interface = PsychHIDGetDeviceInterfacePtrFromIndex(deviceIndex);
        if(interface) error=(*interface)->setReport(interface,reportType-1,reportID,reportBuffer,(psych_uint32) reportSize,50,NULL,NULL,NULL);
        else PrintfExit("PsychHID SetReport: Bad interface.\n");

        // Timestamp:
        if(reportID==0x11) {
            PsychGetPrecisionTimerSeconds(&AInScanStart);
        }

        // Error condition?
        PsychHIDErrors(NULL, error, &name, &description);
#else
        // reportID zero is a special case. We must not touch/overwrite the 1st byte
        // of the usercode provided buffer, but still prefix the buffer passed to HIDLIB
        // with a zero byte:
        if (reportID == 0) {
            memcpy(&scratchBuffer[1], reportBuffer, reportSize);
            scratchBuffer[0] = (unsigned char) 0;
            reportBuffer =  &scratchBuffer[0];
            reportSize++;
        }
        
        // Which report type?
        if (reportType == 2) {
            // Output report:
            
            // Write it: error == -1 would mean error, otherwise it is number of bytes written.
            error = hid_write((hid_device*) device->interface, reportBuffer, (size_t) reportSize);
        }
        else {
            // reportType == 3, aka Feature report:

            // Send it: error == -1 would mean error, otherwise it is number of bytes written.
            error = hid_send_feature_report((hid_device*) device->interface, reportBuffer, (size_t) reportSize);
        }
        
        // Timestamp:
        if(reportID==0x11) {
            PsychGetPrecisionTimerSeconds(&AInScanStart);
        }

        // Positive value means no error:
        if (error >= 0) error = 0;

        // Error condition?
        PsychHIDErrors((hid_device*) device->interface, error, &name, &description);
#endif
	}

    *outErr=mxCreateStructMatrix(1,1,3,fieldNames);
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
    
    return(PsychError_none);	
}
