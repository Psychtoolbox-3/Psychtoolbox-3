/*
 *    PsychToolbox/Source/Common/PsychHID/PsychHIDGiveMeReports.c
 *
 *    PROJECTS: PsychHID
 *
 *    PLATFORMS:  All
 *
 *    AUTHORS:
 *
 *    denis.pelli@nyu.edu           dgp
 *    mario.kleiner.de@gmail.com    mk
 *
 *    HISTORY:
 *
 *    4/7/05  dgp    Wrote it, based on PsychHIDGetReport.c
 *
 */

#include "PsychHID.h"

static char useString[] = "[reports,err]=PsychHID('GiveMeReports',deviceNumber,[reportBytes])";
static char synopsisString[] =
    "Return, as an output argument, all the saved reports from the connected USB HID device.\n"
    "\"deviceNumber\" specifies which device.\n"
    "If supplied, the optional \"reportBytes\" argument imposes a maximum length on each report; "
    "if necessary, reports will be shortened, but not lengthened. "
    "(This feature allows you to receive a report containing just the data you requested when "
    "the firmware insists on providing a fixed length report that may be longer than the valid data.\n"
    "\"reports\" is a struct array, with a struct for each report. "
    "\"reports(i).report\" is a uint8 vector with received report data. "
    "If your device uses reportID then the first byte of the report is the reportID, the following "
    "bytes contain the actual received report data. Otherwise (reportID==0), the received data "
    "starts already in the first byte of reports(i).report.\n"
    "\"reports(i).device\" is the device number of the device.\n"
    "\"reports(i).time\" is the GetSecs time at which it was received from the system. This is *not* the "
    "time when the hardware itself received the report, therefore this value is of limited use and should "
    "be considered unreliable.\n"
    "The returned value \"err.n\" is zero upon success and a nonzero error code upon failure, "
    "as spelled out by \"err.name\" and \"err.description\". ";

static char seeAlsoString[] = "SetReport, GetReport, ReceiveReports, ReceiveReportsStop, GiveMeReports.";

PsychError GiveMeReports(int deviceIndex,int reportBytes); // PsychHIDReceiveReports.c

PsychError PSYCHHIDGiveMeReports(void)
{
    PsychGenericScriptType *outErr;
    const char *fieldNames[] = {"n", "name", "description"};

    char *name = "", *description = "";
    long error = 0;
    int deviceIndex;
    int reportBytes = 1024;

    PsychPushHelp(useString,synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) { PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumOutputArgs(2));
    PsychErrorExit(PsychCapNumInputArgs(2));

    PsychCopyInIntegerArg(1,TRUE,&deviceIndex);
    PsychCopyInIntegerArg(2,false,&reportBytes);

    PsychHIDVerifyInit();

    // Returns 1st return argument 'reports':
    error = GiveMeReports(deviceIndex,reportBytes); // PsychHIDReceiveReports.c

    // Return 2nd return argument 'err' struct:
    PsychHIDErrors(NULL, error, &name, &description);
    PsychAllocOutStructArray(2, kPsychArgOptional, -1, 3, fieldNames, &outErr);
    PsychSetStructArrayStringElement("name", 0, name, outErr);
    PsychSetStructArrayStringElement("description", 0, description, outErr);
    PsychSetStructArrayDoubleElement("n", 0, (double) error, outErr);

    return(PsychError_none);
}
