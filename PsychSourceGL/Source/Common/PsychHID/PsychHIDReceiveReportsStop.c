/*
 *    PsychToolbox/Source/Common/PsychHID/PsychHIDReceiveReportsStop.c
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
 *    4/7/05  dgp    Wrote it, based on PsychHIDGetReport.c
 *
 */

#include "PsychHID.h"

static char useString[] = "err = PsychHID('ReceiveReportsStop', deviceNumber)";
static char synopsisString[] =
    "Stop receiving and saving reports from the specified USB HID device.\n"
    "Calling ReceiveReports enables callbacks (forever) for the incoming reports from that device; "
    "call ReceiveReportsStop to halt acquisition of further reports for this device; "
    "you can resume acquisition by calling ReceiveReports again. "
    "Call GiveMeReports to get all the received reports and empty PsychHID's internal store for that device. "
    "\"deviceNumber\" specifies which device. "
    "The returned value \"err.n\" is zero upon success and a nonzero error code upon failure, "
    "as spelled out by \"err.name\" and \"err.description\". ";

static char seeAlsoString[] = "SetReport, ReceiveReports, GiveMeReports";

PsychError PSYCHHIDReceiveReportsStop(void)
{
    PsychGenericScriptType *outErr;
    const char *fieldNames[] = {"n", "name", "description"};
    char *name = "", *description = "";
    long error = 0;
    int deviceIndex;

    PsychPushHelp(useString,synopsisString,seeAlsoString);
    if (PsychIsGiveHelp()) { PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(1));

    PsychCopyInIntegerArg(1, TRUE, &deviceIndex);

    error = ReceiveReportsStop(deviceIndex);
    PsychHIDErrors(NULL, error, &name, &description);

    PsychAllocOutStructArray(1, kPsychArgOptional, -1, 3, fieldNames, &outErr);
    PsychSetStructArrayStringElement("name", 0, name, outErr);
    PsychSetStructArrayStringElement("description", 0, description, outErr);
    PsychSetStructArrayDoubleElement("n", 0, (double) error, outErr);

    return(PsychError_none);
}
