/*
	/osxptb/trunk/PsychSourceGL/Source/Common/Eyelink/EyelinkGetFloatData.c
 
	PROJECTS: Eyelink
 
	AUTHORS:
		f.w.cornelissen@rug.nl			fwc
		E.Peters@ai.rug.nl				emp
		e_flister@yahoo.com				edf
 
	PLATFORMS:	All
 
	HISTORY:
 
		29/05/2001	emp 		created it (OS9 version)
		30/10/06	fwc			Adapted from early alpha version.
		19/02/09    edf         added GetFloatDataRaw
		23/03/09    edf         adapted to handle LOST_DATA_EVENT and added eye argument to GetFloatDataRaw
 
	TARGET LOCATION:
 
		Eyelink.mexmac resides in:
			PsychHardware/EyelinkToolbox
 
    WARNING returns only samples at the moment ??? doesn't look like that 
*/


#include "PsychEyelink.h"

#define ERR_BUFF_LEN 1000

static char useString[] = "item = Eyelink('GetFloatData', type)";

static char synopsisString[] =
" makes copy of last item."
" You have to supply the type (returned from Eyelink('getnextdatatype'))."
" Returns the item.";

static char seeAlsoString[] = "";

/*
ROUTINE: EyelinkGetFloatData.c
PURPOSE:
   uses INT16 CALLTYPE eyelink_get_float_data(void FARTYPE *buf);
   makes copy of last item from getnextdatatype
   You have to supply the type (returned from eyelink('getnextdatatype')).
   Returns the item.*/

PsychError EyelinkGetFloatData(void)
{
	FSAMPLE   fs;
	FEVENT    fe;
	int type = 0;
	mxArray **mxpp;
	
	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//check to see if the user supplied superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(1));
	PsychErrorExit(PsychRequireNumInputArgs(1));
	PsychErrorExit(PsychCapNumOutputArgs(1));
	
	// Verify eyelink is up and running
	EyelinkSystemIsConnected();
	EyelinkSystemIsInitialized();
	
	PsychCopyInIntegerArg(1, TRUE, &type);
	
	mxpp = PsychGetOutArgMxPtr(1);
	
	switch(type)
{
		case SAMPLE_TYPE:
			if (eyelink_get_float_data(&fs) != type) {
				PsychErrorExitMsg(PsychError_user, "Eyelink: GetFloatData: eyelink_get_float_data did not return sample type as user said it would.");
			}
			(*mxpp) = CreateMXFSample(&fs);
			break;
		case LOST_DATA_EVENT: // queue overflowed, we are not supposed to call eyelink_get_float_data on this
			(*mxpp) = (mxArray*) mxCreateDoubleMatrix(1,1,mxREAL); //no data
			mxGetPr((*mxpp))[0] = 0;
			break;
		case 0:
			(*mxpp) = (mxArray*) mxCreateDoubleMatrix(1,1,mxREAL); //no data
			mxGetPr((*mxpp))[0] = 0;
			break;
		default: // else it is an event
			if (eyelink_get_float_data(&fe) != type) {
				PsychErrorExitMsg(PsychError_user, "Eyelink: GetFloatDataRaw: eyelink_get_float_data did not return event type as user said it would.");
			}
			(*mxpp) = CreateMXFEvent(&fe);
	}
	
	return(PsychError_none);
}


static char useStringRaw[] = "[item, raw] = Eyelink('GetFloatDataRaw', type [, eye])";
static char synopsisStringRaw[] =
"Same as Eyelink('GetFloatData') but additionally returns raw sample data.\n"
"If type is not SAMPLE_TYPE, identical to GetFloatData (raw will be returned as 0).\n"
"Specify the eye for which raw values are desired as LEFT_EYE or RIGHT_EYE as returned by EyelinkInitDefaults.  Omission of this argument is deprecated.\n"
"Normal samples do not contain the corneal reflection data, but some (non-saccade-based) calibration methods require this information.\n"
"The Eyelink 1000 can be configured to send this information as part of 'raw' samples.\n"
"Sol Simpson at SR-Research emphasizes that this is not officially supported or guaranteed.\n"
"1. You may need to install Eyelink.dll from the latest software developer kit at https://www.sr-support.com/forums/showthread.php?t=6 (windows) or https://www.sr-support.com/forums/showthread.php?t=15 (osx)\n"
"2. Issue the commands:\n"
"\t   Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,AREA,GAZERES,HREF,PUPIL,STATUS,INPUT,HMARKER');\n"
"\t   Eyelink('command','inputword_is_window = ON');\n"
"More info:\n"
"\t   HMARKER (originally for Eyelink2's infrared head tracking markers) and INPUT (originally for the TTL lines) are jury-rigged to hold the extra data.\n"
"\t   You can also set file_sample_data to collect raw samples in the .edf file.\n"
"Raw structure fields:\n"
"\t   raw_pupil           raw x, y sensor position of the pupil\n"
"\t   raw_cr              raw x, y sensor position of the corneal reflection\n"
"\t   pupil_area          raw pupil area\n"
"\t   cr_area             raw cornela reflection area\n"
"\t   pupil_dimension     width, height of raw pupil\n"
"\t   cr_dimension        width, height of raw cr\n"
"\t   window_position     x, y position of tracking window on sensor\n"
"\t   pupil_cr            calculated pupil-cr from the raw_pupil and raw_cr fields\n"
"\t   cr_area2            raw area of 2nd corneal reflection candidate\n"
"\t   raw_cr2             raw x, y sensor position of 2nd corneal reflection candidate\n"
"CAUTION: It may or may not work on your setup with your tracker.\n\n";

/*
	ROUTINE: EyelinkGetFloatDataRaw.c
	PURPOSE: As EyelinkGetFloatData, but additionally returns raw sample data.
 */

//forgive the repeated code here, feel free to write a helper
PsychError EyelinkGetFloatDataRaw(void)
{
	FSAMPLE      fs;
	FSAMPLE_RAW  fr;
	FEVENT       fe;
	int type =   0;
	mxArray **mxpp[2];
	int eye;
	int err;
	char errmsg[ERR_BUFF_LEN]="";
	
	//all sub functions should have these two lines
	PsychPushHelp(useStringRaw, synopsisStringRaw, seeAlsoString);
	if (PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };
	
	// Check to see if the user supplied superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(2));
	PsychErrorExit(PsychRequireNumInputArgs(1));
	PsychErrorExit(PsychCapNumOutputArgs(2));
	
	// Verify eyelink is up and running
	EyelinkSystemIsConnected();
	EyelinkSystemIsInitialized();
	
	PsychCopyInIntegerArg(1, TRUE, &type);
	
	mxpp[0] = PsychGetOutArgMxPtr(1);
	mxpp[1] = PsychGetOutArgMxPtr(2);
	
	switch(type){
		case SAMPLE_TYPE:
			if (eyelink_get_float_data(&fs) != type) {
				PsychErrorExitMsg(PsychError_user, "Eyelink: GetFloatDataRaw: eyelink_get_float_data did not return sample type as user said it would.");
			}
			(*mxpp)[0] = CreateMXFSample(&fs);
			
			if (PsychCopyInIntegerArg(2, kPsychArgOptional, &eye)) {
				if (eye!=LEFT_EYE && eye!=RIGHT_EYE) {
					PsychErrorExitMsg(PsychErorr_argumentValueOutOfRange, "EyeLink: GetFloatDataRaw:  eye argument must be LEFT_EYE or RIGHT_EYE as returned by EyelinkInitDefaults\n");
				}
				TrackerOKForRawValues();
				memset(&fr, 0, sizeof(fr));
				if((err = eyelink_get_extra_raw_values_v2(&fs, eye, &fr))){
					sprintf(errmsg, "Eyelink: GetFloatDataRaw: eyelink_get_extra_raw_values_v2 returned error code %d: %s", err, eyelink_get_error(err,"eyelink_get_extra_raw_values_v2")); //no snprintf in msvs?  bug: buff overflow
					PsychErrorExitMsg(PsychError_internal, errmsg);
				}
			} else {
				mexPrintf("EYELINK: WARNING! Omission of the eye argument to GetFloatDataRaw is deprecated.\n");
				if((err = eyelink_get_extra_raw_values(&fs, &fr))){ //deprecated as of Dec 1, 2006 (see eyelink sdk core_expt.h)
					sprintf(errmsg, "Eyelink: GetFloatDataRaw: eyelink_get_extra_raw_values returned error code %d: %s", err, eyelink_get_error(err,"eyelink_get_extra_raw_values")); //no snprintf in msvs?  bug: buff overflow
					PsychErrorExitMsg(PsychError_internal, errmsg);
				}
			}
			
			(*mxpp)[1] = (mxArray *)CreateMXFSampleRaw(&fr); //assumes there is always a raw sample...
			break;
			
		case LOST_DATA_EVENT: // queue overflowed, we are not supposed to call eyelink_get_float_data on this
			(*mxpp)[0] = (mxArray*) mxCreateDoubleMatrix(1,1,mxREAL); //no data
			mxGetPr((*mxpp)[0])[0] = 0;
			(*mxpp)[1] = (mxArray*) mxCreateDoubleMatrix(1,1,mxREAL); //no raw
			mxGetPr((*mxpp)[0])[1] = 0;
			break;
			
		case 0:
			(*mxpp)[0] = (mxArray*) mxCreateDoubleMatrix(1,1,mxREAL); //no data
			mxGetPr((*mxpp)[0])[0] = 0;
			(*mxpp)[1] = (mxArray*) mxCreateDoubleMatrix(1,1,mxREAL); //no raw
			mxGetPr((*mxpp)[0])[1] = 0;
			break;
			
		default: // else it is an event
			if (eyelink_get_float_data(&fe) != type) {
				PsychErrorExitMsg(PsychError_user, "Eyelink: GetFloatDataRaw: eyelink_get_float_data did not return event type as user said it would.");
			}

			(*mxpp)[0] = CreateMXFEvent(&fe);
			(*mxpp)[1] = (mxArray*) mxCreateDoubleMatrix(1,1,mxREAL); //no raw
			mxGetPr((*mxpp)[0])[1] = 0;
	}
	
	return(PsychError_none);
}
