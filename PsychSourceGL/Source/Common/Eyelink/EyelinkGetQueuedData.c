/*
	/osxptb/trunk/PsychSourceGL/Source/Common/Eyelink/EyelinkGetQueuedData.c
 
	PROJECTS: Eyelink
 
	AUTHORS:
		e_flister@yahoo.com			edf
 
	PLATFORMS:	all
 
	HISTORY:
 
		21/03/2009	edf 		created it
 
	TARGET LOCATION:
 
		Eyelink.mexmac resides in:
			PsychHardware/EyelinkToolbox
 */

#include "PsychEyelink.h"

static char useString[] = "[samples, events, drained] = Eyelink('GetQueuedData' [, eye])";

static char synopsisString[] =
"dequeues all samples and events from the link.\n"
"returns double matrices where columns are items and rows are fields from eyelink sample structs.\n"
"return flag 'drained' indicates whether queue was emptied or if this function needs to be called again.\n"
"if you include the eye argument (as LEFT_EYE or RIGHT_EYE as returned by EyelinkInitDefaults), samples will include raw fields for that eye.\n"
"if you don't remove items from the queue often enough, the oldest items will be replaced by a LOSTDATAEVENT, which will appear in the sample records at the location where items were dropped (fields other than type will be set to MISSING_DATA).\n\n"

"sample rows are as follows: \n"
"\t 1: time of sample (when camera imaged eye, in milliseconds since tracker was activated)\n"
"\t 2: type (always SAMPLE_TYPE or LOSTDATAEVENT as returned by EyelinkInitDefaults)\n"
"\t 3: flags (bits indicating what types of data are present, and for which eye(s) - see eye_data.h)\n"
"\t 4: left pupil center x (camera coordinates)\n"
"\t 5: right pupil center x (camera coordinates)\n"
"\t 6: left pupil center y (camera coordinates)\n"
"\t 7: right pupil center y (camera coordinates)\n"
"\t 8: left HEADREF x (angular gaze coordinates)\n"
"\t 9: right HEADREF x (angular gaze coordinates)\n"
"\t 10: left HEADREF y (angular gaze coordinates)\n"
"\t 11: right HEADREF y (angular gaze coordinates)\n"
"\t 12: left pupil size (arbitrary units, area or diameter as selected by pupil_size_diameter command)\n"
"\t 13: right pupil size (arbitrary units, area or diameter as selected by pupil_size_diameter command)\n"
"\t 14: left gaze position x (in pixel coordinates set by screen_pixel_coords command)\n"
"\t 15: right gaze position x (in pixel coordinates set by screen_pixel_coords command)\n"
"\t 16: left gaze position y (in pixel coordinates set by screen_pixel_coords command)\n"
"\t 17: right gaze position y (in pixel coordinates set by screen_pixel_coords command)\n"
"\t 18: angular resolution x (at gaze position in screen pixels per visual degree)\n"
"\t 19: angular resolution y (at gaze position in screen pixels per visual degree)\n"
"\t 20: status (error and status flags (only useful for EyeLink II and EyeLink1000, report CR status and tracking error). see eye_data.h.)\n"
"\t 21: input (data from input port(s))\n"
"\t 22: button input data: high 8 bits indicate changes from last sample, low 8 bits indicate current state of buttons 8 (MSB) to 1 (LSB)\n"
"\t 23: head-tracker data type (0=none)\n"
"\t 24: head-tracker data (not prescaled) 1\n"
"\t 25: head-tracker data (not prescaled) 2\n"
"\t 26: head-tracker data (not prescaled) 3\n"
"\t 27: head-tracker data (not prescaled) 4\n"
"\t 28: head-tracker data (not prescaled) 5\n"
"\t 29: head-tracker data (not prescaled) 6\n"
"\t 30: head-tracker data (not prescaled) 7\n"
"\t 31: head-tracker data (not prescaled) 8\n\n"

"if you request the raw sample fields, they will appear in the following additional rows:\n"
"\t 32: raw x sensor position of the pupil\n"
"\t 33: raw y sensor position of the pupil\n"
"\t 34: raw x sensor position of the corneal reflection\n"
"\t 35: raw y sensor position of the corneal reflection\n"
"\t 36: raw pupil area\n"
"\t 37: raw corneal reflection area\n"
"\t 38: raw width of pupil\n"
"\t 39: raw height of pupil\n"
"\t 40: raw width of corneal reflection\n"
"\t 41: raw height of corneal reflection\n"
"\t 42: raw x position of tracking window on sensor\n"
"\t 43: raw y position of tracking window on sensor\n"
"\t 44: (raw pupil x) - (raw corneal reflection x)\n"
"\t 45: (raw pupil y) - (raw corneal reflection y)\n"
"\t 46: raw area of 2nd corneal reflection candidate\n"
"\t 47: raw x sensor position of the 2nd corneal reflection candidate\n"
"\t 48: raw y sensor position of the 2nd corneal reflection candidate\n\n"

"About raw fields:\n"
"Normal samples do not contain the corneal reflection data, but some (non-saccade-based) calibration methods require this information.\n"
"The Eyelink 1000 can be configured to send this information as part of 'raw' samples, and this is required before use of this function.\n"
"Sol Simpson at SR-Research emphasizes that this is not officially supported or guaranteed.\n"
"1. You need to install Eyelink.dll from the latest software developer kit at https://www.sr-support.com/forums/showthread.php?t=6 (windows) or https://www.sr-support.com/forums/showthread.php?t=15 (osx)\n"
"\t   and the latest tracker host software at https://www.sr-support.com/forums/showthread.php?t=179\n"
"2. Issue the commands:\n"
"\t   Eyelink('command','link_sample_data = LEFT,RIGHT,GAZE,AREA,GAZERES,HREF,PUPIL,STATUS,INPUT,HMARKER');\n"
"\t   Eyelink('command','inputword_is_window = ON');\n"
"More info:\n"
"\t   HMARKER (originally for Eyelink2's infrared head tracking markers) and INPUT (originally for the TTL lines) are jury-rigged to hold the extra data.\n"
"\t   You can also set file_sample_data to collect raw samples in the .edf file.\n"
"CAUTION: It may or may not work on your setup with your tracker.\n\n"
" event rows are as follows: \n"
"\t 1: effective time of event\n"
"\t 2: event type\n"
"\t 3: read (bits indicating which data fields contain valid data - see eye_data.h.)\n"
"\t 4: eye\n"
"\t 5: start time\n"
"\t 6: end time\n"
"\t 7: HEADREF gaze position starting point x\n"
"\t 8: HEADREF gaze position starting point y\n"
"\t 9: display gaze position starting point x (in pixel coordinates set by screen_pixel_coords command)\n"
"\t 10: display gaze position starting point y (in pixel coordinates set by screen_pixel_coords command)\n"
"\t 11: starting pupil size (arbitrary units, area or diameter as selected by pupil_size_diameter command)\n"
"\t 12: HEADREF gaze position ending point x\n"
"\t 13: HEADREF gaze position ending point y\n"
"\t 14: display gaze position ending point x (in pixel coordinates set by screen_pixel_coords command)\n"
"\t 15: display gaze position ending point y (in pixel coordinates set by screen_pixel_coords command)\n"
"\t 16: ending pupil size (arbitrary units, area or diameter as selected by pupil_size_diameter command)\n"
"\t 17: HEADREF gaze position average x\n"
"\t 18: HEADREF gaze position average y\n"
"\t 19: display gaze position average x (in pixel coordinates set by screen_pixel_coords command)\n"
"\t 20: display gaze position average y (in pixel coordinates set by screen_pixel_coords command)\n"
"\t 21: average pupil size (arbitrary units, area or diameter as selected by pupil_size_diameter command)\n"
"\t 22: average gaze velocity magnitude (absolute value) in visual degrees per second\n"
"\t 23: peak gaze velocity magnitude (absolute value) in visual degrees per second\n"
"\t 24: starting gaze velocity in visual degrees per second\n"
"\t 25: ending gaze velocity in visual degrees per second\n"
"\t 26: starting angular resolution x in screen pixels per visual degree\n"
"\t 27: ending angular resolution x in screen pixels per visual degree\n"
"\t 28: starting angular resolution y in screen pixels per visual degree\n"
"\t 29: ending angular resolution y in screen pixels per visual degree\n"
"\t 30: status (collected error and status flags from all samples in the event (only useful for EyeLink II and EyeLink1000, report CR status and tracking error). see eye_data.h.)\n\n";

static char seeAlsoString[] = "";

#define NUM_SAMPLE_FIELDS 31
#define NUM_RAW_SAMPLE_FIELDS 17
#define NUM_EVENT_FIELDS 30
#define FUDGE_FACTOR 2 // how much more space we allocate beyond the reported queue length to account for additional items arriving as we are dequeueing
#define ERR_BUFF_LEN 1000

/*
ROUTINE: EyelinkGetQueuedData
PURPOSE:
	matlab is slow at dealing with structs and looping over eyelink_get_float_data to drain the queue, so we take care of this for the client.
	also eliminates usage error of supplying incorrect type from Eyelink('GetNextDataType') to Eyelink('GetFloatData').
 
 TODO:
    -enable BINOCULAR raws
    -request sr research to let us read inputword_is_window
    -option for output in struct format?  would allow deprecating EyelinkGetFloatData, EyelinkGetFloatDataRaw, and EyelinkGetNextDataType
    -would NaN be a better choice (PsychGetNanValue) than MISSING_DATA for LOST_DATA_EVENT?  right now, MISSING_DATA would not fit into the uint type of many of the fields, if we were to offer them in correctly-typed fields of a struct
    -(also in GetFloatData and ImageModeDisplay) replace all sprintf's with snprintf and _snprintf ala (note mario is adding a macro for this, probably to be called snprintf):
#if PSYCH_SYSTEM != PSYCH_WINDOWS
	snprintf(blah, sizeof(blah) - 1, "str",...);
#else
	_snprintf(blah, sizeof(blah) - 1, "str",...);
#endif
 */

PsychError EyelinkGetQueuedData(void)
{
	FSAMPLE      fs;
	FSAMPLE_RAW  fr;
	FEVENT       fe;
	int numSamples = 0, numEvents = 0, maxSamples, maxEvents, type, eye, fieldNum, index, numSampleFields, err;
	double *samples, *events;
	psych_bool useEye=FALSE;
	PsychNativeBooleanType drained=(PsychNativeBooleanType)FALSE;
	char errmsg[ERR_BUFF_LEN]="";
	
	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//check to see if the user supplied superfluous arguments
	PsychErrorExit(PsychCapNumInputArgs(1));
	PsychErrorExit(PsychRequireNumInputArgs(0));
	PsychErrorExit(PsychCapNumOutputArgs(3));
	
	// Verify eyelink is up and running
	EyelinkSystemIsConnected();
	EyelinkSystemIsInitialized();
	
	if (Verbosity() > 6) mexPrintf("Eyelink: GetQueuedData: boilerplate done\n");
	
	if (PsychCopyInIntegerArg(1, kPsychArgOptional, &eye)) {
		if (eye!=LEFT_EYE && eye!=RIGHT_EYE) {
			PsychErrorExitMsg(PsychErorr_argumentValueOutOfRange, "EyeLink: GetQueuedData:  eye argument must be LEFT_EYE or RIGHT_EYE as returned by EyelinkInitDefaults\n");
		}
		
		TrackerOKForRawValues();
		
		useEye=TRUE;
		numSampleFields=NUM_SAMPLE_FIELDS+NUM_RAW_SAMPLE_FIELDS;
	} else {
		numSampleFields=NUM_SAMPLE_FIELDS;
		useEye=FALSE;
	}
	
	if (Verbosity() > 6) mexPrintf("Eyelink: GetQueuedData: eye/raw chosen\n");
	
	maxSamples = FUDGE_FACTOR*eyelink_data_count(1,0) + 1; //need to allocate at least 1 in case we get a eyelink_get_next_data not predicted by eyelink_data_count
	maxEvents = FUDGE_FACTOR*eyelink_data_count(0,1) + 1;

	samples = (double *)PsychMallocTemp(maxSamples*numSampleFields*sizeof(double)); // according to mario if OOM, ultimately calls to mxCreateNumericArray/mxMalloc will error inside matlab rather than return NULL
	events = (double *)PsychMallocTemp(maxEvents*NUM_EVENT_FIELDS*sizeof(double));
	
	if (Verbosity() > 6) mexPrintf("Eyelink: GetQueuedData: memory allocated for %d samples and %d events\n",maxSamples,maxEvents);
	
	while(!drained && numSamples<maxSamples && numEvents<maxEvents){
		type = eyelink_get_next_data(NULL);
		if (Verbosity() > 6) mexPrintf("Eyelink: GetQueuedData: doing item of type %d\n",type);
		switch(type) {
			case SAMPLE_TYPE:
				if (eyelink_get_float_data(&fs) != type) {
					PsychErrorExitMsg(PsychError_internal, "Eyelink: GetQueuedData: eyelink_get_float_data did not return same sample type as eyelink_get_next_data.");
				}
				if (Verbosity() > 6) mexPrintf("Eyelink: GetQueuedData: get_float called on sample\n");
				if (useEye) {
					if (Verbosity() > 6) mexPrintf("Eyelink: GetQueuedData: calling get_raw\n");
					memset(&fr, 0, sizeof(fr));
					if((err = eyelink_get_extra_raw_values_v2(&fs, eye, &fr))){
						//snprintf(errmsg, ERR_BUFF_LEN, "Eyelink: GetQueuedData: eyelink_get_extra_raw_values_v2 returned error code %d: %s", err, eyelink_get_error(err,"eyelink_get_extra_raw_values_v2"));
						if (Verbosity() > 6) mexPrintf("Eyelink: GetQueuedData: raw value error\n");
						sprintf(errmsg, "Eyelink: GetQueuedData: eyelink_get_extra_raw_values_v2 returned error code %d: %s", err, eyelink_get_error(err,"eyelink_get_extra_raw_values_v2")); //no snprintf in msvs?  bug: buff overflow
						PsychErrorExitMsg(PsychError_internal, errmsg);
					}
				}
				
				index=PsychIndexElementFrom2DArray(numSampleFields, maxSamples, 0, numSamples++);
				samples[index++]=(double)(FLOAT_TIME(&fs)); // 1
				samples[index++]=(double)(fs.type); // 2
				samples[index++]=(double)(fs.flags); // 3
				samples[index++]=(double)(fs.px[0]); // 4
				samples[index++]=(double)(fs.px[1]); // 5
				samples[index++]=(double)(fs.py[0]); // 6
				samples[index++]=(double)(fs.py[1]); // 7
				samples[index++]=(double)(fs.hx[0]); // 8
				samples[index++]=(double)(fs.hx[1]); // 9
				samples[index++]=(double)(fs.hy[0]); // 10
				samples[index++]=(double)(fs.hy[1]); // 11
				samples[index++]=(double)(fs.pa[0]); // 12
				samples[index++]=(double)(fs.pa[1]); // 13
				samples[index++]=(double)(fs.gx[0]); // 14
				samples[index++]=(double)(fs.gx[1]); // 15
				samples[index++]=(double)(fs.gy[0]); // 16
				samples[index++]=(double)(fs.gy[1]); // 17
				samples[index++]=(double)(fs.rx); // 18
				samples[index++]=(double)(fs.ry); // 19
				samples[index++]=(double)(fs.status); // 20
				samples[index++]=(double)(fs.input); // 21
				samples[index++]=(double)(fs.buttons); // 22
				samples[index++]=(double)(fs.htype); // 23
				samples[index++]=(double)(fs.hdata[0]); // 24
				samples[index++]=(double)(fs.hdata[1]); // 25
				samples[index++]=(double)(fs.hdata[2]); // 26
				samples[index++]=(double)(fs.hdata[3]); // 27
				samples[index++]=(double)(fs.hdata[4]); // 28
				samples[index++]=(double)(fs.hdata[5]); // 29
				samples[index++]=(double)(fs.hdata[6]); // 30
				samples[index++]=(double)(fs.hdata[7]); // 31
				
				if (useEye) {
					samples[index++]=(double)(fr.raw_pupil[0]); // 32
					samples[index++]=(double)(fr.raw_pupil[1]); // 33
					samples[index++]=(double)(fr.raw_cr[0]); // 34
					samples[index++]=(double)(fr.raw_cr[1]); // 35
					samples[index++]=(double)(fr.pupil_area); // 36
					samples[index++]=(double)(fr.cr_area); // 37
					samples[index++]=(double)(fr.pupil_dimension[0]); // 38
					samples[index++]=(double)(fr.pupil_dimension[1]); // 39
					samples[index++]=(double)(fr.cr_dimension[0]); // 40
					samples[index++]=(double)(fr.cr_dimension[1]); // 41
					samples[index++]=(double)(fr.window_position[0]); // 42
					samples[index++]=(double)(fr.window_position[1]); // 43
					samples[index++]=(double)(fr.pupil_cr[0]); // 44
					samples[index++]=(double)(fr.pupil_cr[1]); // 45
					samples[index++]=(double)(fr.cr_area2); // 46
					samples[index++]=(double)(fr.raw_cr2[0]); // 47
					samples[index++]=(double)(fr.raw_cr2[1]); // 48
				}
				
				if (Verbosity() > 6) mexPrintf("Eyelink: GetQueuedData: sample copied\n");
				break;
				
			case LOST_DATA_EVENT: // queue overflowed, we are not supposed to call eyelink_get_float_data on this
				index=PsychIndexElementFrom2DArray(numSampleFields, maxSamples, 0, numSamples++);
				for(fieldNum=0; fieldNum<numSampleFields; fieldNum++){
					samples[index++]= (double)((fieldNum==1) ? LOST_DATA_EVENT : MISSING_DATA);
				}
				if (Verbosity() > 6) mexPrintf("Eyelink: GetQueuedData: did lost_data\n");
				break;
				
			case 0: // queue empty
				drained=(PsychNativeBooleanType)TRUE;
				if (Verbosity() > 6) mexPrintf("Eyelink: GetQueuedData: found end of queue\n");
				break;
				
			default: // it is an event
				if (eyelink_get_float_data(&fe) != type) {
					PsychErrorExitMsg(PsychError_internal, "Eyelink: GetQueuedData: eyelink_get_float_data did not return same event type as eyelink_get_next_data.");
				}
				if (Verbosity() > 6) mexPrintf("Eyelink: GetQueuedData: get_float called on event\n");
				index=PsychIndexElementFrom2DArray(NUM_EVENT_FIELDS, maxEvents, 0, numEvents++);
				events[index++]=(double)(fe.time); // 1 %FLOAT_TIME currently a noop on events
				events[index++]=(double)(fe.type); // 2
				events[index++]=(double)(fe.read); // 3
				events[index++]=(double)(fe.eye); // 4
				events[index++]=(double)(fe.sttime); // 5
				events[index++]=(double)(fe.entime); // 6
				events[index++]=(double)(fe.hstx); // 7
				events[index++]=(double)(fe.hsty); // 8
				events[index++]=(double)(fe.gstx); // 9
				events[index++]=(double)(fe.gsty); // 10
				events[index++]=(double)(fe.sta); // 11
				events[index++]=(double)(fe.henx); // 12
				events[index++]=(double)(fe.heny); // 13
				events[index++]=(double)(fe.genx); // 14
				events[index++]=(double)(fe.geny); // 15
				events[index++]=(double)(fe.ena); // 16
				events[index++]=(double)(fe.havx); // 17
				events[index++]=(double)(fe.havy); // 18
				events[index++]=(double)(fe.gavx); // 19
				events[index++]=(double)(fe.gavy); // 20
				events[index++]=(double)(fe.ava); // 21
				events[index++]=(double)(fe.avel); // 22
				events[index++]=(double)(fe.pvel); // 23
				events[index++]=(double)(fe.svel); // 24
				events[index++]=(double)(fe.evel); // 25
				events[index++]=(double)(fe.supd_x); // 26
				events[index++]=(double)(fe.eupd_x); // 27
				events[index++]=(double)(fe.supd_y); // 28
				events[index++]=(double)(fe.eupd_y); // 29
				events[index++]=(double)(fe.status); // 30
		}
	}

	if (Verbosity() > 6) mexPrintf("Eyelink: GetQueuedData: done after dequeueing %d samples (of %d) and %d events (of %d) (drained is %d) (clauses: %d %d %d)\n",numSamples,maxSamples,numEvents,maxEvents,drained,!drained, numSamples<maxSamples, numEvents<maxEvents);
	
	//bug: who frees memory allocated by PsychAllocOutDoubleMatArg if there is no output arg?  according to matlab doc one should not rely on matlab to do it:
	//   "It is more efficient to perform this cleanup in the source MEX-file than to rely on the automatic mechanism."
	//mario reasons that mathworks would have to have done something stupid for this to be true
	
	//these only work for cutting off the unused parts of the buffers because matlab storage is columnwise
	//this is a little fragile - would be nice if PsychCopy* took the dims of the source and handled noncontiguous cases
	//mario says this is purposesly discouraged and inconvenient to encourage columnwise storage and prevent inefficient fragged memcpys
	PsychCopyOutDoubleMatArg(1, kPsychArgOptional, numSampleFields, numSamples, 1, samples);
	PsychCopyOutDoubleMatArg(2, kPsychArgOptional, NUM_EVENT_FIELDS, numEvents, 1, events);
	
	PsychCopyOutBooleanArg(3, kPsychArgOptional, drained);
	
	if (Verbosity() > 6) mexPrintf("Eyelink: GetQueuedData: done with outputs\n");
	
	// PsychFreeTemp(samples); apparently not paradigmatic to do this oneself
	// PsychFreeTemp(events);
	
	return(PsychError_none);
}


psych_bool TrackerOKForRawValues(void) {
	static psych_bool OK=FALSE;
	
	int err,ver[4];
	char errmsg[ERR_BUFF_LEN]="",buf[ERR_BUFF_LEN]="";
	psych_bool tryAgain=TRUE;
	
	if(OK){ //check static value so we don't do the expensive stuff below more than once
		return OK;
	}
	
	while(tryAgain){
		if ((err=eyelink_read_request("link_sample_data"))){
			sprintf(errmsg, "Eyelink: eyelink_read_request returned error code '%d': '%s'", err, eyelink_get_error(err,"eyelink_read_request")); //no snprintf in msvs?  bug: buff overflow
			PsychErrorExitMsg(PsychError_internal, errmsg);
		}
		
		err=NO_REPLY;
		while(err==NO_REPLY){
			err=eyelink_read_reply(buf);
		}
		
		if(err == OK_RESULT){
			if(strlen(buf)>0){
				if (Verbosity() > 6) mexPrintf("Eyelink: TrackerOKForRawValues: link_sample_data is: '%s'\n",buf);
				// currently getting "LEFT, RIGHT, , RAW, HREF, GAZE, GAZERES, AREA, STATUS, INPUT, HMARKER,"
				// what happened to PUPIL?
				
				if(!strstr(buf,"INPUT") || !strstr(buf,"HMARKER")){
					sprintf(errmsg, "Eyelink: your current link_sample_data is '%s', but you must include INPUT and HMARKER to get raw values.  You must also set inputword_is_window = ON.  Use Eyelink('command','link_sample_data = INPUT,HMARKER') and Eyelink('command','inputword_is_window = ON') to do this.",buf); //no snprintf in msvs?  bug: buff overflow
					PsychErrorExitMsg(PsychError_user, errmsg);
				}
				tryAgain=FALSE;
			} else {
				// sent mail to suganthan: eyelink_read_reply(buf) is frequently returning OK_RESULT without setting buf
				// he says send him an eye.log from the tracker host when this happens.
				mexPrintf("Eyelink: TrackerOKForRawValues: eyelink_read_request(\"link_sample_data\") needs to try again\n");
			}
			
		} else {
			sprintf(errmsg, "Eyelink: eyelink_read_reply returned error code '%d': '%s'", err, eyelink_get_error(err,"eyelink_read_reply")); //no snprintf in msvs?  bug: buff overflow
			PsychErrorExitMsg(PsychError_internal, errmsg);
		}
	}
	
    /* disabled this check cuz reading inputword_is_window isn't supported
	if (err=eyelink_read_request("inputword_is_window")){
		sprintf(errmsg, "Eyelink: eyelink_read_request returned error code %d: '%s'", err, eyelink_get_error(err,"eyelink_read_request")); //no snprintf in msvs?  bug: buff overflow
		PsychErrorExitMsg(PsychError_internal, errmsg);
	}
	
	err=NO_REPLY;
	while(err==NO_REPLY){
		err=eyelink_read_reply(buf);
	}
	
	if(err == OK_RESULT){
		if (Verbosity() > 6) mexPrintf("Eyelink: TrackerOKForRawValues: inputword_is_window is: '%s'\n",buf); //must be Eyelink('command','inputword_is_window = ON')
		// currently getting "Variable read not supported"
	 
	 	if(!strstr(buf,"ON")){
			sprintf(errmsg, "Eyelink: your current inputword_is_window is '%s', but must be ON to get raw values.  Use Eyelink('command','inputword_is_window = ON') to do this.",buf); //no snprintf in msvs?  bug: buff overflow
			PsychErrorExitMsg(PsychError_user, errmsg);
		}
	 
	} else {
		sprintf(errmsg, "Eyelink: eyelink_read_reply returned error code '%d': '%s'", err, eyelink_get_error(err,"eyelink_read_reply")); //no snprintf in msvs?  bug: buff overflow
		PsychErrorExitMsg(PsychError_internal, errmsg);
	}
	 */
	
	eyelink_dll_version(buf);
	if (Verbosity() > 6) mexPrintf("Eyelink: TrackerOKForRawValues: eyelink dll version: '%s'\n",buf); //currently "1,8,1,0"
	if(sscanf(buf,"%d,%d,%d,%d",&(ver[0]),&(ver[1]),&(ver[2]),&(ver[3]))==4){
		if(ver[0]<1 || (ver[0]==1 && ver[1]<8) || (ver[0]==1 && ver[1]==8 && ver[2]<1)){
			sprintf(errmsg, "EyeLink: eyelink_dll_version reports '%s' - collecting raw values requires that you install 1.8.1.0 or higher (https://www.sr-support.com/forums/showthread.php?t=6 (windows) or https://www.sr-support.com/forums/showthread.php?t=15 (osx)).",buf); //no snprintf in msvs?  bug: buff overflow
			PsychErrorExitMsg(PsychError_unimplemented, errmsg);
		}
	} else {
		sprintf(errmsg, "Eyelink: couldn't read output of eyelink_dll_version: '%s'",buf);
		PsychErrorExitMsg(PsychError_internal, errmsg);
	}
	
	if(eyelink_get_tracker_version(buf)<3){
		PsychErrorExitMsg(PsychError_unimplemented, "Eyelink: can't get raw values without eyelink 1000 or better");
	} else {
		if (Verbosity() > 6) mexPrintf("Eyelink: TrackerOKForRawValues: tracker version: '%s'\n",buf); //currently "EYELINK CL 4.31"
	}
	if(sscanf(buf,"EYELINK CL %d.%d",&(ver[0]),&(ver[1]))==2){
		if(ver[0]<4 || (ver[0]==4 && ver[1]<31)){
			sprintf(errmsg, "EyeLink: eyelink_get_tracker_version reports '%s' - collecting raw values requires that you install host software 4.31 or higher (https://www.sr-support.com/forums/showthread.php?t=179).",buf); //no snprintf in msvs?  bug: buff overflow
			PsychErrorExitMsg(PsychError_unimplemented,errmsg);
		}
	} else {
		sprintf(errmsg, "Eyelink: couldn't read output of eyelink_get_tracker_version: '%s'",buf);
		PsychErrorExitMsg(PsychError_internal, errmsg);
	}
	
	OK=TRUE;
	return OK;
}
