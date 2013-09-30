/*
 
 Psychtoolbox3/Source/Common/SCREENStartVideoCapture.c		
 
 AUTHORS:
 
 mario.kleiner at tuebingen.mpg.de   mk
 
 PLATFORMS:	
 
 This file should build on any platform. 
 
 HISTORY:
 2/7/06  mk		Created. 
 
 DESCRIPTION:
 
 Start a previously opened video capture device.
 
 TO DO:
 
 */

#include "Screen.h"

static char useString[] = "[fps starttime] = Screen('StartVideoCapture', capturePtr [, captureRateFPS=25] [, dropframes=0] [, startAt]);";
static char synopsisString[] =  "Start video capture device specified by 'capturePtr'. If 'captureRateFPS' "
                                "is provided, the device is requested to capture at that rate. Otherwise it "
                                "is requested to operate at 25 Hz. The real capture rate may differ from the "
                                "requested one, depending on the capabilities of your capture device. For professional "
                                "firewire cameras when used with the IIDC-1394 videocapture engine, PTB tries "
                                "to select the closest supported framerate that is at least as high as the requested "
                                "one. For other videocapture engines, e.g., GStreamer, the result is "
                                "camera dependent. It may choose the closest supported rate, or the maximum supported "
                                "framerate. The real capture rate is returned as return argument 'fps'.\n"
                                "If you pass the special value realmax, then PTB will try to select the highest "
                                "framerate that your capture device supports at the given video resolution and other settings.\n"
                                "Some capture devices may only support a small set of framerates and may fail to start "
                                "completely if you or PTB on your behalf request a unsupported 'captureRateFPS'. Some "
                                "devices ignore the requested framerate completely and run at whatever rate they like, e.g., "
                                "Apples iSight cameras adjust their framerate depending on lighting conditions. You may "
                                "need to tinker a bit here, the auto-selection of framerates is not fool-proof.\n"
                                "If 'dropframes' is provided and set to 1, "
                                "the device is requested to always deliver the most recently acquired frame, "
                                "dropping previously captured but not delivered frames if necessary. The default "
                                "is to queue up as many frames as possible. If you want to do video recordings, you "
                                "want to have the default of zero. If you want to do interactive realtime processing "
                                "of video data (e.g, video feedback for action-perception studies or build your own "
                                "low-cost eyetracker), then you want to use dropframes=1 for lowest possible latency. "
                                "'startAt' This optional argument provides a requested start time in seconds system time, "
                                "e.g., values returned by GetSecs() or Screen('Flip'), when capture should really start. "
                                "Psychtoolbox will wait until that point in time has elapsed before really triggering start "
                                "of capture. This provides a means of soft-synchronizing start of capture with external events. "
                                "The optional return value 'starttime' contains the best estimate (in system time) of when "
                                "video capture was really started. Accuracy of 'startAt' and 'starttime' is highly dependent "
                                "on operating system and capture device. It is assumed to be very accurate on Linux with the "
                                "dc1394 capture engine and professional class IIDC 1394 firewire cameras. No such guarantees "
                                "are possible on other operating systems or with other hardware. Specifically standard consumer "
                                "DV cameras or USB based webcams or built in cameras in Laptops (iSight et al.) may have rather "
                                "unpredictable and variable latency and start timing. Choose a pro operating system and pro "
                                "equipment if you really need timing precision!\n";

static char seeAlsoString[] = "CloseVideoCapture StartVideoCapture StopVideoCapture GetCapturedImage";
	 
PsychError SCREENStartVideoCapture(void) 
{
    int capturehandle = -1;
    double captureFPS = 25;
    int dropframes = 0;
    double starttime = 0;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};
    PsychErrorExit(PsychCapNumInputArgs(4));            // Max. 4 input args.
    PsychErrorExit(PsychRequireNumInputArgs(1));        // Min. 1 input args required.
    PsychErrorExit(PsychCapNumOutputArgs(2));           // Max. 2 output args.
    
    // Get the handle:
    PsychCopyInIntegerArg(1, TRUE, &capturehandle);
    if (capturehandle==-1) {
        PsychErrorExitMsg(PsychError_user, "StartVideoCapture called without valid handle to a capture object.");
    }
    
    PsychCopyInDoubleArg(2, FALSE, &captureFPS);
    if (captureFPS<=0) {
        PsychErrorExitMsg(PsychError_user, "StartVideoCapture called with a negative capture rate.");
    }

    PsychCopyInIntegerArg(3, FALSE, &dropframes);
    if (dropframes<0) {
        PsychErrorExitMsg(PsychError_user, "StartVideoCapture called with invalid (negative) dropframes - argument.");
    }

    PsychCopyInDoubleArg(4, FALSE, &starttime);
    if (starttime<0) {
        PsychErrorExitMsg(PsychError_user, "StartVideoCapture called with invalid (negative) startAt - argument.");
    }

    // Try to start capture:
    captureFPS = (double) PsychVideoCaptureRate(capturehandle, captureFPS, dropframes, &starttime);

    PsychCopyOutDoubleArg(1, FALSE, captureFPS);
    PsychCopyOutDoubleArg(2, FALSE, starttime);

    // Ready!    
    return(PsychError_none);
}
