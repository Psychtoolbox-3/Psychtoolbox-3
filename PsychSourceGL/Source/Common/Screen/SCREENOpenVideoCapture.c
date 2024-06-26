/*
  Psychtoolbox3/Source/Common/SCREENOpenVideoCapture.c		
  
  AUTHORS:
  
  mario.kleiner at tuebingen.mpg.de   mk
  
  PLATFORMS:	
  
  This file should build on any platform. 

  HISTORY:
  
  2/7/06	mk		Created. 
  12/23/07	mk		Extended to allow switching between capture engines on a per-device basis.
 
  DESCRIPTION:
  
  Open a video capture source (framegrabber, webcam, ...), create and initialize a corresponding capture object
  and return a handle to it.
 
  On all systems, the preferred and supported video capture engine is GStreamer.
 
  On Linux and OS/X we support Firewire and USB machine vision cameras which comply to the IIDC-1.x specification,
  via use of the open-source free software library libdc1394 V2. The library itself is most powerful and well
  tested/tuned for GNU/Linux systems, but also well working on OS/X. It has experimental limited support for
  MS-Windows as well, but we don't implement it for Windows yet.
  
  TO DO:
  
*/

#include "Screen.h"

static char useString[] = "videoPtr = Screen('OpenVideoCapture', windowPtr [, deviceIndex][, roirectangle][, pixeldepth][, numbuffers][, allowfallback][, targetmoviename][, recordingflags][, captureEngineType][, bitdepth=8]);";
static char synopsisString[] = 
"Try to open the video source 'deviceIndex' for video capture into onscreen window 'windowPtr' and "
"return a handle 'videoPtr' on success. If 'deviceIndex' is left out, it defaults to zero - use the "
"first capture device attached to your machine. You can get a list of all available video capture "
"devices on your system via a call to Screen('VideoCaptureDevices'). The positive 'deviceIndex' "
"values mentioned there can be used to select a specific video capture device. Certain videocapture "
"devices can't get auto-detected and are therefore not addressable this way. For such devices you "
"must find out the type and unique identifier of the device. Then pass in a negative 'deviceIndex' "
"to tell Screen you want to select a hidden device by its name or id, and pass in the name or id of "
"that hidden device as string in the argument 'targetmoviename', or via the 'SetNextCaptureBinSpec=' "
"option in Screen('SetVideoCaptureParameter'). The special deviceIndex -9 will create a new GStreamer "
"videocapture bin out of a given gst-launch style 'SetNextCaptureBinSpec=' string and attach that as "
"a video source. This allows to create arbitrary GStreamer video capture pipelines.\n"
"'roirectangle' if specified, defines the requested size "
"for captured images. The default is to return the maximum size image provided by the capture device. "
"A 'roirectangle' setting of [0 0 width height] will not define a region of interest, but instead request a "
"video capture resolution of width x height pixels, instead of the default maximum resolution. Settings "
"of [left top right bottom] will leave video capture at maximum resolution, but crop the images to the "
"rectangular subregion as defined by the given left, right, top and bottom borders. Such ROI's are "
"only applied to returned video images, not to recorded video by default. See the 'recordingflags' "
"settings below on how to adjust this behaviour to your needs.\n"
"The real ROI (region of interest) may differ from the requested one, depending on the capabilities of "
"your capture device. 'pixeldepth' if provided, asks for the number of layers that captured textures "
"should have: 1=Luminance image, 2=Luminance+Alpha image, 3=RGB image, 4=RGB+Alpha, 5=YCBCR, 6=I420 image. Default is "
"to take whatever the capture device provides by default. Different devices support different formats "
"so some of these settings may be ignored. Some combinations of video capture devices and graphics "
"cards may support a setting of 5=YCBCR encoding. If they do, then this is an especially efficient way "
"to handle color images, which may result in lower cpu load and higher framerates. A format of 6=YUV-I420 "
"should be supported by all modern graphics cards and may provide some performance benefits, but your mileage "
"may vary. If you need very fast color image capture, try formats 4, 5 and 6 and see which one gives the "
"best performance for your setup.\n"
"'numbuffers' if provided, specifies the number of internal "
"video buffers to use. It defaults to a value that is optimal for your specific hardware for common use. "
"'allowfallback' if set to 1, will allow Psychtoolbox to use a less efficient mode of operation for video "
"capture if your specific hardware or operating system setup doesn't allow to use the high-performance "
"mode. 'allowfallback' defaults to 1 = Allow fallback path.\n"
"'targetmoviename' If you provide a filename for this argument, PTB will record the captured video to the "
"specified  movie file on your filesystem. PTB will use a default video codec for encoding the "
"video stream. If you want to use a specific codec, you can extend the targetmoviename by a string of "
"format :CodecType=xxx , where xxx is the numeric type id or name of the codec.\n"
"Please read 'help VideoRecording' for many more options for tweaking the video recording "
"process via the 'targetmoviename' parameter.\n"
"'recordingflags' specify the behaviour of harddisc-recording and some other capture operations. Please "
"note that not all flags are supported on all capture engines, cameras and operating systems. Unsupported "
"flags will be silently ignored. Most flags are only supported with GStreamer at the moment: "
"0 (default) = Only record video. 2 = Record audio track as well. The value 1 (or 1+2) asks "
"PTB to first record into system memory, and only write the movie file after capture has been stopped. "
"This allows for higher capture framerates, but is limited in recording time by installed memory. Also, "
"this mode currently can sometimes cause hangs and crashes of PTB for unknown reasons, so better avoid!\n"
"A setting of 4 will only enable recording, but no return of captured data, i.e., just record to disk. "
"A setting of 8 will avoid some calls that are supposed to provide better realtime behaviour, but may "
"cause some problems with some video codecs when recording to disk. A setting of 16 will perform most "
"of the heavy work on a separate parallel background thread, utilizing multi-core machines better.\n"
"A setting of 32 will try to select the highest quality codec for texture creation from captured video, "
"instead of the normal quality codec. A setting of 64 will return capture timestamps in the time base of "
"the video engine (e.g., elapsed time since start of capture, or recording time in movie), instead of "
"the default time base, which is regular GetSecs() time.\n"
"A setting of 128 will force use of a videorate converter in pure live capture mode. By default the "
"videorate converter is only used if video recording is active, unless flag 8192 is specified. The "
"converter makes sure that video is "
"recorded (or delivered) at exactly the requested capture framerate, even if the system isn't really "
"capable of maintaining that framerate: If the video source (camera) delivers frames at a too low "
"framerate, the converter will insert duplicated frames to boost up effective framerate. If the source "
"delivers more frames than the engine can handle (e.g., system overload or video encoding too slow) "
"the converter will drop frames to reduce effective framerate. Slight fluctuations are compensated by "
"adjusting the capture timestamps. This mechanism guarantees a constant framerate in recorded video "
"as well as the best possible audio-video sync and smoothness of video, given system constraints. "
"The downside may be that the recorded content and returned timestamps don't reflect the true timing "
"of capture, but a beautified version. In pure live capture, rate conversion is off by default to avoid "
"such potential confounds in the timestamps. Choose this options carefully.\n"
"A setting of 256 in combined video live capture and video recording mode will restrict video framerate "
"conversion to the recorded videostream, but provide mostly untampered true timing to the live capture. "
"By default, framerate conversion applies to recording and live feedback if video recording is enabled.\n"
"A setting of 512 requests that ROI's as defined by the 'roirectangle' parameter get also applied to recorded video. "
"Without this setting, ROI's only apply to live video as returned by Screen('GetCapturedImage',...);\n"
"A setting of 1024 disables application of ROI's to live video as returned by Screen('GetCapturedImage',...);\n"
"A setting of 2048 requests immediate conversion of video textures into a format suitable as offscreen window, "
"for use with Screen('TransformTexture') or for drawing with custom user provided GLSL shaders. Normally this "
"happens automatically on first use, asking for it explicitely may have performance or convenience benefits.\n"
"A setting of 4096 requests to apply some performance optimizations (the setting of filter-caps). This can hurt "
"if a videocapture device refuses to work, with some error message about ''check your filtered caps, if any.''. "
"By default, if the flag is omitted, some performance loss will be present, but capture will be more robust "
"with problematic cameras. "
"A setting of 8192 requests to avoid any use of videorate converters, not even for recording (see above).\n"
"\n"
"'captureEngineType' This optional parameter allows selection of the video capture engine to use for this "
"video source. Allowable values are currently 1 and 3. "
"A value of 1 selects Firewire video capture via the free software library libdc1394-V2. "
"That engine only supports high performance IIDC machine vision cameras that are compliant with the "
"IIDC-1.x standard and are connected via a Firewire (IEEE-1394) bus system or via USB-3 USBVision standard. "
"Use of the engine with such cams allows for much higher flexibility and performance than use of video capture "
"via GStreamer, however one restriction is that sound recording isn't yet supported with this capture engine. "
"This IIDC capture engine is supported on Linux and to a lesser degree on MacOSX, but not on MS-Windows.\n\n"
"A value of 3 selects the GStreamer video capture engine. This engine is supported on all operating systems "
"and allows for video and sound recording of captured video and audio streams. Type 'help GStreamer' for "
"installation and setup instructions for the required GStreamer runtime libraries.\n\n"
"If you don't specify 'captureEngineType', the global "
"setting from Screen('Preference', 'DefaultVideoCaptureEngine') will be used. If you don't specify that either "
"then engine selection will default to GStreamer.\n\n"
"To summarize: \n"
"GStreamer: Is the engine of choice for all operating systems and with most standard applications.\n"
"Firewire engine: Supports only Firewire machine vision cameras, but allows free selection among all connected "
"cameras, simultaneous operation of many cameras, low latency, high framerates and reliability, precise timestamping "
"and low level access to many special features of such cameras, e.g., gain-, shutter-, exposure-, trigger controls etc.\n\n"
"'bitdepth' Optional parameter to ask for video capture in a certain color or luminance resolution, a certain "
"number of bits per color or luminance component, also known as bpc. Defaults to 8 bpc if omitted, ie., 8 bits "
"or 1 Byte resolution per luminance or color channel for classic 256 levels of intensity. Lower values are unsupported "
"and will get rounded up to 8 bpc. Higher values may be supported by some higher end professional class cameras. If "
"you ask for an unsupported value, the engine will try to get the lowest supported value that matches or exceeds "
"what you want. Currently the libdc1394 engine for pro-class IIDC compliant firewire or USB machine vision cameras "
"supports bitdepth > 8 bpc on capable cameras. The GStreamer engine always supports 8 bpc and it may support 16 bpc on "
"capable cameras. Please note that requesting > 8 bpc will cause a substantial increase in both video bus "
"bandwidth and memory consumption: Twice the bus bandwidth and two to four times the amount of memory per video frame, "
"so tread carefully.\n";

static char seeAlsoString[] = "CloseVideoCapture StartVideoCapture StopVideoCapture GetCapturedImage";

PsychError SCREENOpenVideoCapture(void) 
{
	PsychWindowRecordType                   *windowRecord;
	int                                     deviceIndex;
	int                                     capturehandle = -1;
	PsychRectType                           roirectangle;
	psych_bool                              roiassigned;
	int                                     reqdepth = 0;
	int                                     num_dmabuffers = 0;
	int                                     allow_lowperf_fallback = 1;
	char*                                   moviename;
	int                                     recordingflags;
	int                                     engineId;
    int                                     bitdepth = 8;
	
	// All sub functions should have these two lines
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

	PsychErrorExit(PsychCapNumInputArgs(10));           // Max. 10 input args.
	PsychErrorExit(PsychRequireNumInputArgs(1));        // Min. 1 input args required.
	PsychErrorExit(PsychCapNumOutputArgs(1));           // Max. 1 output args.
	
	// Get the window record from the window record argument and get info from the window record
	PsychAllocInWindowRecordArg(kPsychUseDefaultArgPosition, TRUE, &windowRecord);
	// Only onscreen windows allowed:
	if(!PsychIsOnscreenWindow(windowRecord)) {
		PsychErrorExitMsg(PsychError_user, "OpenVideoCapture called on something else than an onscreen window.");
	}
	
	// Get the device index. We default to the first device if none is given:
	deviceIndex=0;
	PsychCopyInIntegerArg(2, FALSE, &deviceIndex);
	
	// Get the optional roi rectangle:
	roiassigned = PsychCopyInRectArg(3, FALSE, roirectangle);
	if (roiassigned && IsPsychRectEmpty(roirectangle)) PsychErrorExitMsg(PsychError_user, "OpenVideoCapture called with invalid (empty) roirectangle argument.");
	
	// Query (optional) output texture pixel depth: By default, we take whatever we get from the capture device:
	PsychCopyInIntegerArg(4, FALSE, &reqdepth);
	if (reqdepth<0) PsychErrorExitMsg(PsychError_user, "OpenVideoCapture called with invalid (negative) pixeldepth argument.");
	// Query number of ringbuffers to use. Our default is coded in the OS dependent subroutines.
	PsychCopyInIntegerArg(5, FALSE, &num_dmabuffers);
	if (num_dmabuffers<0) PsychErrorExitMsg(PsychError_user, "OpenVideoCapture called with invalid (negative) numbuffers argument.");
	// Query, if use of low-performance fallback code is allowed if high-perf engine fails:
	PsychCopyInIntegerArg(6, FALSE, &allow_lowperf_fallback);
	if (allow_lowperf_fallback<0) PsychErrorExitMsg(PsychError_user, "OpenVideoCapture called with invalid (negative) allowfallback flag.");
	
	// Query optional moviename for recording the grabbed video into a movie file:
	moviename = NULL;
	PsychAllocInCharArg(7, FALSE, &moviename);
	
	// Query optional movie recording flags:
	// 0 = Record video, stream to disk immediately (slower, but unlimited recording duration).
	// 1 = Record video, stream to memory, then at end of recording to disk (limited duration by RAM size, but faster).
	// 2 = Record audio as well.
	// 4 = Do not return capture data via Screen('GetCapturedImage') during video recording to disc.
	// 8 = Avoid some performance optimizations which may cause trouble with some codecs.
	// 16= Use multi-threading for automatic background processing and cpu offloading.
	// 32= Return high quality textures via 'GetCapturedImage' if recording in parallel --> Quality tradeoff live feed vs. recording.
	// 64= Return timestamps in engine time instead of GetSecs() time.
    // 128 = Use videorate converter even in pure live-capture mode.
    // 256 = Restrict use of videorate converter to recorded video, do not apply to live feed.
    // 512 = Apply ROI to recorded video as well, not only to live video.
    // 1024 = Do not apply ROI to live video feed.
    // 2048 = Do normalize texture format and orientation of video texture at creation time.
	recordingflags = 0;
	PsychCopyInIntegerArg(8, FALSE, &recordingflags);
	
	// Copy in optional id of video capture engine to use. We default to the one set via the Screen('Preference', 'DefaultVideocaptureEngine');
	// setting, which by itself defaults to 3 (GStreamer) for now.
	engineId = PsychPrefStateGet_VideoCaptureEngine();
	PsychCopyInIntegerArg(9, FALSE, &engineId);
	if (engineId < 0 || engineId > 3) PsychErrorExitMsg(PsychError_user, "OpenVideoCapture called with invalid 'captureEngineType'. Valid are 0,1,2,3.");

	if (engineId == 0 || engineId == 2) {
		printf("\n\n");
		printf("PTB-INFO: Your script explicitely requests use of video capture engine type %i - the %s video capture engine.\n", engineId, (engineId == 0) ? "Quicktime" : "ARVideo");
		printf("PTB-INFO: This engine has been permanently disabled and removed from Psychtoolbox.\n");
		printf("PTB-INFO: We recommend use of the GStreamer video capture engine (engine type 3) as a technically superior replacement.\n");
		printf("PTB-INFO: In most cases, the selected replacement should work without need for any further changes to your code.\n\n");
	}

	// Copy in optional bitdepth, aka bpc setting for capture resolution: Defaults to 8 bpc.
	PsychCopyInIntegerArg(10, FALSE, &bitdepth);
	
	// Try to open the capture device and create & initialize a corresponding capture object.
	// A MATLAB handle to the video capture object is returned upon successfull operation.
	if (roiassigned) {
        PsychOpenVideoCaptureDevice(engineId, windowRecord, deviceIndex, &capturehandle, roirectangle, reqdepth, num_dmabuffers, allow_lowperf_fallback, moviename, recordingflags, bitdepth);
	}
	else {
        PsychOpenVideoCaptureDevice(engineId, windowRecord, deviceIndex, &capturehandle, NULL, reqdepth, num_dmabuffers, allow_lowperf_fallback, moviename, recordingflags, bitdepth);
	}
	
	// Upon sucessfull completion, we'll have a valid handle in 'capturehandle'. Return it to Matlab-world:
	PsychCopyOutDoubleArg(1, TRUE, (double) capturehandle);
	
	// Ready!
	return(PsychError_none);
}

static char useString2[] = "devices = Screen('VideoCaptureDevices' [, engineId]);";
static char synopsisString2[] = 
"Enumerate all available video devices for a given videocapture engine id 'engineId', or "
"for the default engine, if none is given. Returns a struct array with one slot per available "
"device. Entries of each struct in the array are specific to the selected capture engine, "
"except for the entry 'DeviceIndex' which provides a handle that you could pass to "
"Screen('OpenVideoCapture', windowPtr, deviceIndex, ...); as the 'deviceIndex' argument "
"to select the video device corresponding to a given slot.\n"
"The function may return an empty array if no video capture devices could be detected.\n";

PsychError SCREENVideoCaptureDevices(void) 
{
	int engineId;
	
	// All sub functions should have these two lines
	PsychPushHelp(useString2, synopsisString2, seeAlsoString);
	if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none);};

	PsychErrorExit(PsychCapNumInputArgs(1));            // Max. 1 input args.
	PsychErrorExit(PsychRequireNumInputArgs(0));        // Min. 1 input args required.
	PsychErrorExit(PsychCapNumOutputArgs(1));           // Max. 1 output args.
	
	engineId = PsychPrefStateGet_VideoCaptureEngine();
	PsychCopyInIntegerArg(1, FALSE, &engineId);
	
	// Do actual enumeration for the given engineId:
	PsychEnumerateVideoSources(engineId, 1);

	// Ready!
	return(PsychError_none);	
}
