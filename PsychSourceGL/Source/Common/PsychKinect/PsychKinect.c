/*
	PsychToolbox3/Source/Common/PsychKinect/PsychKinect.h
 
	PROJECTS: PsychKinect only.
 
	AUTHORS:
 
	mario.kleiner@tuebingen.mpg.de	mk
 
	PLATFORMS:	All.
 
	HISTORY:
 
	24.11.2010  mk		Created.
	03.04.2011  mk		Make 64-bit clean.
	14.02.2012  mk		Make Linux & OS/X version compatible to libfreenect 0.1.2

	DESCRIPTION:
 
	A Psychtoolbox driver for Microsoft's Kinect 3D-camera, based
	on the free software code of the http://openkinect.org project.
 
 */

#include "PsychKinect.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>

// Include header of libfreeenect:
#include "libfreenect.h"

#if PSYCH_SYSTEM != PSYCH_WINDOWS
typedef void freenect_depth;
typedef void freenect_pixel;
#define freenect_set_rgb_callback freenect_set_video_callback
#define freenect_start_rgb freenect_start_video
#define freenect_stop_rgb freenect_stop_video
#endif

// For Windows build:
#ifndef FREENECT_VIDEO_IR_8BIT
#define FREENECT_VIDEO_IR_8BIT 2
#endif

// Number of maximum simultaneously open kinect devices:
#define MAX_PSYCH_KINECT_DEVS 10
#define MAX_SYNOPSIS_STRINGS 40  

//declare variables local to this file.  
static const char *synopsisSYNOPSIS[MAX_SYNOPSIS_STRINGS];

// Temporary buffer for conversion of depth images into color coded RGB images:
double*		zmap = NULL;
unsigned char	gl_depth_back[640 * 480 * 4];
unsigned short	t_gamma[2048];

// Cached values for depth conversion, to speed up per-pixel computation:
double          depthBaseAndOffset[2];

typedef double (*CALCPROC)(int); 
CALCPROC calcz = NULL;

freenect_context	*f_ctx = NULL;
freenect_device		*f_dev;

typedef struct PsychKNBuffer {
	unsigned short	depth[640 * 480];		// 16-bit disparity pixels.
	unsigned char*	color;				// RGB8 or raw non-bayer-filtered 1 Byte color pixels.
	double		cts;				// Cached captured sample onset time from paCallback.
	int		cwidth, cheight;		// color buffer width x height.
	int		dwidth, dheight;		// color buffer width x height.
} PsychKNBuffer;

// Our device record:
typedef struct PsychKNDevice {
	psych_mutex	mutex;				// Mutex lock for the PsychKNDevice struct.
	psych_condition changeSignal;		        // Condition variable or event object for change signalling (see above).
	psych_thread captureThread;			// Processing thread.
	freenect_device *dev;				// Handle to USB device representing the Kinect.
	PsychKNBuffer* buffers;				// Bufferchain with capture buffers.
	int numbuffers;					// Number of allocated buffers.
	int dwidth, dheight, dsize;			// Width x Height and size in bytes of depth buffer.
	int cwidth, cheight, csize;			// Width x Height and size in bytes of color buffer.
	int dropframes;					// dropframes = 1: Stall capture if FIFO full, 0 = Overwrite oldest frames.
	int bayerFilterMode;				// 0 = Don't: Fetch unfiltered sensor data. 1 = Let libfreenect do Bayer filtering into RGB. 2 = ...
	double	captureStartTime;			// Time of start of capture.
	volatile unsigned int state;		        // Current state: 0=Stopped, 1=Running.
	volatile unsigned int reqstate;		        // Requested state of the stream, as opposed to current 'state'. Written by main-thread.
	unsigned int recposition;			// Current record position in frames since start of capture.
	unsigned int readposition;			// Last read-out frame since start of capture.
	unsigned int frame_valid;			// Current readposition valid?
	unsigned int xruns;				// Number of over-/underflows.
	unsigned int paCalls;				// Number of callback invocations.
	double   R[3][3];                               // Extrinsic (R)otation, (T)ranslation of video camera wrt. depths camera. 
	double   T[3];
	double fx_d;                                    // Depths camera intrinsic parameters. Focal length x/y fx_, fy_.
	double fy_d;
	double cx_d;                                    // Optical center cx_, cy_
	double cy_d;
	double fx_rgb;                                  // RGB video camera intrinsic parameters.
	double fy_rgb;
	double cx_rgb;
	double cy_rgb;
	double undistort_d[5];                          // Optical distortion coefficients depth: k1, k2, p1, p2, k3.
	double undistort_rgb[5];                        // Optical distortion coefficients rgb  : k1, k2, p1, p2, k3.
	double depthBaseAndOffset[2];                   // Base and Offset parameter for mapping of raw depth sensor data to physical distance units.
} PsychKNDevice;

PsychKNDevice kinectdevices[MAX_PSYCH_KINECT_DEVS];
unsigned int  devicecount = 0;
unsigned int  verbosity = 3;
psych_bool initialized = FALSE;

void InitializeSynopsis(void)
{
	int i=0;
	const char **synopsis = synopsisSYNOPSIS;  //abbreviate the long name
	
	synopsis[i++] = "PsychKinectCore - A Psychtoolbox driver for the Microsoft Kinect box.\n";
	synopsis[i++] = "This driver allows to control the box and grab color images and depth";
	synopsis[i++] = "images from the Kinect depth camera.\n";
	synopsis[i++] = "It uses and requires the free software drivers and libraries from the OpenKinect";
	synopsis[i++] = "project (Thanks!): http://openkinect.org\n";
	synopsis[i++] = "Libfreenect is Copyright (C) 2010 - 2012 individual OpenKinect contributors.";
	synopsis[i++] = "Libfreenect requires libusb-1.0, which is licensed under LGPL v2 or later.";
	synopsis[i++] = "See 'help InstallKinect' for more detailed license information.\n";
	synopsis[i++] = "The PsychKinectCore driver is licensed to you under the terms of the MIT license,";
	synopsis[i++] = "with some restrictions imposed when used with Matlab. See 'help License.txt' in";
	synopsis[i++] = "the Psychtoolbox root folder for more details.\n";
	synopsis[i++] = "The driver also uses bits of math inspired by Nicolas Burrus Kinect work:";
	synopsis[i++] = "http://nicolas.burrus.name/index.php/Research/KinectCalibration ";
	synopsis[i++] = "\n";
	
	synopsis[i++] = "Usage:";
	
	// synopsis[i++] = "devices = Screen('VideoCaptureDevices' [, engineId]);";
	synopsis[i++] = "kinectPtr = PsychKinect('Open' [, deviceIndex=0][, numbuffers=2]][, bayerFilterMode=1]);";
	synopsis[i++] = "PsychKinect('Close', kinectPtr);";
	synopsis[i++] = "PsychKinect('SetAngle', kinectPtr, angle);";
	synopsis[i++] = "[starttime, fps, cwidth, cheight, dwidth, dheight] = PsychKinect('Start', kinectPtr [, dropframes=0]);";
	synopsis[i++] = "PsychKinect('Stop', kinectPtr);";
	synopsis[i++] = "status = PsychKinect('GetStatus', kinectPtr);";
	synopsis[i++] = "[...old parameter settings in order of inputs... ] = PsychKinect('SetBaseCalibration', kinectPtr, depthsIntrinsics, rgbIntrinsics, rgbRotation, rgbTranslation, depthsUndistort, rgbUndistort, depthBaseAndOffset);";
	synopsis[i++] = "[result, cts, age] = PsychKinect('GrabFrame', kinectPtr [, waitMode=1][, mostrecent=0]);";
	synopsis[i++] = "PsychKinect('ReleaseFrame', kinectPtr);";
	synopsis[i++] = "[imageOrPtr, width, height, channels, extType, extFormat] = PsychKinect('GetImage', kinectPtr [, imtype=0][, returnTexturePtr=0]);";
	synopsis[i++] =	"[imageOrPtr, width, height, extType, extFormat] = PsychKinect('GetDepthImage', kinectPtr [, format=0][, returnTexturePtr=0]);";
	synopsis[i++] = NULL;  //this tells PsychKinectDisplaySynopsis where to stop
	
	if (i > MAX_SYNOPSIS_STRINGS) {
		PrintfExit("%s: increase dimension of synopsis[] from %ld to at least %ld and recompile.", __FILE__, (long)MAX_SYNOPSIS_STRINGS,(long)i);
	}
}

PsychError PsychKinectDisplaySynopsis(void)
{
	int i;
	
	for (i = 0; synopsisSYNOPSIS[i] != NULL; i++)
		printf("%s\n",synopsisSYNOPSIS[i]);
	
	return(PsychError_none);
}

PsychKNBuffer* PsychGetKNBuffer(PsychKNDevice* kinect, int index)
{
	return(&kinect->buffers[index % kinect->numbuffers]);
}

void PsychDepthCB(freenect_device *dev, freenect_depth *depth, uint32_t timestamp)
{
	PsychKNDevice *kinect = (PsychKNDevice*) freenect_get_user(dev);
	PsychKNBuffer *buffer = PsychGetKNBuffer(kinect, kinect->recposition);

	if (kinect->state == 0) return;
	if (kinect->paCalls & 0x1) return;
	
	memcpy(buffer->depth, depth, kinect->dsize);
	buffer->dwidth = kinect->dwidth;
	buffer->dheight = kinect->dheight;
	
	// Mark depth frame as done:
	kinect->paCalls |= 0x1;
}

void PsychRGBCB(freenect_device *dev, freenect_pixel *rgb, uint32_t timestamp)
{
	PsychKNDevice *kinect = (PsychKNDevice*) freenect_get_user(dev);
	PsychKNBuffer *buffer = PsychGetKNBuffer(kinect, kinect->recposition);

	if (kinect->state == 0) return;
	if (kinect->paCalls & 0x2) return;

	memcpy(buffer->color, rgb, kinect->csize);
	buffer->cwidth = kinect->cwidth;
	buffer->cheight = kinect->cheight;
	
	// Timestamp the buffer:
	buffer->cts = (double) timestamp;
	
	// Mark RGB frame available:
	kinect->paCalls |= 0x2;
}

void* PsychKinectThreadMain(volatile void* deviceToCast)
{
	PsychKNDevice *kinect = (PsychKNDevice*) deviceToCast;
	int abort = FALSE;
	int headroom;

	// Child protection:
	if ((NULL == kinect) || (NULL == kinect->dev)) return(NULL);
	
	// Start kinect's iso streaming:
	freenect_start_depth(kinect->dev);
	freenect_start_rgb(kinect->dev);
	
	// Main processing loop:
	while(!abort) {
		// Check for thread cancellation:
		PsychTestCancelThread(&kinect->captureThread);
		
		// Check again for cancellation, this time the "soft way":
		PsychLockMutex(&kinect->mutex);
		if (kinect->reqstate != 1) {
			abort = TRUE;
			PsychUnlockMutex(&kinect->mutex);
			continue;
		};
		
		headroom = kinect->numbuffers - (kinect->recposition - kinect->readposition);
		
		PsychUnlockMutex(&kinect->mutex);
		
		// Do something?
		if (!abort && ((headroom > 0) || (kinect->dropframes))) {
			// xrun due to oldest frame gets overwritten?
			if (headroom <= 0) kinect->xruns++;
			
			// Kick off USB event handling and our callbacks:
			if (freenect_process_events(f_ctx) < 0) {
				// Error condition! We better abort!
				abort = TRUE;
				printf("PTB-CRITICAL: Error during USB data receive from Kinect! Aborting. Prepare for trouble!\n");
				continue;
			}
			
			// Check if we have a completely new video & depth frame:
			if (kinect->paCalls >= 2+1) {
				// Yes. Reset counter:
				kinect->paCalls = 0;
				
				// Increment our write pointer to finalize this frame and prepare capture of a new one:
				PsychLockMutex(&kinect->mutex);
				kinect->recposition++;
				
				// Signal good news to main thread:
				PsychSignalCondition(&kinect->changeSignal);
				PsychUnlockMutex(&kinect->mutex);
			}
		}
		else {
			// FIFO full and we should retain the oldest frames. Just idle
			// until capacity becomes available:
			PsychYieldIntervalSeconds(0.002);
			
			// xrun due to current frame gets dropped:
			kinect->xruns++;
		}
		
		// Repeat processing loop:
	}
	
	// Stop capture: Unless we're on Linux or OSX where libfreenect-0.1.2
	// (or its specific use of libusb-1.0) has some bug and calling this
    // would cause a hang in the device close routine later...
	#if (PSYCH_SYSTEM != PSYCH_LINUX) && (PSYCH_SYSTEM != PSYCH_OSX)
	freenect_stop_depth(kinect->dev);
	freenect_stop_rgb(kinect->dev);
	#endif

	// Stop of capture:
	PsychLockMutex(&kinect->mutex);
	kinect->reqstate = 0;
	kinect->state = 0;
	kinect->paCalls = 0;
	PsychUnlockMutex(&kinect->mutex);
	
	return(NULL);
}

PsychKNDevice* PsychGetKinect(int handle, psych_bool dontfail)
{
	if (handle < 0 || handle >= MAX_PSYCH_KINECT_DEVS || kinectdevices[handle].dev == NULL) {
		if (!dontfail) {
			printf("PTB-ERROR: Invalid Kinect device handle %i passed. No such device open.\n");
			PsychErrorExitMsg(PsychError_user, "Invalid kinect handle.");
		}
		
		return(NULL);
	}
	
	return(&(kinectdevices[handle]));
}

void PsychKNStop(int handle)
{
	PsychKNDevice* kinect;
	kinect = PsychGetKinect(handle, TRUE);
	if (NULL == kinect) return;
	
	PsychLockMutex(&kinect->mutex);
	
	// Stop capture if it is still running:
	if (kinect->state == 1) {
		// Send abort request to thread:
		kinect->reqstate = 0;
		
		PsychUnlockMutex(&kinect->mutex);
		PsychYieldIntervalSeconds(0.1);
		
		PsychAbortThread(&kinect->captureThread);
		
		// Wait for it to die peacefully:
		PsychDeleteThread(&kinect->captureThread);
		kinect->captureThread = 0;
		
		// Mark as dead:
		kinect->state = 0;
	}
	else {
		PsychUnlockMutex(&kinect->mutex);
	}
	
	return;
}

void PsychKNClose(int handle)
{
	int i;
	PsychKNDevice* kinect;
	kinect = PsychGetKinect(handle, TRUE);
	if (NULL == kinect) return;
	
	// Stop device:
	PsychKNStop(handle);
	
	// Capture is stopped. Release mutexes and other resources:
	PsychDestroyCondition(&kinect->changeSignal);
	PsychDestroyMutex(&kinect->mutex);
	
	if (kinect->buffers) {
		for (i = 0; i < kinect->numbuffers; i++) free(kinect->buffers[i].color);
		free(kinect->buffers);
		kinect->buffers = NULL;
	}
	
	// Close usb connection:
	freenect_close_device(kinect->dev);
	kinect->dev = NULL;
	
	// Done with this device:
	devicecount--;
	
	// Last device closed? If so, shutdown driver:
	if (devicecount <= 0) {
		devicecount = 0;
		initialized = FALSE;
		freenect_shutdown(f_ctx);
		f_ctx = NULL;
	}
}

void PsychKNInit(void) {
	int handle;
	int i;
	float v;
	
	for (handle = 0 ; handle < MAX_PSYCH_KINECT_DEVS; handle++) kinectdevices[handle].dev = NULL;
	devicecount = 0;
	initialized = FALSE;
	
	// Build depths gamma table:
	for (i=0; i<2048; i++) {
		v = (float) i/2048.0f;
		v = powf(v, 3.0f)* 6.0f;
		t_gamma[i] = (unsigned short) (v*6*256);
	}
}

PsychError PsychKNShutdown(void) {
	int handle;
	
	if (initialized) {
		for (handle = 0 ; handle < MAX_PSYCH_KINECT_DEVS; handle++) PsychKNClose(handle);
		free(zmap);
		zmap = NULL;
	}
	
	initialized = FALSE;

	return(PsychError_none);
}

PsychError PSYCHKINECTOpen(void)
{
    static char useString[] = "kinectPtr = PsychKinect('Open' [, deviceIndex=0][, numbuffers=2][, bayerFilterMode=1]);";
    //
    static char synopsisString[] = 
		"Open connection to Kinect box, return a 'kinectPtr' handle to it.\n\n"
		"The call tries to open the Kinect box with index 'deviceIndex', or the "
		"first detected box, if 'deviceIndex' is omitted. It configures the box "
		"for video and depth capture, using a n-buffered pool for internal storage "
		"of up to 'numbuffers' frames. 'numbuffers' defaults to 2 for double-buffering "
		"if omitted.\n"
		"'bayerFilterMode' is optional and selects mode of color reconstruction from "
		"the visual cameras sensor data: A setting of 1 (default) will use a slow and "
		"unoptimized bayer filter algorithm built into the libfreenect driver to return "
		"RGB images via the 'GetImage' function. A setting of 0 will disable filtering and return "
		"a 2D only image with sensor raw data, which you'll need to post-process and filter yourself.\n"
		"A setting of 2 will return infrared image data from the depth camera instead of RGB data.\n"
		"A setting of > 2 would select other filtering methods, but these are not yet implemented. "
		"You would want to use the default value if you need convenient access to RGB images without "
		"any effort on your side. Unfiltered raw data has the advantage that cpu load of kinect "
		"operation is much lower, so you may get higher framerates. Also, raw data only requires "
		"a third of the storage memory if this is of concern, e.g., for high settings of 'numbuffers'.\n"
		"The returned handle can be passed to the other subfunctions to operate the device.\n";
	
	static char seeAlsoString[] = "";	
	
	double R[3][3] = {{  9.9984628826577793e-01, 1.2635359098409581e-03, -1.7487233004436643e-02 },
			  { -1.4779096108364480e-03, 9.9992385683542895e-01, -1.2251380107679535e-02 },
			  {  1.7470421412464927e-02, 1.2275341476520762e-02,  9.9977202419716948e-01 }};

	double T[3] = { 1.9985242312092553e-02, -7.4423738761617583e-04, -1.0916736334336222e-02 };
	double depthIntrinsics[5] = { -2.6386489753128833e-01, 9.9966832163729757e-01, -7.6275862143610667e-04, 5.0350940090814270e-03, -1.3053628089976321e+00 };
	double rgbIntrinsics[5]   = {  2.6451622333009589e-01, -8.3990749424620825e-01, -1.9922302173693159e-03, 1.4371995932897616e-03, 9.1192465078713847e-01 };
	double depthBaseAndOffset[2] = { 0.0 , 0.0 };

	PsychKNDevice* kinect;
	int i, j;
	int deviceIndex = 0;
	int handle = 0;
	freenect_device *dev = NULL;
	int numbuffers = 2;
	int bayerFilterMode = 1;

	#if PSYCH_SYSTEM != PSYCH_WINDOWS
	freenect_frame_mode fmode;
	#endif
	
	// All sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

	//check to see if the user supplied superfluous arguments
	PsychErrorExit(PsychCapNumOutputArgs(1));
	PsychErrorExit(PsychCapNumInputArgs(3));
	
	if (devicecount >= MAX_PSYCH_KINECT_DEVS) PsychErrorExitMsg(PsychError_user, "Maximum number of simultaneously open kinect devices reached.");
	
	// Find a free device slot:
	for (handle = 0; (handle < MAX_PSYCH_KINECT_DEVS) && kinectdevices[handle].dev; handle++);
	if ((handle >= MAX_PSYCH_KINECT_DEVS) || kinectdevices[handle].dev) PsychErrorExitMsg(PsychError_internal, "Maximum number of simultaneously open kinect devices reached.");
	
	// Get optional kinect device index:
	PsychCopyInIntegerArg(1, FALSE, &deviceIndex);
	
	// Don't support anything than a single "default" kinect yet:
	if (deviceIndex < 0) PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' provided. Must be greater or equal to zero!");
	
	// Get optional numbuffers count:
	PsychCopyInIntegerArg(2, FALSE, &numbuffers);
	if (numbuffers < 1) PsychErrorExitMsg(PsychError_user, "Invalid value for 'numbuffers' provided. Must be greater than 1!");
	
	// Get optional bayerFilterMode:
	PsychCopyInIntegerArg(3, FALSE, &bayerFilterMode);
	if (bayerFilterMode < 0) PsychErrorExitMsg(PsychError_user, "Invalid value for 'bayerFilterMode' provided. Must be >= 0!");
	
	// Open and initialize the device:
	if (!initialized) {
		// First time init:
		initialized = TRUE;
		
		// Initialize libusb:
		if (freenect_init(&f_ctx, NULL) < 0) PsychErrorExitMsg(PsychError_system, "Driver initialization of libfreenect failed!");

		// MK FIXME TODO: Adapt to bigger limits for future Kinect devices?
		zmap = malloc(640 * 480 * 6 * sizeof(double));
	}
	
	// Zero init device structure:
	memset(&kinectdevices[handle], 0, sizeof(PsychKNDevice));
	
	kinect = &kinectdevices[handle];

	// Connect with it, get usb connection handle:
	if (freenect_open_device(f_ctx, &dev, deviceIndex) < 0) {
		printf("PsychKinect: ERROR! Failed to connect to kinect with deviceIndex %i. This could mean that the device\n", deviceIndex);
		printf("PsychKinect: ERROR: is already in use by another application or driver. On Linux it could mean it is\n");
		printf("PsychKinect: ERROR: claimed by the Kinect video camera driver. See 'help InstallKinect' for how to resolve this.\n");
		PsychErrorExitMsg(PsychError_user, "Could not connect to kinect device with given 'deviceIndex'! [libusb_open_device failed]");
	}
	
	kinectdevices[handle].dev = dev;
	freenect_set_user(kinectdevices[handle].dev, (void*) &kinectdevices[handle]);

	// Attach callbacks:
	freenect_set_depth_callback(kinect->dev, PsychDepthCB);
	freenect_set_rgb_callback(kinect->dev, PsychRGBCB);

	// Retrieve frame properties of current mode:
	#if PSYCH_SYSTEM != PSYCH_WINDOWS
		// Set video format and resolution:
		if (bayerFilterMode < 2) {
			// RGB video feed:
			fmode = freenect_find_video_mode(FREENECT_RESOLUTION_MEDIUM, ((bayerFilterMode == 1) ? FREENECT_VIDEO_RGB : FREENECT_VIDEO_BAYER));
		} else {
			// IR depth cam video feed:
			fmode = freenect_find_video_mode(FREENECT_RESOLUTION_MEDIUM, FREENECT_VIDEO_IR_8BIT);
		}
		freenect_set_video_mode(kinect->dev, fmode);

		// Set depth format and resolution:
		fmode = freenect_find_depth_mode(FREENECT_RESOLUTION_MEDIUM, FREENECT_DEPTH_11BIT);
		freenect_set_depth_mode(kinect->dev, fmode);

		// Query mode properties for video & depth:
		fmode = freenect_get_current_video_mode(kinect->dev);
		kinect->cwidth = fmode.width;
		kinect->cheight = fmode.height;
		kinect->csize = fmode.bytes;
		fmode = freenect_get_current_depth_mode(kinect->dev);
		kinect->dwidth = fmode.width;
		kinect->dheight = fmode.height;
		kinect->dsize = fmode.bytes;
	#else
		if (bayerFilterMode < 2) {
			// RGB video feed:
			freenect_set_rgb_format(kinect->dev, ((bayerFilterMode == 1) ? FREENECT_FORMAT_RGB : FREENECT_FORMAT_BAYER));
		} else {
			// IR depth cam video feed:
			freenect_set_rgb_format(kinect->dev, FREENECT_VIDEO_IR_8BIT);
		}
		freenect_set_depth_format(kinect->dev, FREENECT_FORMAT_11_BIT);

		kinect->cwidth = FREENECT_FRAME_W;
		kinect->cheight = FREENECT_FRAME_H;
		kinect->csize = (bayerFilterMode == 1) ? FREENECT_RGB_SIZE : FREENECT_BAYER_SIZE;
		kinect->dwidth = FREENECT_FRAME_W;
		kinect->dheight = FREENECT_FRAME_H;
		kinect->dsize = FREENECT_FRAME_W * FREENECT_FRAME_H * 2;
	#endif
		
	// Allocate buffers:
	kinectdevices[handle].buffers = (PsychKNBuffer*) calloc(numbuffers, sizeof(PsychKNBuffer));
	if (NULL == kinectdevices[handle].buffers) {
		printf("PTB-ERROR: Could not create requested %i kinect image buffers. Prepare for trouble!\n", numbuffers);
		PsychErrorExitMsg(PsychError_outofMemory, "Buffer creation failed!");
	}
	
	for (i = 0; i < numbuffers; i++) {
		kinectdevices[handle].buffers[i].color = (unsigned char*) calloc(1, kinect->csize);
		if (NULL == kinectdevices[handle].buffers[i].color) {
			for (j = 0; j < i; j++) free(kinectdevices[handle].buffers[j].color);
			free(kinectdevices[handle].buffers);
			kinectdevices[handle].buffers = NULL;
			printf("PTB-ERROR: Could not create requested %i kinect image buffers - Out of memory at %i'th buffer. Prepare for trouble!\n", numbuffers, i);
			PsychErrorExitMsg(PsychError_outofMemory, "Buffer creation failed!");
		}
	}
	
	kinectdevices[handle].numbuffers = numbuffers;
	kinectdevices[handle].bayerFilterMode = bayerFilterMode;
	
	// Have connection open and ready. Initialize mutexes, condition variables and processing thread:
	if (PsychInitMutex(&(kinectdevices[handle].mutex))) {
		printf("PsychKinect: CRITICAL! Failed to initialize Mutex object for handle %i! Prepare for trouble!\n", handle);
		PsychErrorExitMsg(PsychError_system, "Device mutex creation failed!");
	}
	
	// If we use locking, this will create & init the associated event variable:
	PsychInitCondition(&(kinectdevices[handle].changeSignal), NULL);

	// Preinit kinect camera parameters with ok values:
	// These are a bit wrong for any Kinect except the one
	// they were taken from, but at least they produce an
	// ok initial calibration for testing:
	kinectdevices[handle].fx_d = 5.9421434211923247e+02;
	kinectdevices[handle].fy_d = 5.9104053696870778e+02;
	kinectdevices[handle].cx_d = 3.3930780975300314e+02;
	kinectdevices[handle].cy_d = 2.4273913761751615e+02;

	kinectdevices[handle].fx_rgb = 5.2921508098293293e+02;
	kinectdevices[handle].fy_rgb = 5.2556393630057437e+02;
	kinectdevices[handle].cx_rgb = 3.2894272028759258e+02;
	kinectdevices[handle].cy_rgb = 2.6748068171871557e+02;

	memcpy(&kinectdevices[handle].R, R, sizeof(double) * 3 * 3);
	memcpy(&kinectdevices[handle].T, T, sizeof(double) * 3 * 1);

	memcpy(&kinectdevices[handle].undistort_d, depthIntrinsics, sizeof(double) * 5 * 1);
	memcpy(&kinectdevices[handle].undistort_rgb, rgbIntrinsics, sizeof(double) * 5 * 1);
	memcpy(&kinectdevices[handle].depthBaseAndOffset, depthBaseAndOffset, sizeof(double) * 2 *1);

	// Increment count of open devices:
	devicecount++;
	
	// Return device handle:
	PsychCopyOutDoubleArg(1, FALSE, handle);

	return(PsychError_none);	
}

PsychError PSYCHKINECTClose(void)
{
    static char useString[] = "PsychKinect('Close', kinectPtr);";
    //
    static char synopsisString[] = 
		"Close connection to Kinect box 'kinectPtr'. Do nothing if no such box is open.\n\n";

    static char seeAlsoString[] = "";	

	int handle;

    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(0));
    PsychErrorExit(PsychCapNumInputArgs(1));
    
	// Get device handle:
	PsychCopyInIntegerArg(1, TRUE, &handle);
	
	// Close device:
	PsychKNClose(handle);
	
    return(PsychError_none);	
}

PsychError PSYCHKINECTStart(void)
{
    static char useString[] = "[starttime, fps, cwidth, cheight, dwidth, dheight] = PsychKinect('Start', kinectPtr [, dropframes=0]);";
    //							1		   2	3		4		 5		 6								 1			  2
    static char synopsisString[] = 
		"Start video and depth capture operation of box 'kinectPtr'.\n\n"
		"Starts the capture operation into the internal bufferqueue. Color- and "
		"depths frames can be polled/waited for/dequeued from the queue via a "
		"call to PsychKinect('GrabFrame',...); and then retrieved via the 'GetXXX' "
		"functions. Data is stored in an internal ringbuffer, whose size is "
		"set in the 'Open' function. The optional 'dropframes' flag decides if "
		"the oldest frames should be dropped if the buffer is full (=1) or if "
		"the most recent frame should be dropped if the buffer is full (=0). "
		"\n\n"
		"The function returns the start time of capture and the nominal framerate "
		"'fps', color image width 'cwidth' and height 'cheight', and depth image "
		"'dwidth' and 'dheight'.\n\n";
	
    static char seeAlsoString[] = "";	
	
	int handle, dropframes;
	int rc;
	
	PsychKNDevice *kinect;
	
    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(6));
    PsychErrorExit(PsychCapNumInputArgs(2));
    
	// Get device handle:
	PsychCopyInIntegerArg(1, TRUE, &handle);
	kinect = PsychGetKinect(handle, FALSE);
	
	dropframes = 0;
	PsychCopyInIntegerArg(2, FALSE, &dropframes);
	
	PsychLockMutex(&kinect->mutex);
	if (kinect->state > 0) {
		PsychUnlockMutex(&kinect->mutex);
		PsychErrorExitMsg(PsychError_user, "Tried to start an already started capture operation!");
	}
	
	// Reset device state:
	kinect->recposition = 0;
	kinect->readposition = 0;
	kinect->reqstate = 1;
	kinect->state = 1;
	kinect->xruns = 0;
	kinect->paCalls = 0;
	kinect->frame_valid = 0;
	kinect->dropframes = (dropframes > 0) ? 1 : 0;
	
	// Create and startup thread:
	if ((rc=PsychCreateThread(&(kinectdevices[handle].captureThread), NULL, (void*) PsychKinectThreadMain, (void*) &(kinectdevices[handle])))) {
		printf("PTB-ERROR: Could not create kinect processing thread [%s]. Prepare for trouble!\n", strerror(rc));
		PsychErrorExitMsg(PsychError_system, "Thread creation failed!");
	}
	
	PsychUnlockMutex(&kinect->mutex);
	
	// TODO: Compute sensible starttime:
	PsychGetAdjustedPrecisionTimerSeconds(&(kinect->captureStartTime));
	
	PsychCopyOutDoubleArg(1, FALSE, kinect->captureStartTime);
	PsychCopyOutDoubleArg(2, FALSE, 30);
	PsychCopyOutDoubleArg(3, FALSE, kinect->cwidth);
	PsychCopyOutDoubleArg(4, FALSE, kinect->cheight);
	PsychCopyOutDoubleArg(5, FALSE, kinect->dwidth);
	PsychCopyOutDoubleArg(6, FALSE, kinect->dheight);
	
    return(PsychError_none);	
}

PsychError PSYCHKINECTStop(void)
{
    static char useString[] = "PsychKinect('Stop', kinectPtr);";
    //
    static char synopsisString[] = 
		"Stop capture on Kinect box 'kinectPtr'.\n\n";
	
    static char seeAlsoString[] = "";	
	
	int handle;
	
    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(0));
    PsychErrorExit(PsychCapNumInputArgs(1));
    
	PsychCopyInIntegerArg(1, TRUE, &handle);
	
	// Stop device:
	PsychKNStop(handle);
	
    return(PsychError_none);	
}


PsychError PSYCHKINECTSetAngle(void)
{
	static char useString[] = "PsychKinect('SetAngle', kinectPtr, angle);";
	//
	static char synopsisString[] = 
					"Set tilt angle of box 'kinectPtr' to 'angle' degrees.\n\n"
					"Allowable range for angle is -30 to +30 degrees.\n";

	static char seeAlsoString[] = "";	

	int handle, angle;
	PsychKNDevice *kinect;

	// All sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

	//check to see if the user supplied superfluous arguments
	PsychErrorExit(PsychCapNumOutputArgs(0));
	PsychErrorExit(PsychCapNumInputArgs(2));

	PsychCopyInIntegerArg(1, TRUE, &handle);	
	kinect = PsychGetKinect(handle, FALSE);

	// Get angle and limit it to valid range:
	PsychCopyInIntegerArg(2, TRUE, &angle);
	if (angle < -30) angle = -30;
	if (angle > +30) angle = +30;

	// Set actual tilt angle on device:
	freenect_set_tilt_degs(kinect->dev, angle);	

	return(PsychError_none);
}

PsychError PSYCHKINECTGetStatus(void)
{
    static char useString[] = "status = PsychKinect('GetStatus', kinectPtr);";
    //
    static char synopsisString[] = 
		"Retrieve struct with current status of Kinect box 'kinectPtr'.\n\n"
		"NOT FULLY IMPLEMENTED. JUST RETURNS DEBUG OUTPUT TO CONSOLE!\n";
	
    static char seeAlsoString[] = "";	
	
	int handle;
	PsychKNDevice *kinect;
	
	int navail;
	
    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(1));
    
	PsychCopyInIntegerArg(1, TRUE, &handle);	
	kinect = PsychGetKinect(handle, FALSE);
	
	PsychLockMutex(&kinect->mutex);
	
	navail = kinect->recposition - kinect->readposition;	
	
	printf("\n");
	printf("Running: %i\n", kinect->state);
	printf("Read: %i\n", kinect->readposition);
	printf("Write: %i\n", kinect->recposition);
	printf("navail: %i\n", navail);
	printf("\n");
	
	PsychUnlockMutex(&kinect->mutex);
	
	// TODO: Return sensible status struct:
	PsychCopyOutDoubleArg(1, FALSE, 123456);
	
    return(PsychError_none);	
}

PsychError PSYCHKINECTGrabFrame(void)
{
    static char useString[] = "[result, cts, age] = PsychKinect('GrabFrame', kinectPtr [, waitMode=1][, mostrecent=0]);";
    //							1		2	 3								 1			  2				3
    static char synopsisString[] = 
		"Poll or wait for arrival of new captured frames, release old ones.\n\n"
		"This function checks box 'kinectPtr' for the availability of frames. "
		"If 'waitMode' is 1, the function waits/blocks until data is available "
		"if no new frames are available at time of call. Otherwise it just "
		"polls in a non-blocking fashion.\n"
		"If 'mostrecent' is 1, then the most recent frame is returned and old "
		"frames are discarded - this for low-latency live capture.\n"
		"\n\n"
		"'result' returns the status: 0 = No frame returned, 1 = Frame returned, 2 = Error.\n"
		"'cts' returns a capture time stamp of when the frame was captured.\n"
		"'age' returns the age of the frame -- the number of pending frames in the buffer queue.\n"
		"\n\n";
	
    static char seeAlsoString[] = "";	
	
	int handle;
	PsychKNDevice *kinect;
	PsychKNBuffer *buffer;
	int waitMode, mostrecent;
	int navail;
	double cts;
	
    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(3));
    PsychErrorExit(PsychCapNumInputArgs(3));
    
	PsychCopyInIntegerArg(1, TRUE, &handle);	
	kinect = PsychGetKinect(handle, FALSE);
	
	waitMode = 1;
	PsychCopyInIntegerArg(2, FALSE, &waitMode);	
	
	mostrecent = 0;
	PsychCopyInIntegerArg(3, FALSE, &mostrecent);	
	
	if (kinect->frame_valid) PsychErrorExitMsg(PsychError_user, "Must 'ReleaseFrame' old frame first, before grabbing a new one!");
	if ((kinect->state == 0) && (waitMode!=0)) PsychErrorExitMsg(PsychError_user, "Must 'Start' before first call to 'GrabFrame'!");
	PsychLockMutex(&kinect->mutex);
	
	// New data available?
	navail = kinect->recposition - kinect->readposition;	
	
	if (verbosity > 3) {
		printf("\n");
		printf("Running: %i\n", kinect->state);
		printf("Read: %i\n", kinect->readposition);
		printf("Write: %i\n", kinect->recposition);
		printf("navail: %i\n", navail);
		printf("\n");
	}
	
	while (navail <= 0) {
		// No. Shall we wait?
		if (waitMode == 0) {
			// No, was just polling: Return no data:
			PsychUnlockMutex(&kinect->mutex);
			PsychCopyOutDoubleArg(1, FALSE, 0);
			PsychCopyOutDoubleArg(2, FALSE, -1);
			PsychCopyOutDoubleArg(3, FALSE, 0);
			return(PsychError_none);
		}
		
		// Yes. Wait for change notification from thread:
		PsychTimedWaitCondition(&kinect->changeSignal, &kinect->mutex, 0.5);
		if (kinect->state <= 0) {
			// Something wen't wrong, probably a enect driver error condition.
			// Clear the waitMode flag and continue, so we break out as if we
			// failed to poll for new data:
			navail = 0;
			waitMode = 0;
			if (verbosity > 1) printf("PTB-WARNING: In 'GrabFrame': Kinect driver timed out or emergency stopped!\nPrepare for trouble!\n");
			continue;
		}
		
		// Recheck:
		navail = kinect->recposition - kinect->readposition;		
		
		if (verbosity > 3) {
			printf("\n");
			printf("Running: %i\n", kinect->state);
			printf("Read: %i\n", kinect->readposition);
			printf("Write: %i\n", kinect->recposition);
			printf("navail: %i\n", navail);
			printf("\n");
		}
	}
	
	// Ok, at this point we have at least one frame of data available.
	if ((navail > 1) && (mostrecent == 1)) {
		// Want most recent. Increment readpos to last avail frame:
		kinect->readposition+= navail - 1;
	}
	
	// Recompute available amount:
	navail = kinect->recposition - kinect->readposition;		
	
	// Retrieve bufferptr:
	buffer = PsychGetKNBuffer(kinect, kinect->readposition);
	cts = buffer->cts;
	
	PsychUnlockMutex(&kinect->mutex);
	
	PsychCopyOutDoubleArg(1, FALSE, 1);
	PsychCopyOutDoubleArg(2, FALSE, cts);
	PsychCopyOutDoubleArg(3, FALSE, navail);
	
	kinect->frame_valid = TRUE;
	
    return(PsychError_none);	
}

PsychError PSYCHKINECTReleaseFrame(void)
{
    static char useString[] = "PsychKinect('ReleaseFrame', kinectPtr);";
    static char synopsisString[] = 
		"Release current frame back to bufferqueue.\n\n";
    static char seeAlsoString[] = "";	
	
	int handle;
	PsychKNDevice *kinect;
	
	int navail;
	
    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(3));
    PsychErrorExit(PsychCapNumInputArgs(3));
    
	PsychCopyInIntegerArg(1, TRUE, &handle);	
	kinect = PsychGetKinect(handle, FALSE);
	if (!kinect->frame_valid) PsychErrorExitMsg(PsychError_user, "Must 'GrabFrame' a frame first!");
	
	PsychLockMutex(&kinect->mutex);
	kinect->readposition++;
	navail = kinect->recposition - kinect->readposition;	
	PsychUnlockMutex(&kinect->mutex);
	
	kinect->frame_valid = FALSE;
	
    return(PsychError_none);	
}

PsychError PSYCHKINECTGetImage(void)
{
    static char useString[] = "[imageOrPtr, width, height, channels, extType, extFormat] = PsychKinect('GetImage', kinectPtr [, imtype=0][, returnTexturePtr=0]);";
    //							1			2	   3	   4		 5		  6									   1			2			3
    static char synopsisString[] = 
		"Return the color image data for the frame fetched via 'GrabFrame'.\n\n"
		"If 'returnTexturePtr' is zero (default), a uint8 3D image matrix is returned, which could "
		"be directly fed to image processing toolbox, imwrite() or Screen('MakeTexture').\n"
		"'imtype' if set to 1 returns a color-coded depth image instead of the RGB camera image.\n"
		"If 'returnTexturePtr' is one, a double-encoded memory pointer to a RGBA8 rectangle texture buffer "
		"is returned, for use with Screen('SetOpenGLTextureFromMemPointer') for injection into "
		"PTB's texturing system.\n";
	
    static char seeAlsoString[] = "";	
	
	int handle;
	PsychKNDevice *kinect;
	PsychKNBuffer* buffer;
	unsigned char* outimg;
	int returnTexturePtr, imtype;
	int i, pval, lb;
	int channels = 3;
	
    // All sub functions should have these two lines
    PsychPushHelp(useString, synopsisString,seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
    //check to see if the user supplied superfluous arguments
    PsychErrorExit(PsychCapNumOutputArgs(6));
    PsychErrorExit(PsychCapNumInputArgs(3));
	
	PsychCopyInIntegerArg(1, TRUE, &handle);	
	kinect = PsychGetKinect(handle, FALSE);
	if (!kinect->frame_valid) PsychErrorExitMsg(PsychError_user, "Must 'GrabFrame' a frame first!");
	
	imtype=0;
	PsychCopyInIntegerArg(2, FALSE, &imtype);	
	
	returnTexturePtr=0;
	PsychCopyInIntegerArg(3, FALSE, &returnTexturePtr);	
	
	// Retrieve bufferptr:
	buffer = PsychGetKNBuffer(kinect, kinect->readposition);
	
	if (imtype == 1) {
		// Return color encoded depth image:
		
		// Need to convert into false-color image:		
		for (i=0; i<640*480; i++) {
			pval = t_gamma[buffer->depth[i]];
			lb = pval & 0xff;
			switch (pval>>8) {
				case 0:
					gl_depth_back[3*i+0] = 255;
					gl_depth_back[3*i+1] = 255-lb;
					gl_depth_back[3*i+2] = 255-lb;
					break;
				case 1:
					gl_depth_back[3*i+0] = 255;
					gl_depth_back[3*i+1] = lb;
					gl_depth_back[3*i+2] = 0;
					break;
				case 2:
					gl_depth_back[3*i+0] = 255-lb;
					gl_depth_back[3*i+1] = 255;
					gl_depth_back[3*i+2] = 0;
					break;
				case 3:
					gl_depth_back[3*i+0] = 0;
					gl_depth_back[3*i+1] = 255;
					gl_depth_back[3*i+2] = lb;
					break;
				case 4:
					gl_depth_back[3*i+0] = 0;
					gl_depth_back[3*i+1] = 255-lb;
					gl_depth_back[3*i+2] = 255;
					break;
				case 5:
					gl_depth_back[3*i+0] = 0;
					gl_depth_back[3*i+1] = 0;
					gl_depth_back[3*i+2] = 255-lb;
					break;
				default:
					gl_depth_back[3*i+0] = 0;
					gl_depth_back[3*i+1] = 0;
					gl_depth_back[3*i+2] = 0;
					break;
			}
		}
		
		// Return image data:
		if (returnTexturePtr) {
			// Just return a memory pointer to the colorbuffer:
			PsychCopyOutDoubleArg(1, FALSE, PsychPtrToDouble((void*) gl_depth_back));
		} else {
			PsychAllocOutUnsignedByteMatArg(1, FALSE, 3, buffer->dheight, buffer->dwidth, &outimg);
			memcpy(outimg, gl_depth_back, buffer->dheight * buffer->dwidth * 3);
		}
		
		// Fixed for current Kinect:
		PsychCopyOutDoubleArg(2, FALSE, buffer->dwidth);
		PsychCopyOutDoubleArg(3, FALSE, buffer->dheight);		
	}
	else {
		// Return color image from visual camera:
		
		// If these are not bayer filtered by libfreenect, we only have 1 byte per
		// pixel for a 2D raw sensor data matrix:
		channels = (kinect->bayerFilterMode == 1) ? 3 : 1;
		
		if (returnTexturePtr) {
			// Just return a memory pointer to the colorbuffer:
			PsychCopyOutDoubleArg(1, FALSE, PsychPtrToDouble((void*) buffer->color));
		} else {
			PsychAllocOutUnsignedByteMatArg(1, FALSE, channels, buffer->cheight, buffer->cwidth, &outimg);
			memcpy(outimg, buffer->color, buffer->cheight * buffer->cwidth * channels);
		}
		
		// Fixed for current Kinect:
		PsychCopyOutDoubleArg(2, FALSE, buffer->cwidth);
		PsychCopyOutDoubleArg(3, FALSE, buffer->cheight);		
	}
	
	// Number of image channels:
	PsychCopyOutDoubleArg(4, FALSE, channels);
	
	// GL_RGB8 or GL_LUMINANCE8
	PsychCopyOutDoubleArg(5, FALSE, (channels == 3) ? 32849 : 0x8040);
	
	// GL_UNSIGNED_BYTE:
	PsychCopyOutDoubleArg(6, FALSE, 5121);
	
    return(PsychError_none);	
}

// Fetch zinbuf raw disparity value at location (x,y):
unsigned short getz(unsigned short *zinbuf, int x, int y)
{
	return(zinbuf[y * 640 + x]);
}

// Convert raw disparity value into z distance meters):
// Using Magic formula from www.openkinect.org Wiki,
// "Imaging parameter" section:
static double calclinearz(int raw_depth)
{
	if (raw_depth < 2047) {
		// Valid measurement - Convert:
		return(1.0 / ((double) raw_depth * -0.0030711016 + 3.3309495161));
	}

	// Invalid datapoint - Mark as such:
	return 0;
}

static double calcbaselinez(int raw_depth)
{
	if (raw_depth < 2047) {
		// Valid measurement - Convert:
		return(540.0 * 8.0 * depthBaseAndOffset[0] / (depthBaseAndOffset[1] - ((double) raw_depth)));
	}

	// Invalid datapoint - Mark as such:
	return 0;
}

static double calctanhz(int raw_depth)
{
	const double k1 = 1.1863;
	const double k2 = 2842.5;
	const double k3 = 0.1236;

	if (raw_depth < 2047) {
		// Valid measurement - Convert:
		return(k3 * tan(((double) raw_depth) / k2 + k1));
	}

	// Invalid datapoint - Mark as such:
	return 0;
}


PsychError PSYCHKINECTGetDepthImage(void)
{
	static char useString[] = "[imageOrPtr, width, height, components, extFormat] = PsychKinect('GetDepthImage', kinectPtr [, format=0][, returnTexturePtr=0]);";
	static char synopsisString[] = 
		"Return the depth image data for the frame fetched via 'GrabFrame'.\n\n"
		"If 'returnTexturePtr' is zero (default), a matrix is returned for processing in Matlab/Octave.\n\n"
		"If 'returnTexturePtr' is one, a memory pointer to a buffer is returned.\n\n"
		"'format' defines the type of returned data:\n"
		"0 = Return raw disparity image as 2D double matrix with integral values.\n"
		"1 = Return depths z-image as 2D double matrix with z-distance in meters.\n"
		"2 = Return a vertex buffer with (x,y,z) vertices that define a 3D surface mesh\n"
		"    and (R,G,B) color values for each vertex. -> [x,y,z,r,g,b] per element.\n"
		"3 = Return a vertex buffer with (x,y,z) vertices that define a 3D surface mesh.\n"
		"    and (tx,ty) texture coordinates for vertices -> [x,y,z,tx,ty] per element.\n"
		"4/5 = Return a vertex buffer with (xi,yi,z) vertices that define sensor pixel\n"
		"      position (xi,yi) of depths sensor and reconstructed z value. This needs\n"
		"      to be post-processed in a vertex shader for speedups.\n"
		"6/7 = Return a vertex buffer with vertex id's uniquely identifying each sensor \n"
		"      position of depths sensor and raw sensor value at that location. The whole\n"
		"      3D reconstructin is done on the GPU in a vertex shader for maximum speed.\n"
		"8   = Return a uint16 buffer with a transposed copy of the raw depth buffer.\n"
		"      This is the most compact and efficient way to return raw data to you. The\n"
		"      transposed format is again for efficiency reasons. You need to transpose()\n"
		"      the returned 2D data matrix yourself.\n"
		"      Alternatively return a memory pointer to the 16 bit unsigned integer raw\n"
		"      depth buffer. CAUTION: The pointer becomes invalid as soon as the\n"
		"      current buffer is released via 'ReleaseFrame'! This is the fast-path.\n"
		"      "
		"\n\n";

	static char seeAlsoString[] = "";	
	
	int handle;
	PsychKNDevice *kinect;
	PsychKNBuffer* buffer;
	double* outzmat;
	float* fmap;
	short* inbufs;

	int returnTexturePtr, format, texcoords;
	int i, x, y, m, n, components;
	double P[3], Pr[3], Pt[2];
	double s, zp;
	double fx_d;
	double fy_d;
	double cx_d;
	double cy_d;
	double fx_rgb;
	double fy_rgb;
	double cx_rgb;
	double cy_rgb;
	double R[3][3];
	double T[3];

	const double w = 640;
	const double h = 480;

	// All sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
 
	//check to see if the user supplied superfluous arguments
	PsychErrorExit(PsychCapNumOutputArgs(5));
	PsychErrorExit(PsychCapNumInputArgs(3));

	PsychCopyInIntegerArg(1, TRUE, &handle);	
	kinect = PsychGetKinect(handle, FALSE);
	if (!kinect->frame_valid) PsychErrorExitMsg(PsychError_user, "Must 'GrabFrame' a frame first!");

	// Copy cam parameters to local storage to reduce number
	// of address indirections in perf-critical path:
	fx_d = kinect->fx_d;
	fy_d = kinect->fy_d;
	cx_d = kinect->cx_d;
	cy_d = kinect->cy_d;

	fx_rgb = kinect->fx_rgb;
	fy_rgb = kinect->fy_rgb;
	cx_rgb = kinect->cx_rgb;
	cy_rgb = kinect->cy_rgb;
	
	memcpy(&R, &kinect->R, sizeof(double) * 3 * 3);
	memcpy(&T, &kinect->T, sizeof(double) * 3 * 1);
	memcpy(&depthBaseAndOffset, &kinect->depthBaseAndOffset, sizeof(double) * 2 * 1);

	// Calibration provided?
	if ((depthBaseAndOffset[0] != 0) && (depthBaseAndOffset[1] != 0)) {
		// Yes: Use it with baseline reconstruction method:
		calcz = &calcbaselinez;
	} else {
		// No: Special tanh method requested?
		if (depthBaseAndOffset[1] == 1) {
			// Yes: Use tanh method:
			calcz = &calctanhz;
		} else {
			// No: Use hard-coded linear reconstruction method:
			calcz = &calclinearz;
		}
	}

	format = 0;
	texcoords = 0;
	PsychCopyInIntegerArg(2, FALSE, &format);	
	
	returnTexturePtr=0;
	PsychCopyInIntegerArg(3, FALSE, &returnTexturePtr);	
	
	// Retrieve bufferptr:
	buffer = PsychGetKNBuffer(kinect, kinect->readposition);
	
	switch(format) {
		case 0:
			// Return raw depth image:
			i=0;
			for (x=0; x < 640; x++)
				for (y=0; y < 480; y++) {
					zmap[i++] = (double) getz(buffer->depth, x, y);
				}
					
			// Return image data:
			if (returnTexturePtr) {
				// Just return a memory pointer to the depthbuffer:
				PsychCopyOutPointerArg(1, FALSE, (void*) zmap);
			} else {
				PsychAllocOutDoubleMatArg(1, FALSE, buffer->dheight, buffer->dwidth, 1, &outzmat);
				memcpy(outzmat, zmap, buffer->dheight * buffer->dwidth * 1 * sizeof(double));
			}

			components = 1;
		break;
				
		case 1:
			// Return encoded depth image:
			i=0;
			for (x=0; x < 640; x++)
				for (y=0; y < 480; y++) {
					// Calc z distance in meters:
					zmap[i++] = calcz(getz(buffer->depth, x, y));
				}
					
			// Return image data:
			if (returnTexturePtr) {
				// Just return a memory pointer to the depthbuffer:
				PsychCopyOutPointerArg(1, FALSE, (void*) zmap);
			} else {
				PsychAllocOutDoubleMatArg(1, FALSE, buffer->dheight, buffer->dwidth, 1, &outzmat);
				memcpy(outzmat, zmap, buffer->dheight * buffer->dwidth * 1 * sizeof(double));
			}

			components = 1;
		break;

		case 2:
		case 3:
			// Return (x,y,z) vertex buffer mesh:
			// Mapping equations from
			// http://nicolas.burrus.name/index.php/Research/KinectCalibration
			// Format 3 = Generate texture coords instead of vertex colors:
			if (format == 3) texcoords = 1;

			i=0;
			for (x=0; x < 640; x++) {
				for (y=0; y < 480; y++) {
					// Compute cartesian 3D (x,y,z) vertex coordinates,
					// one vertex per 3D point aka depth camera pixel:
					zp = calcz(getz(buffer->depth, x, y));
					P[0] = ((double) x - cx_d) * zp / fx_d;
					P[1] = ((double) y - cy_d) * zp / fy_d;
					P[2] = zp;

					// Write out vertex 3D (x,y,z):
					zmap[i++] = P[0];
					zmap[i++] = P[1];
					zmap[i++] = P[2];
	
					// Compute 2D texture coordinates for each 3D
					// point, aka (x,y) image coordinates in color image:

					// P3D' = R.P3D + T  --> Project from depth cams
					// reference frame to color cams reference frame,
					// aka apply extrinsic rotation/translation parameters
					// of color cam wrt. depths cam:
					for (m = 0; m < 3; m++) {
						s = 0.0;
						for (n = 0; n < 3; n++) s += R[m][n] * P[n];
						Pr[m] = s + T[m];
					}

					// Apply intrinsic parameters of color cam to project
					// into color cams 2D sensor plane, aka color image
					// coordinates, aka texture coordinates:
					Pt[0] = (Pr[0] * fx_rgb / Pr[2]) + cx_rgb;
					Pt[1] = (Pr[1] * fy_rgb / Pr[2]) + cy_rgb;

					// Return 3-component RGB8 vertex colors, or 2D texture
					// coordinates?
					if (!texcoords) {
						// Perform "manual" texture lookup and return RGB
						// pixel:

						// Clamp:
						if (Pt[0] < 0) Pt[0] = 0;
						if (Pt[1] < 0) Pt[1] = 0;
						if (Pt[0] >= w) Pt[0] = w-1;
						if (Pt[1] >= h) Pt[1] = h-1;

						// "Nearest neighbour texture lookup":
						zmap[i++] = ((double) buffer->color[(int) (3 * ((((int) Pt[1]) * w + ((int) Pt[0]))) + 0)]) / 255.0; // R
						zmap[i++] = ((double) buffer->color[(int) (3 * ((((int) Pt[1]) * w + ((int) Pt[0]))) + 1)]) / 255.0; // G
						zmap[i++] = ((double) buffer->color[(int) (3 * ((((int) Pt[1]) * w + ((int) Pt[0]))) + 2)]) / 255.0; // B
					} else {
						// Output direct texture coordinates into a
						// GL_TEXTURE_RECTANGLE texture:
						zmap[i++] = Pt[0];
						zmap[i++] = Pt[1];
					}
				}
			}

			// Return image data:
			components = 3 + ((texcoords) ? 2 : 3);

			if (returnTexturePtr) {
				// Just return a memory pointer to the depthbuffer:
				PsychCopyOutPointerArg(1, FALSE, (void*) zmap);
			} else {
				PsychAllocOutDoubleMatArg(1, FALSE, components , buffer->dheight, buffer->dwidth, &outzmat);
				memcpy(outzmat, zmap, components * buffer->dheight * buffer->dwidth * sizeof(double));
			}

		break;

		case 4:
		case 5:
			// Return encoded depth image:
			if (format == 4) {
				i = 0;
				for (y=0; y < 480; y++) {
					for (x=0; x < 640; x++) {
						// Calc z distance in meters:
						zmap[i++] = (double) x;
						zmap[i++] = (double) y;
						zmap[i++] = calcz(getz(buffer->depth, x, y));
					}
				}
			} else {
				i = 2;
				for (y=0; y < 480; y++) {
					for (x=0; x < 640; x++) {
						// Calc z distance in meters:
						zmap[i] = calcz(getz(buffer->depth, x, y));
						i+=3;
					}
				}
			}

			// Return image data:
			if (returnTexturePtr) {
				// Just return a memory pointer to the depthbuffer:
				PsychCopyOutPointerArg(1, FALSE, (void*) zmap);
			} else {
				PsychAllocOutDoubleMatArg(1, FALSE, 3, buffer->dheight, buffer->dwidth, &outzmat);
				memcpy(outzmat, zmap, 3 * buffer->dheight * buffer->dwidth * sizeof(double));
			}

			components = 3;
		break;
				
		case 6:
		case 7:
			// Return encoded depth image:
			fmap = (float*) &zmap[0];
			inbufs = (short*) buffer->depth;

			if (format == 6) {
				i = 0;
				while (i < 640 * 480 * 2) {
					fmap[i] = (float) (i / 2);
					i++;
					fmap[i++] = (float) *(inbufs++);
				}
			} else {
				i = 0;
				while (i < 640 * 480 * 2) {
					i++;
					fmap[i++] = (float) *(inbufs++);
				}
			}

			// Return image data:
			if (returnTexturePtr) {
				// Just return a memory pointer to the depthbuffer:
				PsychCopyOutPointerArg(1, FALSE, (void*) zmap);
			} else {
				PsychAllocOutDoubleMatArg(1, FALSE, 1, buffer->dheight, buffer->dwidth, &outzmat);
				memcpy(outzmat, zmap, 1 * buffer->dheight * buffer->dwidth * sizeof(double));
			}

			components = 1;
		break;

		case 8:
			// Return 16 bit unsigned integer raw depths data:
			if (returnTexturePtr) {
				// Just return a memory pointer to the depthbuffer:
				PsychCopyOutPointerArg(1, FALSE, (void*) (buffer->depth));
			} else {
				// Return 16 bit unsigned integer raw depths data in untransposed raw format,
				// (or transposed, from the perspective of the runtime):
				PsychCopyOutUnsignedInt16MatArg(1, FALSE, buffer->dwidth, buffer->dheight, 1, (psych_uint16*) (buffer->depth));
			}

			components = 1;
		break;


		default:
			components = 0;
			PsychErrorExitMsg(PsychError_user, "Invalid 'format' parameter provided!");
	}
	
	// Fixed for current Kinect:
	PsychCopyOutDoubleArg(2, FALSE, buffer->dwidth);
	PsychCopyOutDoubleArg(3, FALSE, buffer->dheight);	
	PsychCopyOutDoubleArg(4, FALSE, components);
	
	// Data format is GL_DOUBLE = 5130 (float would be GL_FLOAT = 5126, GL_INT = 5124):
	// Note: GL_INT was also tried, as it is theoretically more efficient to handle on cpu,
	// but was way slower, as apparently many GPU's have trouble handling it efficiently. Likely
	// we hit some driver fallback or internal data conversion:
	if (format < 8) {
		PsychCopyOutDoubleArg(5, FALSE, (components == 1 && (format == 6 || format == 7)) ? 5126 : 5130);
	} else {
		// GL_UNSIGNED_SHORT format for texture buffer, useable in vertex/geometry shader for
		// vertex-texture-fetch of raw depths images:
		PsychCopyOutDoubleArg(5, FALSE, 5123);
	}

	return(PsychError_none);	
}

PsychError PSYCHKINECTSetBaseCalibration(void)
{
    static char useString[] = "[...old parameter settings in order of inputs... ] = PsychKinect('SetBaseCalibration', kinectPtr, depthsIntrinsics, rgbIntrinsics, rgbRotation, rgbTranslation, depthsUndistort, rgbUndistort, depthBaseAndOffset);";
    //                                                                                                                1          2                 3              4            5               6                7             8
    static char synopsisString[] = 
		"Set and/or Get camera calibration matrices of box 'kinectPtr'.\n\n"
		"This function assigns various parameters needed for 3D reconstruction "
		"of the 3D scene from Kinect raw sensor data.\n"
		"You will need to use an external 3rd party camera calibration tool "
		"to determine the specific parameters of your Kinect. The driver uses "
		"some builtin defaults of a reference Kinect box if you don't provide "
		"updated parameters via this function. This allows to get at least "
		"ok results for quick testing and experimenting.\n"
		"The optical undistortion parameters (arguments 6 and 7) may not get "
		"used by early versions of the driver, or in fast 3D preview modes. "
		"They are accepted but silently ignored on such configurations.\n"
		"'depthsIntrinsics' vector with intrinsic parameters of depth cam:\n "
		"[fx, fy, cx, cy].\n"
		"'rgbIntrinsics' vector with intrinsic parameters of the RGB video camera:\n"
		"[fx, fy, cx, cy].\n"
		"'rgbRotation' a 3x3 rotation matrix which encodes rotation of the "
		"color camera with respect to the depths camera.\n"
		"'rgbTranslation' a 3 component translation vector which encodes "
		"translation of the color camera with respect to the depths camera.\n"
		"'depthsUndistort' undistoration parameters for depths camera.\n"
		"[k1,k2,p1,p2,k3].\n"
		"'rgbUndistort' undistoration parameters for RGB video camera.\n"
		"[k1,k2,p1,p2,k3].\n"
		"'depthBaseAndOffset' a 2 element vector [kinectDepthBase, kinectDephOffset].\n";

    static char seeAlsoString[] = "";	
	
	int handle;	
	PsychKNDevice *kinect;
	double *depthIntrinsics, *rgbIntrinsics;
	double *rgbRotation, *rgbTranslation;
	double *depthBaseAndOffset;
	int m, n, p;

	// All sub functions should have these two lines
	PsychPushHelp(useString, synopsisString,seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};
	
	//check to see if the user supplied superfluous arguments
	PsychErrorExit(PsychCapNumOutputArgs(7));
	PsychErrorExit(PsychCapNumInputArgs(8));

	// Get device handle:
	PsychCopyInIntegerArg(1, TRUE, &handle);
	kinect = PsychGetKinect(handle, FALSE);

	if (PsychAllocOutDoubleMatArg(1, FALSE, 1, 4, 1, &depthIntrinsics)) {
		depthIntrinsics[0] = kinect->fx_d;
		depthIntrinsics[1] = kinect->fy_d;
		depthIntrinsics[2] = kinect->cx_d;
		depthIntrinsics[3] = kinect->cy_d;
	}

	if (PsychAllocOutDoubleMatArg(2, FALSE, 1, 4, 1, &rgbIntrinsics)) {
		rgbIntrinsics[0] = kinect->fx_rgb;
		rgbIntrinsics[1] = kinect->fy_rgb;
		rgbIntrinsics[2] = kinect->cx_rgb;
		rgbIntrinsics[3] = kinect->cy_rgb;
	}

	if (PsychAllocOutDoubleMatArg(3, FALSE, 3, 3, 1, &rgbRotation)) {
		memcpy(rgbRotation, &kinect->R, sizeof(double) * 3 * 3);
	}

	if (PsychAllocOutDoubleMatArg(4, FALSE, 1, 3, 1, &rgbTranslation)) {
		memcpy(rgbTranslation, &kinect->T, sizeof(double) * 3 * 1);
	}

	if (PsychAllocOutDoubleMatArg(5, FALSE, 1, 5, 1, &depthIntrinsics)) {
		memcpy(depthIntrinsics, &kinect->undistort_d, sizeof(double) * 5 * 1);
	}

	if (PsychAllocOutDoubleMatArg(6, FALSE, 1, 5, 1, &rgbIntrinsics)) {
		memcpy(rgbIntrinsics, &kinect->undistort_rgb, sizeof(double) * 5 * 1);
	}

	if (PsychAllocOutDoubleMatArg(7, FALSE, 1, 2, 1, &depthBaseAndOffset)) {
		memcpy(depthBaseAndOffset, &kinect->depthBaseAndOffset, sizeof(double) * 2 * 1);
	}

	if (PsychAllocInDoubleMatArg(2, FALSE, &m, &n, &p, &depthIntrinsics)) {
		if (m * n * p != 4) PsychErrorExitMsg(PsychError_user, "Number of 'depthIntrinsic' elements not 4 as required!");
		kinect->fx_d = depthIntrinsics[0];
		kinect->fy_d = depthIntrinsics[1];
		kinect->cx_d = depthIntrinsics[2];
		kinect->cy_d = depthIntrinsics[3];
	}

	if (PsychAllocInDoubleMatArg(3, FALSE, &m, &n, &p, &rgbIntrinsics)) {
		if (m * n * p != 4) PsychErrorExitMsg(PsychError_user, "Number of 'rgbIntrinsic' elements not 4 as required!");
		kinect->fx_rgb = rgbIntrinsics[0];
		kinect->fy_rgb = rgbIntrinsics[1];
		kinect->cx_rgb = rgbIntrinsics[2];
		kinect->cy_rgb = rgbIntrinsics[3];
	}

	if (PsychAllocInDoubleMatArg(4, FALSE, &m, &n, &p, &rgbRotation)) {
		if ((m!=3) || (n!=3) || (p!= 1)) PsychErrorExitMsg(PsychError_user, "'rgbRotation' matrix isn't a 3-by-3 matrix as required!");
		memcpy(&kinect->R, rgbRotation, sizeof(double) * 3 * 3);
	}

	if (PsychAllocInDoubleMatArg(5, FALSE, &m, &n, &p, &rgbTranslation)) {
		if (m * n * p != 3) PsychErrorExitMsg(PsychError_user, "'rgbTranslation' vector isn't a 3-by-1 vector as required!");
		memcpy(&kinect->T, rgbTranslation, sizeof(double) * 3 * 1);
	}

	if (PsychAllocInDoubleMatArg(6, FALSE, &m, &n, &p, &depthIntrinsics)) {
		if (m * n * p != 5) PsychErrorExitMsg(PsychError_user, "Number of 'depthUndistort' elements not 5 as required!");
		memcpy(&kinect->undistort_d, depthIntrinsics, sizeof(double) * 5 * 1);
	}

	if (PsychAllocInDoubleMatArg(7, FALSE, &m, &n, &p, &rgbIntrinsics)) {
		if (m * n * p != 5) PsychErrorExitMsg(PsychError_user, "Number of 'rgbUndistort' elements not 5 as required!");
		memcpy(&kinect->undistort_rgb, rgbIntrinsics, sizeof(double) * 5 * 1);
	}

	if (PsychAllocInDoubleMatArg(8, FALSE, &m, &n, &p, &depthBaseAndOffset)) {
		if (m * n * p != 2) PsychErrorExitMsg(PsychError_user, "'depthBaseAndOffset' vector isn't a 2-by-1 vector as required!");
		memcpy(&kinect->depthBaseAndOffset, depthBaseAndOffset, sizeof(double) * 2 * 1);
	}

	return(PsychError_none);
}
