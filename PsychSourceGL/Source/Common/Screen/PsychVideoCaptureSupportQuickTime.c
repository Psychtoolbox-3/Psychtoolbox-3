/*
	Common/Screen/PsychVideoCaptureSupportQuickTime.c
	
	PLATFORMS:	
	
		MacOS/X and MS-Windows
		
		This is the videocapture- and videorecording engine based on the
		Apple QuickTime SequenceGrabber API. It works on Apple MacOS/X and MS-Windows.

	AUTHORS:
	
		Mario Kleiner           mk              mario.kleiner@tuebingen.mpg.de

	HISTORY:
	
	DESCRIPTION:
	
		Psychtoolbox functions for dealing with video capture devices. This implementation uses
		Quicktimes Sequence Grabber API, therefore it works on Operating systems with Quicktime support,
		currently OS/X and Windows.

	NOTES:

*/


#include "Screen.h"
#include <float.h>

// No Quicktime Sequence Grabber support for GNU/Linux:
#if (PSYCH_SYSTEM != PSYCH_LINUX) && defined(PTBVIDEOCAPTURE_QT)

#if PSYCH_SYSTEM == PSYCH_OSX
#include <Quicktime/QuickTimeComponents.h>
#endif

#if PSYCH_SYSTEM == PSYCH_WINDOWS
#include <QTML.h>
#include <QuickTimeComponents.h>
#endif

// Forward declaration of internal helper function:
void PsychQTDeleteAllCaptureDevices(void);
OSErr PsychQTSelectVideoSource(SeqGrabComponent seqGrab, SGChannel* sgchanptr, int deviceIndex, psych_bool showOutput, char* retDeviceName);
void* PsychQTVideoCaptureThreadMain(void* vidcapRecordToCast);


// Record which defines all state for a capture device:
typedef struct {
	psych_mutex			mutex;			// Access mutex to coordinate access between worker thread and masterthread, if any.
	unsigned int		recordingflags; // Flags specified at device open time.
    GWorldPtr           gworld;       // Offscreen GWorld into which captured frame is decompressed.
    SeqGrabComponent 	seqGrab;	     // Sequence grabber handle.
    SGChannel           sgchanVideo;  // Handle for video channel of sequence grabber.
	SGChannel			sgchanAudio;  // Handle for audio channel of sequence grabber.
    ImageSequence 	decomSeq;	     // unique identifier for our video decompression sequence
	ComponentInstance   vdig;		  // Actual VDIG component.
    int nrframes;                     // Total count of decompressed images.
    double fps;                       // Acquisition framerate of capture device.
    int width;                        // Width x height of captured images.
    int height;
	int reqpixeldepth;                // Requested depth of single pixel in output texture.
    double last_pts;                  // Capture timestamp of previous frame.
    double current_pts;               // Capture timestamp of current frame
    int nr_droppedframes;             // Counter for dropped frames.
    int frame_ready;                  // Signals availability of new frames for conversion into GL-Texture.
    int grabber_active;               // Grabber running?
    Rect roirect;                     // Region of interest rectangle - denotes subarea of full video capture area.
    double avg_decompresstime;        // Average time spent in Quicktime/Sequence Grabber decompressor.
    double avg_gfxtime;               // Average time spent in GWorld --> OpenGL texture conversion and statistics.
    int nrgfxframes;                  // Count of fetched textures.
	char capDeviceName[1000];		  // Camera name string.
	MatrixRecordPtr	scaleMatrixPtr;	  // Pointer to scaling matrix, if any. NULL otherwise.
	psych_thread		worker;		  // Handle of worker thread, if any. NULL otherwise.
} PsychVidcapRecordType;

static PsychVidcapRecordType vidcapRecordBANK[PSYCH_MAX_CAPTUREDEVICES];
static int numCaptureRecords = 0;
static psych_bool firsttime = TRUE;

/** PsychVideoCaptureDataProc
 *  This callback is called by the SequenceGrabber subsystem whenever a new frame arrives from
 *  the framegrabber hardware. It performs first-time setup of decompression sequence/codecs and
 *  decompression/conversion of grabbed raw data into a bitmap image -> GWorld.
 */
OSErr PsychVideoCaptureDataProc(SGChannel c, Ptr p, long len, long *offset, long chRefCon, TimeValue time, short writeType, long refCon)
{
    int handle;
    OSErr err = noErr;
    CodecFlags	ignore;
    TimeScale 	timeScale;
    double tstart, tend;
	short exph, expw;				
	Fixed scaleX, scaleY;

    PsychGetAdjustedPrecisionTimerSeconds(&tstart);
    
    // Retrieve handle to our capture data structure:
    handle = (int) refCon;

	// Check if we're called on a video channel. If we're called on a sound channel,
	// we simply return -- nothing to do in that case.
	if (vidcapRecordBANK[handle].sgchanVideo != c) {
		// printf("AUDIOCALLBACK AT %lf ...\n", tstart);
		return(noErr);
	}
	
    // Compute capture timestamp:
    err = SGGetChannelTimeScale(c, &timeScale);
    vidcapRecordBANK[handle].current_pts = (double) time / (double) timeScale;
    
    // GWorld for offscreen rendering available?
    if (vidcapRecordBANK[handle].gworld && vidcapRecordBANK[handle].grabber_active) {
        // Yes.
        
        // First time invocation for this sequence grabber?
        if (vidcapRecordBANK[handle].decomSeq == 0) {
            // Need to do one-time setup of decompression sequence:
            ImageDescriptionHandle imageDesc = (ImageDescriptionHandle) NewHandle(0);

            // retrieve a channelÕs current sample description, the channel returns a sample description that is
            // appropriate to the type of data being captured
            err = SGGetChannelSampleDescription(c, (Handle)imageDesc);
			if (PsychPrefStateGet_Verbosity() > 4) {
				printf("PTB-INFO: Video source %i - In proc: Video size is %i x %i pixels.\n", handle, (*imageDesc)->width, (*imageDesc)->height);
				printf("PTB-INFO: Video source %i - In proc: Video res is  %i x %i pixels.\n", handle, Fix2Long((*imageDesc)->hRes), Fix2Long((*imageDesc)->vRes));
			}

			// Matching roirect (expected size of image and dimension of target GWorld) and real image size?
			if ((*imageDesc)->height != (vidcapRecordBANK[handle].roirect.bottom - vidcapRecordBANK[handle].roirect.top) ||
				(*imageDesc)->width != (vidcapRecordBANK[handle].roirect.right - vidcapRecordBANK[handle].roirect.left)) {
				exph = (vidcapRecordBANK[handle].roirect.bottom - vidcapRecordBANK[handle].roirect.top);
				expw = (vidcapRecordBANK[handle].roirect.right - vidcapRecordBANK[handle].roirect.left);

				if (PsychPrefStateGet_Verbosity() > 1) {	
					printf("PTB-WARNING: Real size of video image from camera %i ['%s'] %i x %i doesn't match expected size %i x %i!\n", handle, vidcapRecordBANK[handle].capDeviceName, (*imageDesc)->width, (*imageDesc)->height, expw, exph);
					printf("PTB-WARNING: Will perform some software scaling to compensate for this. However this will incur significant performance loss!\n");
					printf("PTB-WARNING: It is strongly recommended that you specify a roirect of [0, 0, %i, %i] in the Screen('OpenVideocapture')\n", (*imageDesc)->width, (*imageDesc)->height);
					printf("PTB-WARNING: call to set a proper image resolution, so this compute intense scaling can be avoided in future sessions!\n");
				}
				
				// Build scaling matrix to compensate for mismatch in size of input image and
				// destination GWorld / roirect:
				vidcapRecordBANK[handle].scaleMatrixPtr = malloc(sizeof(MatrixRecord));
				SetIdentityMatrix(vidcapRecordBANK[handle].scaleMatrixPtr);
				scaleX = FixRatio(expw, (*imageDesc)->width);
				scaleY = FixRatio(exph, (*imageDesc)->height);
				ScaleMatrix(vidcapRecordBANK[handle].scaleMatrixPtr, scaleX, scaleY, 0, 0);
			}
			else {
				// No scaling matrix needed: Assign NULL matrix:
				vidcapRecordBANK[handle].scaleMatrixPtr = NULL;
			}
			
            // begin the process of decompressing a sequence of frames
            // this is a set-up call and is only called once for the sequence - the ICM will interrogate different codecs
            // and construct a suitable decompression chain, as this is a time consuming process we don't want to do this
            // once per frame (eg. by using DecompressImage)
            // for more information see Ice Floe #8 http://developer.apple.com/quicktime/icefloe/dispatch008.html
            // the destination is specified as the GWorld
            err = DecompressSequenceBegin(&(vidcapRecordBANK[handle].decomSeq),	// pointer to field to receive unique ID for sequence
                                          imageDesc,                            // handle to image description structure
                                          vidcapRecordBANK[handle].gworld,      // port for the DESTINATION image
                                          NULL,					// graphics device handle, if port is set, set to NULL
                                          &(vidcapRecordBANK[handle].roirect),	// source rectangle defining the portion of the image to decompress 
                                          vidcapRecordBANK[handle].scaleMatrixPtr,	// transformation matrix
                                          srcCopy,				// transfer mode specifier
                                          (RgnHandle)NULL,                      // clipping region in dest. coordinate system to use as a mask
                                          0,					// flags
                                          ((vidcapRecordBANK[handle].recordingflags & 32) ? codecMaxQuality : codecNormalQuality),  // accuracy in decompression
                                          ((vidcapRecordBANK[handle].recordingflags & 32) ? bestFidelityCodec : bestSpeedCodec));   // compressor identifier or special identifiers ie. bestSpeedCodec
            if (err!=noErr) {
                printf("PTB-ERROR: Error in Video capture callback!!!\n");
                fflush(NULL);
            }
            
            DisposeHandle((Handle)imageDesc);         
            // printf("DECOMPRESS-ONE_TIME!\n"); fflush(NULL);
            // One-Time setup of decompression engine done.
        }
        
        // Decompress new frame into our offscreen GWorld:
        err = DecompressSequenceFrameS(vidcapRecordBANK[handle].decomSeq,	// sequence ID returned by DecompressSequenceBegin
                                       p,					// pointer to compressed image data
                                       len,					// size of the buffer
                                       0,					// in flags
                                       &ignore,				// out flags
                                       NULL);				// async completion proc
        
        if (err!=noErr) {
            printf("PTB-ERROR: Error in Video capture callback code %i!!!\n", (int) err);
            fflush(NULL);
        }
        
        // Now we should have the required texture data in our GWorld...
        // Increment the newimage - flag:
        vidcapRecordBANK[handle].frame_ready++;
    
        PsychGetAdjustedPrecisionTimerSeconds(&tend);
        vidcapRecordBANK[handle].avg_decompresstime+=(tend - tstart);
        // Update framecounter:
        vidcapRecordBANK[handle].nrframes++;
    }

    return(noErr);
}

/* Main routine of video capture background worker thread: */
void* PsychQTVideoCaptureThreadMain(void* vidcapRecordToCast)
{
	double idlesecs;
	int error;
	
	// Cast to our video capture record struct:
	PsychVidcapRecordType *dev = (PsychVidcapRecordType*) vidcapRecordToCast;
	
	#if PSYCH_SYSTEM != PSYCH_WINDOWS
	// Signal to QuickTime that this is a separate thread.
	if ((error = EnterMoviesOnThread(0)) != noErr) {
		printf("PTB-ERROR: In PsychQTVideoCaptureThreadMain(): EnterMoviesOnThread() returns error code %i\n", (int) error);

		// Emergency exit:
		return (NULL);
	}
	#endif

	// Is grabbing logically active? Otherwise we're done:
	while (dev->grabber_active > 0) {
		// Grabber active. Wait a quarter capture cycle duration. This is save, even
		// if grabber should become inactive during this period, as we will recheck below.
		// We wait imprecisely, as a msecs more or less won't matter here:
		idlesecs = (1.0 / dev->fps) / 4;
		PsychYieldIntervalSeconds(idlesecs);
		
		// Time for work. Lock the mutex:
		PsychLockMutex(&(dev->mutex));
		
		// Do the SGIdle() call if device still active:
		if (dev->grabber_active > 0) {
			if ((error = (int) SGIdle(dev->seqGrab))!= (int) noErr) {
				// We don't abort on non noErr case, but only (optionally) report it. This is because when in harddisc
				// movie recording mode with sound recording, we can get intermittent errors in SGIdle() which are
				// meaningless so they are best ignored and must not lead to interruption of PTB operation.
				if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: In PsychQTVideoCaptureThreadMain(): SGIdle() returns error code %i\n", (int) error);
			}
		}

		// Done. Unlock mutex:
		PsychUnlockMutex(&(dev->mutex));

		// Next iteration...
	}

	#if PSYCH_SYSTEM != PSYCH_WINDOWS
	if ((error = ExitMoviesOnThread()) != noErr) {
		printf("PTB-ERROR: In PsychQTVideoCaptureThreadMain(): ExitMoviesOnThread() returns error code %i\n", (int) error);
	}
	#endif
	
	// Done. Return a NULL result:
	return(NULL);
}

/*
 *     PsychVideoCaptureInit() -- Initialize video capture subsystem.
 *     This routine is called by Screen's RegisterProject.c PsychModuleInit()
 *     routine at Screen load-time. It clears out the vidcapRecordBANK to
 *     bring the subsystem into a clean initial state.
 */
void PsychQTVideoCaptureInit(void)
{
    // Initialize vidcapRecordBANK with NULL-entries:
    int i;
    for (i=0; i < PSYCH_MAX_CAPTUREDEVICES; i++) {
        vidcapRecordBANK[i].gworld = (GWorldPtr) NULL;
        vidcapRecordBANK[i].seqGrab = (SeqGrabComponent) NULL;
        vidcapRecordBANK[i].decomSeq = 0;
        vidcapRecordBANK[i].grabber_active = 0;
		vidcapRecordBANK[i].worker = NULL;
    }    
    numCaptureRecords = 0;
    
    return;
}

void PsychQTEnumerateVideoSources(int outPos)
{
	PsychGenericScriptType 	*devs;
	const char *FieldNames[]={"DeviceIndex", "ClassIndex", "InputIndex", "ClassName", "InputName"};

    SeqGrabComponent	seqGrab;
	SGChannel			sgchannel;
    SGChannel			*sgchanptr;
	OSErr				error;
	SGDeviceList		sgdeviceList;
    int					i, j, selectedIndex;
	char				port_str[64];
	char				class_str[64];
    Str63				devPStr;
	Str255				outDeviceName;
	Str255				outInputName;
	Str255				devName;
    char				msgerr[10000];
	int					deviceClass, deviceInput;
	int					n;
	
    // We startup the Quicktime subsystem only on first invocation.
    if (firsttime) {
#if PSYCH_SYSTEM == PSYCH_WINDOWS
        // Initialize Quicktime for Windows compatibility layer: This will fail if
        // QT isn't installed on the Windows machine...
        error = InitializeQTML(0);
        if (error!=noErr) {
            PsychErrorExitMsg(PsychError_internal, "Quicktime Media Layer initialization failed: Quicktime not properly installed?!?");
        }
#endif

        // Initialize Quicktime-Subsystem:
        error = EnterMovies();
        if (error!=noErr) {
            PsychErrorExitMsg(PsychError_internal, "Quicktime EnterMovies() failed!!!");
        }
        firsttime = FALSE;
    }
	
	// Ok, rather ugly. We create a sequence grabber and associated video channel
	// for the sole purpose of enumeration:
    seqGrab = OpenDefaultComponent(SeqGrabComponentType, 0);
    if (seqGrab == NULL) {
		printf("PTB-ERROR: Failed to open sequence grabber video capture component for enumeration!");
        PsychErrorExitMsg(PsychError_internal, "Failed to open requested sequence grabber for video device enumeration!");
    }

    // Initialize the sequence grabber component:
    error = SGInitialize(seqGrab);
    if (error != noErr) {
        if (seqGrab) CloseComponent(seqGrab);
        PsychErrorExitMsg(PsychError_internal, "SGInitialize() for capture device enumeration failed!"); 
    }
	
	sgchanptr = &sgchannel;
	error = SGNewChannel(seqGrab, VideoMediaType, sgchanptr);
    if (error == noErr) {
		// Enumerate all available devices:
		error = SGGetChannelDeviceList(*sgchanptr, sgDeviceListIncludeInputs, &sgdeviceList);
		if (error == noErr) {
			// First scan: Find out how many devices there in total:
			n = 0;
			for (i = 0; i< (*sgdeviceList)->count; i++) {				
				if(NULL != (((*sgdeviceList)->entry[i]).inputs)) {
					n += (*(((*sgdeviceList)->entry[i]).inputs))->count;

					// One extra slot for default input in this class, if any:
					if ((*(((*sgdeviceList)->entry[i]).inputs))->count > 0) n++;
				}
			}

			// One extra slot for the default device, if any:
			if (n > 0) {
				n++;
			}

			// Create output struct array with n output slots:
			PsychAllocOutStructArray(outPos, TRUE, n, 5, FieldNames, &devs);
			n = 0;

			// Second scan: Actually create the output entries:
			for (i = 0; i< (*sgdeviceList)->count; i++) {
				p2cstrcpy(class_str, (((*sgdeviceList)->entry[i]).name));
				// printf("Device Nr. %i is: %s\n\n", i, class_str);
				
				if(NULL != (((*sgdeviceList)->entry[i]).inputs)) {
					for (j = 0; j < (*(((*sgdeviceList)->entry[i]).inputs))->count; j++) {
						p2cstrcpy(port_str, (*(((*sgdeviceList)->entry[i]).inputs))->entry[j].name);
						// printf("Device [%i -> %i]: %s\n", i, j, port_str);
						
						deviceClass = i + 1;
						deviceInput = j + 1;
						PsychSetStructArrayDoubleElement("DeviceIndex", n, deviceClass * 10000 + deviceInput, devs);
						PsychSetStructArrayDoubleElement("ClassIndex", n, deviceClass, devs);
						PsychSetStructArrayDoubleElement("InputIndex", n, deviceInput, devs);
						PsychSetStructArrayStringElement("ClassName", n, class_str, devs);
						PsychSetStructArrayStringElement("InputName", n, port_str, devs);

						// Next slot:
						n++;
						
						if (j == (*(((*sgdeviceList)->entry[i]).inputs))->selectedIndex) {
							// Special case: Default device input 0 for this class:
							PsychSetStructArrayDoubleElement("DeviceIndex", n, deviceClass * 10000, devs);
							PsychSetStructArrayDoubleElement("ClassIndex", n, deviceClass, devs);
							PsychSetStructArrayDoubleElement("InputIndex", n, deviceInput, devs);
							PsychSetStructArrayStringElement("ClassName", n, class_str, devs);
							PsychSetStructArrayStringElement("InputName", n, port_str, devs);
							
							// Next slot:
							n++;
						}
					}
				}
			}
			
			// If at least one device available, then add the default 0 slot for the default device:
			if (n > 0) {
				i = (*sgdeviceList)->selectedIndex;
				deviceClass = i + 1;
				p2cstrcpy(class_str, (((*sgdeviceList)->entry[i]).name));
				j = (*(((*sgdeviceList)->entry[i]).inputs))->selectedIndex;
				deviceInput = j + 1;
				p2cstrcpy(port_str, (*(((*sgdeviceList)->entry[i]).inputs))->entry[j].name);

				PsychSetStructArrayDoubleElement("DeviceIndex", n, 0, devs);
				PsychSetStructArrayDoubleElement("ClassIndex", n, deviceClass, devs);
				PsychSetStructArrayDoubleElement("InputIndex", n, deviceInput, devs);
				PsychSetStructArrayStringElement("ClassName", n, class_str, devs);
				PsychSetStructArrayStringElement("InputName", n, port_str, devs);
				
				// Next slot:
				n++;
			}

			SGDisposeDeviceList(seqGrab, sgdeviceList);
		}
	}

releaseQTEnum:
	// Done with enumeration: Release ressources:
	SGDisposeChannel(seqGrab, *sgchanptr);
    CloseComponent(seqGrab);

	return;
}

OSErr PsychQTSelectVideoSource(SeqGrabComponent seqGrab, SGChannel* sgchanptr, int deviceIndex, psych_bool showOutput, char* retDeviceName)
{
	ComponentInstance	vdCompInst;
	OSErr				error;
	SGDeviceList		sgdeviceList;
    int					i, j, selectedIndex;
	char				port_str[64];
    Str63				devPStr;
	Str255				outDeviceName;
	Str255				outInputName;
	Str255				devName;
    char				msgerr[10000];
	int					deviceClass, deviceInput;

	// Child protection:
	if (deviceIndex < 0) PsychErrorExitMsg(PsychError_user, "Invalid negative 'deviceIndex' provided. Must be zero or a positive integral number!");

	// Default to no Error:
	error = noErr;

	// Special case zero? This means "default device" or take whatever's setup by default,
	// i.e., change absolutely nothing.
	if (deviceIndex > 0) {
		// Non-Zero deviceIndex. Split up into device class and subdevice (input):
		deviceClass = deviceIndex / 10000;
		deviceInput = deviceIndex % 10000;
		
		// Enumerate all available devices:
		error = SGGetChannelDeviceList(*sgchanptr, sgDeviceListIncludeInputs, &sgdeviceList);
		if (error == noErr) {
			if (showOutput && PsychPrefStateGet_Verbosity() > 4) {
				printf("Number of available input video device classes: %i\n", (*sgdeviceList)->count);
				for (i = 0; i< (*sgdeviceList)->count; i++) {
					p2cstrcpy(port_str, (((*sgdeviceList)->entry[i]).name));
					printf("Device class Nr. %i is: %s\n", i+1, port_str);
				}
				printf("\n\n");
			}

			// deviceClass zero selected? This would mean that we shall accept the default
			// device class selection as made by the OS or external software, but only select
			// the subdevice or input channel within that class:
			if (deviceClass == 0) {
				// Default device class: Just fetch currently selected class:
				selectedIndex = (*sgdeviceList)->selectedIndex;				
			}
			else {
				// Specific device class selected: Validate and try to choose it:
				selectedIndex = deviceClass - 1;
				if (selectedIndex >= (*sgdeviceList)->count) {
					printf("Given deviceIndex %i requests device class %i, but no such device class exists! Failed.\n", deviceIndex, deviceClass);
					PsychErrorExitMsg(PsychError_user, "Invalid 'deviceIndex' for this system setup provided.");
				}
				
				error =	SGSetChannelDevice(*sgchanptr, (StringPtr) &(((*sgdeviceList)->entry[selectedIndex]).name));
				if (error!=noErr) {
					p2cstrcpy(port_str, (((*sgdeviceList)->entry[selectedIndex]).name));
					if (showOutput && PsychPrefStateGet_Verbosity()>2) printf("PTB-ERROR: Failed to select video input device class %i [%s] for some reason!\n", deviceClass, port_str);

					// Revert to default which is still active after failed switch:
					selectedIndex = (*sgdeviceList)->selectedIndex;					
				}
			}

			// At this point we have either the default class selected, or the requested
			// device class. 'selectedIndex' points to the selected class.
			p2cstrcpy(port_str, (((*sgdeviceList)->entry[selectedIndex]).name));
			if (showOutput && PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Selected input video device class now is %i [Quicktime index %i] aka: %s\n", selectedIndex + 1, selectedIndex, port_str);

			// What about the subclass aka input channel? A zero means - Whatever's the default
			// in that class. A non-zero value i means: Selected i'th input. For a zero value,
			// we've got nothing to do:
			if (deviceInput > 0) {
				i = selectedIndex;

				// Any inputs for this class?
				if(NULL != (((*sgdeviceList)->entry[i]).inputs)) {
					if (showOutput && PsychPrefStateGet_Verbosity() > 4) {
						printf("\nMapping before device selection:\n");
						for (j = 0; j < (*(((*sgdeviceList)->entry[i]).inputs))->count; j++) {
							p2cstrcpy(port_str, (*(((*sgdeviceList)->entry[i]).inputs))->entry[j].name);
							printf("Device [%i -> %i]: %s\n", i, j, port_str);
						}
						printf("\n\n");
					}
					
					// Yep. At least deviceInput's input available?
					if (deviceInput - 1 < (*(((*sgdeviceList)->entry[i]).inputs))->count) {
						// Select this input:
						error = SGSetChannelDeviceInput(*sgchanptr, deviceInput - 1);
						if (error!=noErr) {
							// Grabber didn't accept new input video channel:
							if (showOutput && PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Video digitizer component didn't accept input video channel %i as requested by deviceIndex!\n", deviceInput);
						}
						else {
							if (showOutput && PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Video digitizer component accepted input %i!\n", deviceInput);
						}						
					}
					else {
						p2cstrcpy(port_str, (((*sgdeviceList)->entry[selectedIndex]).name));
						if (showOutput && PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Failed to select %i. video input inside device class %i [%s],\nbecause class only has %i inputs!\n", deviceInput, selectedIndex + 1, port_str, (*(((*sgdeviceList)->entry[i]).inputs))->count);
					}

					if (showOutput && PsychPrefStateGet_Verbosity() > 4) {
						printf("\nMapping after device selection:\n");
						for (j = 0; j < (*(((*sgdeviceList)->entry[i]).inputs))->count; j++) {
							p2cstrcpy(port_str, (*(((*sgdeviceList)->entry[i]).inputs))->entry[j].name);
							printf("Device [%i -> %i]: %s\n", i, j, port_str);
						}
						printf("\n\n");
					}
				}
				else {
					// Failed because out of range:
					if (showOutput && PsychPrefStateGet_Verbosity() > 2) printf("PTB-ERROR: Failed to select video input inside device class %i [%s], because class doesn't seem to have inputs available!\n", deviceClass, port_str);
				}
			}
			
			if (0) {
				for (i = 0; i< (*sgdeviceList)->count; i++) {
					p2cstrcpy(port_str, (((*sgdeviceList)->entry[i]).name));
					printf("Device Nr. %i is: %s\n\n", i, port_str);
					
					if(NULL != (((*sgdeviceList)->entry[i]).inputs)) {
						for (j = 0; j < (*(((*sgdeviceList)->entry[i]).inputs))->count; j++) {
							p2cstrcpy(port_str, (*(((*sgdeviceList)->entry[i]).inputs))->entry[j].name);
							printf("Device [%i -> %i]: %s\n", i, j, port_str);
						}
					}
				}
				
				memset(devPStr, 0, sizeof(Str63));
				for (i = 0; i < (*(((*sgdeviceList)->entry[selectedIndex]).inputs))->count; i++) {
					p2cstrcpy(port_str, (*(((*sgdeviceList)->entry[selectedIndex]).inputs))->entry[i].name);
					printf("Device [%i -> %i]: %s\n", selectedIndex, i, port_str);
				}
				
				if (PsychPrefStateGet_Verbosity() > 4) {			
					if ((vdCompInst = SGGetVideoDigitizerComponent(*sgchanptr))) {
						msgerr[0] = 0;
						VDGetInputName(vdCompInst, deviceIndex, devName);
						p2cstrcpy(msgerr, devName);
						printf("PTB-INFO: Selected Video digitizer component has device name %s!\n", msgerr);
					}
				}
				
			}

			// Release device list:
			SGDisposeDeviceList(seqGrab, sgdeviceList);			
		}
		else {
			// Enumeration failure - Game over:
            if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-ERROR: Failed to enumerate available video devices! Will stick to default video device!\n");
		}
	}

	// Retrieve and return capture device name:
	SGGetChannelDeviceAndInputNames(*sgchanptr, outDeviceName, outInputName, nil);
	if (retDeviceName) p2cstrcpy(retDeviceName, outDeviceName);

	// Status output about selection requested?
	if (showOutput && PsychPrefStateGet_Verbosity() > 3) {
		// Display name of selected video device(-class) and input:
		p2cstrcpy(msgerr, outDeviceName);
		printf("PTB-INFO: Selected video input device for deviceIndex %i has outDeviceName: %s\n", deviceIndex, msgerr);
		p2cstrcpy(msgerr, outInputName);
		printf("PTB-INFO: Selected video input device for deviceIndex %i has outInputName: %s\n", deviceIndex, msgerr);
	}

	// Return final status - Should be zero:
	return(error);
}

/*
 *      PsychQTOpenVideoCaptureDevice() -- Create a video capture object.
 *
 *      This function tries to open a Quicktime-Sequencegrabber
 *      and return the associated captureHandle for it.
 *
 *      win = Pointer to window record of associated onscreen window.
 *      deviceIndex = Index of the grabber device. (Currently ignored)
 *      capturehandle = handle to the new capture object.
 *      capturerectangle = If non-NULL a ptr to a PsychRectangle which contains the ROI for capture.
 *      The following arguments are currently ignored on Windows and OS-X:
 *      reqdepth = Number of layers for captured output textures. (0=Don't care, 1=LUMINANCE8, 2=LUMINANCE8_ALPHA8, 3=RGB8, 4=RGBA8)
 *      num_dmabuffers = Number of buffers in the ringbuffer queue (e.g., DMA buffers) - This is OS specific. Zero = Don't care.
 *      allow_lowperf_fallback = If set to 1 then PTB can use a slower, low-performance fallback path to get nasty devices working.
 *		targetmoviefilename = NULL == Only live capture, non-NULL == Pointer to char-string with name of target QT file for video recording.
 *		recordingflags = Only used for recording: Request audio recording, ram recording vs. disk recording and such...
 *		// Query optional movie recording flags:
 *		// 0 = Record video, stream to disk immediately (slower, but unlimited recording duration).
 *		// 1 = Record video, stream to memory, then at end of recording to disk (limited duration by RAM size, but faster).
 *		// 2 = Record audio as well.
 *		// 4 = Only record, but don't provide as live feed. This saves callback overhead if provided for pure recording sessions.
 *		// 8 = Don't call SGPrepare() / SGRelease().
 *		// 16 = Use background worker thread to call SGIdle() for automatic capture/recording.
 *		// 32 = Try to select the highest quality codec for texture creation - codecLosslessQuality instead of codecNormalQuality.
 */
psych_bool PsychQTOpenVideoCaptureDevice(int slotid, PsychWindowRecordType *win, int deviceIndex, int* capturehandle, double* capturerectangle,
				 int reqdepth, int num_dmabuffers, int allow_lowperf_fallback, char* targetmoviefilename, unsigned int recordingflags)
{
    int i;
    OSErr error;
    char msgerr[10000];
    char errdesc[1000];
    Rect movierect, newrect;
    SeqGrabComponent seqGrab = NULL;

    SGChannel *sgchanptr = NULL;
	SGChannel *sgchanaudioptr = NULL;
    ImageDescriptionHandle imageDesc;
    Fixed framerate;
    *capturehandle = -1;
    error=noErr;
#if PSYCH_SYSTEM != PSYCH_WINDOWS
	TimeBase soundTimeBase = NULL;
	TimeBase sgTimeBase = NULL;
	Component clockComponent = NULL;
	ComponentDescription looking;

	CodecNameSpecListPtr	codeclist = NULL;
#endif

    // We startup the Quicktime subsystem only on first invocation.
    if (firsttime) {
#if PSYCH_SYSTEM == PSYCH_WINDOWS
        // Initialize Quicktime for Windows compatibility layer: This will fail if
        // QT isn't installed on the Windows machine...
        error = InitializeQTML(0);
        if (error!=noErr) {
            PsychErrorExitMsg(PsychError_internal, "Quicktime Media Layer initialization failed: Quicktime not properly installed?!?");
        }
#endif

        // Initialize Quicktime-Subsystem:
        error = EnterMovies();
        if (error!=noErr) {
            PsychErrorExitMsg(PsychError_internal, "Quicktime EnterMovies() failed!!!");
        }
        firsttime = FALSE;
    }
    
    // Slot 'slotid' will contain the record for our new capture object:
    
    // Zero-out new record:
    vidcapRecordBANK[slotid].gworld=NULL;
    vidcapRecordBANK[slotid].decomSeq=0;    
    vidcapRecordBANK[slotid].grabber_active = 0;
	vidcapRecordBANK[slotid].scaleMatrixPtr = NULL;

    // Open sequence grabber:
    // ======================

    seqGrab = OpenDefaultComponent(SeqGrabComponentType, 0);
    if (seqGrab == NULL) {
		printf("PTB-ERROR: Failed to open sequence grabber video capture component for deviceIndex %i!", deviceIndex);
        PsychErrorExitMsg(PsychError_internal, "Failed to open requested sequence grabber for video capture!");
    }
    
    // Initialize the sequence grabber component:
    error=noErr;

    error = SGInitialize(seqGrab);
    if (error != noErr) {
        if (seqGrab) CloseComponent(seqGrab);
        PsychErrorExitMsg(PsychError_internal, "SGInitialize() for capture device failed!"); 
    }

    // Specify the destination data reference for a record operation
    // tell it we're not making a movie
    // if the flag seqGrabDontMakeMovie is used, the sequence grabber still calls
    // your data function, but does not write any data to the movie file
    // writeType will always be set to seqGrabWriteAppend
    error = SGSetDataRef(seqGrab, 0, 0, (targetmoviefilename == NULL) ? seqGrabDontMakeMovie : seqGrabToDisk);
    if (error !=noErr) {
        if (seqGrab) CloseComponent(seqGrab);
        PsychErrorExitMsg(PsychError_internal, "SGSetDataRef for capture device failed!");            
    }

	// Set dummy GWorld - we need this to prevent SGNewChannel from crashing on Windoze.
	SGSetGWorld(seqGrab, 0, 0);

    // Create and setup video channel on sequence grabber:
    sgchanptr = &(vidcapRecordBANK[slotid].sgchanVideo);
	error = SGNewChannel(seqGrab, VideoMediaType, sgchanptr);
    if (error == noErr) {
		// Select device and input for this video channel: Returned error code is more
		// informative than important: The function will always select some device if
		// it returns, even if it isn't according to given spec, which would be apparent
		// by  a non-zero error:
		error = PsychQTSelectVideoSource(seqGrab, sgchanptr, deviceIndex, FALSE, (char*) &(vidcapRecordBANK[slotid].capDeviceName));

		// Some low-level debug output, if requested:
        SGGetVideoRect(*sgchanptr, &newrect);
		if (PsychPrefStateGet_Verbosity()>4) printf("PTB-INFO: VideoRect ROI is %i,%i,%i,%i\n", newrect.left, newrect.top, newrect.right, newrect.bottom);

		SGGetChannelBounds(*sgchanptr, &newrect);
		if (PsychPrefStateGet_Verbosity()>4) printf("PTB-INFO: ChannelBounds ROI is %i,%i,%i,%i\n", newrect.left, newrect.top, newrect.right, newrect.bottom);

		// Retrieve size of the capture rectangle - and therefore size of
        // our GWorld for offscreen rendering:
        error = SGGetSrcVideoBounds(*sgchanptr, &movierect);
		if (PsychPrefStateGet_Verbosity()>4) printf("PTB-INFO: [%i] SrcVideoBounds ROI is %i,%i,%i,%i\n", error, movierect.left, movierect.top, movierect.right, movierect.bottom);

		// Reasonable movierect returned?
		if ((noErr != error) || ((movierect.right == 1600) && (movierect.bottom == 1200) && (PSYCH_SYSTEM == PSYCH_OSX))) {
			// Either none returned or the typical 1600 x 1200 default IIDC max rectangle on OS/X.
			// In both cases we can be quite certain that the values are bogus.
			if (NULL == capturerectangle) {
				// Only output warning message about bogus capture area settings and resulting 640 x 480 override
				// if the device namestring doesn't contain "iSight", as we know the iSight behaves like that on
				// every Apple computer, so no point of cluttering the screen with warnings about the obvious.
				// Same goes for the "unibrain Fire-i":
				if ((NULL == strstr(vidcapRecordBANK[slotid].capDeviceName, "iSight")) && (NULL == strstr(vidcapRecordBANK[slotid].capDeviceName, "unibrain Fire-i")) && (PsychPrefStateGet_Verbosity() > 2)) {
					printf("\nPTB-INFO: In Screen('OpenVideoCapture',...) for video capture device '%s' with deviceIndex %i:\n", vidcapRecordBANK[slotid].capDeviceName, deviceIndex);
					printf("PTB-INFO: Your code doesn't specify an explicit 'roirectangle' for requested camera resolution, but wants me to\n");
					printf("PTB-INFO: auto-detect the optimal image size. The operating system recommends a size of 1600 x 1200 pixels for your\n");
					printf("PTB-INFO: camera, which is almost always a wrong and bogus setting! I'll therefore default to the most commonly found\n");
					printf("PTB-INFO: reasonable default of 640 x 480 pixels and hope that this is ok. If this doesn't work or if you have indeed\n");
					printf("PTB-INFO: a high-end camera with 1600 x 1200 pixels resolution connected, then please specify the proper settings via\n");
					printf("PTB-INFO: the 'roirectangle' parameter in Screen('OpenVideoCapture', ...);\n\n");
				}
				
				// Assign more reasonable default: Works with most entry level and intermediate level
				// IIDC cameras, most NTSC video sources and with the Apple iSight:
				movierect.left = 0;
				movierect.top = 0;
				movierect.right = 640;
				movierect.bottom = 480;
			}
		}

        error = noErr;
		
        // Capture-Rectangle (ROI) specified?
        if (capturerectangle) {
            // Yes. Try to set it up...
            // Assign new requested ROI:
            newrect.left=(short) capturerectangle[kPsychLeft];
            newrect.top=(short)  capturerectangle[kPsychTop];
            newrect.right=(short)  capturerectangle[kPsychRight];
            newrect.bottom=(short) capturerectangle[kPsychBottom];
            
            if (PsychPrefStateGet_Verbosity()>3) printf("PTB-INFO: Selected video capture ROI is %i,%i,%i,%i\n", newrect.left, newrect.top, newrect.right, newrect.bottom);
            
            if ((int) capturerectangle[kPsychLeft]<movierect.left || (int) capturerectangle[kPsychTop]<movierect.top ||
                (int) capturerectangle[kPsychRight]>movierect.right || (int) capturerectangle[kPsychBottom]>movierect.bottom) {
                // ROI is not a subrectangle of video digitizers maximum ROI!
                SGDisposeChannel(seqGrab, *sgchanptr);
                *sgchanptr = NULL;
                CloseComponent(seqGrab);
                PsychErrorExitMsg(PsychError_user, "Invalid video region of interest (not inside image) specified!");
            }
            
            // Try to set our own custom video capture rectangle for the digitizer hardware:
			// Disabled: Does more harm than good! Instead just assign our override movierect:
//            error=SGSetVideoRect(*sgchanptr, &newrect);
//            if (error!=noErr) {
//                // Grabber didn't accept new rectangle :(
//                if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Video capture device didn't accept new capture area. Reverting to full hardware capture area,\n");
//                if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Trying to only process specified ROI by restricting conversion to ROI in software...\n");

                movierect.left=(int) capturerectangle[kPsychLeft];
                movierect.top=(int) capturerectangle[kPsychTop];
                movierect.right=(int) capturerectangle[kPsychRight];
                movierect.bottom=(int) capturerectangle[kPsychBottom];

//            }
//            else {
//                // Retrieve new capture rectangle settings:
//                error = SGGetVideoRect(*sgchanptr, &movierect);
//            }
        }
        
        // Store our roirect in structure:
        vidcapRecordBANK[slotid].roirect = movierect;
        movierect.right-=movierect.left;
        movierect.bottom-=movierect.top;
        movierect.left=0;
        movierect.top=0;
		
		// retrieve a channelÕs current sample description, the channel returns a sample description that is
		// appropriate to the type of data being captured
		//ImageDescriptionHandle imageDesc = (ImageDescriptionHandle) NewHandle(0);
		//error = SGGetChannelSampleDescription(*sgchanptr, (Handle)imageDesc);
		//printf("PTB-INFO: Capture rectangle width = %i x height = %i ...\n", (int) ((ImageDescription)(**imageDesc)).width, (int) ((ImageDescription)(**imageDesc)).height);
		//DisposeHandle((Handle) imageDesc);
		
        // Now that we know the movierect of our ROI,
        // destroy the channel, assign a properly sized GWorld
        // and recreate the channel:
        SGDisposeChannel(seqGrab, *sgchanptr);
        *sgchanptr = NULL;

		// Check if the requested pixel depth matches what we can give:
		switch (reqdepth) {
			case 2:
				// A no-go: Instead we use 1 channel luminance8:
				if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Video capture engine doesn't support requested Luminance+Alpha format. Will revert to pure Luminance instead...\n");
			case 1:
				reqdepth = 1;
			break;
			
			case 3:
				// A n-go: Instead we use 4 channels RGBA with a all 255 alpha channel --> More efficient.
				if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Video capture engine doesn't support requested RGB format. Will use slower software conversion of RGBA->RGB...\n");
				break;
			case 4:
			case 0:
				reqdepth = 4;
			break;
			
			default:
				// Unknown format:
				CloseComponent(seqGrab);
				PsychErrorExitMsg(PsychError_user, "You requested an invalid image depths (not one of 0, 1, 2, 3 or 4). Aborted.");
		} 

        // Create GWorld for this grabber object:
        error = QTNewGWorld(&vidcapRecordBANK[slotid].gworld, (reqdepth < 3) ? k8IndexedGrayPixelFormat : 0, &movierect,  NULL, NULL, 0);
        if (error!=noErr) {
          CloseComponent(seqGrab);
          PsychErrorExitMsg(PsychError_internal, "Quicktime GWorld creation for capture device failed!");
        }
    
        // Set grabbers graphics world to our GWorld:
        error = SGSetGWorld(seqGrab, vidcapRecordBANK[slotid].gworld, NULL );
        if (error !=noErr) {
          DisposeGWorld(vidcapRecordBANK[slotid].gworld);
          vidcapRecordBANK[slotid].gworld = NULL;
          CloseComponent(seqGrab);
          PsychErrorExitMsg(PsychError_internal, "Assignment of GWorld to capture device failed!");            
        }

        // Create and setup video channel on sequence grabber:
		error = SGNewChannel(seqGrab, VideoMediaType, sgchanptr);
        if (error !=noErr) {
          DisposeGWorld(vidcapRecordBANK[slotid].gworld);
          vidcapRecordBANK[slotid].gworld = NULL;
          CloseComponent(seqGrab);
          PsychErrorExitMsg(PsychError_internal, "Assignment of GWorld to capture device failed!");            
        }

		// Select device and input for this video channel:
		error = PsychQTSelectVideoSource(seqGrab, sgchanptr, deviceIndex, TRUE, (char*) &(vidcapRecordBANK[slotid].capDeviceName));
        if (error!=noErr) {
            // Grabber didn't accept new rectangle :(
            if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Failed to select video capture device according to requested deviceIndex %i\n", deviceIndex);
        }

		// Assign VDIG handle:
		vidcapRecordBANK[slotid].vdig = SGGetVideoDigitizerComponent(*sgchanptr);

		// More debug output:
		VDGetMaxSrcRect(vidcapRecordBANK[slotid].vdig, currentIn, &newrect);
		if (PsychPrefStateGet_Verbosity()>4) printf("PTB-INFO: VDGetMaxSrcRect ROI is %i,%i,%i,%i\n", newrect.left, newrect.top, newrect.right, newrect.bottom);

		VDGetDigitizerRect(vidcapRecordBANK[slotid].vdig, &newrect);
		if (PsychPrefStateGet_Verbosity()>4) printf("PTB-INFO: VDGetDigitizerRect ROI is %i,%i,%i,%i\n", newrect.left, newrect.top, newrect.right, newrect.bottom);
		
        // Try to set our own custom video capture rectangle:
		error = noErr;

        // Disabled: Does more harm than good: error=SGSetVideoRect(*sgchanptr, &newrect);
        if (error!=noErr) {
            // Grabber didn't accept new rectangle :(
            if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Video capture device didn't accept new capture area. Reverting to full area...\n"); fflush(NULL);
        }

        error = SGSetChannelBounds(*sgchanptr, &movierect);
        if (error == noErr) {
            // set usage for new video channel to avoid playthrough
            // note we don't set seqGrabPlayDuringRecord
            error = SGSetChannelUsage(*sgchanptr, seqGrabRecord | seqGrabLowLatencyCapture | seqGrabAlwaysUseTimeBase);
        }

        if (error != noErr) {
            // clean up on failure
            SGDisposeChannel(seqGrab, *sgchanptr);
            *sgchanptr = NULL;
            DisposeGWorld(vidcapRecordBANK[slotid].gworld);
            vidcapRecordBANK[slotid].gworld = NULL;
            if (seqGrab) CloseComponent(seqGrab);
            PsychErrorExitMsg(PsychError_internal, "SGSetChannelBounds() or SGSetChannelUsage() for capture device failed!");            
        }

        SGGetVideoRect(*sgchanptr, &newrect);
		if (PsychPrefStateGet_Verbosity()>4) printf("PTB-INFO: VideoRect ROI is %i,%i,%i,%i\n", newrect.left, newrect.top, newrect.right, newrect.bottom);
    }
    else {
        // clean up on failure
        *sgchanptr = NULL;
        if (seqGrab) CloseComponent(seqGrab);
        
        if (error == -9400 || error== -9405) PsychErrorExitMsg(PsychError_user, "Couldn't connect to video capture device! Device offline or disconnected?");
        PsychErrorExitMsg(PsychError_internal, "SGNewChannel() for capture device failed!");            
    }

    // Retrieve a channels current sample description, the channel returns a sample description that is
    // appropriate to the type of data being captured
    //    imageDesc = (ImageDescriptionHandle) NewHandle(0);
    //    error = SGGetChannelSampleDescription(vidcapRecordBANK[slotid].sgchanVideo, (Handle)imageDesc);
    //    movierect.left = 0;
    //    movierect.top = 0;
    //    movierect.right = (**imageDesc).width;
    //    movierect.bottom = (**imageDesc).height 
    //    DisposeHandle((Handle)imageDesc);

    // Specify a data callback function: Gets called whenever a new frame is ready.
	// We only spec a callback if this is either pure live video capture and processing,
	// or if this is video (and sound) recording as well and recordingflags don't tell
	// us this is "recording only". The callback is needed for video processing by
	// Matlab+PTB but it comes at a bit of extra overhead:
	if ((targetmoviefilename == NULL) || !(recordingflags & 4)) {
		error = SGSetDataProc(seqGrab, NewSGDataUPP(PsychVideoCaptureDataProc), slotid);
		if (error !=noErr) {
			DisposeGWorld(vidcapRecordBANK[slotid].gworld);
			vidcapRecordBANK[slotid].gworld = NULL;
			SGDisposeChannel(seqGrab, *sgchanptr);
			*sgchanptr = NULL;
			CloseComponent(seqGrab);
			PsychErrorExitMsg(PsychError_internal, "Assignment of capture callback fcn. to capture device failed!");            
		}
	}
	
    // Store a reference to our slotid for this channel. This gets passed to the
    // videocapture callback fcn. so it knows to which capture object to relate to...
    error = SGSetChannelRefCon(vidcapRecordBANK[slotid].sgchanVideo, slotid);
    if (error !=noErr) {
        DisposeGWorld(vidcapRecordBANK[slotid].gworld);
        vidcapRecordBANK[slotid].gworld = NULL;
        SGDisposeChannel(seqGrab, *sgchanptr);
        *sgchanptr = NULL;
        CloseComponent(seqGrab);
        PsychErrorExitMsg(PsychError_internal, "Assignment of Refcon to capture device failed!");            
    }
    
	// Sound channel requested as well?
	if (targetmoviefilename && (recordingflags & 2)) {
		// Create sound grabber component as well:
		sgchanaudioptr = &(vidcapRecordBANK[slotid].sgchanAudio);
		error = SGNewChannel(seqGrab, SoundMediaType, sgchanaudioptr);
		SGSetChannelRefCon(vidcapRecordBANK[slotid].sgchanAudio, slotid);
		// More setup code needed?
		if (error == noErr) error = SGSetChannelUsage(*sgchanaudioptr,  seqGrabRecord | seqGrabLowLatencyCapture | seqGrabAlwaysUseTimeBase);
	}
	
	// Recording to disc requested?
	if (targetmoviefilename==NULL) {
		// No. Just live capture and processing:
		error = SGSetDataOutput(seqGrab, NULL, seqGrabDontMakeMovie);
	}
	else {
		// Yes. Set it up:
        char* codecstr;
		FSSpec recfile;
		CodecType codectypeid;
		CompressorComponent codeccomp;
		codectypeid = 0;
		codecstr = strstr(targetmoviefilename, ":CodecType="); 
		if (codecstr) {
			sscanf(codecstr, ":CodecType= %i", &codectypeid);
			*codecstr = 0;
		}
		
		if (PsychPrefStateGet_Verbosity() > 3) {
			// Query and print list of all supported codecs on this setup:
			#if PSYCH_SYSTEM == PSYCH_OSX
			if ((GetCodecNameList(&codeclist, 1) == noErr) && (codeclist != NULL)) {
				printf("PTB-INFO: The following codecs are supported on your system for video recording:\n");
				for (i=0; i<codeclist->count; i++) {
					printf(":CodecType=%i  - '%s'  \n", (int) codeclist->list[i].cType, codeclist->list[i].typeName);
				}
				printf("--------------------------------------------------------------------------------\n\n");
				DisposeCodecNameList(codeclist); codeclist = NULL;
			}
			#else
				printf("PTB-INFO: Printout of codec list not yet implemented on M$-Windows due to some Microsoft compiler brain-damage.\n");
			#endif
		}
		
		NativePathNameToFSSpec(targetmoviefilename, &recfile, 0);

		if ((recordingflags & 1) && (recordingflags & 2)) {
			// Sound recording and memory recording requested? This is a no go, it doesn't work
			// as of OS/X 10.4.10, so sound overrides this: Switch to disk recording.
			if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Both, sound recording and recording to memory requested. This doesn't work! Switching to direct disk recording instead...\n");
			recordingflags = recordingflags & (~1);
		}
		 
		// If recordingflags & 1, then we request capture to memory with writeout at end of capture operation. Otherwise
		// we request immediate capture to disk. We always append to an existing movie file, instead of overwriting it.
		error = SGSetDataOutput(seqGrab, &recfile, ((recordingflags & 1) ? seqGrabToMemory : seqGrabToDisk));
		
		// This call would select a specific video compressor, if we had any except the default one ;)
		if (error==noErr && codectypeid!=0) {
			// MK: Could do this to change compression settings, but don't:
			//			error = SGSetVideoCompressorType(vidcapRecordBANK[slotid].sgchanVideo, codectypeid);
			//			error = SGGetVideoCompressor(vidcapRecordBANK[slotid].sgchanVideo, NULL, &codeccomp, NULL, NULL, NULL);
			//			error = SGSetVideoCompressor(vidcapRecordBANK[slotid].sgchanVideo, 0, codeccomp, codecHighQuality, 0, 1);
			error = SGSetVideoCompressorType(vidcapRecordBANK[slotid].sgchanVideo, codectypeid);
			if (error != noErr && PsychPrefStateGet_Verbosity() > 1) {
				printf("PTB-WARNING: Video recording engine could not enable requested codec of type id %i: QT Error code %i.\n", (int) codectypeid, (int) error);
				printf("PTB-HINT:    Rerun the script with a 'Verbosity' preference setting of greater than or equal to 4 to print out a list of supported codecs.\n");		
			}
			error=noErr;
		}
	}
	
    if (error !=noErr) {
        DisposeGWorld(vidcapRecordBANK[slotid].gworld);
        vidcapRecordBANK[slotid].gworld = NULL;
        SGDisposeChannel(seqGrab, *sgchanptr);
        *sgchanptr = NULL;
        CloseComponent(seqGrab);
        PsychErrorExitMsg(PsychError_internal, "Assignment of SGSetDataOutput() to capture device failed!");            
    }

    // Get ready! Unless recordingflags set to 8, in which case we won't SGPrepare():
	error = noErr;
    if (!(recordingflags & 8)) error = SGPrepare(seqGrab, false, true);
    if (error !=noErr) {
        DisposeGWorld(vidcapRecordBANK[slotid].gworld);
        vidcapRecordBANK[slotid].gworld = NULL;
        SGDisposeChannel(seqGrab, *sgchanptr);
        *sgchanptr = NULL;
        CloseComponent(seqGrab);
        PsychErrorExitMsg(PsychError_internal, "SGPrepare() for capture device failed!");            
    }

	// Setup of timebase for sequence grabber: This only works on OS/X, not
	// with Quicktime SDK for Windows, at least not with the brain-damaged MS-Compiler. Also,
	// SGPrepare() must have been called beforehand if audio clock is requested, because otherwise
	// audio timebase won't be available yet:

#if PSYCH_SYSTEM != PSYCH_WINDOWS
	
	// Retrieve current timebase of sequence grabber:
	if (noErr != SGGetTimeBase(seqGrab, &sgTimeBase)) {
		if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Video capture engine could not retrieve timebase. Capture timestamps and/or audio-video sync may be inaccurate!\n");
	}
	else {
		// Set error: This will trigger assignment of system clock below if
		// no sound recording requested:
		error = 1;
		
		// Sound recording requested?
		if (recordingflags & 2) {
			// Sound recording: We assign the masterclock of the sound timebase as master for sequence grabber.
			// This should give best possible audio-video sync:
			
			if ((noErr != (error = SGGetChannelTimeBase(vidcapRecordBANK[slotid].sgchanAudio, &soundTimeBase))) || (NULL == soundTimeBase)) {
				if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: (I) Video recording engine could not assign sound timebase as master timebase: Audio-video sync may be inaccurate! [QT-Error %i]\n", error);
			}
			else {
				SetTimeBaseMasterClock(sgTimeBase, (Component) GetTimeBaseMasterClock(soundTimeBase), NULL);
				error = GetMoviesError();
			}
		}
		
		// No sound recording requested (error == 1) or error during audio clock assignment?
		// If so, error != noErr and we assign the system clock as timebase:
		if (noErr != error) {
			// Only video recording/capture: We assign the system clock as masterclock for sequence grabber.
			// This provides capturetimestamps that are in sync with all other PTB clocks:
			looking.componentType = clockComponentType;
			looking.componentSubType = systemMicrosecondClock;
			looking.componentManufacturer = 0;
			looking.componentFlags = 0;
			looking.componentFlagsMask = 0;
			clockComponent = FindNextComponent(NULL, &looking);
			if (clockComponent == NULL) {
				if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Video recording engine could not assign system timebase as master timebase: Capture timestamps may be inaccurate!\n");
			}
			else {
				SetTimeBaseMasterClock(sgTimeBase, clockComponent, NULL);
				error = GetMoviesError();
			}
		}

		if ((noErr != error) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: Video recording engine could not assign new timebase: Capture timestamps may be inaccurate! [QT-Error %i]\n", error);
	}

#endif
    
    // Grabber should be ready now.
    
	// Create and initialize Mutex:
	PsychInitMutex(&(vidcapRecordBANK[slotid].mutex));

    // Assign new record:
    vidcapRecordBANK[slotid].seqGrab=seqGrab;    

    // Assign final handle:
    *capturehandle = slotid;

    // Increase counter:
    numCaptureRecords++;

    // Query capture framerate: MK This doesn't return meaningful results for
    // some reason :(
    error=SGGetFrameRate(vidcapRecordBANK[slotid].sgchanVideo, &framerate);
    if (error==noErr) {
      // Query worked: Assign it.
	   vidcapRecordBANK[slotid].fps = (double) FixedToFloat(framerate);
	 }
    else {
      // Query failed: Assign a dummy value of 25 Hz and output a warning.
      vidcapRecordBANK[slotid].fps = 25;
      if (PsychPrefStateGet_Verbosity()>1) printf("PTB-WARNING: Couldn't determine real capture framerate of grabber device %i. Assigning dummy value of 25 fps.\n", slotid);
    }

    // Determine size of images in movie:
    vidcapRecordBANK[slotid].width = movierect.right - movierect.left;
    vidcapRecordBANK[slotid].height = movierect.bottom - movierect.top;

	// Store pixeldepth in layers:
	vidcapRecordBANK[slotid].reqpixeldepth = reqdepth;

    // Reset framecounter:
    vidcapRecordBANK[slotid].nrframes = 0;

	// Store final recordingflags:
	vidcapRecordBANK[slotid].recordingflags = recordingflags;
	
    if (PsychPrefStateGet_Verbosity()>4) printf("PTB-INFO: Final assigned size of video input image is W x h = %i x  %i.\n", vidcapRecordBANK[slotid].width, vidcapRecordBANK[slotid].height);

    return(TRUE);
}

/*
 *  PsychQTCloseVideoCaptureDevice() -- Close a capture device and release all associated ressources.
 */
void PsychQTCloseVideoCaptureDevice(int capturehandle)
{
	OSErr err = noErr;

    if (capturehandle < 0 || capturehandle >= PSYCH_MAX_CAPTUREDEVICES) {
        PsychErrorExitMsg(PsychError_user, "Invalid capturehandle provided!");
    }
    
    if (vidcapRecordBANK[capturehandle].gworld == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid capturehandle provided. No capture device associated with this handle !!!");
    }

	// Clear grabber_active flag, also as a signal to potential worker thread to terminate:
	vidcapRecordBANK[capturehandle].grabber_active = 0;

	// Worker thread still active?
	if (vidcapRecordBANK[capturehandle].worker) {
		// Destroy it, wait for destruction, clean it up:
		PsychDeleteThread(&vidcapRecordBANK[capturehandle].worker);
	}

    // Stop capture immediately:
    SGStop(vidcapRecordBANK[capturehandle].seqGrab);
    
	// Destroy Mutex:
	PsychDestroyMutex(&(vidcapRecordBANK[capturehandle].mutex));

	// Release SG if previously prepared:
    if (!(vidcapRecordBANK[capturehandle].recordingflags & 8)) SGRelease(vidcapRecordBANK[capturehandle].seqGrab);

	// Release decompression sequence, if any:
	if (vidcapRecordBANK[capturehandle].decomSeq > 0) {
		err = CDSequenceEnd(vidcapRecordBANK[capturehandle].decomSeq);
		if ((noErr != err) && PsychPrefStateGet_Verbosity() > 1) {
			printf("PTB-WARNING: In shutdown of video capture device %i: Release of decompression sequence reports QT error %i. Ressource leakage??\n", capturehandle, (int) err);
		}
		vidcapRecordBANK[capturehandle].decomSeq = 0;
	}

	// Delete scaling matrix, if any:
	if (vidcapRecordBANK[capturehandle].scaleMatrixPtr) {
		free(vidcapRecordBANK[capturehandle].scaleMatrixPtr);
		vidcapRecordBANK[capturehandle].scaleMatrixPtr = NULL;
	}

    // Delete GWorld if any:
    if (vidcapRecordBANK[capturehandle].gworld) DisposeGWorld(vidcapRecordBANK[capturehandle].gworld);
    vidcapRecordBANK[capturehandle].gworld = NULL;

    if (vidcapRecordBANK[capturehandle].sgchanVideo) CloseComponent(vidcapRecordBANK[capturehandle].sgchanVideo);
	vidcapRecordBANK[capturehandle].sgchanVideo = NULL;

    if (vidcapRecordBANK[capturehandle].sgchanAudio) CloseComponent(vidcapRecordBANK[capturehandle].sgchanAudio);
	vidcapRecordBANK[capturehandle].sgchanAudio = NULL;
    
    // Release grabber:
    CloseComponent(vidcapRecordBANK[capturehandle].seqGrab);
    vidcapRecordBANK[capturehandle].seqGrab = NULL;
    vidcapRecordBANK[capturehandle].grabber_active = 0;
    
    // Decrease counter:
    if (numCaptureRecords>0) numCaptureRecords--;
        
    return;
}

/*
 *  PsychDeleteAllCaptureDevices() -- Delete all capture objects and release all associated ressources.
 */
void PsychQTDeleteAllCaptureDevices(void)
{
    int i;
    for (i=0; i<PSYCH_MAX_CAPTUREDEVICES; i++) {
        if (vidcapRecordBANK[i].gworld) PsychCloseVideoCaptureDevice(i);
    }
    return;
}


/*
 *  PsychQTGetTextureFromCapture() -- Create an OpenGL texturemap from a specific videoframe from given capture object.
 *
 *  win = Window pointer of onscreen window for which a OpenGL texture should be created.
 *  capturehandle = Handle to the capture object.
 *  checkForImage = >0 == Just check if new image available, 0 == really retrieve the image, blocking if necessary.
 *                   2 == Check for new image, block inside this function (if possible) if no image available.
 *  timeindex = This parameter is currently ignored and reserved for future use.
 *  out_texture = Pointer to the Psychtoolbox texture-record where the new texture should be stored.
 *  presentation_timestamp = A ptr to a double variable, where the presentation timestamp of the returned frame should be stored.
 *  summed_intensity = An optional ptr to a double variable. If non-NULL, then sum of intensities over all channels is calculated and returned.
 *  outrawbuffer = An optional ptr to a memory buffer of sufficient size. If non-NULL, the buffer will be filled with the captured raw image data, e.g., for use inside Matlab or whatever...
 *  Returns Number of pending or dropped frames after fetch on success (>=0), -1 if no new image available yet, -2 if no new image available and there won't be any in future.
 */
int PsychQTGetTextureFromCapture(PsychWindowRecordType *win, int capturehandle, int checkForImage, double timeindex, PsychWindowRecordType *out_texture, double *presentation_timestamp, double* summed_intensity, rawcapimgdata* outrawbuffer)
{
    OSErr		error = noErr;
    GLuint texid;
    int w, h, padding;
    double targetdelta, realdelta, frames;
    unsigned int intensity = 0;
    unsigned int count, i;
    unsigned char* pixptr;
	unsigned char* outpixptr;
	unsigned char swapbyte;
    psych_bool newframe = FALSE;
    double tstart, tend;
    unsigned int pixval, alphacount;
    int nrdropped;
	int reqdepth;
	
    PsychGetAdjustedPrecisionTimerSeconds(&tstart);

    // Sanity checks:
    if (capturehandle < 0 || capturehandle >= PSYCH_MAX_CAPTUREDEVICES || vidcapRecordBANK[capturehandle].gworld == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid capturehandle provided.");
    }

	// Compute width, height and padding for later creation of textures etc. Need to do this here,
	// so we can return the values for raw data retrieval:
    w=vidcapRecordBANK[capturehandle].width;
    h=vidcapRecordBANK[capturehandle].height;
	reqdepth = vidcapRecordBANK[capturehandle].reqpixeldepth;

    // Hack: Need to extend rect by 4 pixels, because GWorlds are 4 pixels-aligned via
    // image row padding:
#if PSYCH_SYSTEM == PSYCH_OSX
	// RGBA format is 4-Byte aligned...
    padding = 4 + (4 - (w % 4)) % 4;
	// Special case for luminance pixel formats: 16-Byte aligned...
	if (reqdepth < 3) padding = 16 + (16 - (w % 16)) % 16;
#else
    padding= 0;
#endif

	// If a outrawbuffer struct is provided, we fill it with info needed to allocate a
	// sufficient memory buffer for returned raw image data later on:
	if (outrawbuffer) {
		outrawbuffer->w = w;
		outrawbuffer->h = h;
		outrawbuffer->depth = reqdepth;
	}
    
	// Lock device mutex so we have exclusive access:
	PsychLockMutex(&(vidcapRecordBANK[capturehandle].mutex));

	// Actual capture ops driven by parallel background worker thread, so nothing to do for us?
	if (NULL == vidcapRecordBANK[capturehandle].worker) {
		// No, we need to do it: Grant some processing time to the sequence grabber engine:
		if ((error=SGIdle(vidcapRecordBANK[capturehandle].seqGrab))!=noErr) {
			// We don't abort on non noErr case, but only (optionally) report it. This is because when in harddisc
			// movie recording mode with sound recording, we can get intermittent errors in SGIdle() which are
			// meaningless so they are best ignored and must not lead to interruption of PTB operation.
			if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: In PsychGetTextureFromCapture(): SGIdle() returns error code %i\n", (int) error);
			
			// Unlock mutex:
			PsychUnlockMutex(&(vidcapRecordBANK[capturehandle].mutex));
			
			// Return abort error code:
			return(-2);
		}
	}

	// Ist this device in "recording only" mode? If so then it doesn't have a
	// callback proc assigned, therefore all the processing below this line would
	// just cause a deadlock. The only allowed mode is checkForImage == 4 in this
	// case:
	if ((vidcapRecordBANK[capturehandle].recordingflags & 4) && (checkForImage!=4)) {
		// Invalid mode of invocation!

		// Unlock mutex:
		PsychUnlockMutex(&(vidcapRecordBANK[capturehandle].mutex));

		PsychErrorExitMsg(PsychError_user, "Capturedevice was opened in ''disk recording only'' mode: You must specify a 'waitForImage' flag of 4 in your Screen('GetCapturedImage') call!");
	}
	
    // Check if a new captured frame is ready for retrieval...
    newframe = (psych_bool) vidcapRecordBANK[capturehandle].frame_ready;
    // ...and clear out the ready flag immediately:
    vidcapRecordBANK[capturehandle].frame_ready = 0;

    // Presentation timestamp requested?
    if (presentation_timestamp) {
        // Already available? Return it:
        *presentation_timestamp =vidcapRecordBANK[capturehandle].current_pts;
    }
    
    // Should we just check for new image? If so, just return availability status:
    if (checkForImage) {
		// Unlock mutex: No further access to shared variables until return():
		PsychUnlockMutex(&(vidcapRecordBANK[capturehandle].mutex));

        if (vidcapRecordBANK[capturehandle].grabber_active == 0) {
            // Grabber stopped. We'll never get a new image:
            return(-2);
        }

        // Timestamping:
        PsychGetAdjustedPrecisionTimerSeconds(&tend);
        vidcapRecordBANK[capturehandle].nrgfxframes++;
        vidcapRecordBANK[capturehandle].avg_gfxtime+=(tend - tstart);

        // Grabber active. Just return availability status: -1 = none avail. yet, 0 = At least one avail.
        return((newframe) ? 0 : -1);
    }
    
    // This point is only reached if checkForImage == FALSE, which only happens
    // if a new frame is available in our GWorld:

    // Synchronous texture fetch code for GWorld rendering mode:
    
    // Build a standard PTB texture record:    
    
    // Assign texture rectangle:
    if (out_texture) {
        // Oldstyle: No longer needed due to contraction code below: PsychMakeRect(out_texture->rect, 0, 0, w+padding, h);    
        PsychMakeRect(out_texture->rect, 0, 0, w, h);    
        
        // Set NULL - special texture object as part of the PTB texture record:
        out_texture->targetSpecific.QuickTimeGLTexture = NULL;
        
        // Set texture orientation as if it were an inverted Offscreen window: Upside-down.
        out_texture->textureOrientation = 3;
        
        // Setup a pointer to our GWorld as texture data pointer: Settin memsize to zero
        // prevents unwanted free() operation in PsychDeleteTexture...
        out_texture->textureMemorySizeBytes = 0;
    }
    
    // Lock GWorld:
    if(!LockPixels(GetGWorldPixMap(vidcapRecordBANK[capturehandle].gworld))) {
        // Locking surface failed! We abort.
		PsychUnlockMutex(&(vidcapRecordBANK[capturehandle].mutex));
        PsychErrorExitMsg(PsychError_internal, "PsychGetTextureFromCapture(): Locking GWorld pixmap surface failed!!!");
    }

	// Preprocessing requested?
	if (out_texture || summed_intensity || outrawbuffer) {
		// Manual negating required?	
		if (reqdepth < 3) {
			// Luminance formats (1 Byte Luminance) are inverted for weird reasons. Need to "negate" all
			// pixels before any further processing can occur. Also, we need to "contract" the image to
			// remove the 16-Byte padding:
			pixptr = (unsigned char*) GetPixBaseAddr(GetGWorldPixMap(vidcapRecordBANK[capturehandle].gworld));
			outpixptr = pixptr;
			count = ((w+padding) * h);
			for (i=0; i<count; i++) {
				if ((i % (w+padding)) < w) {
					*(outpixptr++) = 255 - *(pixptr++);
				}
				else {
					pixptr++;
				}
			}
		}
		else {
			// Need to "contract" the image to get rid of padding:
			pixptr = (unsigned char*) GetPixBaseAddr(GetGWorldPixMap(vidcapRecordBANK[capturehandle].gworld));
			outpixptr = pixptr;
			for (i=0; i<h; i++) {
				for (count=0; count < w*4; count++) *(outpixptr++) = *(pixptr++);
				for (count=w*4; count < (w+padding)*4; count++) pixptr++;
			}
			
			// At requested depth 3 (RGB), need to contract again from RGBA -> RGB
			if (reqdepth == 3) {
				pixptr = (unsigned char*) GetPixBaseAddr(GetGWorldPixMap(vidcapRecordBANK[capturehandle].gworld));
				outpixptr = pixptr;
				count = w * h * 4;
				for (i=0; i<count; i++) {
					// Each 4th byte needs to be discarded, because it is the "dead" alpha byte:
					if (i % 4 == 3) {
						// Only advance read-ptr to get rid of alpha byte:
						pixptr++;
						// This is ugly: The byte order that QT delivers is not what OpenGL needs for RGB textures,
						// ie., the order of QT is BGR,instead of required RGB. Need to swap R<->B. Is this ugly, or what?!?
						// TODO: Works on Little Endian Intel systems, but is it needed o Big Endian PowerPC as well?
						swapbyte = outpixptr[-3];
						outpixptr[-3] = outpixptr[-1];
						outpixptr[-1] = swapbyte;
					}
					else {
						// Standard copy of a R, G or B byte:
						*(outpixptr++) = *(pixptr++);
					}
				}
			}
		}
	}
	
    if (out_texture) {
        // This will retrieve an OpenGL compatible pointer to the GWorlds pixel data and assign it to our texmemptr:
        out_texture->textureMemory = (GLuint*) GetPixBaseAddr(GetGWorldPixMap(vidcapRecordBANK[capturehandle].gworld));

		// Activate OpenGL context of target window: We'll need it for texture fetch...
		PsychSetGLContext(win);
		
		// Disable client storage on OS-X, if it is enabled for some reason:
		#if PSYCH_SYSTEM == PSYCH_OSX
		glPixelStorei(GL_UNPACK_CLIENT_STORAGE_APPLE, GL_FALSE);
		#endif

		out_texture->depth = 8 * reqdepth;
		out_texture->nrchannels = reqdepth;

        // Let PsychCreateTexture() do the rest of the job of creating, setting up and
        // filling an OpenGL texture with GWorlds content:
		//
		// CAUTION: An error abort inside this function could leave us with a locked mutex!
		// However, error aborts inside this are not likely, and users will hopefully "clear all"
		// and thereby resolve dangling mutex issues after such a fatal error:
        PsychCreateTexture(out_texture);

        // Undo hack from above after texture creation: Now we need the real width of the
        // texture for proper texture coordinate assignments in drawing code et al.
        PsychMakeRect(out_texture->rect, 0, 0, w, h);    
        // Ready to use the texture...
    }
    
    // Sum of pixel intensities requested?
    if(summed_intensity) {
        pixptr = (unsigned char*) GetPixBaseAddr(GetGWorldPixMap(vidcapRecordBANK[capturehandle].gworld));
        count = (w * h * reqdepth);
        alphacount = 0;
        for (i=0; i<count; i++) {
            pixval=(unsigned int) pixptr[i];
            // Is this the alpha-channel component which is fixed to 255? If so, count it.
            if (pixval==255) alphacount++;
            intensity+=pixval;
        }
		
        // Try to discount the w*h*255 alpha channel values, if alpha channel is fixed to 255:
        // Some video digitizers set alpha component correctly to 255, some leave it at the
        // wrong value of zero :(
        if ((reqdepth > 1) && (alphacount >= w*h)) intensity = intensity - (w * h * 255);
        *summed_intensity = ((double) intensity) / w / h / ((reqdepth > 1) ? 3 : 1);
    }

	// Raw data requested?
	if (outrawbuffer) {
		// Copy it out:
		outrawbuffer->w = w;
		outrawbuffer->h = h;
		outrawbuffer->depth = reqdepth;
		count = (outrawbuffer->w * outrawbuffer->h * outrawbuffer->depth);
		memcpy(outrawbuffer->data, (const void *) GetPixBaseAddr(GetGWorldPixMap(vidcapRecordBANK[capturehandle].gworld)), count);
	}

    // Unlock GWorld surface.
    UnlockPixels(GetGWorldPixMap(vidcapRecordBANK[capturehandle].gworld));
    
    // Detection of dropped frames: This is a heuristic. We'll see how well it works out...
    
    // Expected delta between successive presentation timestamps:
    targetdelta = 1.0f / vidcapRecordBANK[capturehandle].fps;
    
    // Compute real delta, given rate and playback direction:
    realdelta = *presentation_timestamp - vidcapRecordBANK[capturehandle].last_pts;
    if (realdelta<0) realdelta = 0;
    frames = realdelta / targetdelta;

    // Dropped frames?
    if (frames > 1 && vidcapRecordBANK[capturehandle].last_pts>=0) {
      nrdropped = (int) (frames - 1 + 0.5);
      vidcapRecordBANK[capturehandle].nr_droppedframes += nrdropped;
    }
    else {
      nrdropped = 0;
    }

    // Record timestamp as reference for next check:    
    vidcapRecordBANK[capturehandle].last_pts = *presentation_timestamp;

	// Unlock mutex:
	PsychUnlockMutex(&(vidcapRecordBANK[capturehandle].mutex));
    
    // Timestamping:
    PsychGetAdjustedPrecisionTimerSeconds(&tend);
    vidcapRecordBANK[capturehandle].nrgfxframes++;
    vidcapRecordBANK[capturehandle].avg_gfxtime+=(tend - tstart);
    
    // We're successfully done. Return number of frames that had to be dropped, or are
    // pending in internal ringbuffers:
    return(nrdropped);
}

/*
 *  PsychQTVideoCaptureRate() - Start- and stop video capture.
 *
 *  capturehandle = Grabber to start-/stop.
 *  playbackrate = zero == Stop capture, non-zero == Capture
 *  dropframes = Currently ignored.
 *  startattime = Deadline (in system time) to wait for before real start of capture.
 *  Returns Number of dropped frames during capture.
 */
int PsychQTVideoCaptureRate(int capturehandle, double capturerate, int dropframes, double* startattime)
{
    int dropped = 0;
    OSErr error = noErr;
    Fixed framerate;
	long usage;
	Fixed maxframerate;
	long milliSecPerFrame, bps;
	int retrycount = 0;
	unsigned long storageRem;
	
    if (capturehandle < 0 || capturehandle >= PSYCH_MAX_CAPTUREDEVICES) {
        PsychErrorExitMsg(PsychError_user, "Invalid capturehandle provided!");
    }
        
    // Fetch references to objects we need:
    if (vidcapRecordBANK[capturehandle].seqGrab == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid capturehandle provided. No movie associated with this handle !!!");
    }
    
    if (capturerate > 0) {
        // Start capture:
        if (vidcapRecordBANK[capturehandle].grabber_active) PsychErrorExitMsg(PsychError_user, "You tried to start video capture, but capture is already started!");

	// Low latency capture disabled?
	if (dropframes == 0) {
		// Yes. Need to clear the lowlat flag from our channel config:
		if (!(vidcapRecordBANK[capturehandle].recordingflags & 8)) error = SGRelease(vidcapRecordBANK[capturehandle].seqGrab);
		if ((noErr != error) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: In highlat: SGRelease() reports error %i in start function.\n", (int) error);
		
		error = SGGetChannelUsage(vidcapRecordBANK[capturehandle].sgchanVideo, &usage);
		usage&=~seqGrabLowLatencyCapture;
		error = SGSetChannelUsage(vidcapRecordBANK[capturehandle].sgchanVideo, usage);
		if ((noErr != error) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: In highlat: Video SGSetChannelUsage() reports error %i in start function.\n", (int) error);
		
		if (vidcapRecordBANK[capturehandle].sgchanAudio) {
			error = SGGetChannelUsage(vidcapRecordBANK[capturehandle].sgchanAudio, &usage);
			usage&=~seqGrabLowLatencyCapture;
			error = SGSetChannelUsage(vidcapRecordBANK[capturehandle].sgchanAudio, usage);
			if ((noErr != error) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: In highlat: Audio SGSetChannelUsage() reports error %i in start function.\n", (int) error);
		}
		
		if (!(vidcapRecordBANK[capturehandle].recordingflags & 8)) error = SGPrepare(vidcapRecordBANK[capturehandle].seqGrab, false, true);
		if ((noErr != error) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: In highlat: SGPrepare() reports error %i in start function.\n", (int) error);
	}
	else {
		// dropframes non-null. lowlat flag set?
		error = SGGetChannelUsage(vidcapRecordBANK[capturehandle].sgchanVideo, &usage);
		if (!(usage & seqGrabLowLatencyCapture)) {
			// Lowlatency requested, but flag not set. Need to set it:			
			if (!(vidcapRecordBANK[capturehandle].recordingflags & 8)) error = SGRelease(vidcapRecordBANK[capturehandle].seqGrab);
			if ((noErr != error) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: In lowlat: SGRelease() reports error %i in start function.\n", (int) error);
			
			error = SGGetChannelUsage(vidcapRecordBANK[capturehandle].sgchanVideo, &usage);
			usage|= seqGrabLowLatencyCapture;
			error = SGSetChannelUsage(vidcapRecordBANK[capturehandle].sgchanVideo, usage);
			if ((noErr != error) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: In lowlat: Video SGSetChannelUsage() reports error %i in start function.\n", (int) error);
			
			if (vidcapRecordBANK[capturehandle].sgchanAudio) {
				error = SGGetChannelUsage(vidcapRecordBANK[capturehandle].sgchanAudio, &usage);
				usage|= seqGrabLowLatencyCapture;
				error = SGSetChannelUsage(vidcapRecordBANK[capturehandle].sgchanAudio, usage);
				if ((noErr != error) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: In lowlat: Audio SGSetChannelUsage() reports error %i in start function.\n", (int) error);
			}
			
			if (!(vidcapRecordBANK[capturehandle].recordingflags & 8)) error = SGPrepare(vidcapRecordBANK[capturehandle].seqGrab, false, true);
			if ((noErr != error) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: In lowlat: SGPrepare() reports error %i in start function.\n", (int) error);
		}
	}

	framerate = FloatToFixed((float) ((capturerate == DBL_MAX) ? 0 : capturerate));
	error = SGSetFrameRate(vidcapRecordBANK[capturehandle].sgchanVideo, framerate);
	if ((noErr != error) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: SGSetFrameRate() reports error %i in start function.\n", (int) error);

	// Wait until start deadline reached:
	if (*startattime != 0) PsychWaitUntilSeconds(*startattime);

	// Start the engine!
	vidcapRecordBANK[capturehandle].grabber_active = 1;
	error = SGStartRecord(vidcapRecordBANK[capturehandle].seqGrab);
	if ((noErr != error) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: SGStartRecord() reports error %i in start function.\n", (int) error);

	// Record real start time:
	PsychGetAdjustedPrecisionTimerSeconds(startattime);

        vidcapRecordBANK[capturehandle].last_pts = -1.0;
        vidcapRecordBANK[capturehandle].nr_droppedframes = 0;
        vidcapRecordBANK[capturehandle].frame_ready = 0;
		framerate = (Fixed) 0;
        error = SGGetFrameRate(vidcapRecordBANK[capturehandle].sgchanVideo, &framerate);
		if ((noErr != error) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: SGGetFrameRate() reports error %i in start function.\n", (int) error);

        vidcapRecordBANK[capturehandle].fps = (double) FixedToFloat(framerate);
		VDGetDataRate(vidcapRecordBANK[capturehandle].vdig, &milliSecPerFrame, &maxframerate, &bps);
		
        if (PsychPrefStateGet_Verbosity()>3) printf("PTB-INFO: Capture framerate is reported as: %lf  [max = %f]\n", vidcapRecordBANK[capturehandle].fps, FixedToFloat(maxframerate));
        if (PsychPrefStateGet_Verbosity()>3) printf("PTB-INFO: Capture latency is reported as: %i msecs\n", milliSecPerFrame);
		if ((FixedToFloat(maxframerate) > 0) && ((double) FixedToFloat(maxframerate) < vidcapRecordBANK[capturehandle].fps)) {
			// Maxframerate valid and lower than reported framerate?!? Clamp to maximum:
			vidcapRecordBANK[capturehandle].fps = (double) FixedToFloat(maxframerate);
			if (PsychPrefStateGet_Verbosity()>3) printf("PTB-INFO: Maximum framerate overrides reported one! Will assume a capture framerate of: %lf fps for further operation.\n", vidcapRecordBANK[capturehandle].fps);
		}
		
		// Start worker thread for video capture, if requested:
		if (vidcapRecordBANK[capturehandle].recordingflags & 16) {
			if ((error = (OSErr) PsychCreateThread(&vidcapRecordBANK[capturehandle].worker, NULL, PsychQTVideoCaptureThreadMain, (void*) &vidcapRecordBANK[capturehandle]))) {
				if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: In start of videocapture for device %i: Failed to start worker thread [%s]!\n", capturehandle, strerror((int) error));
			}
			else {
				if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Videocapture for device %i uses a parallel background capture thread for operation as requested.\n", capturehandle);
			}
		}
    }
    else {
		if (PsychPrefStateGet_Verbosity() > 3) {
			SGGetStorageSpaceRemaining(vidcapRecordBANK[capturehandle].seqGrab, &storageRem);
			printf("PTB-INFO: At stop time, device %i reports %f MB of storage space remaining.\n", capturehandle, ((float) storageRem) / 1024 / 1024);
		}

		// Clear grabber_active flag, also as a signal to potential worker thread to terminate:
        vidcapRecordBANK[capturehandle].grabber_active = 0;

		// Worker thread active?
		if (vidcapRecordBANK[capturehandle].worker) {
			// Destroy it, wait for destruction, clean it up:
			PsychDeleteThread(&vidcapRecordBANK[capturehandle].worker);
		}

        // Stop capture:
		retrycount = 0;
		while ( (retrycount < 4) && ((error = SGStop(vidcapRecordBANK[capturehandle].seqGrab)) != noErr) ) {
			retrycount++;
			if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: %i st. invocation of SGStop() reports error %i in stop function.\n", retrycount, (int) error);
		}

		// Reset frame_ready:
        vidcapRecordBANK[capturehandle].frame_ready = 0;

        // Output count of dropped frames:
        if ((dropped=vidcapRecordBANK[capturehandle].nr_droppedframes) > 0) {
            if (PsychPrefStateGet_Verbosity()>3) printf("PTB-INFO: Video capture dropped %i frames on device %i to keep pipe running.\n", vidcapRecordBANK[capturehandle].nr_droppedframes, capturehandle); 
        }
        if (vidcapRecordBANK[capturehandle].nrframes>0)  vidcapRecordBANK[capturehandle].avg_decompresstime/= (double)vidcapRecordBANK[capturehandle].nrframes;
        if (PsychPrefStateGet_Verbosity()>3) printf("PTB-INFO: Average time spent in video decompressor was %lf milliseconds.\n", vidcapRecordBANK[capturehandle].avg_decompresstime * 1000.0f);
        if (vidcapRecordBANK[capturehandle].nrgfxframes>0)  vidcapRecordBANK[capturehandle].avg_gfxtime/= (double)vidcapRecordBANK[capturehandle].nrgfxframes;
        if (PsychPrefStateGet_Verbosity()>3) printf("PTB-INFO: Average time spent in GetCapturedImage was %lf milliseconds.\n", vidcapRecordBANK[capturehandle].avg_gfxtime * 1000.0f);
    }
    
    // Reset framecounter:
    vidcapRecordBANK[capturehandle].nrframes = 0;
    vidcapRecordBANK[capturehandle].avg_decompresstime = 0;
    vidcapRecordBANK[capturehandle].nrgfxframes = 0;
    vidcapRecordBANK[capturehandle].avg_gfxtime = 0;
    
    // Return either real capture framerate (at start of capture) or count of dropped frames - at end of capture.
    return((capturerate!=0) ? (int) (vidcapRecordBANK[capturehandle].fps + 0.5) : dropped);
}

/* Set capture device specific parameters:
 * On OS-X and Windows (and therefore in this implementation) this is currently a no-op, until
 * we find out how to do this with the Sequence-Grabber API.
 */
double PsychQTVideoCaptureSetParameter(int capturehandle, const char* pname, double value)
{
	PsychRectType roirect;
	QTAtomContainer         atomContainer;
    QTAtom                  featureAtom;
    VDIIDCFeatureSettings   settings;
    VideoDigitizerComponent vd;
    ComponentDescription    desc;
    ComponentResult         result = noErr;
	FSSpec					recfile;
	
	// Valid handle provided?
	if (capturehandle < 0 || capturehandle >= PSYCH_MAX_CAPTUREDEVICES) {
		PsychErrorExitMsg(PsychError_user, "Invalid capturehandle provided!");
	}
	
	// Fetch references to objects we need:
    if (vidcapRecordBANK[capturehandle].seqGrab == NULL) {
        PsychErrorExitMsg(PsychError_user, "Invalid capturehandle provided. No capture device associated with this handle !!!");
    }

	// Set a new target movie name for video recordings:
	if (strstr(pname, "SetNewMoviename=")) {
		// Find start of movie namestring and assign to pname:
		pname = strstr(pname, "=");
		pname++;
		
		// Convert to QT format:
		NativePathNameToFSSpec(pname, &recfile, 0);

		// If recordingflags & 1, then we request capture to memory with writeout at end of capture operation. Otherwise
		// we request immediate capture to disk. We always append to an existing movie file, instead of overwriting it.
		result = SGSetDataOutput(vidcapRecordBANK[capturehandle].seqGrab, &recfile, ((vidcapRecordBANK[capturehandle].recordingflags & 1) ? seqGrabToMemory : seqGrabToDisk));
		if (result != noErr) {
			printf("PTB-ERROR: 'SetVideoCaptureParameter for device %i: While trying to set new output movie name %s. Quicktime returned error code %i! Failed!\n", capturehandle, pname, (int) result);
			PsychErrorExitMsg(PsychError_user, "Selecting a new output movie file for video recording failed!");
		}
		
		if (!(vidcapRecordBANK[capturehandle].recordingflags & 8)) result = SGPrepare(vidcapRecordBANK[capturehandle].seqGrab, false, true);
		if (result !=noErr) {
			if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: 'SetVideoCaptureParameter for device %i: While trying to set new output movie name %s. SGPrepare() for capture device failed with error code %i! Failed!\n", capturehandle, pname, (int) result);

			// Try to get over it and hope for the best...
			return(0);
		}

		// Return success:
		if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: 'SetVideoCaptureParameter for device %i: Set new output movie filename to '%s'.\n", capturehandle, pname);
		return(0);
	}
	
	// Return current framerate:
	if (strcmp(pname, "GetFramerate")==0) {
		PsychCopyOutDoubleArg(1, FALSE, vidcapRecordBANK[capturehandle].fps);
		return(0);
	}
	
	// Return current ROI of camera, as requested (and potentially modified during
	// PsychOpenCaptureDevice(). This is a read-only parameter, as the ROI can
	// only be set during Screen('OpenVideoCapture').
	if (strcmp(pname, "GetROI")==0) {
		PsychMakeRect(roirect, vidcapRecordBANK[capturehandle].roirect.left, vidcapRecordBANK[capturehandle].roirect.top,
					  vidcapRecordBANK[capturehandle].roirect.right, vidcapRecordBANK[capturehandle].roirect.bottom); 
		PsychCopyOutRectArg(1, FALSE, roirect);
		return(0);
	}
	
	// Return vendor name string:
	if (strcmp(pname, "GetVendorname")==0) {
		PsychCopyOutCharArg(1, FALSE, "Unknown Vendor.");
		return(0);
	}
	
	// Return model name string:
	if (strcmp(pname, "GetModelname")==0) {
		PsychCopyOutCharArg(1, FALSE, "Unknown Model.");
		return(0);
	}
	
    if (vidcapRecordBANK[capturehandle].sgchanVideo == NULL) {
        PsychErrorExitMsg(PsychError_user, "Sorry, provided capture device doesn't have a video digitizer attached!");
    }
	
    // Get the video digitizer and make sure it's legit
    vd = SGGetVideoDigitizerComponent(vidcapRecordBANK[capturehandle].sgchanVideo);
    if (vd == NULL) {
        PsychErrorExitMsg(PsychError_user, "Sorry, unable to query video digitizer component!");
    }

	// Check if it is IIDC compliant. We only support IIDC compliant devices (usually Firewire):
    GetComponentInfo((Component)vd, &desc, NULL, NULL, NULL);
    if (vdSubtypeIIDC != desc.componentSubType) {
		if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Tried to query or set video capture parameters on a non-supported (non IIDC) video device. Skipped.\n");
		goto bail;
	}
			printf("XXX\n");
/*	unsigned short             brightness = 65535;
	result = VDSetSaturation(vd, &brightness);
	if (result==noErr) {
		printf("Newval %i\n", (int) brightness);
	}
	else {
		printf("failed. %i\n", (int) brightness);
	}
*/	
	
    // *** now do the real work ***
	// vdIIDCFeatureGain
    // return the gain feature in an atom container
    result = VDIIDCGetFeaturesForSpecifier(vd, vdIIDCFeatureSaturation, &atomContainer);
    if (noErr == result) {
		printf("AAA\n");
        // Find the feature atom:
        featureAtom = QTFindChildByIndex(atomContainer, kParentAtomIsContainer,
                                         vdIIDCAtomTypeFeature, 1, NULL);
        if (0 == featureAtom) { result = cannotFindAtomErr; goto bail; }
		printf("BBB\n");
		
        // find the gain settings from the feature atom and copy the data
        // into our settings
        result = QTCopyAtomDataToPtr(atomContainer,
                                     QTFindChildByID(atomContainer, featureAtom,
													 vdIIDCAtomTypeFeatureSettings,
													 vdIIDCAtomIDFeatureSettings, NULL),
                                     true, sizeof(settings), &settings, NULL);
        if (noErr == result) {
		printf("CCC\n");

            /* When indicating capabilities, the flag being set indicates that the
			feature can be put into the given state.
			When indicating/setting state, the flag represents the current/desired
			state. Note that certain combinations of flags are valid for capabilities
			(i.e. vdIIDCFeatureFlagOn | vdIIDCFeatureFlagOff) but are mutually
			exclusive for state.
			*/
            // is the setting supported?
            if (settings.capabilities.flags & (vdIIDCFeatureFlagOn |
                                               vdIIDCFeatureFlagManual |
                                               vdIIDCFeatureFlagRawControl)) {
                // set state flags
                settings.state.flags = (vdIIDCFeatureFlagOn |
                                        vdIIDCFeatureFlagManual |
                                        vdIIDCFeatureFlagRawControl);
						printf("DDD\n");

                // set value - will either be 500 or the max value supported by
                // the camera represented in a float between 0 and 1.0
                settings.state.value = (1.0 / settings.capabilities.rawMaximum) *
					((settings.capabilities.rawMaximum > 500) ? 500 :
					 settings.capabilities.rawMaximum);
				
                // store the result back in the container
                result = QTSetAtomData(atomContainer,
                                       QTFindChildByID(atomContainer, featureAtom,
													   vdIIDCAtomTypeFeatureSettings,
													   vdIIDCAtomIDFeatureSettings, NULL),
                                       sizeof(settings), &settings);
                if (noErr == result) {
						printf("EEE\n");

                    // set it on the device
                    result = VDIIDCSetFeatures(vd, atomContainer);
                }
            } else {
                // can't do it!
                result = featureUnsupported;
            }
        }
    }
	
// We jump here if something's going wrong with the sequence grabber code:
bail:

	// Unsupported parameter:
	if (result!=noErr && PsychPrefStateGet_Verbosity() > 1) printf("PTB-INFO: Video capture parameter %s not supported on video capture device %i. Skipped.\n", pname, capturehandle);

	// Just return the "not supported" value DBL_MAX:
	return(DBL_MAX);
}

/*
 *  void PsychQTExitVideoCapture() - Shutdown handler.
 *
 *  This routine is called by Screen('CloseAll') and on clear Screen time to
 *  do final cleanup. It deletes all capture objects
 *
 */
void PsychQTExitVideoCapture(void)
{
    // Release all capture devices
    PsychQTDeleteAllCaptureDevices();
    
    // Shutdown Quicktime toolbox: We skip this, because according to Apple its not necessary,
    // and for some reason it reliably hangs Matlab, so one has to force-quit it :-(
    // Don't do this: ExitMovies();
#if PSYCH_SYSTEM == PSYCH_WINDOWS
    // Do not do it even on Windows. It can cause crashes...
    // Shutdown Quicktime core system:
    // ExitMovies();
    
    // Shutdown Quicktime for Windows compatibility layer:
    // TerminateQTML();
#endif

    firsttime = TRUE;
    return;
}

#endif
