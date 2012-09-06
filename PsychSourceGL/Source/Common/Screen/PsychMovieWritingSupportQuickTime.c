/*
	Common/Screen/PsychMovieWritingSupportQuickTime.c
	
	PLATFORMS:	
	
		MacOS/X and MS-Windows, but only if GStreamer support is disabled.
		
		This is the movie editing and writing/creation engine based on the
		Apple QuickTime API. It works on Apple MacOS/X and MS-Windows.

	AUTHORS:
	
		Mario Kleiner           mk              mario.kleiner@tuebingen.mpg.de

	HISTORY:
	
	04/18/10		mk		Wrote it. 
	
	DESCRIPTION:
	
		Psychtoolbox functions for dealing with Quicktime movie editing. This implementation uses
		Quicktimes API's, therefore it works on Operating systems with Quicktime support,
		currently OS/X and Windows.

	NOTES:

*/

#include "Screen.h"

// Surrogates to prevent linker failure if neither GStreamer nor Quicktime available:
#if !defined(PTB_USE_GSTREAMER) && !defined(PSYCHQTAVAIL)
void PsychMovieWritingInit(void) { return; }
void PsychExitMovieWriting(void) { return; }
void PsychDeleteAllMovieWriters(void) { return; }
int PsychCreateNewMovieFile(char* moviefile, int width, int height, double framerate, char* movieoptions)
{
    PsychErrorExitMsg(PsychError_unimplemented, "Sorry, movie writing not supported on this operating system");
    return(0);
}

int PsychFinalizeNewMovieFile(int movieHandle) {
    PsychErrorExitMsg(PsychError_unimplemented, "Sorry, movie writing not supported on this operating system");
    return FALSE;
}
int PsychAddVideoFrameToMovie(int moviehandle, int frameDurationUnits, psych_bool isUpsideDown)
{
    PsychErrorExitMsg(PsychError_unimplemented, "Sorry, movie writing not supported on this operating system");
    return(1);
}

unsigned char*	PsychGetVideoFrameForMoviePtr(int moviehandle, unsigned int* twidth, unsigned int* theight)
{
    PsychErrorExitMsg(PsychError_unimplemented, "Sorry, movie writing not supported on this operating system");
    return(NULL);
}

psych_bool PsychAddAudioBufferToMovie(int moviehandle, unsigned int nrChannels, unsigned int nrSamples, double* buffer)
{
    PsychErrorExitMsg(PsychError_unimplemented, "Sorry, movie writing not supported on this operating system");
    return FALSE;
}

// End of surrogate routines.
#endif

// No Quicktime support for GNU/Linux:
#if PSYCH_SYSTEM != PSYCH_LINUX

// Shall we use the QuickTime implementation? It takes precedence over
// GStreamer on the only system config which still supports QT == 32-Bit OSX:
#if defined(PSYCHQTAVAIL)

#if PSYCH_SYSTEM == PSYCH_OSX
#include <Quicktime/QuickTimeComponents.h>
#endif

#if PSYCH_SYSTEM == PSYCH_WINDOWS
#include <QTML.h>
#include <QuickTimeComponents.h>
#endif

// Forward declaration of internal helper function:


// Record which defines all state for a capture device:
typedef struct {
	Movie						Movie;
	Track						Track;
	Media						Media;
	FSSpec						File;
	short						ResRefNum;
	short						ResID;
	GWorldPtr					GWorld;
	PixMapHandle				PixMap;
	CodecType					CodecType;
	CodecQ						CodecQuality;
	long						MaxComprSize;
	Handle						ComprDataHdl;
	Ptr							ComprDataPtr;
	ImageDescriptionHandle		ImageDesc;
	CGrafPtr 					SavedPort;
	GDHandle					SavedDevice;
	Rect						Rect;
	int							height;
	int							width;
	int							padding;
} PsychMovieWriterRecordType;

static PsychMovieWriterRecordType moviewriterRecordBANK[PSYCH_MAX_MOVIEWRITERDEVICES];
static int moviewritercount = 0;
static psych_bool firsttime = TRUE;

void PsychMovieWritingInit(void)
{
	int i;
	
	for (i = 0; i < PSYCH_MAX_MOVIEWRITERDEVICES; i++) {
		memset(&(moviewriterRecordBANK[i]), 0, sizeof(PsychMovieWriterRecordType));
	}
	
	moviewritercount = 0;
	
	return;
}

void PsychDeleteAllMovieWriters(void)
{
	int i;
		
	for (i = 0; i < PSYCH_MAX_MOVIEWRITERDEVICES; i++) {
		if (moviewriterRecordBANK[i].Movie) PsychFinalizeNewMovieFile(i);
	}
}

void PsychExitMovieWriting(void)
{
	PsychDeleteAllMovieWriters();
    firsttime = TRUE;
	return;
}

PsychMovieWriterRecordType* PsychGetMovieWriter(int moviehandle, psych_bool unsafe)
{
	if (moviehandle < 0 || moviehandle >= PSYCH_MAX_MOVIEWRITERDEVICES) PsychErrorExitMsg(PsychError_user, "Invalid handle for moviewriter provided!");
	if (!unsafe && (NULL == moviewriterRecordBANK[moviehandle].Movie)) PsychErrorExitMsg(PsychError_user, "Invalid handle for moviewriter provided! No such writer open.");
	return(&(moviewriterRecordBANK[moviehandle]));
}

unsigned char*	PsychGetVideoFrameForMoviePtr(int moviehandle, unsigned int* twidth, unsigned int* theight)
{
	PsychMovieWriterRecordType* pwriterRec = PsychGetMovieWriter(moviehandle, FALSE);
	if (NULL == pwriterRec->PixMap) return(NULL);

	*twidth = pwriterRec->width;
	*theight = pwriterRec->height - pwriterRec->padding; // Hack hack pwriterRec->padding !
	return((unsigned char*) GetPixBaseAddr(pwriterRec->PixMap));
}

int PsychAddVideoFrameToMovie(int moviehandle, int frameDurationUnits, psych_bool isUpsideDown)
{
	PsychMovieWriterRecordType* pwriterRec = PsychGetMovieWriter(moviehandle, FALSE);

	OSErr				myErr = noErr;
	int					x, y, w, h;
	unsigned char*		pixptr   = (unsigned char*) GetPixBaseAddr(pwriterRec->PixMap);
	unsigned int*		wordptr  = (unsigned int*) GetPixBaseAddr(pwriterRec->PixMap);
	unsigned int		*wordptr2, *wordptr1;
	unsigned int		dummy;

	if (NULL == pwriterRec->Media) return(0);

	if ((frameDurationUnits < 1) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING:In AddFrameToMovie: Negative or zero 'frameduration' %i units for moviehandle %i provided! Sounds like trouble ahead.\n", frameDurationUnits, moviehandle);

	// Draw testpattern: Disabled at compile-time by default:
	if (FALSE) {
		for (y = 0; y < pwriterRec->height; y++) {
			for (x = 0; x < pwriterRec->width; x++) {
				*(pixptr++) = (unsigned char) 255; // alpha
				*(pixptr++) = (unsigned char) y; // Red
				*(pixptr++) = (unsigned char) x; // Green
				*(pixptr++) = (unsigned char) 0; // Blue
			}
		}
	}
	
	// Imagebuffer is upside-down: Need to flip it vertically:
	if (isUpsideDown) {
		h = pwriterRec->height - pwriterRec->padding; // Hack pwriterRec->padding !
		w = pwriterRec->width;
		wordptr1 = wordptr;
		for (y = 0; y < h/2; y++) {
			wordptr2 = wordptr;
			wordptr2 += ((h - 1 - y) * w);
			for (x = 0; x < w; x++) {
				dummy = *wordptr1;
				*(wordptr1++) = *wordptr2;
				*(wordptr2++) = dummy;				 
			}
		}
	}
	
	// Apply codec to compress image:
	myErr = CompressImage(	pwriterRec->PixMap, 
							&pwriterRec->Rect, 
							pwriterRec->CodecQuality,
							pwriterRec->CodecType,
							pwriterRec->ImageDesc, 
							pwriterRec->ComprDataPtr);
	if (myErr != noErr) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In AddFrameToMovie: Video compression on current frame for moviehandle %i failed [CompessImage() returned QT error code %i]!\n", moviehandle, (int) myErr);
		goto bail;
	}
	
	// Add encoded buffer to movie:
	myErr = AddMediaSample(	pwriterRec->Media, 
							pwriterRec->ComprDataHdl,
							0,								// no offset in data
							(**(pwriterRec->ImageDesc)).dataSize, 
							frameDurationUnits,			// frame duration
							(SampleDescriptionHandle) pwriterRec->ImageDesc, 
							1,								// one sample
							0,								// self-contained samples
							NULL);
	if (myErr != noErr) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In AddFrameToMovie: Adding current frame to moviehandle %i failed [AddMediaSample() returned QT error code %i]!\n", moviehandle, (int) myErr);
		goto bail;
	}

	if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG:In AddFrameToMovie: Added new videoframe with %i units duration and upsidedown = %i to moviehandle %i.\n", frameDurationUnits, (int) isUpsideDown, moviehandle);

bail:
	return(myErr);
}

int PsychCreateNewMovieFile(char* moviefile, int width, int height, double framerate, char* movieoptions)
{
	PsychMovieWriterRecordType* pwriterRec = NULL;
	int						moviehandle = 0;
	long					myFlags = createMovieFileDeleteCurFile | createMovieFileDontCreateResFile;
	OSErr					myErr = noErr;
	char*					poption;
	int						dummyInt;
	float					dummyFloat;
	char					myfourcc[4];

	// Still capacity left?
	if (moviewritercount >= PSYCH_MAX_MOVIEWRITERDEVICES) PsychErrorExitMsg(PsychError_user, "Maximum number of movie writers exceeded. Please close some first!");

	// Find first free (i.e., NULL) slot and assign moviehandle:
	while ((pwriterRec = PsychGetMovieWriter(moviehandle, TRUE)) && pwriterRec->Movie) moviehandle++;

	if (firsttime) {
#if PSYCH_SYSTEM == PSYCH_WINDOWS
        // Initialize Quicktime for Windows compatibility layer: This will fail if
        // QT isn't installed on the Windows machine...
        myErr = InitializeQTML(0);
        if (myErr!=noErr) {
            printf("PTB-ERROR: Quicktime Media Layer initialization failed with error code %i in call to InitializeQTML(0): Quicktime not properly installed?!?", (int) myErr);
            PsychErrorExitMsg(PsychError_system, "Quicktime Media Layer initialization failed. Quicktime not properly installed or not installed at all?!?");
        }
#endif

        // Initialize Quicktime-Subsystem:
        myErr = EnterMovies();
        if (myErr!=noErr) {
            printf("PTB-ERROR: Quicktime initialization failed with error code %i in call to EnterMovies().", (int) myErr);
            PsychErrorExitMsg(PsychError_system, "Quicktime initialization in EnterMovies() failed!!!");
        }
		firsttime = FALSE;
	}

	// Translate char string with movie name to FSSpec:
	NativePathNameToFSSpec(moviefile, &pwriterRec->File, 0);

	pwriterRec->CodecType = kH264CodecType;
	pwriterRec->CodecQuality = codecNormalQuality; 
	pwriterRec->ResID = movieInDataForkResID;	
	pwriterRec->height = height;
	pwriterRec->width = width;
	pwriterRec->padding = (PSYCH_SYSTEM == PSYCH_OSX) ? 1 : 0;
	
	// Assign numeric 32-bit FOURCC equivalent code to select codec:
	// This is optional. We default to kH264CodecType:
	if ((poption = strstr(movieoptions, "CodecFOURCCId="))) {
		if (sscanf(poption, "CodecFOURCCId=%i", &dummyInt) == 1) {
			pwriterRec->CodecType = (CodecType) dummyInt;
			if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Codec with FOURCC numeric id %i selected for encoding of movie %i [%s].\n", dummyInt, moviehandle, moviefile);
		}
		else PsychErrorExitMsg(PsychError_user, "Invalid CodecFOURCCId= parameter provided in movieoptions parameter. Parse error!");
	}

	// Assign 4 character string FOURCC code to select codec:
	if ((poption = strstr(movieoptions, "CodecFOURCC="))) {
		if (sscanf(poption, "CodecFOURCC=%c%c%c%c", &myfourcc[3], &myfourcc[2], &myfourcc[1], &myfourcc[0]) == 4) {
			pwriterRec->CodecType = (CodecType) (*((unsigned int*) (&myfourcc[0])));
			if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Codec with FOURCC '%c%c%c%c' = numeric id %i selected for encoding of movie %i [%s].\n", myfourcc[3], myfourcc[2], myfourcc[1], myfourcc[0], (int) pwriterRec->CodecType, moviehandle, moviefile);
		}
		else PsychErrorExitMsg(PsychError_user, "Invalid CodecFOURCC= parameter provided in movieoptions parameter. Must be exactly 4 characters! Parse error!");
	}

	// Assign numeric encoding quality level:
	// This is optional. We default to "normal quality":
	if ((poption = strstr(movieoptions, "EncodingQuality="))) {
		if (sscanf(poption, "EncodingQuality=%f", &dummyFloat) == 1) {
			// Map floating point quality level between 0.0 and 1.0 to 10 discrete levels:
			dummyInt = (int)(10.0 * dummyFloat + 0.5);
			dummyInt = (dummyInt < 0)  ?  0 : dummyInt;
			dummyInt = (dummyInt > 10) ? 10 : dummyInt;
			
			// Assign one of Quicktime's 6 code quality setting to the 10 levels:
			switch(dummyInt) {
				case 0:
					pwriterRec->CodecQuality = codecMinQuality;
				break;

				case 1:
					pwriterRec->CodecQuality = codecMinQuality;
				break;

				case 2:
					pwriterRec->CodecQuality = codecLowQuality;
				break;

				case 3:
					pwriterRec->CodecQuality = codecLowQuality;
				break;

				case 4:
					pwriterRec->CodecQuality = codecNormalQuality;
				break;

				case 5:
					pwriterRec->CodecQuality = codecNormalQuality;
				break;

				case 6:
					pwriterRec->CodecQuality = codecNormalQuality;
				break;

				case 7:
					pwriterRec->CodecQuality = codecHighQuality;
				break;

				case 8:
					pwriterRec->CodecQuality = codecHighQuality;
				break;

				case 9:
					pwriterRec->CodecQuality = codecMaxQuality;
				break;

				case 10:
					pwriterRec->CodecQuality = codecLosslessQuality;
				break;
				
			}
			
			if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Encoding quality level %i selected for encoding of movie %i [%s].\n", dummyInt, moviehandle, moviefile);
		}
		else PsychErrorExitMsg(PsychError_user, "Invalid EncodingQuality= parameter provided in movieoptions parameter. Parse error!");
	}

	// Check for valid parameters. Also warn if some parameters are borderline for certain codecs:
	if ((framerate < 1) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING:In CreateMovie: Negative or zero 'framerate' %f units for moviehandle %i provided! Sounds like trouble ahead.\n", (float) framerate, moviehandle);
	if (width < 1) PsychErrorExitMsg(PsychError_user, "In CreateMovie: Invalid zero or negative 'width' for video frame size provided!");
	if ((width < 4) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING:In CreateMovie: 'width' of %i pixels for moviehandle %i provided! Some video codecs may malfunction with such a small width.\n", width, moviehandle);
	if ((width % 4 != 0) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING:In CreateMovie: 'width' of %i pixels for moviehandle %i provided! Some video codecs may malfunction with a width which is not a multiple of 4 or 16.\n", width, moviehandle);
	if (height < 1) PsychErrorExitMsg(PsychError_user, "In CreateMovie: Invalid zero or negative 'height' for video frame size provided!");
	if ((height < 4) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING:In CreateMovie: 'height' of %i pixels for moviehandle %i provided! Some video codecs may malfunction with such a small height.\n", height, moviehandle);

	// Create a movie file for the destination movie:
	myErr = CreateMovieFile(&pwriterRec->File, FOUR_CHAR_CODE('TVOD'), smCurrentScript, myFlags, &pwriterRec->ResRefNum, &pwriterRec->Movie);
	if (myErr != noErr) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In CreateMovie: Creating movie file with handle %i [%s] failed [CreateMovieFile() returned QT error code %i]!\n", moviehandle, moviefile, (int) myErr);
		if ((myErr == -47) && (PsychPrefStateGet_Verbosity() > 0)) printf("PTB-ERROR:In CreateMovie: Do you have this file open or highlighted in some other application like Quicktime player or Finder?!?\n");
		goto bail;
	}

	// Create the movie track and media:
	pwriterRec->Track = NewMovieTrack(pwriterRec->Movie, FixRatio(width, 1), FixRatio(height, 1), kNoVolume);
	myErr = GetMoviesError();
	if (myErr != noErr) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In CreateMovie: Creating track with videoframes of size %i x %i pixels for movie file %i [%s] failed [NewMovieTrack() returned QT error code %i]!\n", width, height, moviehandle, moviefile, (int) myErr);
		goto bail;
	}
		
	pwriterRec->Media = NewTrackMedia(pwriterRec->Track, VideoMediaType, (int) framerate, NULL, 0);
	myErr = GetMoviesError();
	if (myErr != noErr) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In CreateMovie: Creating video media stream with framerate %i for file with handle %i [%s] failed [NewTrackMedia() returned QT error code %i]!\n", (int) framerate, moviehandle, moviefile, (int) myErr);
		goto bail;
	}

	// Create the media samples:
	myErr = BeginMediaEdits(pwriterRec->Media);
	if (myErr != noErr) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In CreateMovie: Creating movie file with handle %i [%s] failed [BeginMediaEdits() returned QT error code %i]!\n", moviehandle, moviefile, (int) myErr);
		goto bail;
	}

	// Create a GWorld as target with a pixedepth of 32 bpp:
	// Hack hack pwriterRec->padding!!
	MacSetRect(&pwriterRec->Rect, 0, 0, width - pwriterRec->padding, height - pwriterRec->padding);
	myErr = NewGWorld(&pwriterRec->GWorld, 32, &pwriterRec->Rect, NULL, NULL, (GWorldFlags)0);
	if (myErr != noErr) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In CreateMovie: Creating movie file with handle %i [%s] failed [NewGWorld(rectsize = %i,%i and padding %i) returned QT error code %i]!\n", moviehandle, moviefile, width - pwriterRec->padding, height - pwriterRec->padding, pwriterRec->padding, (int) myErr);
		goto bail;
	}

	pwriterRec->PixMap = GetGWorldPixMap(pwriterRec->GWorld);
	if (pwriterRec->PixMap == NULL) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In CreateMovie: Creating movie file with handle %i [%s] failed [GetGWorldPixMap() returned NULL-Ptr]!\n", moviehandle, moviefile);
		goto bail;
	}

	LockPixels(pwriterRec->PixMap);
	myErr = GetMaxCompressionSize(	pwriterRec->PixMap,
									&pwriterRec->Rect, 
									0,							// let ICM choose depth
									pwriterRec->CodecQuality, 
									pwriterRec->CodecType, 
									(CompressorComponent) anyCodec,
									&pwriterRec->MaxComprSize);
	if (myErr != noErr) {
		if (PsychPrefStateGet_Verbosity() > 0) {
			printf("PTB-ERROR:In CreateMovie: Creating movie file with handle %i [%s] failed [GetMaxCompressionSize() returned QT error code %i]!\n", moviehandle, moviefile, (int) myErr);
			if (myErr == noCodecErr) printf("PTB-ERROR:In CreateMovie: The system doesn't know or support the requested type of video movie codec.\n");
		}
		goto bail;
	}

	pwriterRec->ComprDataHdl = NewHandle(pwriterRec->MaxComprSize);
	if (pwriterRec->ComprDataHdl == NULL) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In CreateMovie: Creating movie file with handle %i [%s] failed [NewHandle(for compressor) returned NULL-Ptr]!\n", moviehandle, moviefile);
		goto bail;
	}

	HLockHi(pwriterRec->ComprDataHdl);
	pwriterRec->ComprDataPtr = *(pwriterRec->ComprDataHdl);

	pwriterRec->ImageDesc = (ImageDescriptionHandle)NewHandle(4);
	if (pwriterRec->ImageDesc == NULL) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In CreateMovie: Creating movie file with handle %i [%s] failed [NewHandle(4) returned NULL-Ptr]!\n", moviehandle, moviefile);
		goto bail;
	}

	SetGWorld(pwriterRec->GWorld, NULL);
	
	// Increment count of open movie writers:
	moviewritercount++;
	
	if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Moviehandle %i successfully opened for movie writing into file '%s'.\n", moviehandle, moviefile);

	// Return new handle:
	return(moviehandle);

bail:

	if (pwriterRec->ComprDataHdl != NULL)
		DisposeHandle(pwriterRec->ComprDataHdl);

	if (pwriterRec->GWorld != NULL)
		DisposeGWorld(pwriterRec->GWorld);
		
	if (pwriterRec->ResRefNum != 0)
		CloseMovieFile(pwriterRec->ResRefNum);

	if (pwriterRec->Movie != NULL)
		DisposeMovie(pwriterRec->Movie);

	pwriterRec->Movie = NULL;
	pwriterRec->Track = NULL;
	pwriterRec->Media = NULL;
	pwriterRec->ResRefNum = 0;
	pwriterRec->GWorld = NULL;
	pwriterRec->PixMap = NULL;
	pwriterRec->MaxComprSize = 0L;
	pwriterRec->ComprDataHdl = NULL;
	pwriterRec->ComprDataPtr = NULL;
	pwriterRec->ImageDesc = NULL;
	pwriterRec->SavedPort = NULL;
	pwriterRec->SavedDevice = NULL;

	// Return failure:
	return(-1);
}

int PsychFinalizeNewMovieFile(int movieHandle)
{
	OSErr					myErr = noErr;

	PsychMovieWriterRecordType* pwriterRec = PsychGetMovieWriter(movieHandle, FALSE);

	if (NULL == pwriterRec->Media) return(0);

	myErr = EndMediaEdits(pwriterRec->Media);
	if (myErr != noErr) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In FinalizeMovie: Finalizing/Closing movie file with handle %i failed [EndMediaEdits() returned QT error code %i]!\n", movieHandle, (int) myErr);
		goto bail;
	}
	
	// add the media to the track
	myErr = InsertMediaIntoTrack(pwriterRec->Track, 0, 0, GetMediaDuration(pwriterRec->Media), fixed1);
	if (myErr != noErr) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In FinalizeMovie: Finalizing/Closing movie file with handle %i failed [InsertMediaIntoTrack() returned QT error code %i]!\n", movieHandle, (int) myErr);
		goto bail;
	}
	
	// add the movie atom to the movie file
	myErr = AddMovieResource(pwriterRec->Movie, pwriterRec->ResRefNum, &pwriterRec->ResID, NULL);
	if (myErr != noErr) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In FinalizeMovie: Finalizing/Closing movie file with handle %i failed [AddMovieResource() returned QT error code %i]!\n", movieHandle, (int) myErr);
		goto bail;
	}

bail:
	if (pwriterRec->ImageDesc != NULL)
		DisposeHandle((Handle)pwriterRec->ImageDesc);

	if (pwriterRec->ComprDataHdl != NULL)
		DisposeHandle(pwriterRec->ComprDataHdl);

	if (pwriterRec->GWorld != NULL)
		DisposeGWorld(pwriterRec->GWorld);

	if (pwriterRec->ResRefNum != 0)
		CloseMovieFile(pwriterRec->ResRefNum);

	if (pwriterRec->Movie != NULL)
		DisposeMovie(pwriterRec->Movie);

	myErr = GetMoviesError();
	if (myErr != noErr) {
		if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR:In FinalizeMovie: Finalizing/Closing movie file with handle %i failed [CloseMovieFile() or DisposeMovie() returned QT error code %i]!\n", movieHandle, (int) myErr);
	}

	pwriterRec->Movie = NULL;
	pwriterRec->Track = NULL;
	pwriterRec->Media = NULL;
	pwriterRec->ResRefNum = 0;
	pwriterRec->GWorld = NULL;
	pwriterRec->PixMap = NULL;
	pwriterRec->MaxComprSize = 0L;
	pwriterRec->ComprDataHdl = NULL;
	pwriterRec->ComprDataPtr = NULL;
	pwriterRec->ImageDesc = NULL;
	pwriterRec->SavedPort = NULL;
	pwriterRec->SavedDevice = NULL;
	
	// Decrement count of active writers:
	moviewritercount--;

	// Return success/fail status:
	if ((myErr == noErr) && (PsychPrefStateGet_Verbosity() > 3)) printf("PTB-INFO: Moviehandle %i successfully closed and movie written to filesystem.\n", movieHandle);
	return(myErr == 0);
}

psych_bool PsychAddAudioBufferToMovie(int moviehandle, unsigned int nrChannels, unsigned int nrSamples, double* buffer)
{
    PsychErrorExitMsg(PsychError_unimplemented, "Sorry, storing audio tracks in movies is not supported by the Quicktime based movie writing functions.");
    return(0);
}

// End of routines.
#endif // QT-Implementation.
#endif // non-Linux.
