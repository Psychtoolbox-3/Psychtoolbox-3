/*
 *        PsychToolbox3/Source/Common/PsychPortAudio/PsychPortAudio.c
 *
 *        PLATFORMS:        All
 *
 *        AUTHORS:
 *
 *        Mario Kleiner     mk        mario.kleiner.de@gmail.com
 *
 *        HISTORY:
 *
 *        21.03.2007        mk        wrote it.
 *        03.04.2011        mk        Make 64 bit clean. Allow 64-bit sized operations and float matrices.
 *        03.04.2011        mk        License changed to MIT with some restrictions.
 *
 *        DESCRIPTION:
 *
 *        Low level Psychtoolbox sound i/o driver. Useful for audio capture, playback and
 *        feedback with well controlled timing and low latency. Uses the free software
 *        PortAudio library, API Version 19 (http://www.portaudio.com), which has a MIT style license.
 *
 *        This seems to link statically against libportaudio.a on OS/X. However, it doesn't!
 *        Due to some restrictions of the OS/X linker we can't link statically against portaudio,
 *        so we have to have a .dylib version of the library installed in a system library
 *        search path, both for compiling and using PsychPortAudio. The current universal
 *        binary version of libportaudio.a for OS/X can be found in the...
 *        PsychSourceGL/Cohorts/PortAudio/ subfolder, which also contains versions for
 *        Linux and Windows.
 *
 */

#include "PsychPortAudio.h"

#if PSYCH_SYSTEM == PSYCH_OSX
#include "pa_mac_core.h"
#endif

#ifdef PTB_USE_ASIO
#include "pa_asio.h"
#endif

#if PSYCH_SYSTEM == PSYCH_WINDOWS
#include "pa_win_wasapi.h"
#endif

#if PSYCH_SYSTEM == PSYCH_LINUX
#include "pa_linux_alsa.h"
#include <alsa/asoundlib.h>
#include <dlfcn.h>

void (*myjack_set_error_function)(void(*)(const char *)) = NULL;

// Dummy error handler to swallow pointless ALSA debug/warning/errr messages, if handler is attached:
static void ALSAErrorHandler(const char *file, int line, const char *function, int err, const char *fmt, ...)
{
}

// Pseudo-Defines of Portaudio internal structs, mimicking the memory layout
// of the true Portaudio internal structs just enough so that we can cast a PaStream*
// into a PaAlsaStream*, then access a stream-internal variable of the buffer-processor
// and fudge with it in the most hackish and puke inducing way to work around a bug
// present in the audio-capture only mode setup code (ie. no playback) inside the ALSA
// backend since 24th May 2009! This over 9 years old bug causes massive distortion if
// one only uses audio-capture, but no playback or full-duplex mode.
//
// A clean way of fixing this will be to get a fix into upstream, but we need to support
// all the existing libportaudio implementations in all Linux distributions from 2009 to
// 2019. The only option for shipping distros therefore would be to go back to static
// linking against our own bug-fixed variant of libportaudio, but i haven't spent the
// better part of a week to get us to finally be able to drop our own private libportaudio
// builds, just to reintroduce them now! Therefore ugly hack it is...

typedef enum {
    paUtilFixedHostBufferSize,
    paUtilBoundedHostBufferSize,
    paUtilUnknownHostBufferSize,
    paUtilVariableHostBufferSizePartialUsageAllowed
} PaUtilHostBufferSizeMode;

typedef struct {
    unsigned long framesPerUserBuffer;
    unsigned long framesPerHostBuffer;

    PaUtilHostBufferSizeMode hostBufferSizeMode;
    int useNonAdaptingProcess;
    int userOutputSampleFormatIsEqualToHost;
    int userInputSampleFormatIsEqualToHost;
    unsigned long framesPerTempBuffer;

    unsigned int inputChannelCount;
    unsigned int bytesPerHostInputSample;
    unsigned int bytesPerUserInputSample;
    int userInputIsInterleaved;
    // More stuff follows in real data struct...
} PseudoBufferProcessor;

typedef struct {
    double samplingPeriod;
    double measurementStartTime;
    double averageLoad;
} PseudoCpuLoadMeasurer;

#define PA_STREAM_MAGIC (0x18273645)

typedef struct PseudoStreamRepresentation {
    unsigned long magic;    /**< set to PA_STREAM_MAGIC */
    struct PseudoStreamRepresentation *nextOpenStream; /**< field used by multi-api code */
    void *streamInterface;
    void *streamCallback;
    void *streamFinishedCallback;
    void *userData;
    PaStreamInfo streamInfo;
} PseudoStreamRepresentation;

typedef struct PseudoStreamRepresentationExt {
    unsigned long magic;    /**< set to PA_STREAM_MAGIC */
    struct PseudoStreamRepresentationExt *nextOpenStream; /**< field used by multi-api code */
    void *streamInterface;
    void *streamCallback;
    void *streamFinishedCallback;
    void *userData;
    PaStreamInfo streamInfo;
    PaHostApiTypeId hostApiType;
} PseudoStreamRepresentationExt;

typedef struct PseudoAlsaStream {
    PseudoStreamRepresentation streamRepresentation;
    PseudoCpuLoadMeasurer cpuLoadMeasurer;
    PseudoBufferProcessor bufferProcessor;
    // More stuff follows in real data struct...
} PseudoAlsaStream;

typedef struct PseudoAlsaStreamExt {
    PseudoStreamRepresentationExt streamRepresentation;
    PseudoCpuLoadMeasurer cpuLoadMeasurer;
    PseudoBufferProcessor bufferProcessor;
    // More stuff follows in real data struct...
} PseudoAlsaStreamExt;

#endif

// Need to define these as they aren't defined in portaudio.h
// for some mysterious reason:
typedef void (*PaUtilLogCallback ) (const char *log);
void PaUtil_SetDebugPrintFunction(PaUtilLogCallback  cb);

#if PSYCH_SYSTEM == PSYCH_LINUX
void (*myPaUtil_SetDebugPrintFunction)(PaUtilLogCallback  cb) = NULL;

// Wrapper implementation, as many libportaudio.so implementations seem to lack this function :(:
void PaUtil_SetDebugPrintFunction(PaUtilLogCallback  cb)
{
    // Try to get function dynamically:
    myPaUtil_SetDebugPrintFunction = dlsym(RTLD_NEXT, "PaUtil_SetDebugPrintFunction");

    if (myPaUtil_SetDebugPrintFunction)
        myPaUtil_SetDebugPrintFunction(cb);

    return;
}
#endif

// Forward define of prototype of our own custom new PortAudio extension function for Zero latency direct input monitoring:
PaError Pa_DirectInputMonitoring(PaStream *stream, int enable, int inputChannel, int outputChannel, double gain, double pan);

#define MAX_SYNOPSIS_STRINGS 50
static const char *synopsisSYNOPSIS[MAX_SYNOPSIS_STRINGS];

#define kPortAudioPlayBack   1
#define kPortAudioCapture    2
#define kPortAudioFullDuplex 3
#define kPortAudioMonitoring 4
#define kPortAudioIsMaster   8
#define kPortAudioIsSlave    16
#define kPortAudioIsAMModulator         32
#define kPortAudioIsOutputCapture       64
#define kPortAudioIsAMModulatorForSlave 128

// Maximum number of audio devices we handle:
// This consumes around 200 Bytes static memory per potential device, so
// we waste about 200 kb or RAM here, which is acceptable nowadays. However: If a device
// is actually opened, it will consume additional memory ressources in PortAudio,
// the operating systems sound subsystem and kernel, and in the audio hardware device
// driver. It will also consume (potentially limited) audio hardware ressources,
// and most importantly, it will consume cpu time, bus cycles and hw/interrupts for the running
// audio processing threads inside PortAudio, the OS sound system and device drivers!
// --> Too many simultaneously open devices may hurt performance and general timing, even
// if the devices are mostly idle!!
// --> Anybody that manages to even come close to this limit does something horribly wrong and certainly non-portable!
#define MAX_PSYCH_AUDIO_DEVS 1024

// Maximum number of audio channels we support per open device:
#define MAX_PSYCH_AUDIO_CHANNELS_PER_DEVICE 256

// Maximum number of attached slave devices we support per open master device:
#define MAX_PSYCH_AUDIO_SLAVES_PER_DEVICE 1024

// Initial size (and increment) of audio buffer list. List will grow by that
// many slots whenever it needs to grow:
#define PSYCH_AUDIO_BUFFERLIST_INCREMENT 1024

// PA_ANTICLAMPGAIN is premultiplied onto any sample provided by usercode, reducing
// signal amplitude by a tiny fraction. This is a workaround for a bug in the
// sampleformat converters in Portaudio for float -> 32 bit int and float -> 24 bit int.
// They cause integer wraparound and thereby audio artifacts if the signal has an
// amplitude of almost +1.0f. This attenuates conformant signals in -1 to +1 range just
// enough to "fix it". Note: The single precision float math inside our driver and inside
// PortAudio limits the effective precision of sound signals to 23 bits + sign = 24 bits,
// so we just manage to properly feed a 24 bit audio DAC, but there's zero headroom left
// for higher precision and our true precision when using gain/volume modulation or mixing
// may be a bit lower than 24 bits due to accumulated numeric roundoff errors. Theoretically
// we could switch to double precision math inside our driver to at least retain 24 bits
// in the mixing/modulation stage.
#define PA_ANTICLAMPGAIN 0.9999999

// Uncomment this define MUTEX_LOCK_TIME_STATS to enable tracing of
// mutex lock hold times for low-level debugging and tuning:
//#define MUTEX_LOCK_TIME_STATS 1

typedef struct PsychPASchedule {
    unsigned int    mode;                   // Mode of schedule slot: 0 = Invalid slot, > 0 valid slot, where different bits in the int mean something...
    double          repetitions;            // Number of repetitions for the playloop defined in this slot.
    psych_int64     loopStartFrame;         // Start of playloop in frames.
    psych_int64     loopEndFrame;           // End of playloop in frames.
    int             bufferhandle;           // Handle of the playout buffer to use. Zero is the standard playbuffer as set by 'FillBuffer'. Negative handles
                                            // may have special meaning in future implementations.
    double          tWhen;                  // Time in seconds, either absolute or relative spec, depending on command.
    unsigned int    command;                // Command code: 0 = Normal playback buffer. 1 = Pause & Restart playback, 2 = Schedule end of playback, ..
} PsychPASchedule;

// Our device record:
typedef struct PsychPADevice {
    psych_mutex             mutex;          // Mutex lock for the PsychPADevice struct.
    psych_condition         changeSignal;   // Condition variable or event object for change signalling (see above).
    int                     opmode;         // Mode of operation: Playback, capture or full duplex? Master, Slave or standalone?
    int                     runMode;        // Runmode: 0 = Stop engine at end of playback, 1 = Keep engine running in hot-standby, ...
    PaStream *stream;                       // Pointer to associated portaudio stream.
    const PaStreamInfo*     streaminfo;     // Pointer to stream info structure, provided by PortAudio.
    PaHostApiTypeId         hostAPI;        // Type of host API.
    int                     indeviceidx;    // Device index of capture device. -1 if none open.
    int                     outdeviceidx;   // Device index of output device. -1 if none open.
    volatile double         reqStartTime;   // Requested start time in system time (secs).
    volatile double         startTime;      // Real start time in system time (secs). Returns real start time after start.
    // The real start time is the time when the first sample hit the speaker in playback or full-duplex mode.
    // Its the time when the first sample was captured in pure capture mode - or when the first sample should
    // be captured in pure capture mode with scheduled start. Whenever playback is active, startTime only
    // refers to the output stage.
    volatile double         captureStartTime;       // Time when first captured sample entered the sound hardware in system time (secs). This information is
    // redundant in pure capture mode - its identical to startTime. In full duplex mode, this is an estimate of
    // when the first sample was captured, whereas startTime is an estimate of when the first sample was output.
    volatile double reqStopTime;    // Requested stop time in system time (secs). Set to DBL_MAX if none requested.
    volatile double estStopTime;    // Estimated sound offset time after stop of playback.
    volatile double currentTime;    // Current playout time of the last sound sample submitted to the engine. Will be wrong in case of playback abort!
    volatile unsigned int state;    // Current state of the stream: 0=Stopped, 1=Hot Standby, 2=Playing, 3=Aborting playback. Mostly written/updated by paCallback.
    volatile unsigned int reqstate; // Requested state of the stream, as opposed to current 'state'. Written by main-thread, read & processed by paCallback.
    double     repeatCount;         // Number of repetitions: -1 = Loop forever, 1 = Once, n = n repetitions.
    float*     outputbuffer;        // Pointer to float memory buffer with sound output data.
    psych_int64 outputbuffersize;   // Size of output buffer in bytes.
    psych_int64 loopStartFrame;     // Start of current playloop in frames.
    psych_int64 loopEndFrame;       // End of current playloop in frames.
    psych_int64 playposition;       // Current playposition in samples since start of playback for current buffer and playloop (not frames, not bytes!)
    psych_int64 writeposition;      // Current writeposition in samples since start of playback (for incremental filling).
    psych_int64 totalplaycount;     // Total running count of samples since start of playback, accumulated over all buffers and playloop(not frames, not bytes!)
    float*     inputbuffer;         // Pointer to float memory buffer with sound input data (captured sound data).
    psych_int64 inputbuffersize;    // Size of input buffer in bytes.
    psych_int64 recposition;        // Current record position in samples since start of capture.
    psych_int64 readposition;       // Last read-out sample since start of capture.
    psych_int64 outchannels;        // Number of output channels.
    psych_int64 inchannels;         // Number of input channels.
    unsigned int xruns;             // Number of over-/underflows of input-/output channel for this stream.
    unsigned int paCalls;           // Number of callback invocations.
    unsigned int noTime;            // Number of timestamp malfunction - Should not happen anymore.
    psych_int64 batchsize;          // Maximum number of frames requested during callback invokation: Estimate of real buffersize.
    double     predictedLatency;    // Latency that PortAudio predicts for current callbackinvocation. We will compensate for that when starting audio.
    double   latencyBias;           // A bias value to add to the value that PortAudio reports for total buffer->Speaker latency.
    // This value defaults to zero, but can be set up automatically on OS/X or manually on other OSes to compensate
    // for slight mistakes in PA's estimate.

    // Audio schedule related:
    PsychPASchedule* schedule;      // Pointer to start of array with playback schedule, or a NULL pointer if none defined.
    volatile unsigned int schedule_size;    // Size of schedule array in slots.
    volatile unsigned int schedule_pos;     // Current position in schedule (in slots).
    unsigned int schedule_writepos;         // Current position in schedule (in slots).

    // Master-Slave virtual device related:
    int*    outputmappings;         // Mapping array of output slave channels to associated master channels for mix and merge. NULL on master devices.
    int*    inputmappings;          // Mapping array of input slave channels to associated master channels for distribution. NULL on master devices.
    int    slaveCount;              // Number of attached slave devices. Zero on slave devices.
    int*    slaves;                 // Array of pahandle's of all attached slave devices, ie., an array with slaveCount valid (non -1) entries. NULL on slaves.
    int    pamaster;                // pahandle of master device for a slave. -1 on master devices.
    int    slaveDirty;              // Flag: 0 means that a slave didn't do anything, so no mixdown/merge by master required. 1 means: Do mixdown/merge.
    float*    slaveOutBuffer;       // Temporary output buffer for slaves to store their output data. Used as input for output mix/merge. NULL on non-masters.
    float*    slaveInBuffer;        // Temporary input buffer for slaves to receive their input data. Used as output from distributor. NULL on non-masters.
    float*    slaveGainBuffer;      // Temporary output buffer for AM modulator slaves to store their gain output data. NULL on non AMModulators for slaves.
    int    modulatorSlave;          // pahandle of a slave device that acts as a modulator for this device. -1 if none assigned.
    double    firstsampleonset;     // Cached sample onset time from paCallback.
    double    cst;                  // Cached captured sample onset time from paCallback.
    double    now;                  // Cached invocation time from paCallback.

    // Mixer volume related:
    float*    outChannelVolumes;    // Array of per-outputchannel volume settings on slave devices, NULL and not used on non-slave devices.
    float    masterVolume;          // Master volume setting for all non-slave audio devices, i.e., masters and regular devices. Unused on slaves.
} PsychPADevice;

PsychPADevice audiodevices[MAX_PSYCH_AUDIO_DEVS];
unsigned int  audiodevicecount = 0;
unsigned int  verbosity = 4;
double        yieldInterval = 0.001;            // How long to wait in calls to PsychYieldIntervalSeconds().
psych_bool    uselocking = TRUE;                // Use Mutex locking and signalling code for thread synchronization?
psych_bool    lockToCore1 = TRUE;               // NO LONGER USED: Lock all engine threads to run on cpu core 1 on Windows to work around broken TSC sync on multi-cores?
psych_bool    pulseaudio_autosuspend = TRUE;    // Should we try to suspend the Pulseaudio sound server on Linux while we're active?
psych_bool    pulseaudio_isSuspended = FALSE;   // Is PulseAudio suspended by us?

double debugdummy1, debugdummy2;

psych_bool pa_initialized = FALSE;

// Definition of an audio buffer:
struct PsychPABuffer_Struct {
    unsigned int locked;            // locked: >= 1 = Buffer in use by some active audio device. 0 = Buffer unused.
    float*     outputbuffer;        // Pointer to float memory buffer with sound output data.
    psych_int64 outputbuffersize;   // Size of output buffer in bytes.
    psych_int64 outchannels;        // Number of channels.
};

typedef struct PsychPABuffer_Struct PsychPABuffer;

psych_mutex    bufferListmutex;            // Mutex lock for the audio bufferList.
PsychPABuffer*  bufferList;                // Pointer to start of audio bufferList.
int    bufferListCount;                    // Number of slots allocated in bufferList.

// Return the first unused/closed device handle:
unsigned int PsychPANextHandle(void)
{
    int i;
    for (i = 0; i < MAX_PSYCH_AUDIO_DEVS; i++) {
        if (NULL == audiodevices[i].stream)
            break;
    }

    return(i);
}

// Scan all schedules of all active and open audio devices to check if
// given audiobuffer is referenced. Invalidate reference, if so:
// The special handle == -1 invalidates all references except the ones to special buffer zero.
psych_bool PsychPAInvalidateBufferReferences(int handle)
{
    unsigned int i, j;
    psych_bool anylocked = FALSE;

    // Scan all open audio devices:
    for(i = 0; i < MAX_PSYCH_AUDIO_DEVS; i++) {
        // Device open?
        if (audiodevices[i].stream) {
            // Schedule attached?
            if (audiodevices[i].schedule) {
                // Scan it and mark all references to our buffer(s) as invalid:
                for (j = 0; j < audiodevices[i].schedule_size; j++) {
                    // Slot active and with relevant bufferhandle?
                    if ((audiodevices[i].schedule[j].bufferhandle == handle) || ((audiodevices[i].schedule[j].bufferhandle !=0) && (handle == -1)) ) {
                        // Invalidate this reference:
                        audiodevices[i].schedule[j].mode = 0;
                        audiodevices[i].schedule[j].bufferhandle = 0;
                        anylocked = TRUE;
                    }
                }
            }
        }
    }

    return(anylocked);
}

// Create a new audiobuffer for 'outchannels' audio channels and 'nrFrames' samples
// per channel. Init header, allocate zero-filled memory, enqeue in bufferList.
// Resize/Grow bufferList if neccessary. Return handle to buffer.
int PsychPACreateAudioBuffer(psych_int64 outchannels, psych_int64 nrFrames)
{
    PsychPABuffer* tmpptr;
    int i, handle;

    // Does a bufferList exist?
    if ((bufferListCount <= 0) || (bufferList == NULL)) {
        // First call. Allocate and zero-fill initial bufferList:
        bufferList = (PsychPABuffer*) calloc(PSYCH_AUDIO_BUFFERLIST_INCREMENT, sizeof(PsychPABuffer));
        if (NULL == bufferList) PsychErrorExitMsg(PsychError_outofMemory, "Insufficient free memory for allocating new audio buffers when trying to create internal bufferlist!");

        bufferListCount = PSYCH_AUDIO_BUFFERLIST_INCREMENT;
    }

    // Search a free slot in bufferList: We start at slot 1, ie., we skip slot 0.
    // This because we don't want to ever return a handle of zero, as zero denotes the
    // special per-audiodevice playback buffer.
    i = 1; while ((i < bufferListCount) && (NULL != bufferList[i].outputbuffer)) i++;

    // Success?
    if ((i >= bufferListCount)) {
        // Nope. Could not find free slot. Need to resize the bufferList with more capacity.

        // Need to lock bufferList lock to do this:
        PsychLockMutex(&bufferListmutex);

        // Reallocate bufferList: This may relocate the bufferList:
        tmpptr = (PsychPABuffer*) realloc((void*) bufferList, (bufferListCount + PSYCH_AUDIO_BUFFERLIST_INCREMENT) * sizeof(PsychPABuffer) );
        if (NULL == tmpptr) {
            // Failed! Unlock mutex:
            PsychUnlockMutex(&bufferListmutex);

            // Error out. The old allocation and parameters are still valid:
            PsychErrorExitMsg(PsychError_outofMemory, "Insufficient free memory for allocating new audio buffers when trying to grow internal bufferlist!");
        }

        // Assign new pointer and size:
        bufferList = tmpptr;
        tmpptr = &(bufferList[i]);
        bufferListCount += PSYCH_AUDIO_BUFFERLIST_INCREMENT;

        // Initialize new segment of list to zero:
        memset((void*) tmpptr, 0, PSYCH_AUDIO_BUFFERLIST_INCREMENT * sizeof(PsychPABuffer));

        // Done resizing bufferlist. Unlock mutex:
        PsychUnlockMutex(&bufferListmutex);

        // Ready. 'i' now points to first free slot in new segment of extended bufferList.
    }

    // Assign slotid of bufferList slot in handle:
    handle = i;

    // Invalidate all potential stale references to the new 'handle' in all schedules:
    PsychPAInvalidateBufferReferences(handle);

    // Allocate actual data buffer:
    bufferList[handle].outputbuffersize = outchannels * nrFrames * sizeof(float);
    bufferList[handle].outchannels = outchannels;

    if (NULL == ( bufferList[handle].outputbuffer = (float*) calloc(1, (size_t) bufferList[handle].outputbuffersize) )) {
        // Out of memory: Release bufferList header and error out:
        PsychErrorExitMsg(PsychError_outofMemory, "Insufficient free memory for allocating new audio buffer when trying to allocate actual buffer!");
    }

    // Ok, we're ready with an empty, silence filled audiobuffer. Return its handle:
    return(handle);
}

// Delete all audio buffers and bufferList itself: Called during shutdown.
void PsychPADeleteAllAudioBuffers(void)
{
    int i;

    if (bufferListCount > 0) {
        // Lock list:
        PsychLockMutex(&bufferListmutex);

        // Invalidate all referencing slots in all schedules:
        PsychPAInvalidateBufferReferences(-1);

        // Free all audio buffers:
        for (i = 0; i < bufferListCount; i++) {
            if (NULL != bufferList[i].outputbuffer) free(bufferList[i].outputbuffer);
        }

        // Release memory for bufferheader array itself:
        free(bufferList);
        bufferList = NULL;
        bufferListCount = 0;

        // Unlock list:
        PsychUnlockMutex(&bufferListmutex);
    }

    return;
}

PsychPABuffer* PsychPAGetAudioBuffer(int handle)
{
    // Does buffer with given handle exist?
    if ((handle < 0) || (handle >= bufferListCount) || (bufferList[handle].outputbuffer == NULL)) {
        PsychErrorExitMsg(PsychError_user, "Invalid audio bufferhandle provided! The handle doesn't correspond to an existing audiobuffer.");
    }

    return( &(bufferList[handle]) );
}

// Scan all schedules of all active and open audio devices to check which
// audiobuffers are active and lock them:
psych_bool PsychPAUpdateBufferReferences(void)
{
    int i;
    unsigned int j;
    psych_bool anylocked = FALSE;

    // First we reset all locked flags of all buffers:
    for (i = 0; i < bufferListCount; i++) bufferList[i].locked = 0;

    // Scan all open audio devices:
    for(i = 0; i < MAX_PSYCH_AUDIO_DEVS; i++) {
        // Device open?
        if (audiodevices[i].stream) {
            // Schedule attached and device active?
            if ((audiodevices[i].schedule) && ((audiodevices[i].state > 0) && Pa_IsStreamActive(audiodevices[i].stream))) {
                // Active schedule. Scan it and mark all referenced buffers as locked:
                for (j = 0; j < audiodevices[i].schedule_size; j++) {
                    // Slot active and with valid bufferhandle?
                    if ((audiodevices[i].schedule[j].mode & 2) && (audiodevices[i].schedule[j].bufferhandle > 0)) {
                        // Mark used and active audiobuffer as locked:
                        bufferList[ audiodevices[i].schedule[j].bufferhandle ].locked = 1;
                        anylocked = TRUE;
                    }
                }
            }
        }
    }

    return(anylocked);
}

// Delete audiobuffer 'handle' if this is possible. If it isn't possible
// at the moment, 'waitmode' will determine the strategy:
int PsychPADeleteAudioBuffer(int handle, int waitmode)
{
    // Retrieve buffer:
    PsychPABuffer* buffer = PsychPAGetAudioBuffer(handle);

    // Make sure all buffer locked flags are up to date:
    PsychPAUpdateBufferReferences();

    // Buffer locked?
    if (buffer->locked) {
        // Yes :-( In 'waitmode' zero we fail:
        if (waitmode == 0) return(0);

        // In waitmode 1, we retry spin-waiting until buffer available:
        while (buffer->locked) {
            PsychYieldIntervalSeconds(yieldInterval);
            PsychPAUpdateBufferReferences();
        }
    }

    // Delete buffer:
    if (NULL != buffer->outputbuffer) free(buffer->outputbuffer);
    memset(buffer, 0, sizeof(PsychPABuffer));

    // Success:
    return(1);
}

static void PsychPALockDeviceMutex(PsychPADevice* dev)
{
    #ifdef MUTEX_LOCK_TIME_STATS
    // Compute effective mutex lock hold time. Typical duration on a 2nd Rev. MacBookPro DualCore
    // under typical playback load in the paCallback() is less than 100 microseconds, usually around 50 usecs, so
    // nothing to worry about yet.
    PsychGetAdjustedPrecisionTimerSeconds(&debugdummy1);
    #endif

    if (uselocking) {
        PsychLockMutex(&(dev->mutex));
    }
}

static void PsychPAUnlockDeviceMutex(PsychPADevice* dev)
{
    if (uselocking) {
        PsychUnlockMutex(&(dev->mutex));
    }

    #ifdef MUTEX_LOCK_TIME_STATS
    PsychGetAdjustedPrecisionTimerSeconds(&debugdummy2);
    dev->predictedLatency = debugdummy2 - debugdummy1;
    #endif
}

static void PsychPACreateSignal(PsychPADevice* dev)
{
    if (uselocking) {
        PsychInitCondition(&(dev->changeSignal), NULL);
    }
}


static void PsychPADestroySignal(PsychPADevice* dev)
{
    if (uselocking) {
        PsychDestroyCondition(&(dev->changeSignal));
    }
}

static void PsychPASignalChange(PsychPADevice* dev)
{
    if (uselocking) {
        PsychSignalCondition(&(dev->changeSignal));
    }
}

static void PsychPAWaitForChange(PsychPADevice* dev)
{
    if (uselocking) {
        // Locking and signalling: We have to wait for a signal, and we
        // enter here with the device mutex held. Waiting is operating system dependent:
        PsychWaitCondition(&(dev->changeSignal), &(dev->mutex));
    }
    else {
        // No locking and signalling: Just yield for a bit, then retry...
        PsychYieldIntervalSeconds(yieldInterval);
    }
}

// Callback function which gets called when a portaudio stream (aka our engine) goes idle for any reason:
// This will reset the device state to "idle/stopped" aka 0, reset pending stop requests and signal
// the master thread if it is waiting for this to happen:
void PAStreamFinishedCallback(void *userData)
{
    PsychPADevice* dev = (PsychPADevice*) userData;

    // Lock device struct:
    PsychPALockDeviceMutex(dev);

    // Reset state to zero aka idle / stopped and reset pending requests:
    dev->reqstate = 255;

    // Update "true" state to inactive:
    dev->state = 0;

    // If estimated stop time is still undefined at this point, it won't
    // get computed anymore because the engine is stopped. We choose the
    // last submitted samples playout time as best guess of the real
    // stop time. On a regular stop (where hardware plays out all pending
    // audio buffers) this is probably spot-on. On a fast abort however,
    // this may be a too late estimate if the hardware really stopped
    // immediately and dropped pending audio buffers, but we don't know
    // how much was really played out as this is highly hardware dependent:
    if (dev->estStopTime == 0) dev->estStopTime = dev->currentTime;

    // Signal state change:
    PsychPASignalChange(dev);

    // Unlock device struct:
    PsychPAUnlockDeviceMutex(dev);

    // Ready.
    return;
}

/* Logger callback function to output PortAudio debug messages at 'verbosity' > 5. */
void PALogger(const char* msg)
{
    if (verbosity > 5) printf("PTB-DEBUG: PortAudio says: %s\n", msg);
    return;
}


// Called exclusively from paCallback, with device-mutex held.
// Check if a schedule is defined. If not, return repetition, playloop and bufferparameters
// from the device struct, ie., old behaviour. If yes, check if an update of the schedule is
// needed (ie., progressing to the next slot) and do so if needed. Return parameters from
// current slot, or an abort signal if end of schedule reached.
//
// Return codes:
// 0 = Valid buffer and playback parameters assigned. Proceed with audio playback as scheduled.
// 1 = Regular end of playback or schedule reached. Finish playback.
// 2 = Error abort, something is totally frelled.
// 4 = Abort bufferfill operation for this host audio buffer via zerofill, but don't switch to idle mode / don't stop engine.
//     Instead switch back to hot-standby so playback can be picked up again at a later point in time.
//       This is used to reschedule start of playback for a following slot at a later time.
int PsychPAProcessSchedule(PsychPADevice* dev, psych_int64 *playposition, float** ret_playoutbuffer, psych_int64* ret_outsbsize, psych_int64* ret_outsboffset, double* ret_repeatCount, psych_int64* ret_playpositionlimit)
{
    psych_int64     loopStartFrame, loopEndFrame;
    psych_int64     outsbsize, outsboffset;
    psych_int64     outchannels = dev->outchannels;
    unsigned int    slotid, cmd;
    double          repeatCount;
    double          reqTime = 0;
    psych_int64     playpositionlimit;

    // NULL-Schedule?
    if (dev->schedule == NULL) {
        // Yes: Assign settings from dev-struct:
        *ret_playoutbuffer = dev->outputbuffer;
        outsbsize = dev->outputbuffersize / sizeof(float);

        // Fetch boundaries of playback loop:
        loopStartFrame = dev->loopStartFrame;
        loopEndFrame = dev->loopEndFrame;
        repeatCount = dev->repeatCount;

        // Revalidate boundaries of playback loop:
        if (loopStartFrame * outchannels >= outsbsize) loopStartFrame = (outsbsize / outchannels) - 1;
        if (loopStartFrame < 0) loopStartFrame = 0;
        if (loopEndFrame * outchannels >= outsbsize) loopEndFrame = (outsbsize / outchannels) - 1;
        if (loopEndFrame < 0) loopEndFrame = 0;
        if (loopEndFrame < loopStartFrame) loopEndFrame = loopStartFrame;

        // Remap defined playback loop to "corrected" outsbsize and offset for later copy-op:
        outsbsize   = (loopEndFrame - loopStartFrame + 1) * outchannels;
        outsboffset = loopStartFrame * outchannels;

        // Compute playpositionlimit, the upper limit of played out samples from loop duration and repeatCount...
        playpositionlimit = ((psych_int64) (repeatCount * outsbsize));
        // ...and make sure it ends on integral sample frame boundaries:
        playpositionlimit -= playpositionlimit % outchannels;

        // Check if loop and repetition constraints are still valid:
        if ( !((repeatCount == -1) || (*playposition < playpositionlimit)) || (NULL == *ret_playoutbuffer)) {
            // No. Signal regular end of playback:
            return(1);
        }
    }
    else {
        // No: Real schedule:

        do {
            // Find current slot (with wraparound):
            slotid = dev->schedule_pos % dev->schedule_size;

            // Current slot valid and pending?
            if ((dev->schedule[slotid].mode & 2) == 0) {
                // No: End of schedule reached - Signal regular abort request and that's it:
                return(1);
            }

            // Current slot is valid: Assign it:
            cmd = dev->schedule[slotid].command;
            if (cmd > 0) {
                // Special command buffer: Doesn't contain sound, but some special
                // control commands. Process it, then advance to next slot...

                // This makes sure we repeat the loop, advancing to the next slot:
                *ret_playoutbuffer = NULL;
                outsbsize = 0;

                // Compute absolute deadline from given tWhen timespec and type of timespec:
                if (cmd & 4)  reqTime = dev->schedule[slotid].tWhen;                        // Absolute system time specified.
                // Relative to last requested start time. We use last true start time as fallback if the requested start time is undefined:
                if (cmd & 8)  reqTime = ((dev->reqStartTime > 0.0) ? dev->reqStartTime : dev->startTime) + dev->schedule[slotid].tWhen;
                if (cmd & 16) reqTime = dev->startTime + dev->schedule[slotid].tWhen;        // Relative to last true start time.
                // TODO: The following two "end time" related ones are pretty broken - Can't
                // work the way i want it to work, as relevant information is not available at
                // the time we'd need it, ie., at this point in execution flow...
                // Need to think about this, or scrap the idea of end-time related timeoffsets...
                if (cmd & 32) reqTime = dev->reqStopTime + dev->schedule[slotid].tWhen;        // Relative to last requested end time.
                if (cmd & 64) reqTime = dev->estStopTime + dev->schedule[slotid].tWhen;        // Relative to last true end time.

                // Pause-playback-and-restart command?
                if (cmd & 1) {
                    // Yes. We should abort playback at this point, ie. fill the remainder
                    // of the audiobuffer with silence. Then we shall switch back to hot-standby
                    // and reschedule ourselves for restart of playback at a given target reqTime:
                    dev->reqStartTime = reqTime;

                    // Manually invalidate this slot and advance schedule to next one:
                    *playposition = 0;
                    // Only disable if the flag 4 aka "don't auto-disable" isn't set:
                    if (!(dev->schedule[slotid].mode & 4)) dev->schedule[slotid].mode &= ~2;
                    dev->schedule_pos++;

                    // Return with special code 4 to reschedule:
                    return(4);
                }

                // Reschedule end-of-playback command?
                if (cmd & 2) {
                    // Yes. We should end playback at given target reqTime, then go idle:
                    dev->reqStopTime = reqTime;

                    // Nothing further to do here. Will auto-advance to next schedule slot after
                    // this one, possibly outputting further sound...
                }

                // End of command buffer processing.
            } // Regular audio buffer: First assign outbuffer size and pointer...
            else if (dev->schedule[slotid].bufferhandle <= 0) {
                // Default device playoutbuffer:
                *ret_playoutbuffer = dev->outputbuffer;
                outsbsize = dev->outputbuffersize / sizeof(float);
            }
            else
            {
                // Dynamic buffer: Dereference bufferhandle and fetch buffer data for later use:

                // Need to lock bufferList lock to do this:
                PsychLockMutex(&bufferListmutex);

                if (bufferList && (dev->schedule[slotid].bufferhandle < bufferListCount)) {
                    // Fetch pointer to actual audio data buffer:
                    *ret_playoutbuffer = bufferList[dev->schedule[slotid].bufferhandle].outputbuffer;

                    // Retrieve buffersize in samples:
                    outsbsize = bufferList[dev->schedule[slotid].bufferhandle].outputbuffersize / sizeof(float);

                    // Another child protection:
                    if (outchannels != bufferList[dev->schedule[slotid].bufferhandle].outchannels) {
                        *ret_playoutbuffer = NULL;
                        outsbsize = 0;
                    }
                }
                else {
                    *ret_playoutbuffer = NULL;
                    outsbsize = 0;
                }

                // Unlock bufferList mutex again:
                PsychUnlockMutex(&bufferListmutex);
            }

            // ... then loop and repeat parameters:
            loopStartFrame = dev->schedule[slotid].loopStartFrame;
            loopEndFrame   = dev->schedule[slotid].loopEndFrame;
            repeatCount    = dev->schedule[slotid].repetitions;

            // Revalidate boundaries of playback loop:
            if (loopStartFrame * outchannels >= outsbsize) loopStartFrame = (outsbsize / outchannels) - 1;
            if (loopStartFrame < 0) loopStartFrame = 0;
            if (loopEndFrame * outchannels >= outsbsize) loopEndFrame = (outsbsize / outchannels) - 1;
            if (loopEndFrame < 0) loopEndFrame = 0;
            if (loopEndFrame < loopStartFrame) loopEndFrame = loopStartFrame;

            // Remap defined playback loop to "corrected" outsbsize and offset for later copy-op:
            outsbsize = (psych_int64) ((loopEndFrame - loopStartFrame + 1) * outchannels);
            outsboffset = (psych_int64) (loopStartFrame * outchannels);

            // Compute playpositionlimit, the upper limit of played out samples from loop duration and repeatCount...
            playpositionlimit = ((psych_int64) (repeatCount * outsbsize));
            // ...and make sure it ends on integral sample frame boundaries:
            playpositionlimit -= playpositionlimit % outchannels;

            // Check if loop and repetition constraints as well as actual audio buffer for this slot are still valid:
            if ( !((repeatCount == -1) || (*playposition < playpositionlimit)) || (NULL == *ret_playoutbuffer) ) {
                // Constraints violated. This slot is used up: Reset playposition and advance to next slot:
                *playposition = 0;
                // Only disable if the flag 4 aka "don't auto-disable" isn't set:
                if (!(dev->schedule[slotid].mode & 4)) dev->schedule[slotid].mode &= ~2;
                dev->schedule_pos++;
            }
            else {
                // Constraints ok, break out of do-while slot advance loop:
                break;
            }
        } while(TRUE);
    }

    *ret_outsbsize = outsbsize;
    *ret_outsboffset = outsboffset;
    *ret_repeatCount = repeatCount;
    *ret_playpositionlimit = playpositionlimit;

    // Safety check: If playoutbuffer is NULL at this point, then somethings screwed
    // and we request abort of playback:
    if (NULL == *ret_playoutbuffer) return(2);

    // Return 0 exit to signal a valid update:
    return(0);
}

/* paCallback: PortAudo I/O processing callback.
 *
 * This callback is called by PortAudios playback/capture engine whenever
 * it needs new data for playback or has new data from capture. We are expected
 * to take the inputBuffer's content and store it in our own recording buffers,
 * and push data from our playback buffer queue into the outputBuffer.
 *
 * timeInfo tells us useful timing information, so we can estimate latencies,
 * compensate for them, and so on...
 *
 * This callback is part of a realtime/interrupt/system context so don't do
 * things like calling PortAudio functions, allocating memory, file i/o or
 * other unbounded operations!
 */
static int paCallback( const void *inputBuffer, void *outputBuffer,
                       unsigned long framesPerBuffer,
                       const PaStreamCallbackTimeInfo* timeInfo,
                       PaStreamCallbackFlags statusFlags,
                       void *userData )
{
    // Assign all variables, especially our dev device structure
    // with info about this stream:
    PsychPADevice* dev = (PsychPADevice*) userData;
    float *out = (float*) outputBuffer;
    float *in = (float*) inputBuffer;
    float *playoutbuffer;
    float *tmpBuffer, *mixBuffer;
    float masterVolume, neutralValue;
    psych_int64  j, k;
    psych_int64 i, silenceframes, committedFrames, max_i;
    psych_int64 inchannels, outchannels;
    psych_int64  playposition, outsbsize, insbsize, recposition;
    psych_int64  outsboffset;
    unsigned int reqstate;
    double now, firstsampleonset, onsetDelta, offsetDelta, captureStartTime;
    double repeatCount;
    psych_int64 playpositionlimit;
    PaHostApiTypeId hA;
    psych_bool stopEngine;
    psych_bool isMaster, isSlave;
    int slaveId, modulatorSlave, parc, numSlavesHandled;

    // Device struct attached to stream? If no device struct
    // is attached, we can't continue and tell the engine to abort
    // processing of this stream:
    if (dev == NULL) return(paAbort);

    // Define our role:
    isMaster = (dev->opmode & kPortAudioIsMaster) ? TRUE : FALSE;
    isSlave = (dev->opmode & kPortAudioIsSlave) ? TRUE : FALSE;

    // Query host API: Done without mutex held, as it doesn't change during device lifetime:
    hA=dev->hostAPI;

    // Only compute timestamps from raw data if we're not a slave:
    if (!isSlave) {
        // Buffer timestamp computation code:
        //
        // This is executed at each callback iteration to provide PTB timebase onset times for
        // the very first sample in this output buffer.
        //
        // The code is time-critical: To be executed as soon as possible after paCallback invocation,
        // therefore at the beginning of control-flow before trying to lock (and potentially stall at) the device mutex!
        //
        // Executing without the mutex locked is safe: It only modifies/reads local variables, and only
        // reads a few device struct variables which are all guaranteed to remain constant while the
        // engine is running.

        // Retrieve current system time:
        PsychGetAdjustedPrecisionTimerSeconds(&now);

        // FIXME: PortAudio stable sets timeInfo->currentTime == 0 --> Breakage!!!
        // That's why we currently have our own PortAudio version.

        #if PSYCH_SYSTEM == PSYCH_LINUX
        if (hA == paALSA) {
            // ALSA on Linux can return timestamps in CLOCK_MONOTONIC time
            // instead of our standard GetSecs() [gettimeofday] timesbase.
            // It is up to the audio driver and its config, which one is returned.
            //
            // This effect was first noticed during testing on Ubuntu 10.04 LTS,
            // so this will be part of all PTB releases after 8.5.2010.
            //
            // Therefore we need to check if CLOCK_MONOTONIC time is returned
            // and remap such timestamps to our standard GetSecs timebase before
            // further processing.

            // Get current CLOCK_MONOTONIC time:
            double tMonotonic = PsychOSGetLinuxMonotonicTime();

            // Returned current time timestamp closer to tMonotonic than to GetSecs time?
            if (fabs(timeInfo->currentTime - tMonotonic) < fabs(timeInfo->currentTime - now)) {
                // Timestamps are in monotonic time! Need to remap.
                // tMonotonic shall be the offset between GetSecs and monotonic time,
                // i.e., the offset that needs to be added to monotonic timestamps to
                // remap them to GetSecs time:
                tMonotonic = now - tMonotonic;

                // Correct all PortAudio timestamps by adding corrective offset:
                ((PaStreamCallbackTimeInfo*) timeInfo)->currentTime += tMonotonic;
                ((PaStreamCallbackTimeInfo*) timeInfo)->outputBufferDacTime += tMonotonic;
                ((PaStreamCallbackTimeInfo*) timeInfo)->inputBufferAdcTime += tMonotonic;
            }
        }
        #endif

        if (hA==paCoreAudio || hA==paDirectSound || hA==paMME || hA==paALSA) {
            // On these systems, DAC-time is already returned in the system timebase,
            // at least with our modified version of PortAudio, so a simple
            // query will return the onset time of the first sample. Well,
            // looks as if we need to add the device inherent latency, because
            // this describes everything up to the point where DMA transfer is
            // initiated, but not device inherent latency. This additional latency
            // is added via latencyBias, which is initialized in the Open-Function
            // via a low-level driver query to CoreAudio:
            if (dev->opmode & kPortAudioPlayBack) {
                // Playback enabled: Use DAC time as basis for timing:
                firstsampleonset = (double) timeInfo->outputBufferDacTime + dev->latencyBias;
            }
            else {
                // Recording (only): Use ADC time as basis for timing:
                firstsampleonset = (double) timeInfo->inputBufferAdcTime + dev->latencyBias;
            }

            // Compute estimated capturetime in captureStartTime. This is only important in
            // full-duplex mode, redundant in pure half-duplex capture mode:
            captureStartTime = (double) timeInfo->inputBufferAdcTime;
        }
        else {
            // Not yet verified how these other audio APIs behave. Play safe
            // and perform timebase remapping: This also needs our special fixed
            // PortAudio version where currentTime actually has a value:
            if (dev->opmode & kPortAudioPlayBack) {
                // Playback enabled: Use DAC time as basis for timing:
                // Assign predicted (remapped to our time system) audio onset time for this buffer:
                firstsampleonset = now + ((double) (timeInfo->outputBufferDacTime - timeInfo->currentTime)) + dev->latencyBias;
            }
            else {
                // Recording (only): Use ADC time as basis for timing:
                // Assign predicted (remapped to our time system) audio onset time for this buffer:
                firstsampleonset = now + ((double) (timeInfo->inputBufferAdcTime - timeInfo->currentTime)) + dev->latencyBias;
            }

            // Compute estimated capturetime in captureStartTime. This is only important in
            // full-duplex mode, redundant in pure half-duplex capture mode:
            captureStartTime = now + ((double) (timeInfo->inputBufferAdcTime - timeInfo->currentTime));
        }

        if (FALSE) {
            // Debug code to compare our two timebases against each other: On OS/X,
            // luckily both timebases are identical, ie. our UpTime() timebase used
            // everywhere in PTB and the CoreAudio AudioHostClock() are identical.
            psych_uint64 ticks;
            double  tickssec;
            PsychGetPrecisionTimerTicksPerSecond(&tickssec);
            PsychGetPrecisionTimerTicks(&ticks);
            printf("AudioHostClock: %lf   vs. System clock: %lf\n", ((double) (psych_int64) ticks) / tickssec, now);
        }

        // End of timestamp computation:
        // If audio capture is active, then captureStartTime contains the running estimate of when the first input sample
        // in inbuffer was hitting the audio input - This will be used in state 1 aka hotstandby when deciding to actually start "real" capture.
        //
        // now contains system time when entering this callback - Time in PTB timebase.
        //
        // Most importantly: firstsampleonset contains the estimated onset time (in PTB timebase) of the first sample that
        // we will store in the outputbuffer during this callback invocation. This timestamp is used in multiple places below,
        // e.g., in hotstandby mode (state == 1) to decide when to start actual playback by switching to state 2 and emitting
        // samples to the outputbuffer. Also used if a specific reqEndTime is selected, ie., a sound offset at a scheduled
        // offset time to compute when to stop. It's also used for checking for skipped buffers and other problems...
    }
    else {
        // We're a slave device: Just fetch precooked timestamps from our master:
        firstsampleonset = audiodevices[dev->pamaster].firstsampleonset;
        captureStartTime = audiodevices[dev->pamaster].cst;
        now = audiodevices[dev->pamaster].now;
    }

    // Cache cooked timestamps:
    dev->firstsampleonset = firstsampleonset;
    dev->cst = captureStartTime;
    dev->now = now;

    // Acquire device lock: We'll likely hold it until exit from paCallback:
    PsychPALockDeviceMutex(dev);

    // Cache requested state:
    reqstate = dev->reqstate;

    // Count total number of calls:
    dev->paCalls++;

    // Call number of timestamp failures:
    if (timeInfo->currentTime == 0) dev->noTime++;

    // Keep track of maximum number of frames requested:
    if (dev->batchsize < (psych_int64) framesPerBuffer) dev->batchsize = (psych_int64) framesPerBuffer;

    // Keep track of buffer over-/underflows:
    if (statusFlags & (paInputOverflow | paInputUnderflow | paOutputOverflow | paOutputUnderflow)) dev->xruns++;

    // Reset number of already committed sample frames for this buffer fill iteration to zero:
    // This is a running count of how much of the current output buffer has been filled with
    // sound content (ie., not simply zero silence padding frames):
    committedFrames = 0;

    // Reset number of silence frames so far:
    silenceframes = 0;

    // NULL-out pointer to buffer with sound data to play. It will get initialized later on
    // in PsychPAProcessSchedule():
    playoutbuffer = NULL;

    // Query number of output channels:
    outchannels = (psych_int64) dev->outchannels;

    // Query number of output channels:
    inchannels = (psych_int64) dev->inchannels;

    // Query number of repetitions:
    repeatCount = dev->repeatCount;

    // Get our current playback position in samples (not frames or bytes!!):
    playposition = dev->playposition;
    recposition = dev->recposition;

    // Compute size of soundbuffers in samples:
    outsbsize = dev->outputbuffersize / sizeof(float);
    insbsize = dev->inputbuffersize / sizeof(float);

    // Cache masterVolume locally for higher efficiency in integration loops:
    masterVolume = dev->masterVolume;

    // Assign 1 as neutral value for AM modulator slaves, as 1 is the neutral
    // element for gain-modulation (multiplication), as opposed to zero as neutral
    // element for mixing and for outputting "silence":
    neutralValue = (dev->opmode & kPortAudioIsAMModulator) ? 1.0f : 0.0f;

    // Requested logical playback state is "stopped" or "aborting" ? If so, abort.
    if ((reqstate == 0) || (reqstate == 3)) {
        // Set estimated stop time to last committed sample time, unless it is already set from
        // a previous callback invocation:
        if (dev->estStopTime == 0) dev->estStopTime = dev->currentTime;

        // Acknowledge request by resetting it:
        dev->reqstate = 255;

        // Update "true" state to inactive:
        dev->state = 0;

        // Signal state change:
        PsychPASignalChange(dev);

        // Release mutex here, because dev->runMode never changes below us, and
        // all other ops are on local variables:
        PsychPAUnlockDeviceMutex(dev);

        // Prime the outputbuffer with silence, so playback is effectively stopped:
        if (outputBuffer && !isSlave) memset(outputBuffer, 0, (size_t) (framesPerBuffer * outchannels * sizeof(float)));

        if (dev->runMode == 0) {
            // Runmode 0: We shall really stop the engine:

            // Either paComplete gracefully, playing out pending buffers, or
            // request a hard paAbort if abortion is requested:
            return((reqstate == 0) ? paComplete : paAbort);
        }
        else {
            // Runmode 1: We just "idle", ie., set our state to "idle" and return with
            // no further operations. Future invocations of the callback will also turn
            // into no-ops due to the nominal idle state:
            return(paContinue);
        }
    }

    // Are we in a nominally "idle" / inactive state?
    if (dev->state == 0) {
        // We are effectively idle, but the engine shall keep running. This is usually
        // the case in runMode > 0. Here we "simulate" idle state by simply outputting
        // silence and then returning control to PortAudio:

        // Set estimated stop time to last committed sample time, unless it is already set from
        // a previous callback invocation:
        if (dev->estStopTime == 0) dev->estStopTime = dev->currentTime;

        // Release mutex here, as memset() only operates on "local" data:
        PsychPAUnlockDeviceMutex(dev);

        // Prime the outputbuffer with silence to simulate a stopped audio device:
        if (outputBuffer && !isSlave) memset(outputBuffer, 0, (size_t) (framesPerBuffer * outchannels * sizeof(float)));

        // Done:
        return(paContinue);
    }

    // This point is only reached in hot-standby or active playback/capture/feedback modes:

    // Are we already playing back and/or capturing real audio data,
    // or are we still on hot-standby? PsychPortAudio tries to start
    // playback/capture of a sound exactly at a requested point in
    // system time (i.e. the first sound sample should hit the speaker
    // as closely as possible to that point in time) by use of the
    // following trick: When the user script executes PsychPortAudio('Start'),
    // our routine immediately starts processing of the portaudio stream,
    // starting up the audio hardwares DACs/ADCs and Portaudios engine.
    // After a short latency, our paCallback() (this routine) gets called
    // by the realtime audio scheduler, requesting or providing audio sample
    // data. In the provided timeInfo - variable, we are provided with the
    // current playback time of the audio device (in seconds since stream start)
    // and an estimate of when our first provided sample will hit the speakers.
    // We convert these timestamps into system time, so we'll know how far the user
    // provided onset deadline is away. If the deadline is far away, so the samples
    // for this callback iteration would hit the speaker too early, we simply return
    // a zero filled outputbuffer --> We output silence. If the deadline is somewhere
    // in the middle of this outputbuffer, we fill the appropriate amount of bufferspace
    // with zeros, then copy in our first real samples into the remaining buffer.
    // After that, we switch to real playback, all future calls will provide PA
    // with real sampledata. Assuming the sound onset estimate provided by PA is
    // correct, this should allow accurate sound onset. It all depends on the
    // latency estimate...

    // Hot standby?
    if (dev->state == 1) {
        // Hot standby: Compare running buffer timestamps 'now', 'firstsampleonset' and 'captureStartTime'
        // to requested sound onset time and decide if and when to start playback and capture:

        // Store our measured/updated PortAudio + HostAPI + Driver + Hardware latency:
        // We'll update the running estimate until transition from hot standby to active:
        dev->predictedLatency = firstsampleonset - now;

        // Ditto for capture starttime:
        if (dev->opmode & kPortAudioCapture) {
            // Store estimated capturetime in captureStartTime. This is only important in
            // full-duplex mode, redundant in pure half-duplex capture mode.
            // If we're not a regular input capture device, but a output capture device,
            // then our effective captureStartTime is the *playback* starttime, as we're
            // capturing what is fed to the outputs/speakers, not what comes from the inputs!
            dev->captureStartTime = (dev->opmode & kPortAudioIsOutputCapture) ? firstsampleonset : captureStartTime;
        }

        // Compute difference between requested onset time and presentation time
        // of the first sample of this callbacks returned buffer:
        onsetDelta = dev->reqStartTime - firstsampleonset;

        // Time left until onset and in playback mode?
        if ((onsetDelta > 0) && (dev->opmode & kPortAudioPlayBack)) {
            // Some time left: A full buffer duration?
            if (onsetDelta >= ((double) framesPerBuffer / (double) dev->streaminfo->sampleRate)) {
                // At least one buffer away...

                // Release mutex, as remainder only operates on locals:
                PsychPAUnlockDeviceMutex(dev);

                // At least one buffer away. Fill our buffer with zeros, aka silence:
                if (!isSlave) memset(outputBuffer, 0, (size_t) (framesPerBuffer * outchannels * sizeof(float)));

                // Ready. Tell engine to continue stream processing, i.e., call us again...
                return(paContinue);
            }
            else {
                // A bit time left, but less than a full buffer. Need to pad the head of
                // this buffer with zeros, aka silence, then fill the rest with real data:
                silenceframes = (psych_int64) (onsetDelta * ((double) dev->streaminfo->sampleRate));

                // Fill in some silence:
                if (neutralValue == 0) {
                    // Fast-path: Zerofill for silence...
                    memset(outputBuffer, 0, (size_t) (silenceframes * outchannels * sizeof(float)));
                    out+= (silenceframes * outchannels);
                }
                else {
                    // Slow-path: Usually a 1.0 fill for AM modulator mode:
                    for (i = 0; i < silenceframes * outchannels; i++) *(out++) = neutralValue;
                }

                // Decrement remaining real audio data count:
                framesPerBuffer-= (unsigned long) silenceframes;

                // Advance silenceframes into the buffer:
                committedFrames+=silenceframes;

                // dev->startTime now exactly corresponds to onset of first non-silence sample at dev->reqStartTime.
                // At least if everything works properly.
                dev->startTime = dev->reqStartTime;

                // Mark us as running:
                dev->state = 2;

                // Signal state change:
                PsychPASignalChange(dev);
            }
        }
        else {
            // Ooops! We are late! Store real estimated onset in the startTime field and
            // then hurry up! If we are in "capture only" mode, we disregard the 'when'
            // deadline and always start immediately.
            dev->startTime = firstsampleonset;

            // Mark us as running:
            dev->state = 2;

            // Signal state change:
            PsychPASignalChange(dev);
        }
    }

    // Master device implementation. Dispatch into all attached slave device callbacks,
    // distribute captured data to the slaves, collect and merge or mix output data from
    // the slaves:
    if (isMaster) {
        // Scratch buffers for slave callbacks already allocated?
        if ((dev->opmode & kPortAudioCapture) && (dev->slaveInBuffer == NULL)) {
            // Allocate input distribution buffer:
            dev->slaveInBuffer = (float*) malloc(sizeof(float) * (size_t) (dev->batchsize * inchannels));
            if (NULL == dev->slaveInBuffer) {
                // Out of memory! Perform an emergency abort:
                dev->reqstate = 255;
                dev->state = 0;
                PsychPASignalChange(dev);
                PsychPAUnlockDeviceMutex(dev);
                return(paAbort);
            }
        }

        if ((dev->opmode & kPortAudioPlayBack) && (dev->slaveOutBuffer == NULL)) {
            // Allocate output receive buffer and output gain receive buffer:
            dev->slaveOutBuffer = (float*) malloc(sizeof(float) * (size_t) (dev->batchsize * outchannels));
            dev->slaveGainBuffer = (float*) malloc(sizeof(float) * (size_t) (dev->batchsize * outchannels));
            if ((NULL == dev->slaveOutBuffer) || (NULL == dev->slaveGainBuffer)) {
                // Out of memory! Perform an emergency abort:
                dev->reqstate = 255;
                dev->state = 0;
                PsychPASignalChange(dev);
                PsychPAUnlockDeviceMutex(dev);
                return(paAbort);
            }
        }

        if (NULL != outputBuffer) {
            // Have scratch buffers ready. Clear output intermix buffer:
            memset(outputBuffer, 0, (size_t) (dev->batchsize * outchannels * sizeof(float)));
        }

        // Iterate over all slave device callbacks: Or at least until all registered slaves are handled.
        numSlavesHandled = 0;
        for (i = 0; (i < MAX_PSYCH_AUDIO_SLAVES_PER_DEVICE) && (numSlavesHandled < dev->slaveCount); i++) {
            // Valid slave slot?
            slaveId = dev->slaves[i];

            // We skip invalid slots and output capturer slaves:
            if ((slaveId > -1) && !(audiodevices[slaveId].opmode & kPortAudioIsOutputCapture)) {
                // Valid slave:

                // Is this device an AM modulator attached to a slave?
                if (audiodevices[slaveId].opmode & kPortAudioIsAMModulatorForSlave) {
                    // Yes. Mark it as handled, then skip it. It will be called as part
                    // of processing of its parent slave:
                    numSlavesHandled++;
                    continue;
                }

                // This is a "real" audio slave, not a modulator or such:

                // Gain modulator slave for this real slave attached, valid and active?
                // If this is the case, we need to unconditionally execute it here, regardless
                // of what the actual 'slaveId' device is up to. Otherwise we can run into
                // time sync issues and ugly deadlocks in the calling code:
                modulatorSlave = audiodevices[slaveId].modulatorSlave;
                if ((modulatorSlave > -1) && (audiodevices[modulatorSlave].stream) &&
                    (audiodevices[modulatorSlave].opmode & kPortAudioIsAMModulatorForSlave) && (audiodevices[modulatorSlave].state > 0)) {
                    // Yes. Execute it:
                    audiodevices[modulatorSlave].slaveDirty = 0;

                // Prefill buffer with neutral 1.0:
                tmpBuffer = dev->slaveGainBuffer;
                for (j = 0; j < dev->batchsize * audiodevices[modulatorSlave].outchannels; j++) *(tmpBuffer++) = 1.0;

                // This will potentially fill the slaveGainBuffer with gain modulation values.
                // The passed slaveInBuffer is meaningless for a modulator slave and only contains random junk...
                paCallback( (const void*) dev->slaveInBuffer, (void*) dev->slaveGainBuffer, (unsigned long) dev->batchsize, timeInfo, statusFlags, (void*) &(audiodevices[modulatorSlave]));
                    }
                    else {
                        // No. Either no modulator slave or slave not currently active. Signal this
                        // by setting modulatorSlave to a -1 value:
                        modulatorSlave = -1;
                    }

                    // Skip actual slaves processing if its state is zero == completely inactive.
                    if (audiodevices[slaveId].state > 0) {
                        // Slave is active, need to process it:

                        // Reset dirty flag for this slave:
                        audiodevices[slaveId].slaveDirty = 0;

                        // Is this a playback slave?
                        if (audiodevices[slaveId].opmode & kPortAudioPlayBack) {
                            // Prefill slaves output buffer with 1.0, a neutral gain value for playback slaves
                            // without a AM modulator attached. The same prefill is needed with AM modulator,
                            // this time to make the modulator itself happy:
                            tmpBuffer = dev->slaveOutBuffer;
                            for (j = 0; j < dev->batchsize * audiodevices[slaveId].outchannels; j++) *(tmpBuffer++) = 1.0;

                            // Ok, the outbuffer is filled with a neutral 1.0 gain value. This will work
                            // even if no per-slave gain modulation is provided by a modulator slave.

                            // Is a modulator slave active and did it write any gain AM values?
                            if ((modulatorSlave > -1) && (audiodevices[modulatorSlave].slaveDirty)) {
                                // Yes. Need to distribute them to proper channels in slaveOutBuffer:
                                tmpBuffer = dev->slaveGainBuffer;
                                mixBuffer = dev->slaveOutBuffer;
                                for (j = 0; j < dev->batchsize; j++) {
                                    // Iterate over all target channels in the slave device outputbuffer:
                                    for (k = 0; k < audiodevices[modulatorSlave].outchannels; k++) {
                                        // Modulate current sample in intermixbuffer via multiplication:
                                        mixBuffer[(j * audiodevices[slaveId].outchannels) + audiodevices[modulatorSlave].outputmappings[k]] = *(tmpBuffer++) * audiodevices[modulatorSlave].outChannelVolumes[k];
                                    }
                                }
                            }
                        }    // Ok, the slaveOutBuffer for this playback slave is prefilled with valid gain modulation data to apply to the actual sound output.

                        // Capture enabled on slave? If so, we need to distribute our captured audio data to it:
                        if (audiodevices[slaveId].opmode & kPortAudioCapture) {
                            tmpBuffer = dev->slaveInBuffer;
                            // For each sampleFrame in the input buffer:
                            for (j = 0; j < dev->batchsize; j++) {
                                // Iterate over all target channels in the slave devices inputbuffer:
                                for (k = 0; k < audiodevices[slaveId].inchannels; k++) {
                                    // And fetch from corrsponding source channel of our device:
                                    *(tmpBuffer++) = in[(j * inchannels) + audiodevices[slaveId].inputmappings[k]];
                                }
                            }
                        }

                        // Temporary input buffer is filled for slave callback: Execute it.
                        paCallback( (const void*) dev->slaveInBuffer, (void*) dev->slaveOutBuffer, (unsigned long) dev->batchsize, timeInfo, statusFlags, (void*) &(audiodevices[slaveId]));

                        // Check if the paCallback actually filled anything into the dev->slaveOutBuffer:
                        if ((audiodevices[slaveId].opmode & kPortAudioPlayBack) && audiodevices[slaveId].slaveDirty) {
                            // Slave has written meaningful data to its output buffer. Merge & mix it:

                            // Process from first non-silence sample slot (after silenceframes prefix) until end of buffer:
                            tmpBuffer = &(dev->slaveOutBuffer[committedFrames * audiodevices[slaveId].outchannels]);
                            mixBuffer = (float*) outputBuffer;

                            // Special AM-Modulator slave?
                            if (audiodevices[slaveId].opmode & kPortAudioIsAMModulator) {
                                // Yes: This slave doesn't provide audio data for mixing, but instead
                                // a time-series of gain modulation samples for amplitude modulation.
                                // Multiply the master channels samples with the slaves "gain samples"
                                // to apply AM modulation:
                                for (j = committedFrames; j < dev->batchsize; j++) {
                                    // Iterate over all target channels in the slave device outputbuffer:
                                    for (k = 0; k < audiodevices[slaveId].outchannels; k++) {
                                        // Modulate current sample in intermixbuffer via multiplication:
                                        mixBuffer[(j * outchannels) + audiodevices[slaveId].outputmappings[k]] *= *(tmpBuffer++) * audiodevices[slaveId].outChannelVolumes[k];
                                    }
                                }
                            }
                            else {
                                // Regular mix: Mix all output channels of the slave into the proper target channels
                                // of the master by simple addition. Apply per-channel volume settings of the slave
                                // during mix:
                                for (j = committedFrames; j < dev->batchsize; j++) {
                                    // Iterate over all target channels in the slave device outputbuffer:
                                    for (k = 0; k < audiodevices[slaveId].outchannels; k++) {
                                        // Mix current sample via addition:
                                        mixBuffer[(j * outchannels) + audiodevices[slaveId].outputmappings[k]] += *(tmpBuffer++) * audiodevices[slaveId].outChannelVolumes[k];
                                    }
                                }
                            }
                        }
                    }

                    // One more slave handled:
                    numSlavesHandled++;
            }
        }    // Next slave...

        // Done merging sound data from slaves. Mastercode can now process special output capture slaves
        // and other special post-mix slaves:
        for (i = 0; (i < MAX_PSYCH_AUDIO_SLAVES_PER_DEVICE) && (numSlavesHandled < dev->slaveCount); i++) {
            // Valid slave slot?
            slaveId = dev->slaves[i];

            // We skip invalid slots and anything except output capturer slaves:
            if ((slaveId > -1) && (audiodevices[slaveId].opmode & kPortAudioIsOutputCapture)) {
                // Valid slave: Skip its processing if its state is zero == completely inactive.
                if (audiodevices[slaveId].state > 0) {
                    // Reset dirty flag for this slave: Not strictly needed for output capture slaves...
                    audiodevices[slaveId].slaveDirty = 0;

                    // Output capture enabled on slave? If so, we need to distribute our output audio data to it:
                    if ((audiodevices[slaveId].opmode & kPortAudioCapture) && (audiodevices[slaveId].opmode & kPortAudioIsOutputCapture)) {
                        // Our target buffer is the slaveOutBuffer here, because it is guaranteed to exist and
                        // have sufficient capacity:
                        tmpBuffer = dev->slaveOutBuffer;

                        // Our input is the mixBuffer from previous mixes:
                        mixBuffer = (float*) outputBuffer;

                        // For each sampleFrame in the mixBuffer:
                        for (j = 0; j < dev->batchsize; j++) {
                            // Iterate over all target channels in the slave devices inputbuffer:
                            for (k = 0; k < audiodevices[slaveId].inchannels; k++) {
                                // And fetch from corrsponding mixBuffer channel of our device, applying the same
                                // masterVolume setting that the master output device will apply later:
                                *(tmpBuffer++) = masterVolume * mixBuffer[(j * outchannels) + audiodevices[slaveId].inputmappings[k]];
                            }
                        }
                    }

                    // Temporary input buffer is filled for slave callback: dev->slaveOutBuffer acts as the input
                    // buffer for the callback, as it contains the properly formatted data from our mixBuffer. It
                    // also acts as output buffer from the slave, so that slaves with playback enabled are happy.
                    // However, this is a pure dummy-sink, the data written by the slave isn't used for anything,
                    // but simply discarded, as output capture slaves have no meaningful sink for their output.
                    paCallback( (const void*) dev->slaveOutBuffer, (void*) dev->slaveOutBuffer, (unsigned long) dev->batchsize, timeInfo, statusFlags, (void*) &(audiodevices[slaveId]));
                }

                // One more slave handled:
                numSlavesHandled++;
            }
        }    // Next outputcapture-slave...
        // Done processing slaves. Mastercode can now continue just as slave or regular device code would do.
    }

    // This code only executes on non-masters in live monitoring mode - a mode where all
    // sound data is fed back immediately with shortest possible latency
    // from input to output, without any involvement of Matlab/Octave code:
    // Note: "Non-master" also means: A standard sound device without any slaves!
    if (!isMaster && (dev->opmode & kPortAudioMonitoring)) {
        // Copy input buffer to output buffer:
        // We need to offset our copy by committedFrames into the "in"putbuffer.
        // Reason: committedFrames output buffer frames have been filled with
        // silence already, and framesPerBuffer has been decremented by these
        // committedFrames - We can only copy this reduced amount from the input
        // buffer to the output buffer. We choose the most recent 'framesPerBuffer'
        // frames from the inputbuffer, ie., at offset committedFrames, because the
        // freshest frames are towards the end of the buffer, not at the beginning.
        // (Yes this may twist your mind, but it makes sense, given the way the buffers
        // are filled during capture and emptied during playback, believe me!)
        memcpy(out, &(in[committedFrames * inchannels]), (size_t) (framesPerBuffer * outchannels * sizeof(float)));

        // Store updated positions in device structure:
        dev->playposition = playposition + (framesPerBuffer * outchannels);
        dev->recposition = dev->playposition;

        // Compute output time of last fed back sample from this iteration:
        committedFrames += framesPerBuffer;
        dev->currentTime = firstsampleonset + ((double) committedFrames / (double) dev->streaminfo->sampleRate);

        // Mark our output buffer as dirty:
        dev->slaveDirty = 1;

        // Return from callback:
        PsychPAUnlockDeviceMutex(dev);
        return(paContinue);
    }

    // This code retrieves and stores captured sound data on non-master devices, if any:
    if (!isMaster && (dev->opmode & kPortAudioCapture)) {
        // Check if required input buffer is there.
        if (dev->inputbuffer == NULL) {
            // Ouch! Perform emergency shutdown:
            dev->reqstate = 255;
            dev->state = 0;

            PsychPASignalChange(dev);
            PsychPAUnlockDeviceMutex(dev);

            return(paAbort);
        }

        // This is the simple case (compared to playback processing).
        // Just copy all available data to our internal buffer:
        for (i=0; (i < dev->batchsize * inchannels); i++) {
            dev->inputbuffer[recposition % insbsize] = (float) *in++;
            recposition++;
        }

        // Store updated recording position in device structure:
        dev->recposition = recposition;
    }

    // This code emits actual sound data to the engine:
    if (dev->opmode & kPortAudioPlayBack) {
        // Set stopEngine request to false, ie. continue...
        stopEngine = FALSE;
        parc = 0;

        // Mark our output buffer as dirty:
        dev->slaveDirty = 1;

        // Last chance to honor a potential playback abort request. Check once more...
        reqstate = dev->reqstate;

        // Compute time delta between requested sound stop time (sound offset time) and
        // time of next sample in to-be-filled buffer, taking potential committedFrames for
        // zero padding at head of buffer into account:
        offsetDelta = dev->reqStopTime - (firstsampleonset + ((double) committedFrames / (double) dev->streaminfo->sampleRate));

        // Clamp to at most 10 seconds ahead, because that is more than enough even
        // for the largest conceivable hostbuffersizes, and it prevents numeric overflow
        // in the math below when converting to psych_int64 ints:
        offsetDelta = (offsetDelta > 10.0) ? 10.0 : offsetDelta;

        // Convert remaining time until requested stop time into sample frames until stop:
        offsetDelta = offsetDelta * (double)(dev->streaminfo->sampleRate);

        // Convert into samples: max_i is the maximum allowable value for 'i'
        // in order to satisfy the dev->reqStopTime:
        max_i  = (psych_int64) (offsetDelta * (double) outchannels);
        max_i -= max_i % outchannels;

        // Count of outputted frames in this part of the code:
        i=0;

        // Stoptime already reached or abort request from master thread received? If so, stop the engine:
        if (reqstate == 0 || reqstate == 3 || (offsetDelta <= 0) ) stopEngine = TRUE;

        // Repeat until stopEngine condition, or this callbacks host output buffer is full,
        // or max_i timeout reached for end of processing, or no more valid slots available
        // in current schedule. Assign all relevant parameters from schedule:
        while (!stopEngine && (i < framesPerBuffer * outchannels) && (i < max_i) &&
            ((parc = PsychPAProcessSchedule(dev, &playposition, &playoutbuffer, &outsbsize, &outsboffset, &repeatCount, &playpositionlimit)) == 0)) {
            // Process this slot:

            if (!isMaster && !isSlave) {
                // Non-master, non-slave device: This is a regular sound device.
                // Copy requested number of samples for each channel into the output buffer: Take the case of
                // "loop forever" and "loop repeatCount" times into account, as well as stop times:
                for (; (i < framesPerBuffer * outchannels) && (i < max_i) && ((repeatCount == -1) || (playposition < playpositionlimit)); i++) {
                    *(out++) = playoutbuffer[outsboffset + ( playposition % outsbsize )] * masterVolume;
                    playposition++;
                }
            }
            else if (!isMaster) {
                // Non-master device: This is a slave.
                // Copy requested number of samples for each channel into the output buffer: Take the case of
                // "loop forever" and "loop repeatCount" times into account, as well as stop times:
                for (; (i < framesPerBuffer * outchannels) && (i < max_i) && ((repeatCount == -1) || (playposition < playpositionlimit)); i++) {
                    // We multiply in order to apply possible per-channel, per-sample gain values as
                    // defined by the master - i.e., by an AM modulator that is attached to us:
                    *(out++) *= playoutbuffer[outsboffset + ( playposition % outsbsize )] * masterVolume;
                    playposition++;
                }
            }
            else {
                // Master device: We don't output our own audio data. Just apply the masterVolume
                // gain setting common to all output channels of the device:
                for (; (i < framesPerBuffer * outchannels) && (i < max_i) && ((repeatCount == -1) || (playposition < playpositionlimit)); i++) {
                    *(out++) *= masterVolume;
                    playposition++;
                }
            }

            // Store updated playposition in device structure:
            dev->playposition = playposition;

            // Abort condition?
            if ((i >= max_i) || ((parc = PsychPAProcessSchedule(dev, &playposition, &playoutbuffer, &outsbsize, &outsboffset, &repeatCount, &playpositionlimit)) > 0)) stopEngine = TRUE;
            }

            // Store updated playposition in device structure:
            dev->playposition = playposition;

            // Update total count of emitted sample frames during this callback by number of non-silence frames:
            committedFrames += i / outchannels;

            // Compute output time of last outputted sample from this iteration:
            dev->currentTime = firstsampleonset + ((double) committedFrames / (double) dev->streaminfo->sampleRate);

            // Update total count of emitted samples since start of playback:
            dev->totalplaycount+= (committedFrames - silenceframes) * outchannels;

            // Another end-of-playback check:
            if (parc > 0) {
                stopEngine = TRUE;
            }

            // End of playback reached due to maximum number of possible output samples reached or abortcondition satisfied?
            if ((i < framesPerBuffer * outchannels) || stopEngine) {
                // Premature stop of buffer filling because abortcondition satisfied: Need to go idle or stop the whole engine:

                // We need to zero-fill the remainder of the buffer and tell the engine
                // to finish playback:
                while(i < framesPerBuffer * outchannels) {
                    *out++ = neutralValue;
                    i++;
                }

                // Signal that engine is stopped/will stop very soonish:
                // Unless parc == 4 request a rescheduled restart, ie., switching to hot-standby
                // instead of idle:
                dev->state = (parc == 4) ? 1 : 0;
                dev->reqstate = 255;

                // Set estimated stop time to last committed sample time:
                if ((dev->estStopTime == 0) && (dev->state == 0)) dev->estStopTime = dev->currentTime;

                // Signal state change:
                PsychPASignalChange(dev);

                // Safe to unlock, as dev->runMode never changes below us:
                PsychPAUnlockDeviceMutex(dev);

                if ((dev->runMode == 0) && (dev->state == 0)) {
                    // Either paComplete gracefully, playing out pending buffers, or
                    // request a hard paAbort if abortion is requested:
                    return((reqstate == 3) ? paAbort : paComplete);
                }
                else {
                    // Runmode 1 or greater: Have set our state to "idle" already. Just
                    // continue:
                    return(paContinue);
                }
            }
    }

    // Tell engine to continue stream processing, i.e., call us again...
    PsychPAUnlockDeviceMutex(dev);
    return(paContinue);
}

void PsychPACloseStream(int id)
{
    int pamaster, i;
    PaStream* stream = audiodevices[id].stream;

    // Valid and active device?
    if (stream) {
        // Need different destruction procedures for normals vs. masters vs. slaves:
        if (audiodevices[id].opmode & kPortAudioIsSlave) {
            // Virtual slave device.
            pamaster = audiodevices[id].pamaster;

            // First we need to get the lock of the masterdevice to prevent it
            // from operating on its list of slaves or from calling the to-be-deleted
            // slave:
            PsychPALockDeviceMutex(&(audiodevices[pamaster]));

            // At this point we have the mutex and neither the master nor our slave are operating.
            // We also know that nobody is waiting on any finalization or state-change signals from
            // the slave device, because the only thread that could wait on such events is the
            // main interpreter-thread on which we are executing at this moment -- obviously we're
            // not blocked. The slave is not really operating on hardware that needs to be shut down
            // and the slave is not operating on its own thread, so no thread cleanup to do either.
            //
            // Simply put, we can simply detach ourselves from the master, then skip destruction steps
            // for real audio devices and continue with release operations for our data structures,
            // buffers and sync primitives.

            // Find our slot in the master:
            for (i=0; (i < MAX_PSYCH_AUDIO_SLAVES_PER_DEVICE) && (audiodevices[pamaster].slaves[i] != id); i++);

            // Detach us from master queue:
            audiodevices[pamaster].slaves[i] = -1;
            audiodevices[pamaster].slaveCount--;

            // Detach master from us:
            audiodevices[id].pamaster = -1;

            // If we're an AM modulator for another slave, we need to detach
            // from that slave:
            if (audiodevices[id].opmode & kPortAudioIsAMModulatorForSlave) {
                for (i = 0; i < MAX_PSYCH_AUDIO_DEVS; i++) if (audiodevices[i].modulatorSlave == id) { audiodevices[i].modulatorSlave = -1; }
            }

            // Detached :-) -- Master can continue with whatever...
            PsychPAUnlockDeviceMutex(&audiodevices[pamaster]);

            // Continue with data structure destruction etc.
        }
        else {
            // Regular or master audio device: Real hardware with real need for
            // Portaudio shutdown.

            // Stop, shutdown and release audio stream:
            Pa_StopStream(stream);

            // Unregister the stream finished callback:
            Pa_SetStreamFinishedCallback(stream, NULL);

            // Our device thread, callbacks and hardware are stopped, all mutexes are unlocked,
            // all our potential slaves are inactive as well. We can safely destroy our slaves,
            // if any, ie., if we are a master:
            if ((audiodevices[id].slaveCount > 0) && (audiodevices[id].slaves)) {
                // Iterate over our slaves and close them:
                for (i=0; i < MAX_PSYCH_AUDIO_SLAVES_PER_DEVICE; i++) {
                    if (audiodevices[id].slaves[i] > -1) PsychPACloseStream(audiodevices[id].slaves[i]);
                }

                // Sanity-check:
                if (audiodevices[id].slaveCount > 0) printf("PsychPortAudio: CRITICAL! To be deleted master device %i has non-zero slaveCount %i after destroying slaves! BUG!!\n", id, audiodevices[id].slaveCount);

                // No slaves anymore. Done with master specific destruction.
            }

            // Destruction for both master- and regular audio devices:

            // Close and destroy the hardware portaudio stream:
            Pa_CloseStream(stream);
        }

        // Common destruct path for all types of devices:

        // Release stream reference to now dead stream:
        audiodevices[id].stream = NULL;

        // Free associated sound outputbuffer:
        if(audiodevices[id].outputbuffer) {
            free(audiodevices[id].outputbuffer);
            audiodevices[id].outputbuffer = NULL;
            audiodevices[id].outputbuffersize = 0;
        }

        // Free associated sound inputbuffer:
        if(audiodevices[id].inputbuffer) {
            free(audiodevices[id].inputbuffer);
            audiodevices[id].inputbuffer = NULL;
            audiodevices[id].inputbuffersize = 0;
        }

        // Free associated schedule, if any:
        if(audiodevices[id].schedule) {
            free(audiodevices[id].schedule);
            audiodevices[id].schedule = NULL;
            audiodevices[id].schedule_size = 0;
        }

        // Free associated sound intermixbuffers:
        if(audiodevices[id].slaveOutBuffer) {
            free(audiodevices[id].slaveOutBuffer);
            audiodevices[id].slaveOutBuffer = NULL;
        }

        if(audiodevices[id].slaveGainBuffer) {
            free(audiodevices[id].slaveGainBuffer);
            audiodevices[id].slaveGainBuffer = NULL;
        }

        if(audiodevices[id].slaveInBuffer) {
            free(audiodevices[id].slaveInBuffer);
            audiodevices[id].slaveInBuffer = NULL;
        }

        // Free slave array:
        if(audiodevices[id].slaves) {
            free(audiodevices[id].slaves);
            audiodevices[id].slaves = NULL;
        }

        // Free in-/outputmapping:
        if(audiodevices[id].inputmappings) {
            free(audiodevices[id].inputmappings);
            audiodevices[id].inputmappings = NULL;
        }

        if(audiodevices[id].outputmappings) {
            free(audiodevices[id].outputmappings);
            audiodevices[id].outputmappings = NULL;
        }

        // Free vector of outChannelVolumes:
        if(audiodevices[id].outChannelVolumes) {
            free(audiodevices[id].outChannelVolumes);
            audiodevices[id].outChannelVolumes = NULL;
        }

        // If we use locking, we need to destroy the per-device mutex:
        if (uselocking && PsychDestroyMutex(&(audiodevices[id].mutex))) printf("PsychPortAudio: CRITICAL! Failed to release Mutex object for pahandle %i! Prepare for trouble!\n", id);

        // If we use locking, this will destroy the associated event variable:
        PsychPADestroySignal(&(audiodevices[id]));

        // One device less:
        if (verbosity > 4) {
            printf("PTB-INFO: Closing handle %i.\n", id);
        }
        audiodevicecount--;
    }

    return;
}

const char** InitializeSynopsis(void)
{
    int i=0;
    const char **synopsis = synopsisSYNOPSIS;  //abbreviate the long name

    synopsis[i++] = "PsychPortAudio - A sound driver built around the PortAudio sound library:\n";
    synopsis[i++] = "\nGeneral information:\n";
    synopsis[i++] = "version = PsychPortAudio('Version');";
    synopsis[i++] = "oldlevel = PsychPortAudio('Verbosity' [,level]);";
    synopsis[i++] = "count = PsychPortAudio('GetOpenDeviceCount');";
    synopsis[i++] = "devices = PsychPortAudio('GetDevices' [,devicetype] [, deviceIndex]);";
    synopsis[i++] = "\nGeneral settings:\n";
    synopsis[i++] = "[oldyieldInterval, oldMutexEnable, lockToCore1, audioserver_autosuspend] = PsychPortAudio('EngineTunables' [, yieldInterval] [, MutexEnable] [, lockToCore1] [, audioserver_autosuspend]);";
    synopsis[i++] = "oldRunMode = PsychPortAudio('RunMode', pahandle [,runMode]);";
    synopsis[i++] = "\n\nDevice setup and shutdown:\n";
    synopsis[i++] = "pahandle = PsychPortAudio('Open' [, deviceid][, mode][, reqlatencyclass][, freq][, channels][, buffersize][, suggestedLatency][, selectchannels][, specialFlags=0]);";
    synopsis[i++] = "pahandle = PsychPortAudio('OpenSlave', pamaster [, mode][, channels][, selectchannels]);";
    synopsis[i++] = "PsychPortAudio('Close' [, pahandle]);";
    synopsis[i++] = "oldOpMode = PsychPortAudio('SetOpMode', pahandle [, opModeOverride]);";
    synopsis[i++] = "oldbias = PsychPortAudio('LatencyBias', pahandle [,biasSecs]);";
    synopsis[i++] = "[oldMasterVolume, oldChannelVolumes] = PsychPortAudio('Volume', pahandle [, masterVolume][, channelVolumes]);";
    #if (PSYCH_SYSTEM == PSYCH_OSX) && !defined(paMacCoreChangeDeviceParameters)
    synopsis[i++] = "enable = PsychPortAudio('DirectInputMonitoring', pahandle, enable [, inputChannel = -1][, outputChannel = 0][, gainLevel = 0.0][, stereoPan = 0.5]);";
    #endif
    synopsis[i++] = "[underflow, nextSampleStartIndex, nextSampleETASecs] = PsychPortAudio('FillBuffer', pahandle, bufferdata [, streamingrefill=0][, startIndex=Append]);";
    synopsis[i++] = "bufferhandle = PsychPortAudio('CreateBuffer' [, pahandle], bufferdata);";
    synopsis[i++] = "PsychPortAudio('DeleteBuffer'[, bufferhandle] [, waitmode]);";
    synopsis[i++] = "PsychPortAudio('RefillBuffer', pahandle [, bufferhandle=0], bufferdata [, startIndex=0]);";
    synopsis[i++] = "PsychPortAudio('SetLoop', pahandle[, startSample=0][, endSample=max][, UnitIsSeconds=0]);";
    synopsis[i++] = "startTime = PsychPortAudio('Start', pahandle [, repetitions=1] [, when=0] [, waitForStart=0] [, stopTime=inf] [, resume=0]);";
    synopsis[i++] = "startTime = PsychPortAudio('RescheduleStart', pahandle, when [, waitForStart=0] [, repetitions] [, stopTime]);";
    synopsis[i++] = "status = PsychPortAudio('GetStatus' pahandle);";
    #if PSYCH_LANGUAGE == PSYCH_MATLAB
    synopsis[i++] = "[audiodata absrecposition overflow cstarttime] = PsychPortAudio('GetAudioData', pahandle [, amountToAllocateSecs][, minimumAmountToReturnSecs][, maximumAmountToReturnSecs][, singleType=0]);";
    #else
    synopsis[i++] = "[audiodata absrecposition overflow cstarttime] = PsychPortAudio('GetAudioData', pahandle [, amountToAllocateSecs][, minimumAmountToReturnSecs][, maximumAmountToReturnSecs][, singleType=1]);";
    #endif
    synopsis[i++] = "[startTime endPositionSecs xruns estStopTime] = PsychPortAudio('Stop', pahandle [,waitForEndOfPlayback=0] [, blockUntilStopped=1] [, repetitions] [, stopTime]);";
    synopsis[i++] = "PsychPortAudio('UseSchedule', pahandle, enableSchedule [, maxSize = 128]);";
    synopsis[i++] = "[success, freeslots] = PsychPortAudio('AddToSchedule', pahandle [, bufferHandle=0][, repetitions=1][, startSample=0][, endSample=max][, UnitIsSeconds=0][, specialFlags=0]);";

    synopsis[i++] = NULL;  //this tells PsychDisplayScreenSynopsis where to stop
    if (i > MAX_SYNOPSIS_STRINGS) {
        PrintfExit("%s: increase dimension of synopsis[] from %ld to at least %ld and recompile.",__FILE__,(long)MAX_SYNOPSIS_STRINGS,(long)i);
    }

    return(synopsisSYNOPSIS);
}

PaHostApiIndex PsychPAGetLowestLatencyHostAPI(void)
{
    PaHostApiIndex ai;

    #if PSYCH_SYSTEM == PSYCH_OSX
    // CoreAudio or nothing ;-)
    return(Pa_HostApiTypeIdToHostApiIndex(paCoreAudio));
    #endif

    #if PSYCH_SYSTEM == PSYCH_LINUX
    // Try ALSA first...
    if (((ai=Pa_HostApiTypeIdToHostApiIndex(paALSA))!=paHostApiNotFound) && (Pa_GetHostApiInfo(ai)->deviceCount > 0)) return(ai);
    // Then JACK...
    if (((ai=Pa_HostApiTypeIdToHostApiIndex(paJACK))!=paHostApiNotFound) && (Pa_GetHostApiInfo(ai)->deviceCount > 0)) return(ai);
    // Then OSS...
    if (((ai=Pa_HostApiTypeIdToHostApiIndex(paOSS))!=paHostApiNotFound) && (Pa_GetHostApiInfo(ai)->deviceCount > 0)) return(ai);
    // then give up!
    printf("PTB-ERROR: Could not find an operational audio subsystem on this Linux machine! Soundcard and driver installed and enabled?!?\n");
    return(paHostApiNotFound);
    #endif

    #if PSYCH_SYSTEM == PSYCH_WINDOWS
    #ifdef PTB_USE_ASIO
    // Try ASIO first. It's supposed to be the lowest latency Windows API on soundcards that suppport it.
    if (((ai=Pa_HostApiTypeIdToHostApiIndex(paASIO))!=paHostApiNotFound) && (Pa_GetHostApiInfo(ai)->deviceCount > 0)) return(ai);
    #endif
    // Then Vistas new WASAPI, which is supposed to be usable since around Windows-7 and pretty good since Windows-10:
    if (((ai=Pa_HostApiTypeIdToHostApiIndex(paWASAPI))!=paHostApiNotFound) && (Pa_GetHostApiInfo(ai)->deviceCount > 0)) return(ai);
    // Then WDM kernel streaming Win2000 and later. This is the best builtin sound system we get on pre Windows-7:
    if (((ai=Pa_HostApiTypeIdToHostApiIndex(paWDMKS))!=paHostApiNotFound) && (Pa_GetHostApiInfo(ai)->deviceCount > 0)) return(ai);
    // Then DirectSound: Bad, but not a complete disaster if the sound card has DS native drivers:
    if (((ai=Pa_HostApiTypeIdToHostApiIndex(paDirectSound))!=paHostApiNotFound) && (Pa_GetHostApiInfo(ai)->deviceCount > 0)) return(ai);
    // Then Windows MME, a complete disaster, but better than silence...?!?
    if (((ai=Pa_HostApiTypeIdToHostApiIndex(paMME))!=paHostApiNotFound) && (Pa_GetHostApiInfo(ai)->deviceCount > 0)) return(ai);
    // then give up!
    printf("PTB-ERROR: Could not find an operational audio subsystem on this Windows machine! Soundcard and driver installed and enabled?!?\n");
    return(paHostApiNotFound);
    #endif

    printf("PTB-FATAL-ERROR: Impossible point in code execution reached! (End of PsychPAGetLowestLatencyHostAPI()\n");
    return(paHostApiNotFound);
}

PsychError PSYCHPORTAUDIODisplaySynopsis(void)
{
    int i;

    for (i = 0; synopsisSYNOPSIS[i] != NULL; i++)
        printf("%s\n",synopsisSYNOPSIS[i]);

    return(PsychError_none);
}

// Module exit function: Stop and close all audio streams, terminate PortAudio...
PsychError PsychPortAudioExit(void)
{
    PaError err;
    int i;

    if (pa_initialized) {
        for(i=0; i<MAX_PSYCH_AUDIO_DEVS; i++) {
            // Close i'th stream, if it is open:
            PsychPACloseStream(i);
        }
        audiodevicecount = 0;

        // Delete all audio buffers and the bufferlist itself:
        PsychPADeleteAllAudioBuffers();

        // Release audiobufferlist mutex lock:
        PsychDestroyMutex(&bufferListmutex);

        // Shutdown PortAudio itself:
        err = Pa_Terminate();
        if (err) {
            printf("PTB-FATAL-ERROR: PsychPortAudio: Shutdown of PortAudio subsystem failed. Depending on the quality\n");
            printf("PTB-FATAL-ERROR: of your operating system, this may leave the sound system of your machine dead or confused.\n");
            printf("PTB-FATAL-ERROR: Exit and restart Matlab/Octave/Python. Windows users additionally may want to reboot...\n");
            printf("PTB-FATAL-ERRRO: PortAudio reported the following error: %s\n\n", Pa_GetErrorText(err));
        }
        else {
            pa_initialized = FALSE;
        }

        // Detach our callback function for low-level debug output:
        PaUtil_SetDebugPrintFunction(NULL);

        #if PSYCH_SYSTEM == PSYCH_LINUX
            // Disable ALSA error handler:
            snd_lib_error_set_handler(NULL);
            if (myjack_set_error_function) {
                myjack_set_error_function(NULL);
                myjack_set_error_function = NULL;
            }
        #endif

        // Restart suspended PulseAudio server if it was suspended by us:
        if (pulseaudio_isSuspended) {
            int rc = 0;

            // Call external "pactl" utility via shell to ask
            // the server to resume all sinks and sources.
            // These calls should fail silently if either the
            // pactl utility isn't installed or no sound
            // server is running:
            if (verbosity > 4) printf("PTB-INFO: Trying to resume potentially suspended PulseAudio server... ");
            rc += system("pactl suspend-sink 0");
            rc += system("pactl suspend-source 0");
            if (verbosity > 4) printf("... status %i\n", rc);

            pulseaudio_isSuspended = FALSE;
        }
    }

    return(PsychError_none);
}

void PsychPortAudioInitialize(void)
{
    PaError err;
    int i;

    // PortAudio already initialized?
    if (!pa_initialized) {
        // No - Initialize:

        // Suspend possibly running PulseAudio sound server if requested.
        // We do this before PortAudio initialization, so the PulseAudio
        // server can't interfere with device enumeration as done during
        // PortAudio init. This is a Linux only thing...
        if (pulseaudio_autosuspend && (PSYCH_SYSTEM == PSYCH_LINUX)) {
            int rc = 0;

            // Call external "pactl" utility via runtime to ask
            // the server to suspend all sinks and sources,
            // hopefully releasing the underlying low-level audio
            // devices for use by use. These calls should fail silently
            // if either the pactl utility isn't installed or no sound
            // server is running:
            if (verbosity > 4) printf("PTB-INFO: Trying to suspend potentially running PulseAudio server... ");
            rc += system("pactl suspend-sink 1");
            rc += system("pactl suspend-source 1");
            if (verbosity > 4) printf("... status %i\n", rc);

            pulseaudio_isSuspended = TRUE;
        }

        #if PSYCH_SYSTEM == PSYCH_WINDOWS
        // Sanity check dynamic portaudio dll loading on Windows:
        if ((NULL == LoadLibrary("portaudio_x86.dll")) && (NULL == LoadLibrary("portaudio_x64.dll"))) {
            // Failed:
            printf("\n\nPTB-ERROR: Tried to initialize PsychPortAudio's PortAudio engine. This didn't work,\n");
            printf("PTB-ERROR: because i couldn't find or load the required portaudio_x86.dll or portaudio_x64.dll library.\n");
            printf("PTB-ERROR: Please make sure to call the InitializePsychSound function before first use of\n");
            printf("PTB-ERROR: PsychPortAudio, otherwise this error will happen.\n\n");
            PsychErrorExitMsg(PsychError_user, "Failed to initialize due to portaudio DLL loading problem. Call InitializePsychSound first! Aborted.");
        }
        #endif

        // Setup callback function for low-level debug output:
        PaUtil_SetDebugPrintFunction(PALogger);

        #if PSYCH_SYSTEM == PSYCH_LINUX
            // Set an error handler for ALSA debug output/errors to stop the spewage of utterly
            // pointless ALSA warning messages to stderr. At verbosity <= 5 we sent ALSA chatter
            // to a dummy error handler. At levels > 5 we disable our error handler, so ALSA
            // chatter goes to stderr. The same is true for the JACK backend, which does only
            // set its own error callback if it makes it through initialization in Pa_Initialize().
            // During the attempt to connect to the Jack server, it doesn't set its own callback,
            // so all error output during jack_client_open() spills into our stderr console -
            // and considerable pointless spillage there is :( -- Try to set our own override
            // just for the time during Pa_Initialize() -- it will get overriden, but at least
            // prevent the ugly spillage during startup:
            myjack_set_error_function = dlsym(RTLD_DEFAULT, "jack_set_error_function");
            if (myjack_set_error_function)
                myjack_set_error_function(PALogger);

            if (verbosity <= 5)
                snd_lib_error_set_handler(ALSAErrorHandler);
            else
                snd_lib_error_set_handler(NULL);
        #endif

        if ((err=Pa_Initialize())!=paNoError) {
            printf("PTB-ERROR: Portaudio initialization failed with following port audio error: %s \n", Pa_GetErrorText(err));
            PaUtil_SetDebugPrintFunction(NULL);
            PsychErrorExitMsg(PsychError_system, "Failed to initialize PortAudio subsystem.");
        }
        else {
            if(verbosity>2) {
                printf("PTB-INFO: Using modified %s\n", Pa_GetVersionText());
            }
        }

        for(i=0; i<MAX_PSYCH_AUDIO_DEVS; i++) {
            audiodevices[i].stream = NULL;
        }

        audiodevicecount=0;

        // Init audio bufferList to empty and Mutex to unlocked:
        bufferListCount = 0;
        bufferList = NULL;
        PsychInitMutex(&bufferListmutex);

        // On Vista systems and later, we assume everything will be fine wrt. to timing and multi-core
        // systems, but still perform consistency checks at each call to PsychGetPrecisionTimerSeconds().
        // Therefore we don't lock our threads to a single core by default. On pre-Vista systems, we
        // lock all threads to core 1 by default:
        lockToCore1 = (PsychIsMSVista()) ? FALSE : TRUE;

        pa_initialized = TRUE;
    }
}

/* PsychPortAudio('Open') - Open and initialize an audio device via PortAudio.
 */
PsychError PSYCHPORTAUDIOOpen(void)
{
    static char useString[] = "pahandle = PsychPortAudio('Open' [, deviceid][, mode][, reqlatencyclass][, freq][, channels][, buffersize][, suggestedLatency][, selectchannels][, specialFlags=0]);";
    //                                                            1             2         3                    4        5            6              7                      8                      9
    static char synopsisString[] =
    "Open a PortAudio audio device and initialize it. Returns a 'pahandle' device handle for the device.\n\n"
    "On most operating systems you can open each physical sound device only once per running session. If "
    "you feel the need to call 'Open' multiple times on the same audio device, read the section about slave "
    "devices and the help 'PsychPortAudio OpenSlave?' instead for a suitable solution.\n"
    "All parameters are optional and have reasonable defaults. 'deviceid' Index to select amongst multiple "
    "logical audio devices supported by PortAudio. Defaults to whatever the systems default sound device is. "
    "Different device id's may select the same physical device, but controlled by a different low-level sound "
    "system. E.g., Windows has about five different sound subsystems. 'mode' Mode of operation. Defaults to "
    "1 == sound playback only. Can be set to 2 == audio capture, or 3 for simultaneous capture and playback of sound. "
    "Note however that mode 3 (full duplex) does not work reliably on all sound hardware. On some hardware this mode "
    "may crash hard! There is also a special monitoring mode == 7, which only works for full duplex devices "
    "when using the same number of input- and outputchannels. This mode allows direct feedback of captured sounds "
    "back to the speakers with minimal latency and without involvement of your script at all, however no sound "
    "can be captured during this time and your code mostly doesn't have any control over timing etc. \n"
    "You can also define a audio device as a master device by adding the value 8 to mode. Master devices themselves "
    "are not directly used to playback or capture sound. Instead one can create (multiple) slave devices that are "
    "attached to a master device. Each slave can be controlled independently to playback or record sound through "
    "a subset of the channels of the master device. This basically allows to virtualize a soundcard. See "
    "help for subfunction 'OpenSlave' for more info.\n"
    "'reqlatencyclass' Allows to select how aggressive PsychPortAudio should be about minimizing sound latency "
    "and getting good deterministic timing, i.e. how to trade off latency vs. system load and playing nicely "
    "with other sound applications on the system. Level 0 means: Don't care about latency or timing precision. This mode works always "
    "and with all settings, plays nicely with other sound applications. Level 1 (the default) means: Try to get "
    "the lowest latency that is possible under the constraint of reliable playback, freedom of choice for all parameters and "
    "interoperability with other applications. Level 2 means: Take full control over the audio device, even if this "
    "causes other sound applications to fail or shutdown. Level 3 means: As level 2, but request the most aggressive "
    "settings for the given device. Level 4: Same as 3, but fail if device can't meet the strictest requirements. "
    "'freq' Requested playback/capture rate in samples per second (Hz). Defaults to a value that depends on the "
    "requested latency mode. 'channels' Number of audio channels to use, defaults to 2 for stereo. If you perform "
    "simultaneous playback and recording, you can provide a 2 element vector for 'channels', specifying different "
    "numbers of output channels and input channels. The first element in such a vector defines the number of playback "
    "channels, the 2nd element defines capture channels. E.g., [2, 1] would define 2 playback channels (stereo) and 1 "
    "recording channel. See the optional 'selectchannels' argument for selection of physical device channels on multi- "
    "channel cards.\n"
    "'buffersize' "
    "requested size and number of internal audio buffers, smaller numbers mean lower latency but higher system load "
    "and some risk of overloading, which would cause audio dropouts. 'suggestedLatency' optional requested latency in "
    "seconds. PortAudio selects internal operating parameters depending on sampleRate, suggestedLatency and buffersize "
    "as well as device internal properties to optimize for low latency output. Best left alone, only here as manual "
    "override in case all the auto-tuning cleverness fails.\n"
    "'selectchannels' optional matrix with mappings of logical channels to device channels: If you only want to use "
    "a subset of the channels present on your sound card, e.g., only 2 playback channels on a 16 channel soundcard, "
    "then you'd set the 'channels' argument to 2. The 'selectchannels' argument allows you to select, e.g.,  which "
    "two of the 16 channels to use for playback. 'selectchannels' is a one row by 'channels' matrix with mappings "
    "for pure playback or pure capture. For full-duplex mode (playback and capture), 'selectchannels' must be a "
    "2 rows by max(channels) column matrix. row 1 will define playback channel mappings, whereas row 2 will then "
    "define capture channel mappings. In any case, the number in the i'th column will define which physical device "
    "channel will be used for playback or capture of the i'th PsychPortAudio channel (the i'th row of your sound "
    "matrix). Numbering of physical device channels starts with zero! Example: Both, playback and simultaneous "
    "recording are requested and 'channels' equals 2, ie, two playback channels and two capture channels. If you'd "
    "specify 'selectchannels' as [0, 6 ; 12, 14], then playback would happen to device channels zero and six, "
    "sound would be captured from device channels 12 and 14. Please note that channel selection is currently "
    "only supported on some sound cards. The parameter is silently ignored on non-capable hardware or driver software.\n\n"
    "'specialFlags' Optional flags: Default to zero, can be or'ed or added together with the following flags/settings:\n"
    "1 = Never prime output stream. By default the output stream is primed. Don't bother if you don't know what this means.\n"
    "2 = Always clamp audio data to the valid -1.0 to 1.0 range. Clamping is enabled by default.\n"
    "4 = Never clamp audio data.\n"
    "8 = Always dither output audio data. By default, dithering is enabled in normal mode, and disabled in low-latency mode. "
    "Dithering adds a stochastic noise pattern to the least significant bit of each output sample to reduce the impact of "
    "audio quantization artifacts. Dithering can improve signal to noise ratio and quality of output sound, but it is more "
    "compute intense and it could change very low-level properties of the audio signal, because what you hear is not exactly "
    "what you specified.\n"
    "16 = Never dither audio data, not even in normal mode.\n\n";

    static char seeAlsoString[] = "Close GetDeviceSettings ";

    int buffersize, latencyclass, mode, deviceid, i, numel, specialFlags;
    double freq;
    int* nrchannels;
    int  mynrchannels[2];
    int  m, n, p;
    double* mychannelmap;
    double suggestedLatency, lowlatency;
    PaHostApiIndex paHostAPI;
    PaStreamParameters outputParameters;
    PaStreamParameters inputParameters;
    const PaDeviceInfo* inputDevInfo, *outputDevInfo, *referenceDevInfo;
    PaStreamFlags sflags;
    PaError err;
    PaStream *stream = NULL;

    #if PSYCH_SYSTEM == PSYCH_OSX
        #ifdef paMacCoreChangeDeviceParameters
            // New style since at least 19.6.0:
            PaMacCoreStreamInfo outhostapisettings, inhostapisettings;
            SInt32 CoreAudioOutChannelMap[MAX_PSYCH_AUDIO_CHANNELS_PER_DEVICE];
            SInt32 CoreAudioInChannelMap[MAX_PSYCH_AUDIO_CHANNELS_PER_DEVICE];
        #else
            // Old style:
            paMacCoreStreamInfo hostapisettings;
        #endif
    #endif

    #ifdef PA_ASIO_H
        int outputmappings[MAX_PSYCH_AUDIO_CHANNELS_PER_DEVICE];
        int inputmappings[MAX_PSYCH_AUDIO_CHANNELS_PER_DEVICE];
        PaAsioStreamInfo inasioapisettings;
        PaAsioStreamInfo outasioapisettings;
    #endif

    #if PSYCH_SYSTEM == PSYCH_WINDOWS
        PaWasapiStreamInfo inwasapiapisettings;
        PaWasapiStreamInfo outwasapiapisettings;
    #endif

    unsigned int id = PsychPANextHandle();

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(9));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(0)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    if (id >= MAX_PSYCH_AUDIO_DEVS) PsychErrorExitMsg(PsychError_user, "Maximum number of simultaneously open audio devices reached.");

    freq = 0;
    buffersize = 0;
    latencyclass = 1;
    mode = kPortAudioPlayBack;
    deviceid = -1;
    specialFlags = 0;

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    // Sanity check: Any hardware found?
    if (Pa_GetDeviceCount() == 0) PsychErrorExitMsg(PsychError_user, "Could not find *any* audio hardware on your system! Either your machine doesn't have audio hardware, or somethings seriously screwed.");

    // We default to generic system settings for host api specific settings:
    outputParameters.hostApiSpecificStreamInfo = NULL;
    inputParameters.hostApiSpecificStreamInfo  = NULL;

    // Request optional deviceid:
    PsychCopyInIntegerArg(1, kPsychArgOptional, &deviceid);
    if (deviceid < -1) PsychErrorExitMsg(PsychError_user, "Invalid deviceid provided. Valid values are -1 to maximum number of devices.");

    // Request optional mode of operation:
    PsychCopyInIntegerArg(2, kPsychArgOptional, &mode);
    if (mode < 1 || mode > 15 || mode & kPortAudioIsAMModulator || mode & kPortAudioIsAMModulatorForSlave || mode & kPortAudioIsOutputCapture || ((mode & kPortAudioMonitoring) && ((mode & kPortAudioFullDuplex) != kPortAudioFullDuplex))) {
        PsychErrorExitMsg(PsychError_user, "Invalid mode for regular- or master-audio device provided: Outside valid range or invalid combination of flags.");
    }

    if (!(mode & (kPortAudioCapture | kPortAudioPlayBack)))
        PsychErrorExitMsg(PsychError_user, "Invalid mode for regular- or master-audio device provided: mode must contain at least playback (1), capture (2) or full-duplex (3).");

    // Request optional latency class:
    PsychCopyInIntegerArg(3, kPsychArgOptional, &latencyclass);
    if (latencyclass < 0 || latencyclass > 4) PsychErrorExitMsg(PsychError_user, "Invalid reqlatencyclass provided. Valid values are 0 to 4.");

    if (deviceid == -1) {
        // Default devices requested:
        if ((latencyclass == 0) && (PSYCH_SYSTEM != PSYCH_LINUX)) {
            // High latency mode on non-Linux. Simply pick system default devices.
            // We don't pick these on Linux, because we'd end up with the ancient
            // OSS, which is almost always a suboptimal choice for our purpose, even
            // in high-latency mode.
            outputParameters.device = Pa_GetDefaultOutputDevice(); /* Default output device. */
            inputParameters.device  = Pa_GetDefaultInputDevice(); /* Default input device. */
        }
        else {
            // Low latency mode: Try to find the host API which is supposed to be the fastest on
            // a platform, then pick its default devices:
            paHostAPI = PsychPAGetLowestLatencyHostAPI();
            outputParameters.device = Pa_GetHostApiInfo(paHostAPI)->defaultOutputDevice;
            inputParameters.device  = Pa_GetHostApiInfo(paHostAPI)->defaultInputDevice;
        }

        // Make sure we don't choose a default audio output device which is likely to
        // send its output to nirvana. If this is the case, try to find a better alternative.
        // Currently blacklisted are HDMI and DisplayPort video outputs of graphics cards
        // with sound output over video.
        // Additionally we blacklist the "default", "sysdefault" and "pulse" devices if
        // we operate under Linux + ALSA and the Pulseaudio sound server got suspended, as
        // those devices most likely output their sound over Pulseaudio, which would end in
        // a hang if PulseAudio is suspended:
        if ((mode & kPortAudioPlayBack) || (mode & kPortAudioMonitoring)) {
            outputDevInfo = Pa_GetDeviceInfo(outputParameters.device);
            if (outputDevInfo && (Pa_GetDeviceCount() > 1) &&
                (strstr(outputDevInfo->name, "HDMI") || strstr(outputDevInfo->name, "hdmi") || strstr(outputDevInfo->name, "isplay") ||
                 ((outputDevInfo->hostApi == Pa_HostApiTypeIdToHostApiIndex(paALSA)) && pulseaudio_isSuspended && (strstr(outputDevInfo->name, "default") || strstr(outputDevInfo->name, "pulse"))))) {
                // Selected output default device seems to be a HDMI or DisplayPort output
                // of a graphics card. Try to find a better default choice.
                paHostAPI = outputDevInfo->hostApi;
                referenceDevInfo = NULL; // Make compiler happy.
                for (deviceid = 0; deviceid < (int) Pa_GetDeviceCount(); deviceid++) {
                    referenceDevInfo = Pa_GetDeviceInfo(deviceid);
                    if (!referenceDevInfo || (referenceDevInfo->hostApi != paHostAPI) ||
                        (referenceDevInfo->maxOutputChannels < 1) ||
                        (strstr(referenceDevInfo->name, "HDMI") || strstr(referenceDevInfo->name, "hdmi") || strstr(referenceDevInfo->name, "isplay") ||
                        ((referenceDevInfo->hostApi == Pa_HostApiTypeIdToHostApiIndex(paALSA)) && pulseaudio_isSuspended && (strstr(referenceDevInfo->name, "default") || strstr(referenceDevInfo->name, "pulse"))))) {
                        // Unsuitable.
                        continue;
                    }

                    // Found it:
                    break;
                }

                // Found something better? Otherwise we stick to the original choice.
                if (deviceid < Pa_GetDeviceCount()) {
                    // Yes.
                    if (verbosity > 2) printf("PTB-INFO: Choosing deviceIndex %i [%s] as default output audio device.\n", deviceid, referenceDevInfo->name);
                    outputParameters.device = (PaDeviceIndex) deviceid;
                }
                else {
                    // No, warn user about possible silence:
                    if (verbosity > 2) {
                        printf("PTB-INFO: Chosen default audio device with deviceIndex %i seems to be a HDMI or DisplayPort\n", (int) outputParameters.device);
                        printf("PTB-INFO: video output of your graphics card [Name = %s], or not a true hardware device.\n", outputDevInfo->name);
                        printf("PTB-INFO: Tried to find an alternative default output device but couldn't find a suitable one.\n");
                        printf("PTB-INFO: If you don't hear any sound, or the software appears to be hanging, then that is likely\n");
                        printf("PTB-INFO: the reason - sound playing out to a connected display device without any speakers.\n");
                        printf("PTB-INFO: See 'PsychPortAudio GetDevices?' for available devices.\n");
                    }
                }

                // Reset our temporaries:
                referenceDevInfo = NULL;
                deviceid = -1;
            }
        }

        // Same filterig applies to capture devices...
        if (mode & kPortAudioCapture) {
            inputDevInfo = Pa_GetDeviceInfo(inputParameters.device);
            if (inputDevInfo && (Pa_GetDeviceCount() > 1) &&
                (strstr(inputDevInfo->name, "HDMI") || strstr(inputDevInfo->name, "hdmi") || strstr(inputDevInfo->name, "isplay") ||
                 ((inputDevInfo->hostApi == Pa_HostApiTypeIdToHostApiIndex(paALSA)) && pulseaudio_isSuspended && (strstr(inputDevInfo->name, "default") || strstr(inputDevInfo->name, "pulse"))))) {
                // Selected input default device seems to be a HDMI or DisplayPort output
                // of a graphics card. Try to find a better default choice.
                paHostAPI = inputDevInfo->hostApi;
                referenceDevInfo = NULL; // Make compiler happy.
                for (deviceid = 0; deviceid < (int) Pa_GetDeviceCount(); deviceid++) {
                    referenceDevInfo = Pa_GetDeviceInfo(deviceid);
                    if (!referenceDevInfo || (referenceDevInfo->hostApi != paHostAPI) ||
                        (referenceDevInfo->maxOutputChannels < 1) ||
                        (strstr(referenceDevInfo->name, "HDMI") || strstr(referenceDevInfo->name, "hdmi") || strstr(referenceDevInfo->name, "isplay") ||
                        ((referenceDevInfo->hostApi == Pa_HostApiTypeIdToHostApiIndex(paALSA)) && pulseaudio_isSuspended && (strstr(referenceDevInfo->name, "default") || strstr(referenceDevInfo->name, "pulse"))))) {
                        // Unsuitable.
                        continue;
                    }

                    // Found it:
                    break;
                }

                // Found something better? Otherwise we stick to the original choice.
                if (deviceid < Pa_GetDeviceCount()) {
                    // Yes.
                    if (verbosity > 2) printf("PTB-INFO: Choosing deviceIndex %i [%s] as default input audio device.\n", deviceid, referenceDevInfo->name);
                    inputParameters.device = (PaDeviceIndex) deviceid;
                }
                else {
                    // No, warn user about possible silence:
                    if (verbosity > 2) {
                        printf("PTB-INFO: Chosen default audio device with deviceIndex %i seems to be a HDMI or DisplayPort\n", (int) inputParameters.device);
                        printf("PTB-INFO: video output of your graphics card [Name = %s], or not a true hardware device.\n", inputDevInfo->name);
                        printf("PTB-INFO: Tried to find an alternative default input device but couldn't find a suitable one.\n");
                        printf("PTB-INFO: If you record any sound, or the software appears to be hanging, then that is likely\n");
                        printf("PTB-INFO: the reason - capturing from a connected display device without any microphone.\n");
                        printf("PTB-INFO: See 'PsychPortAudio GetDevices?' for available devices.\n");
                    }
                }

                // Reset our temporaries:
                referenceDevInfo = NULL;
                deviceid = -1;
            }
        }
    }
    else {
        // Specific device requested: In valid range?
        if (deviceid >= Pa_GetDeviceCount() || deviceid < 0) {
            PsychErrorExitMsg(PsychError_user, "Invalid deviceid provided. Higher than the number of devices - 1 or lower than zero.");
        }

        outputParameters.device = (PaDeviceIndex) deviceid;
        inputParameters.device = (PaDeviceIndex) deviceid;
    }

    // Query properties of selected device(s):
    inputDevInfo  = Pa_GetDeviceInfo(inputParameters.device);
    outputDevInfo = Pa_GetDeviceInfo(outputParameters.device);

    // Select one of them as "reference" info devices: It's properties are used whenever
    // no more specialized info is available. We use the output device (if any) as reference,
    // otherwise fall back to the inputdevice.
    referenceDevInfo = (outputDevInfo) ? outputDevInfo : inputDevInfo;

    // Sanity check: Any hardware found?
    if (referenceDevInfo == NULL) PsychErrorExitMsg(PsychError_user, "Could not find *any* audio hardware on your system - or at least not with the provided deviceid, if any!");

    // Check if current set of selected/available devices is compatible with our playback mode:
    if (((mode & kPortAudioPlayBack) || (mode & kPortAudioMonitoring)) && ((outputDevInfo == NULL) || (outputDevInfo && outputDevInfo->maxOutputChannels <= 0))) {
        PsychErrorExitMsg(PsychError_user, "Audio output requested, but there isn't any audio output device available or you provided a deviceid for something else than an output device!");
    }

    if (((mode & kPortAudioCapture) || (mode & kPortAudioMonitoring)) && ((inputDevInfo == NULL) || (inputDevInfo && inputDevInfo->maxInputChannels <= 0))) {
        PsychErrorExitMsg(PsychError_user, "Audio input requested, but there isn't any audio input device available or you provided a deviceid for something else than an input device!");
    }

    // Request optional frequency:
    PsychCopyInDoubleArg(4, kPsychArgOptional, &freq);
    if (freq < 0) PsychErrorExitMsg(PsychError_user, "Invalid frequency provided. Must be greater than 0 Hz, or 0 for auto-select.");

    // Request optional number of channels:
    numel = 0; nrchannels = NULL;
    PsychAllocInIntegerListArg(5, kPsychArgOptional, &numel, &nrchannels);
    if (numel == 0) {
        // No optional channelcount argument provided: Default to two for playback and recording, unless device is
        // a mono device -- then we default to one. Note: It is important to never assign zero channels, as some
        // code makes assumptions about these always being non-zero, and would crash due to divide-by-zero otherwise.
        // Assigning at least 1 channel, even if the device doesn't have one, is safe by design of the remaining code.
        mynrchannels[0] = (outputDevInfo && outputDevInfo->maxOutputChannels < 2) ? 1 : 2;
        mynrchannels[1] = (inputDevInfo && inputDevInfo->maxInputChannels < 2) ? 1 : 2;
    }
    else if (numel == 1) {
        // One argument provided: Set same count for playback and recording:
        if (*nrchannels < 1 || *nrchannels > MAX_PSYCH_AUDIO_CHANNELS_PER_DEVICE) PsychErrorExitMsg(PsychError_user, "Invalid number of channels provided. Valid values are 1 to device maximum.");
        mynrchannels[0] = *nrchannels;
        mynrchannels[1] = *nrchannels;
    }
    else if (numel == 2) {
        // Separate counts for playback and recording provided: Set'em up.
        if (nrchannels[0] < 1 || nrchannels[0] > MAX_PSYCH_AUDIO_CHANNELS_PER_DEVICE) PsychErrorExitMsg(PsychError_user, "Invalid number of playback channels provided. Valid values are 1 to device maximum.");
        mynrchannels[0] = nrchannels[0];
        if (nrchannels[1] < 1 || nrchannels[1] > MAX_PSYCH_AUDIO_CHANNELS_PER_DEVICE) PsychErrorExitMsg(PsychError_user, "Invalid number of capture channels provided. Valid values are 1 to device maximum.");
        mynrchannels[1] = nrchannels[1];
    }
    else {
        // More than 2 channel counts provided? Impossible.
        mynrchannels[0] = mynrchannels[1] = 0; // Make compiler happy.
        PsychErrorExitMsg(PsychError_user, "You specified a list with more than two 'channels' entries? Can only be max 2 for playback- and capture.");
    }

    // Make sure that number of capture and playback channels is the same for fast monitoring/feedback mode:
    if ((mode & kPortAudioMonitoring) && (mynrchannels[0] != mynrchannels[1])) PsychErrorExitMsg(PsychError_user, "Fast monitoring/feedback mode selected, but number of capture and playback channels differs! They must be the same for this mode!");

    // Request optional buffersize:
    PsychCopyInIntegerArg(6, kPsychArgOptional, &buffersize);
    if (buffersize < 0 || buffersize > 4096) PsychErrorExitMsg(PsychError_user, "Invalid buffersize provided. Valid values are 0 to 4096 samples.");

    // Request optional suggestedLatency:
    suggestedLatency = -1.0;
    PsychCopyInDoubleArg(7, kPsychArgOptional, &suggestedLatency);
    if (suggestedLatency!=-1 && (suggestedLatency < 0.0 || suggestedLatency > 1.0)) PsychErrorExitMsg(PsychError_user, "Invalid suggestedLatency provided. Valid values are 0.0 to 1.0 seconds.");

    #if PSYCH_SYSTEM == PSYCH_OSX
        // Apply OS/X CoreAudio specific optimizations:
        #ifdef paMacCoreChangeDeviceParameters
            // New style since at least 19.6.0:
            unsigned long flags = (latencyclass > 1) ? paMacCoreChangeDeviceParameters : 0;
            if (latencyclass > 3) flags|= paMacCoreFailIfConversionRequired;
            PaMacCore_SetupStreamInfo( &outhostapisettings, flags);
            PaMacCore_SetupStreamInfo( &inhostapisettings, flags);
            outputParameters.hostApiSpecificStreamInfo = (PaMacCoreStreamInfo*) &outhostapisettings;
            inputParameters.hostApiSpecificStreamInfo = (PaMacCoreStreamInfo*) &inhostapisettings;
        #else
            // Old style:
            unsigned long flags = (latencyclass > 1) ? paMacCore_ChangeDeviceParameters : 0;
            if (latencyclass > 3) flags|= paMacCore_FailIfConversionRequired;
            paSetupMacCoreStreamInfo( &hostapisettings, flags);
            outputParameters.hostApiSpecificStreamInfo = (paMacCoreStreamInfo*) &hostapisettings;
            inputParameters.hostApiSpecificStreamInfo = (paMacCoreStreamInfo*) &hostapisettings;
        #endif
    #endif

    #if PSYCH_SYSTEM == PSYCH_WINDOWS
    if (Pa_GetHostApiInfo(referenceDevInfo->hostApi)->type == paWASAPI) {
        inputParameters.hostApiSpecificStreamInfo = (PaWasapiStreamInfo*) &inwasapiapisettings;
        outputParameters.hostApiSpecificStreamInfo = (PaWasapiStreamInfo*) &outwasapiapisettings;

        // Init wasapi stream settings struct to vanilla defaults:
        memset(&inwasapiapisettings, 0, sizeof(inwasapiapisettings));
        memset(&outwasapiapisettings, 0, sizeof(outwasapiapisettings));
        inwasapiapisettings.size = sizeof(PaWasapiStreamInfo);
        outwasapiapisettings.size = sizeof(PaWasapiStreamInfo);
        inwasapiapisettings.hostApiType = paWASAPI;
        outwasapiapisettings.hostApiType = paWASAPI;
        inwasapiapisettings.version = 1;
        outwasapiapisettings.version = 1;

        // Default, but set for clarity - We don't use a special category like voice chat, media, movies etc.:
        inwasapiapisettings.streamCategory = eAudioCategoryOther;
        outwasapiapisettings.streamCategory = eAudioCategoryOther;

        // Disable digital signal processing (equalizers, bass boost, automatic gain control, 3d spatialisation and such)
        // by default on Windows 8.1 and later - unsupported on Windows 8 and earlier:
        inwasapiapisettings.streamOption = eStreamOptionRaw;
        outwasapiapisettings.streamOption = eStreamOptionRaw;
    }
    #endif

    // Get optional channel map:
    mychannelmap = NULL;
    PsychAllocInDoubleMatArg(8, kPsychArgOptional, &m, &n, &p, &mychannelmap);
    if (mychannelmap) {
        // Channelmapping provided: Sanity check it.
        if (m<1 || m>2 || p!=1 || (n!=((mynrchannels[0] > mynrchannels[1]) ? mynrchannels[0] : mynrchannels[1]))) {
            PsychErrorExitMsg(PsychError_user, "Invalid size of 'selectchannels' matrix argument: Must be a one- or two row by max(channels) column matrix!");
        }

        #ifdef PA_ASIO_H
        if (Pa_GetHostApiInfo(referenceDevInfo->hostApi)->type == paASIO) {
            // MS-Windows and connected to an ASIO device. Try to assign channel mapping:
            if (mode & kPortAudioPlayBack) {
                // Playback mappings:
                outputParameters.hostApiSpecificStreamInfo = (PaAsioStreamInfo*) &outasioapisettings;
                outasioapisettings.size = sizeof(PaAsioStreamInfo);
                outasioapisettings.hostApiType = paASIO;
                outasioapisettings.version = 1;
                outasioapisettings.flags = paAsioUseChannelSelectors;
                outasioapisettings.channelSelectors = (int*) &outputmappings[0];
                for (i=0; i<mynrchannels[0]; i++) outputmappings[i] = (int) mychannelmap[i * m];
                if (verbosity > 3) {
                    printf("PTB-INFO: Will try to use the following logical channel -> device channel mappings for sound output to audio stream %i :\n", id);
                    for (i=0; i<mynrchannels[0]; i++) printf("%i --> %i : ", i+1, outputmappings[i]);
                    printf("\n\n");
                }
            }

            if (mode & kPortAudioCapture) {
                // Capture mappings:
                inputParameters.hostApiSpecificStreamInfo = (PaAsioStreamInfo*) &inasioapisettings;
                inasioapisettings.size = sizeof(PaAsioStreamInfo);
                inasioapisettings.hostApiType = paASIO;
                inasioapisettings.version = 1;
                inasioapisettings.flags = paAsioUseChannelSelectors;
                inasioapisettings.channelSelectors = (int*) &inputmappings[0];
                // Index into first row of one-row matrix or 2nd row of two-row matrix:
                for (i=0; i<mynrchannels[1]; i++) inputmappings[i] = (int) mychannelmap[(i * m) + (m-1)];
                if (verbosity > 3) {
                    printf("PTB-INFO: Will try to use the following logical channel -> device channel mappings for sound capture from audio stream %i :\n", id);
                    for (i=0; i<mynrchannels[1]; i++) printf("%i --> %i : ", i+1, inputmappings[i]);
                    printf("\n\n");
                }
            }
            // Mappings setup up. The PortAudio library will sanity check this further against device constraints...
        }
        else
        #endif
        #if PSYCH_SYSTEM == PSYCH_WINDOWS
        if (Pa_GetHostApiInfo(referenceDevInfo->hostApi)->type != paASIO) {
            // Compute and apply channelmap:
            if (mode & kPortAudioPlayBack) {
                PaWinWaveFormatChannelMask channelMask = (PaWinWaveFormatChannelMask) 0;
                for (i = 0; i < mynrchannels[0]; i++)
                    channelMask |= (PaWinWaveFormatChannelMask) (1 << ((int) mychannelmap[i * m]));

                if (verbosity > 3) {
                    printf("PTB-INFO: Will try to use the following logical channel -> device channel mappings for sound output to audio stream %i :\n", id);
                    for (i = 0; i < mynrchannels[0]; i++)
                        printf("%i --> %i : ", i + 1, (int) mychannelmap[i * m]);

                    printf("\n\n");
                }

                switch (Pa_GetHostApiInfo(referenceDevInfo->hostApi)->type) {
                    case paWASAPI:
                        outwasapiapisettings.flags |= paWinWasapiUseChannelMask;
                        outwasapiapisettings.channelMask = channelMask;
                        break;

                    /* TODO:
                    case paWDMKS:
                        outwasapiapisettings.flags |= ;
                        outwasapiapisettings.channelMask = channelMask;
                        break;

                    case paDirectSound:
                        outwasapiapisettings.flags |= ;
                        outwasapiapisettings.channelMask = channelMask;
                        break;

                    case paMME:
                        outwasapiapisettings.flags |= ;
                        outwasapiapisettings.channelMask = channelMask;
                        break;
                    */
                }
            }

            // Windows api's can't do channel mapping on the capture side, i think?
            if ((mode & kPortAudioCapture) && (verbosity > 3)) {
                printf("PTB-INFO: Audio capture enabled, but 'selectchannels' mapping for audio capture is ignored on this hardware + driver combo.\n");
            }
        }
        else
        #endif

        #if PSYCH_SYSTEM == PSYCH_OSX && defined(paMacCoreChangeDeviceParameters)
            if (Pa_GetHostApiInfo(referenceDevInfo->hostApi)->type == paCoreAudio) {
                // macOS CoreAudio. Try to assign channel mapping:
                if (mode & kPortAudioPlayBack) {
                    // Playback mappings:
                    for (i = 0; i < outputDevInfo->maxOutputChannels && i < MAX_PSYCH_AUDIO_CHANNELS_PER_DEVICE; i++)
                        CoreAudioOutChannelMap[i] = -1;

                    for (i = 0; i < mynrchannels[0]; i++) {
                        if (((int) mychannelmap[i * m] >= outputDevInfo->maxOutputChannels) || ((int) mychannelmap[i * m] < 0) ||
                            ((int) mychannelmap[i * m] >= MAX_PSYCH_AUDIO_CHANNELS_PER_DEVICE)) {
                            if (verbosity > 0)
                                printf("PTB-ERROR: Invalid output channel mapping %i -> %i. Output device only has %i channels.\n",
                                       i+1, (int) mychannelmap[i * m], outputDevInfo->maxOutputChannels);
                            PsychErrorExitMsg(PsychError_user, "Invalid 'selectchannels' matrix argument for output mapping.");
                        }

                        CoreAudioOutChannelMap[(int) mychannelmap[i * m]] = i;
                    }

                    if (verbosity > 3) {
                        printf("PTB-INFO: Will try to use the following logical channel -> device channel mappings for sound output to audio stream %i :\n", id);
                        for (i = 0; i < mynrchannels[0]; i++)
                            printf("%i --> %i : ", i+1, (int) mychannelmap[i * m]);
                        printf("\n\n");
                    }

                    PaMacCore_SetupChannelMap(&outhostapisettings,
                                              (const SInt32 * const) CoreAudioOutChannelMap,
                                              (unsigned long) outputDevInfo->maxOutputChannels);
                }

                if (mode & kPortAudioCapture) {
                    // Capture mappings:
                    for (i = 0; i < mynrchannels[1]; i++)
                        CoreAudioInChannelMap[i] = -1;

                    // Index into first row of one-row matrix or 2nd row of two-row matrix:
                    for (i = 0; i < mynrchannels[1]; i++) {
                        // Input channels are mapped differently in CoreAudio from output channels. Target channel on our
                        // side is array slot, source device channel is assigned value:
                        if (((int) mychannelmap[(i * m) + (m - 1)] < 0) ||
                            ((int) mychannelmap[(i * m) + (m - 1)] >= inputDevInfo->maxInputChannels )) {
                            if (verbosity > 0)
                                printf("PTB-ERROR: Invalid input channel mapping %i <- %i. Input device only has %i channels.\n",
                                       i+1, (int) mychannelmap[(i * m) + (m - 1)], inputDevInfo->maxInputChannels);
                            PsychErrorExitMsg(PsychError_user, "Invalid 'selectchannels' matrix argument for input mapping.");
                        }

                        CoreAudioInChannelMap[i] = (int) mychannelmap[(i * m) + (m - 1)];
                    }

                    if (verbosity > 3) {
                        printf("PTB-INFO: Will try to use the following logical channel <- device channel mappings for sound capture from audio stream %i :\n", id);
                        for (i = 0; i < mynrchannels[1]; i++)
                            printf("%i <-- %i : ", i+1, (int) mychannelmap[(i * m) + (m - 1)]);
                        printf("\n\n");
                    }

                    PaMacCore_SetupChannelMap(&inhostapisettings,
                                              (const SInt32 * const) CoreAudioInChannelMap,
                                              (unsigned long) mynrchannels[1]);
                }
                // Mappings setup up. The PortAudio library will sanity check this further against device constraints...
            }
            else
        #endif

            // Backend without channel mapping support.
            if (verbosity > 2) printf("PTB-WARNING: Provided 'selectchannels' channel mapping is ignored because this audio backend does not support it.\n");
    }

    // Copy in optional specialFlags:
    PsychCopyInIntegerArg(9, kPsychArgOptional, &specialFlags);

    // Set channel count:
    outputParameters.channelCount = mynrchannels[0];    // Number of output channels.
    inputParameters.channelCount = mynrchannels[1];     // Number of input channels.

    // Fix sample format to float for now...
    outputParameters.sampleFormat = paFloat32;
    inputParameters.sampleFormat  = paFloat32;

    // Setup buffersize:
    if (buffersize == 0) {
        // No specific buffersize requested: Leave it unspecified to get optimal selection by lower level driver.
        // Especially on Windows and Linux, basically any choice other than paFramesPerBufferUnspecified
        // will just lead to trouble and less optimal results, as the sound subsystems are very good at
        // choosing this parameter optimally if allowed, but have a hard job to cope with any user-enforced choices.
        buffersize = paFramesPerBufferUnspecified;
    }
    else {
        // Extra buffersize validation possible on OSX with upstream portaudio:
        #if (PSYCH_SYSTEM == PSYCH_OSX) && defined(paMacCoreChangeDeviceParameters)
        long minBufferSizeFrames, maxBufferSizeFrames;
        if (paNoError == PaMacCore_GetBufferSizeRange(outputDevInfo ? outputParameters.device : inputParameters.device,
                                                      &minBufferSizeFrames, &maxBufferSizeFrames)) {
            if (verbosity > 3)
                printf("PTB-INFO: Allowable host audiobuffersize range is %i to %i sample frames.\n",
                       (int) minBufferSizeFrames, (int) maxBufferSizeFrames);

            if ((long) buffersize < minBufferSizeFrames) {
                if (verbosity > 2)
                    printf("PTB-WARNING: Requested buffersize %i samples is lower than allowed minimum %i. Clamping to minimum.\n",
                           buffersize, (int) minBufferSizeFrames);
                buffersize = (int) minBufferSizeFrames;
            }

            if ((long) buffersize > maxBufferSizeFrames) {
                if (verbosity > 2)
                    printf("PTB-WARNING: Requested buffersize %i samples is higher than allowed maximum %i. Clamping to maximum.\n",
                           buffersize, (int) maxBufferSizeFrames);
                buffersize = (int) maxBufferSizeFrames;
            }
        }
        #endif
    }

    // Now we have auto-selected buffersize or user provided override...

    // Setup samplerate:
    if (freq == 0) {
        // No specific frequency requested:
        if (latencyclass < 3) {
            // At levels < 3, we select the device specific default.
            freq = referenceDevInfo->defaultSampleRate;
        }
        else {
            freq = 96000; // Go really high...
        }
    }

    // Now we have auto-selected frequency or user provided override...

    // Set requested latency: In class 0 we choose device recommendation for dropout-free operation, in
    // all higher (lowlat) classes we request zero or low latency. PortAudio will
    // clamp this request to something safe internally.
    switch (Pa_GetHostApiInfo(referenceDevInfo->hostApi)->type) {
        case paCoreAudio:    // CoreAudio is sometimes broken and causes sound dropouts at low latencies on 10.14 Mojave.
            lowlatency = (latencyclass > 1) ? 0.0 : 0.010;    // Go easy on it and only ask for 10 msecs latency at default
            break;                                            // latency class, otherwise for minimum, as in the past.
        case paWASAPI:       // WASAPI will automatically clamp to safe minimum.
        case paWDMKS:        // dto. for WDMKS Windows kernel streaming.
            lowlatency = 0.0;
            break;

        case paMME:             // No such a thing as low latency, but easy to kill the machine with too low settings!
            lowlatency = 0.1;   // Play safe, request 100 msecs. Otherwise terrible things may happen!
            break;

        case paDirectSound:     // DirectSound defaults to 120 msecs, which is way too much! It doesn't accept 0.0 msecs.
            lowlatency = 0.05;  // Choose some half-way safe tradeoff: 50 msecs.
            break;

        #ifdef PTB_USE_ASIO
        case paASIO:
            // ASIO: A value of zero would set safe (and high latency!) defaults. Too small values get
            // clamped to a safe minimum by the driver, so we select a very small positive value, say
            // 1 msec to get lowest possible latency for latencyclass of at least 2. In latency class 1
            // we play a bit safer and go for 5 msecs:
            lowlatency = (latencyclass >= 2) ? 0.001 : 0.005;
            break;
        #endif

        case paALSA:
            // For ALSA we choose 10 msecs by default, lowering to 5 msecs if exp. requested. Experience
            // shows that the effective latency will be only a fraction of this, so we are good.
            // Otoh going too low could cause dropouts on the tested 2006 MacbookPro 2nd generation...
            // This will need more lab testing and tweaking - and the user can override anyway...
            lowlatency = (latencyclass > 2) ? 0.005 : 0.010;
            break;

        default:            // Not the safest assumption for non-verified Api's, but we'll see...
            lowlatency = 0.0;
    }

    if (suggestedLatency == -1.0) {
        // None provided: Choose default based on latency mode:
        outputParameters.suggestedLatency = (latencyclass == 0 && outputDevInfo) ? outputDevInfo->defaultHighOutputLatency : lowlatency;
        inputParameters.suggestedLatency  = (latencyclass == 0 && inputDevInfo) ? inputDevInfo->defaultHighInputLatency : lowlatency;

        // Make sure that requested high or default output latency on Apples trainwreck is never lower than 10 msecs.
        // Especially on macOS 10.14 this seems to be neccessary for crackle-free playback: (Forum message #23422)
        if ((latencyclass <= 1) && (PSYCH_SYSTEM == PSYCH_OSX) && (outputParameters.suggestedLatency < 0.010))
            outputParameters.suggestedLatency = 0.010;
    }
    else {
        // Override provided: Use it.
        outputParameters.suggestedLatency = suggestedLatency;
        inputParameters.suggestedLatency  = suggestedLatency;
    }

    #if PSYCH_SYSTEM == PSYCH_WINDOWS
    if (Pa_GetHostApiInfo(referenceDevInfo->hostApi)->type == paWASAPI) {
        if (latencyclass > 1) {
            inwasapiapisettings.flags |= paWinWasapiExclusive;
            outwasapiapisettings.flags |= paWinWasapiExclusive;
        }

        if (latencyclass > 2) {
            // Windows 10 or later?
            if (PsychOSIsMSWin10()) {
                // Force audio engine into using the optimal sample format for us:
                // This is only supported on Windows-10 and later.
                inwasapiapisettings.streamOption = eStreamOptionMatchFormat;
                outwasapiapisettings.streamOption = eStreamOptionMatchFormat;
            }
            else if (latencyclass > 3) {
                // At least class 4, which means to fail if settings can't be applied,
                // and this is a pre-Windows 10 system, where the settings are unsupported,
                // so fail loudly:
                printf("PTB-ERROR: Desired ultralow-latency audio parameters for device %i unsupported by operating system, as this would require Windows-10 or later.\n", deviceid);
                PsychErrorExitMsg(PsychError_system, "Failed to open PortAudio audio device due to unsupported reqlatencyclass 4 on operating system older than Windows 10.");
            }
        }
    }
    #endif

    // Our stream shall be primed initially with our callbacks data, not just with zeroes.
    // In high latency-mode 0, we request sample clipping and dithering, so sound is more
    // high quality on Windows. In low-latency mode, we don't dither by default to save computation time.
    sflags = paPrimeOutputBuffersUsingStreamCallback;
    sflags = (latencyclass <= 0) ? sflags : (sflags | paDitherOff);

    // specialFlags 1: Don't prime outputbuffer via streamcallback:
    if (specialFlags & 1) sflags &= ~paPrimeOutputBuffersUsingStreamCallback;

    // specialFlags 2: Always clip audio data:
    if (specialFlags & 2) sflags &= ~paClipOff;

    // specialFlags 4: Never clip audio data:
    if (specialFlags & 4) sflags |= paClipOff;

    // specialFlags 8: Always dither audio data:
    if (specialFlags & 8) sflags &= ~paDitherOff;

    // specialFlags 16: Never dither audio data:
    if (specialFlags & 16) sflags |= paDitherOff;

    #if PSYCH_SYSTEM == PSYCH_LINUX
        // On ALSA in aggressive low-latency mode, reduce number of periods (aka device buffers) to 2 for double-buffering.
        // The default in Portaudio is 4 periods, so that's what we use in non-aggressive mode
        if (Pa_GetHostApiInfo(referenceDevInfo->hostApi)->type == paALSA)
            PaAlsa_SetNumPeriods((latencyclass > 2) ? 2 : 4);
    #endif

    // Check if the requested sample format and settings are likely supported by Audio API:
    err = Pa_IsFormatSupported(((mode & kPortAudioCapture) ?  &inputParameters : NULL), ((mode & kPortAudioPlayBack) ? &outputParameters : NULL), freq);
    if (err != paNoError && err != paDeviceUnavailable) {
        printf("PTB-ERROR: Desired audio parameters for device %i unsupported by audio device: %s \n", deviceid, Pa_GetErrorText(err));
        if (err == paInvalidSampleRate)
            printf("PTB-ERROR: Seems the requested audio sample rate %lf Hz is not supported by this combo of hardware and sound driver.\n", freq);
        else if (err == paInvalidChannelCount)
            printf("PTB-ERROR: Seems the requested number of audio channels is not supported by this combo of hardware and sound driver.\n");
        else
            printf("PTB-ERROR: This could be, e.g., due to an unsupported combination of audio sample rate, audio channel count/allocation, or audio sample format.\n");

        if (PSYCH_SYSTEM == PSYCH_LINUX)
            printf("PTB-ERROR: On Linux you may be able to use ALSA audio converter plugins to make this work.\n");
        PsychErrorExitMsg(PsychError_user, "Failed to open PortAudio audio device due to unsupported combination of audio parameters.");
    }

    // Try to create & open stream:
    if (err == paNoError)
        err = Pa_OpenStream(
                            &stream,                                                        /* Return stream pointer here on success. */
                            ((mode & kPortAudioCapture) ?  &inputParameters : NULL),        /* Requested input settings, or NULL in pure playback case. */
                            ((mode & kPortAudioPlayBack) ? &outputParameters : NULL),       /* Requested input settings, or NULL in pure playback case. */
                            freq,                                                           /* Requested sampling rate. */
                            buffersize,                                                     /* Requested buffer size. */
                            sflags,                                                         /* Define special stream property flags. */
                            paCallback,                                                     /* Our processing callback. */
                            &audiodevices[id]);                               /* Our own device info structure */

    if (err != paNoError || stream == NULL) {
        printf("PTB-ERROR: Failed to open audio device %i. PortAudio reports this error: %s \n", deviceid, Pa_GetErrorText(err));
        if (err == paDeviceUnavailable) {
            printf("PTB-ERROR: Could not open audio device, most likely because it is already in use by a previous call to\n");
            printf("PTB-ERROR: PsychPortAudio('Open', ...). You can open each audio device only once per session. If you need\n");
            printf("PTB-ERROR: multiple independent devices simulated on one physical audio device, look into use of audio\n");
            printf("PTB-ERROR: slave devices. See help for this by typing 'PsychPortAudio OpenSlave?'.\n");
            PsychErrorExitMsg(PsychError_user, "Audio device unavailable. Most likely tried to open device multiple times.");
        }

        PsychErrorExitMsg(PsychError_system, "Failed to open PortAudio audio device.");
    }

    // Setup our final device structure:
    audiodevices[id].opmode = mode;
    audiodevices[id].runMode = 1; // Keep engine running by default. Minimal extra cpu-load for significant reduction in startup latency.
    audiodevices[id].stream = stream;
    audiodevices[id].streaminfo = Pa_GetStreamInfo(stream);
    audiodevices[id].hostAPI = Pa_GetHostApiInfo(referenceDevInfo->hostApi)->type;
    audiodevices[id].startTime = 0.0;
    audiodevices[id].reqStartTime = 0.0;
    audiodevices[id].reqStopTime = DBL_MAX;
    audiodevices[id].estStopTime = 0;
    audiodevices[id].currentTime = 0;
    audiodevices[id].state = 0;
    audiodevices[id].reqstate = 255;
    audiodevices[id].repeatCount = 1;
    audiodevices[id].outputbuffer = NULL;
    audiodevices[id].outputbuffersize = 0;
    audiodevices[id].inputbuffer = NULL;
    audiodevices[id].inputbuffersize = 0;
    audiodevices[id].outchannels = mynrchannels[0];
    audiodevices[id].inchannels = mynrchannels[1];
    audiodevices[id].latencyBias = 0.0;
    audiodevices[id].schedule = NULL;
    audiodevices[id].schedule_size = 0;
    audiodevices[id].schedule_pos = 0;
    audiodevices[id].schedule_writepos = 0;
    audiodevices[id].outdeviceidx = (audiodevices[id].opmode & kPortAudioPlayBack) ? outputParameters.device : -1;
    audiodevices[id].indeviceidx  = (audiodevices[id].opmode & kPortAudioCapture)  ? inputParameters.device  : -1;
    audiodevices[id].outputmappings = NULL;
    audiodevices[id].inputmappings = NULL;
    audiodevices[id].slaveCount = 0;
    audiodevices[id].slaves = NULL;
    audiodevices[id].pamaster = -1;
    audiodevices[id].modulatorSlave = -1;
    audiodevices[id].slaveOutBuffer = NULL;
    audiodevices[id].slaveGainBuffer = NULL;
    audiodevices[id].slaveInBuffer = NULL;
    audiodevices[id].outChannelVolumes = NULL;
    audiodevices[id].masterVolume = 1.0;
    audiodevices[id].playposition = 0;
    audiodevices[id].totalplaycount = 0;

    // If this is a master, create a slave device list and init it to "empty":
    if (mode & kPortAudioIsMaster) {
        audiodevices[id].slaves = (int*) malloc(sizeof(int) * MAX_PSYCH_AUDIO_SLAVES_PER_DEVICE);
        if (NULL == audiodevices[id].slaves) PsychErrorExitMsg(PsychError_outofMemory, "Insufficient memory during slave devicelist creation!");
        for (i=0; i < MAX_PSYCH_AUDIO_SLAVES_PER_DEVICE; i++) audiodevices[id].slaves[i] = -1;

        if (mode & kPortAudioPlayBack) {
            // Allocate a dummy outputbuffer with one sampleframe:
            audiodevices[id].outputbuffersize = sizeof(float) * audiodevices[id].outchannels * 1;
            audiodevices[id].outputbuffer = (float*) malloc((size_t) audiodevices[id].outputbuffersize);
            if (audiodevices[id].outputbuffer==NULL) PsychErrorExitMsg(PsychError_outofMemory, "Out of system memory when trying to allocate audio buffer.");
        }

        if (mode & kPortAudioCapture) {
            // Allocate a dummy inputbuffer with one sampleframe:
            audiodevices[id].inputbuffersize = sizeof(float) * audiodevices[id].inchannels * 1;
            audiodevices[id].inputbuffer = (float*) calloc(1, (size_t) audiodevices[id].inputbuffersize);
            if (audiodevices[id].inputbuffer == NULL) PsychErrorExitMsg(PsychError_outofMemory, "Free system memory exhausted when trying to allocate audio recording buffer!");
        }
    }

    // If we use locking, we need to initialize the per-device mutex:
    if (uselocking && PsychInitMutex(&(audiodevices[id].mutex))) {
        printf("PsychPortAudio: CRITICAL! Failed to initialize Mutex object for pahandle %i! Prepare for trouble!\n", id);
        PsychErrorExitMsg(PsychError_system, "Audio device mutex creation failed!");
    }

    // If we use locking, this will create & init the associated event variable:
    PsychPACreateSignal(&(audiodevices[id]));

    // Register the stream finished callback:
    Pa_SetStreamFinishedCallback(audiodevices[id].stream, PAStreamFinishedCallback);

    #if PSYCH_SYSTEM == PSYCH_LINUX
        // Enable realtime scheduling for the portaudio audio processing thread on ALSA:
        if (audiodevices[id].hostAPI == paALSA) {
            PaAlsa_EnableRealtimeScheduling(audiodevices[id].stream, 1);

            // Is a hack-fix needed for a bug present in all versions of libportaudio
            // earlier than 19.7? The bug triggers in pure half-duplex audio capture
            // mode (== no playback requested) and consists of erroneously switching
            // the hostBufferSizeMode of the PaBufferProcessor for capture streams from
            // paUtilFixedHostBufferSize to paUtilBoundedHostBufferSize, causing
            // misconversion of captured audio somewhere in the bufferProcessor.
            //
            // The offending (buggy) if () conditional is at the end of the function
            // PaAlsaStream_DetermineFramesPerBuffer() inside pa_linux_alsa.c. It looks
            // like this: if( !self->playback.canMmap || !accurate ) ..., when it should look
            // like this: if((!self->playback.canMmap && self->playback.pcm) || !accurate ) ...
            // in order to prevent the invalid !self->playback.canMmap check, which always
            // evaluates to "true" if playback is not actually used in pure capture mode.
            //
            // We can't prevent the mis-setting of hostBufferSizeMode to paUtilBoundedHostBufferSize,
            // but we can reset it to the correct paUtilFixedHostBufferSize after Pa_OpenStream()
            // and before actual start of capture operations, by hacking ourselves into the internal
            // data structure of the PaBufferProcessor associated with the audio stream.
            //
            // We don't apply the fix if usercode specified a buffersize other than
            // the default paFramesPerBufferUnspecified -- in such a case, switching
            // to a different hostBufferSizeMode may be required...
            if (((mode & kPortAudioFullDuplex) == kPortAudioCapture) &&
                (buffersize == paFramesPerBufferUnspecified) &&
                (Pa_GetVersion() < paMakeVersionNumber(19,6,1))) {
                PseudoBufferProcessor* bp;

                // Cast PaStream* into ALSA backend specific AlsaStream*, using the pseudo-
                // datastructures defined to mimick PortAudio 19.6.0 an earliers true internal
                // structures:
                PseudoAlsaStreamExt* myalsaext = (PseudoAlsaStreamExt*) audiodevices[id].stream;
                if ((myalsaext->streamRepresentation.magic == PA_STREAM_MAGIC) &&
                    (myalsaext->streamRepresentation.hostApiType == paALSA)) {
                    bp = &myalsaext->bufferProcessor;
                    if (verbosity > 5)
                        printf("PTB-DEBUG: Probed PaAlsaStream struct is from patched portaudio with portmixer.patch (Ubuntu style).\n");
                }
                else {
                    PseudoAlsaStream* myalsa = (PseudoAlsaStream*) audiodevices[id].stream;
                    if (myalsa->streamRepresentation.magic == PA_STREAM_MAGIC) {
                        bp = &myalsa->bufferProcessor;
                        if (verbosity > 5)
                            printf("PTB-DEBUG: Probed PaAlsaStream struct is from vanilla upstream portaudio.\n");
                    }
                    else {
                        bp = NULL;
                        if (verbosity > 1)
                            printf("PTB-WARNING: Probing PaAlsaStream struct failed! Can't apply half-duplex capture bug workaround for Portaudio v19.6.0!\n");
                    }
                }

                // Additional sanity check if we got the casting right, and if our fix-hack is applicable:
                if (bp && (bp->hostBufferSizeMode == paUtilBoundedHostBufferSize)) {
                    // Reset to sane mode:
                    bp->hostBufferSizeMode = paUtilFixedHostBufferSize;

                    if (verbosity > 4)
                        printf("PTB-INFO: Applying paUtilFixedHostBufferSize workaround for pure half-duplex capture mode.\n");
                }
            }
        }
    #endif

    if (verbosity > 3) {
        printf("PTB-INFO: New audio device %i with handle %i opened as PortAudio stream:\n", deviceid, id);

        if (audiodevices[id].opmode & kPortAudioPlayBack) {
            printf("PTB-INFO: For %i channels Playback: Audio subsystem is %s, Audio device name is ", (int) audiodevices[id].outchannels, Pa_GetHostApiInfo(outputDevInfo->hostApi)->name);
            printf("%s\n", outputDevInfo->name);
        }

        if (audiodevices[id].opmode & kPortAudioCapture) {
            printf("PTB-INFO: For %i channels Capture: Audio subsystem is %s, Audio device name is ", (int) audiodevices[id].inchannels, Pa_GetHostApiInfo(inputDevInfo->hostApi)->name);
            printf("%s\n", inputDevInfo->name);
        }

        printf("PTB-INFO: Real samplerate %f Hz. Input latency %f msecs, Output latency %f msecs.\n",
               audiodevices[id].streaminfo->sampleRate, audiodevices[id].streaminfo->inputLatency * 1000.0,
               audiodevices[id].streaminfo->outputLatency * 1000.0);
    }

    // Return device handle:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) id);

    // One more audio device...
    audiodevicecount++;

    return(PsychError_none);
}

/* PsychPortAudio('OpenSlave') - Open and initialize a virtual audio slave device.
 */
PsychError PSYCHPORTAUDIOOpenSlave(void)
{
    static char useString[] = "pahandle = PsychPortAudio('OpenSlave', pamaster [, mode][, channels][, selectchannels]);";
    //                                                                  1              2          3              4
    static char synopsisString[] =
    "Open a virtual slave audio device and initialize it. Returns a 'pahandle' device handle for the device.\n"
    "Slave audio devices are almost always attached to a 'pamaster' audio device that you need to open and initialize "
    "first via a call to pamaster = PsychPortAudio('Open', ...); The 'pamaster' corresponds to a real soundcard "
    "present in your system which has been setup for a certain mode of operation, timing mode, frequency etc. "
    "The exception is if you create an AMModulator slave (see below). It can be attached to another playback "
    "slave device by passing the handle of that device as 'pamaster'. "
    "To each real pamaster you can attach multiple slave devices. Each slave represents a subset of the audio channels "
    "of the pamaster device. Each slave can be used as if it is a real independent soundcard for capture or "
    "playback with or without its own set of playback buffers, its own schedule and mode of operation, its own "
    "start and stop times, status etc. Under the hood, the slave device will use the defined subset of channels "
    "of the master device to achieve this illusion of independent soundcards. Typical usage scenarios would be:\n\n"
    "* Mixing of independent soundtracks: Create multiple slave devices, one for each soundtrack. Attach all of "
    "them to the same set of channels on a common master device. Now you can feed each slave with its own soundtrack, "
    "start, stop and schedule playback of each slave/soundtrack independently. Sound from all active channels of "
    "all active slaves will be mixed together and output through the channels of the master device.\n\n"
    "* Independent control of channels: Same as above, but the slaves don't share the same subset of channels "
    "on the master, but each slave chooses a distinctive subset of channels, therefore no mixing will occur."
    "\n\n"
    "'pamaster' is the mandatory pahandle of a master audio device. The master must have been opened via "
    "PsychPortAudio('Open', ..., mode, ...); with the mode flag value '8' included to define it as a master (see 'Open').\n\n"
    "'mode' Mode of operation. See help of 'Open' function for valid settings. However, a mode flag of '8' for "
    "master operation is not allowed for obvious reasons. Also, the given 'mode' flags must be a subset of the "
    "'mode' set on the master device! E.g., you can't specify a audio capture mode flag if the master device isn't "
    "enabled for capture as well. As an exception kPortAudioMonitoring mode is allowed to be present on a slave "
    "while missing on the associated master device.\n"
    "The slave-only mode flag 32 (kPortAudioIsAMModulator) defines a slave device not as a source of audio data, "
    "but as a source of amplitude modulation (AM) envelope data. Its samples don't create sound, but gain-modulate "
    "the sound of other slaves attached to the master to allow precisely timed AM modulation. See the help for "
    "PsychPortAudio('Volume') for more details about AM.\n\n"
    "The slave-only mode flag 64 (kPortAudioIsOutputCapture) defines a slave capture device that captures audio "
    "data from the *output channels* of the master device, i.e., it records the audio stream that is sent to the "
    "speakers. This may be useful for capturing PsychPortAudio's audio output for documentation or debug purposes.\n\n"
    "All slave devices share the same settings for latencymode, timing, sampling frequency "
    "and other low-level tunable parameters, because they operate on the same underlying audio hardware.\n\n"
    "'channels' Define total number of playback and capture channels to use. See help for 'Open?' for explanation. "
    "By default, all channels of the master are used.\n\n"
    "'selectchannels' optional matrix with mappings of logical slave channels to logical master channels: "
    "selectchannels' is a one row by 'channels' matrix with mappings "
    "for pure playback or pure capture. For full-duplex mode (playback and capture), 'selectchannels' must be a "
    "2 rows by max(channels) column matrix. row 1 will define playback channel mappings, whereas row 2 will then "
    "define capture channel mappings. In any case, the number in the i'th column will define which master device "
    "channel will be used for playback or capture of the i'th PsychPortAudio channel (the i'th row of your sound "
    "matrix). Numbering of master device channels starts with one! Example: Both, playback and simultaneous "
    "recording are requested from a master device which represents a 16 channel soundcard with all 16 channels open. "
    "If you'd specify 'selectchannels' as [1, 6 ; 12, 14], then playback would happen to master channels one and six, "
    "sound would be captured from master channels 12 and 14.\n\n";

    static char seeAlsoString[] = "Open Close GetDeviceSettings ";

    int mode, pamaster, paparent, i, numel;
    int* nrchannels;
    int  mynrchannels[2];
    int  m, n, p;
    double* mychannelmap;
    int modeExceptions = 0;
    unsigned int id = PsychPANextHandle();

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(4));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    if (id >= MAX_PSYCH_AUDIO_DEVS) PsychErrorExitMsg(PsychError_user, "Maximum number of simultaneously open audio devices exceeded. You need to close some before you can open new ones.");

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    // Request mandatory pamaster device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &pamaster);

    // Valid device?
    if (pamaster < 0 || pamaster >= MAX_PSYCH_AUDIO_DEVS || audiodevices[pamaster].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio master device handle 'pamaster' provided. No such master device open!");

    // Request optional mode of operation:
    if (PsychCopyInIntegerArg(2, kPsychArgOptional, &mode)) {
        // Mode provided. Validate, but remove potential master flag first:
        mode = mode & (~kPortAudioIsMaster);

        // Being an AM modulator implies being a playback device, as we use the
        // whole audio playback scheduling and buffer facilities:
        if (mode & kPortAudioIsAMModulator) mode |= kPortAudioPlayBack;

        // Being an output capturer implies special rules:
        if (mode & kPortAudioIsOutputCapture) {
            // The associated master must have output that we can capture, ie., it must be configured for playback:
            if (!(audiodevices[pamaster].opmode & kPortAudioPlayBack)) PsychErrorExitMsg(PsychError_user, "Invalid mode: This shall be an output capture device, but master isn't configured for output/playback!");

            // We are a capture device:
            mode |= kPortAudioCapture;

            // Master doesn't need to be a capture device, as we're capturing its output, not its input:
            modeExceptions |= kPortAudioCapture;
        }

        // Is mode a subset of the masters mode as required? One valid exception is kPortAudioMonitoring,
        // which is allowed to be present on the slave but missing on the master, as well as the AMModulator
        // mode and output capturer mode:
        modeExceptions |= kPortAudioMonitoring | kPortAudioIsAMModulator | kPortAudioIsAMModulatorForSlave | kPortAudioIsOutputCapture;
        if (((mode & ~modeExceptions) & audiodevices[pamaster].opmode) != (mode & ~modeExceptions)) PsychErrorExitMsg(PsychError_user, "Invalid mode provided: Mode flags are not a subset of the mode flags of the pamaster device!");
        if ((mode < 1) || ((mode & kPortAudioMonitoring) && ((mode & kPortAudioFullDuplex) != kPortAudioFullDuplex))) {
            PsychErrorExitMsg(PsychError_user, "Invalid mode provided: Outside valid range or invalid combination of flags.");
        }
    }
    else {
        // No mode provided: Use settings from pamaster, minus the master flag of course.
        mode = audiodevices[pamaster].opmode & (~kPortAudioIsMaster);
    }

    // Mode set and validated: Add the slave flag to mark us as slave device:
    mode |= kPortAudioIsSlave;

    // Is our 'pamaster' handle an actual master device?
    if (!(audiodevices[pamaster].opmode & kPortAudioIsMaster)) {
        // No, it is a slave device itself.

        // This is only allowed if we are a modulator device and the target is configured for playback. Error out if this isn't the case:
        if (!(mode & kPortAudioIsAMModulator) || !(audiodevices[pamaster].opmode & kPortAudioPlayBack)) {
            PsychErrorExitMsg(PsychError_user, "Provided device handle 'pamaster' is not configured as a master device. Attaching to a slave is only allowed if it is playback capable and this device is configured as an AM modulator!");
        }

        // We're a modulator for a slave.
        mode |= kPortAudioIsAMModulatorForSlave;
    }

    // Request optional number of channels:
    numel = 0; nrchannels = NULL;
    PsychAllocInIntegerListArg(3, kPsychArgOptional, &numel, &nrchannels);
    if (numel == 0) {
        // No optional channelcount argument provided: Default to pamaster settings:
        mynrchannels[0] = (int) audiodevices[pamaster].outchannels;
        mynrchannels[1] = (int) audiodevices[pamaster].inchannels;
        if (mode & kPortAudioIsOutputCapture) mynrchannels[1] = (int) audiodevices[pamaster].outchannels;
    }
    else if (numel == 1) {
        // One argument provided: Set same count for playback and recording:
        if (*nrchannels < 1 || *nrchannels > MAX_PSYCH_AUDIO_CHANNELS_PER_DEVICE) PsychErrorExitMsg(PsychError_user, "Invalid number of channels provided. Valid values are 1 to device maximum.");
        if (*nrchannels > audiodevices[pamaster].outchannels) PsychErrorExitMsg(PsychError_user, "Invalid number of channels provided. More output channels than pamaster device.");
        if (*nrchannels > audiodevices[pamaster].inchannels) PsychErrorExitMsg(PsychError_user, "Invalid number of channels provided. More input channels than pamaster device.");

        mynrchannels[0] = *nrchannels;
        mynrchannels[1] = *nrchannels;
    }
    else if (numel == 2) {
        // Separate counts for playback and recording provided: Set'em up.
        if (nrchannels[0] < 1 || nrchannels[0] > MAX_PSYCH_AUDIO_CHANNELS_PER_DEVICE || nrchannels[0] > audiodevices[pamaster].outchannels) PsychErrorExitMsg(PsychError_user, "Invalid number of playback channels provided. Valid values are 1 to number of output channels on pamaster device.");
        mynrchannels[0] = nrchannels[0];
        if (mode & kPortAudioIsOutputCapture) {
            // Output capture device:
            if (nrchannels[1] < 1 || nrchannels[1] > MAX_PSYCH_AUDIO_CHANNELS_PER_DEVICE || nrchannels[1] > audiodevices[pamaster].outchannels) PsychErrorExitMsg(PsychError_user, "Invalid number of capture channels provided. Valid values are 1 to number of output channels on pamaster device.");
            mynrchannels[1] = nrchannels[1];
        }
        else {
            // Regular capture device:
            if (nrchannels[1] < 1 || nrchannels[1] > MAX_PSYCH_AUDIO_CHANNELS_PER_DEVICE || nrchannels[1] > audiodevices[pamaster].inchannels) PsychErrorExitMsg(PsychError_user, "Invalid number of capture channels provided. Valid values are 1 to number of capture channels on pamaster device.");
            mynrchannels[1] = nrchannels[1];
        }
    }
    else {
        // More than 2 channel counts provided? Impossible.
        mynrchannels[0] = mynrchannels[1] = 0; // Make compiler happy.
        PsychErrorExitMsg(PsychError_user, "You specified a list with more than two 'channels' entries? Can only be max 2 for playback- and capture.");
    }

    // Make sure that number of capture and playback channels is the same for fast monitoring/feedback mode:
    if ((mode & kPortAudioMonitoring) && (mynrchannels[0] != mynrchannels[1])) PsychErrorExitMsg(PsychError_user, "Fast monitoring/feedback mode selected, but number of capture and playback channels differs! They must be the same for this mode!");

    // Get optional channel map:
    audiodevices[id].outputmappings = NULL;
    audiodevices[id].inputmappings = NULL;
    mychannelmap = NULL;
    PsychAllocInDoubleMatArg(4, kPsychArgOptional, &m, &n, &p, &mychannelmap);

    if (mychannelmap) {
        // Channelmapping provided: Sanity check it.
        if (m<1 || m>2 || p!=1 || (n!=((mynrchannels[0] > mynrchannels[1]) ? mynrchannels[0] : mynrchannels[1]))) {
            PsychErrorExitMsg(PsychError_user, "Invalid size of 'selectchannels' matrix argument: Must be a one- or two row by max(channels) column matrix!");
        }

        // Basic check ok. Build mapping structure:
        if (mode & kPortAudioPlayBack) {
            // Assign playback mappings:
            audiodevices[id].outputmappings = malloc(sizeof(int) * mynrchannels[0]);
            for (i = 0; i < mynrchannels[0]; i++) audiodevices[id].outputmappings[i] = (int) mychannelmap[i * m] - 1;

            if (verbosity > 4) {
                printf("PTB-INFO: Will use the following logical slave-channel -> master-channel mappings for sound output to audio stream %i:\n", id);
                for (i=0; i < mynrchannels[0]; i++) printf("%i --> %i : ", i+1, audiodevices[id].outputmappings[i] + 1);
                printf("\n\n");
            }

            // Revalidate:
            for (i = 0; i < mynrchannels[0]; i++)
                if (audiodevices[id].outputmappings[i] < 0 || audiodevices[id].outputmappings[i] >= audiodevices[pamaster].outchannels) {
                    printf("PTB-ERROR: Slot %i of output channel selection list contains invalid master device channel id %i [Valid between 1 and %i]!\n", i+1, audiodevices[id].outputmappings[i] + 1, (int)  audiodevices[pamaster].outchannels);
                    free(audiodevices[id].outputmappings);
                    audiodevices[id].outputmappings = NULL;
                    PsychErrorExitMsg(PsychError_user, "Invalid items in 'selectchannels' matrix argument!");
                }
        }

        if ((mode & kPortAudioCapture) && !(mode & kPortAudioIsOutputCapture)) {
            // Assign Capture mappings:
            audiodevices[id].inputmappings = malloc(sizeof(int) * mynrchannels[1]);
            // Index into first row of one-row matrix or 2nd row of two-row matrix:
            for (i = 0; i < mynrchannels[1]; i++) audiodevices[id].inputmappings[i] = (int) mychannelmap[(i * m) + (m-1)] - 1;

            if (verbosity > 4) {
                printf("PTB-INFO: Will use the following logical slave-channel -> master-channel mappings for sound capture from audio stream %i:\n", id);
                for (i=0; i < mynrchannels[1]; i++) printf("%i --> %i : ", i+1, audiodevices[id].inputmappings[i] + 1);
                printf("\n\n");
            }

            // Revalidate:
            for (i = 0; i < mynrchannels[1]; i++)
                if (audiodevices[id].inputmappings[i] < 0 || audiodevices[id].inputmappings[i] >= audiodevices[pamaster].inchannels) {
                    printf("PTB-ERROR: Slot %i of capture channel selection list contains invalid master device channel id %i [Valid between 1 and %i]!\n", i+1, audiodevices[id].inputmappings[i] + 1, (int) audiodevices[pamaster].inchannels);
                    free(audiodevices[id].inputmappings);
                    audiodevices[id].inputmappings = NULL;
                    PsychErrorExitMsg(PsychError_user, "Invalid items in 'selectchannels' matrix argument!");
                }
        }

        if (mode & kPortAudioIsOutputCapture) {
            // Assign Output Capture mappings:
            audiodevices[id].inputmappings = malloc(sizeof(int) * mynrchannels[1]);
            // Index into first row of one-row matrix or 2nd row of two-row matrix:
            for (i = 0; i < mynrchannels[1]; i++) audiodevices[id].inputmappings[i] = (int) mychannelmap[(i * m) + (m-1)] - 1;

            if (verbosity > 4) {
                printf("PTB-INFO: Will use the following logical slave-channel -> master-channel mappings for sound capture from audio output stream %i:\n", id);
                for (i=0; i < mynrchannels[1]; i++) printf("%i --> %i : ", i+1, audiodevices[id].inputmappings[i] + 1);
                printf("\n\n");
            }

            // Revalidate:
            for (i = 0; i < mynrchannels[1]; i++)
                if (audiodevices[id].inputmappings[i] < 0 || audiodevices[id].inputmappings[i] >= audiodevices[pamaster].outchannels) {
                    printf("PTB-ERROR: Slot %i of capture channel selection list contains invalid master device channel id %i [Valid between 1 and %i]!\n", i+1, audiodevices[id].inputmappings[i] + 1, (int)  audiodevices[pamaster].outchannels);
                    free(audiodevices[id].inputmappings);
                    audiodevices[id].inputmappings = NULL;
                    PsychErrorExitMsg(PsychError_user, "Invalid items in 'selectchannels' matrix argument!");
                }
        }

    }
    else {
        // No channelmap provided: Assign one that identity-maps all master device channels:
        if (mode & kPortAudioPlayBack) {
            // Assign playback mappings:
            audiodevices[id].outputmappings = malloc(sizeof(int) * mynrchannels[0]);
            for (i = 0; i < mynrchannels[0]; i++) audiodevices[id].outputmappings[i] = i;

            if (verbosity > 4) {
                printf("PTB-INFO: Will use the following logical slave-channel -> master-channel mappings for sound output to audio stream %i:\n", id);
                for (i=0; i < mynrchannels[0]; i++) printf("%i --> %i : ", i+1, audiodevices[id].outputmappings[i]+1);
                printf("\n\n");
            }
        }

        if (mode & kPortAudioCapture) {
            // Assign Capture mappings:
            audiodevices[id].inputmappings = malloc(sizeof(int) * mynrchannels[1]);
            // Index into first row of one-row matrix or 2nd row of two-row matrix:
            for (i = 0; i < mynrchannels[1]; i++) audiodevices[id].inputmappings[i] = i;

            if (verbosity > 4) {
                printf("PTB-INFO: Will use the following logical slave-channel -> master-channel mappings for sound capture from audio stream %i:\n", id);
                for (i=0; i < mynrchannels[1]; i++) printf("%i --> %i : ", i+1, audiodevices[id].inputmappings[i]+1);
                printf("\n\n");
            }
        }
    }

    // Setup our final device structure. Mostly settings from the master device, with
    // a few settings specific to the slave:
    audiodevices[id].opmode = mode;
    audiodevices[id].runMode = 1;
    audiodevices[id].stream = audiodevices[pamaster].stream;
    audiodevices[id].streaminfo = Pa_GetStreamInfo(audiodevices[pamaster].stream);
    audiodevices[id].hostAPI = audiodevices[pamaster].hostAPI;
    audiodevices[id].startTime = 0.0;
    audiodevices[id].reqStartTime = 0.0;
    audiodevices[id].reqStopTime = DBL_MAX;
    audiodevices[id].estStopTime = 0;
    audiodevices[id].currentTime = 0;
    audiodevices[id].state = 0;
    audiodevices[id].reqstate = 255;
    audiodevices[id].repeatCount = 1;
    audiodevices[id].outputbuffer = NULL;
    audiodevices[id].outputbuffersize = 0;
    audiodevices[id].inputbuffer = NULL;
    audiodevices[id].inputbuffersize = 0;
    audiodevices[id].outchannels = mynrchannels[0];
    audiodevices[id].inchannels = mynrchannels[1];
    audiodevices[id].latencyBias = 0.0;
    audiodevices[id].schedule = NULL;
    audiodevices[id].schedule_size = 0;
    audiodevices[id].schedule_pos = 0;
    audiodevices[id].schedule_writepos = 0;
    audiodevices[id].outdeviceidx = audiodevices[pamaster].outdeviceidx;
    audiodevices[id].indeviceidx  = audiodevices[pamaster].indeviceidx;
    audiodevices[id].slaveCount = 0;
    audiodevices[id].slaves = NULL;
    audiodevices[id].pamaster = -1;
    audiodevices[id].modulatorSlave = -1;
    audiodevices[id].slaveOutBuffer = NULL;
    audiodevices[id].slaveGainBuffer = NULL;
    audiodevices[id].slaveInBuffer = NULL;
    audiodevices[id].masterVolume = 1.0;
    audiodevices[id].playposition = 0;
    audiodevices[id].totalplaycount = 0;

    // Setup per-channel output volumes for slave: Each channel starts with a 1.0 setting, ie., max volume:
    if (audiodevices[id].outchannels > 0) {
        audiodevices[id].outChannelVolumes = (float*) malloc(sizeof(float) * (size_t) audiodevices[id].outchannels);
        if (audiodevices[id].outChannelVolumes == NULL) PsychErrorExitMsg(PsychError_outofMemory, "Memory exhausted during audio volume vector allocation.");
        for (i = 0; i < audiodevices[id].outchannels; i++) audiodevices[id].outChannelVolumes[i] = 1.0;
    }
    else {
        audiodevices[id].outChannelVolumes = NULL;
    }

    // If we use locking, we need to initialize the per-device mutex:
    if (uselocking && PsychInitMutex(&(audiodevices[id].mutex))) {
        printf("PsychPortAudio: CRITICAL! Failed to initialize Mutex object for pahandle %i! Prepare for trouble!\n", id);
        PsychErrorExitMsg(PsychError_system, "Audio device mutex creation failed!");
    }

    // If we use locking, this will create & init the associated event variable:
    PsychPACreateSignal(&(audiodevices[id]));

    // Assign parent:
    paparent = pamaster;

    // Remap the meaning of master if we are a modulator for another slave. Our pamaster is actually our parents pamaster:
    if (mode & kPortAudioIsAMModulatorForSlave) pamaster = audiodevices[paparent].pamaster;

    // Attach us to master device: This needs to be done under master mutex protection.
    PsychPALockDeviceMutex(&audiodevices[pamaster]);

    // Maximum number of allowable slaves exceeded?
    if (audiodevices[pamaster].slaveCount >= MAX_PSYCH_AUDIO_SLAVES_PER_DEVICE) {
        // Ouch. Have to abort here:
        PsychPAUnlockDeviceMutex(&audiodevices[pamaster]);
        PsychErrorExitMsg(PsychError_user, "Can't attach slave audio device to specified master. Maximum number of allowable slaves for master exceeded!");
    }

    // Sufficient space for slaves: Find us the first free slot and attach us:
    for (i=0; (i < MAX_PSYCH_AUDIO_SLAVES_PER_DEVICE) && (audiodevices[pamaster].slaves[i] > -1); i++);

    // Attach us to master:
    audiodevices[pamaster].slaves[i] = id;
    audiodevices[pamaster].slaveCount++;

    // Attach master to us:
    audiodevices[id].pamaster = pamaster;

    // If we are a modulator for another slave then attach us as modulator to that slave.
    // No need to mutex-lock the slave, as that protection is already implied by the masters mutex:
    if (mode & kPortAudioIsAMModulatorForSlave) audiodevices[paparent].modulatorSlave = id;

    // Attached :-)
    PsychPAUnlockDeviceMutex(&audiodevices[pamaster]);

    if (verbosity > 4) {
        printf("PTB-INFO: New virtual audio slave device with handle %i opened and attached to parent device handle %i [master %i].\n", id, paparent, pamaster);

        if (audiodevices[id].opmode & kPortAudioIsAMModulator) {
            printf("PTB-INFO: For %i channels amplitude modulation.\n", (int) audiodevices[id].outchannels);
        }
        else if (audiodevices[id].opmode & kPortAudioPlayBack) {
            printf("PTB-INFO: For %i channels playback.\n", (int) audiodevices[id].outchannels);
        }

        if (audiodevices[id].opmode & kPortAudioIsOutputCapture) {
            printf("PTB-INFO: For %i channels capture of master output mix.\n", (int) audiodevices[id].inchannels);
        } else if (audiodevices[id].opmode & kPortAudioCapture) {
            printf("PTB-INFO: For %i channels capture.\n", (int) audiodevices[id].inchannels);
        }

        printf("PTB-INFO: Real samplerate %f Hz. Input latency %f msecs, Output latency %f msecs.\n",
               audiodevices[id].streaminfo->sampleRate, audiodevices[id].streaminfo->inputLatency * 1000.0,
               audiodevices[id].streaminfo->outputLatency * 1000.0);
    }

    // Return device handle:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) id);

    // One more audio device...
    audiodevicecount++;

    return(PsychError_none);
}

/* PsychPortAudio('Close') - Close an audio device via PortAudio.
 */
PsychError PSYCHPORTAUDIOClose(void)
{
    static char useString[] = "PsychPortAudio('Close' [, pahandle]);";
    static char synopsisString[] =
    "Close a PortAudio audio device. The optional 'pahandle' is the handle of the device to close. If pahandle "
    "is omitted, all audio devices will be closed and the driver will shut down.\n";
    static char seeAlsoString[] = "Open GetDeviceSettings ";

    int pahandle= -1;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(1));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(0)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(0));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgOptional, &pahandle);
    if (pahandle == -1) {
        // Full shutdown requested:
        PsychPortAudioExit();
    }
    else {
        // Close one device:
        if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");
        PsychPACloseStream(pahandle);

        // All devices down? Shutdown PortAudio if so:
        if (audiodevicecount == 0) PsychPortAudioExit();
    }

    return(PsychError_none);
}

/* PsychPortAudio('FillBuffer') - Fill audio outputbuffer of a device with data.
 */
PsychError PSYCHPORTAUDIOFillAudioBuffer(void)
{
    static char useString[] = "[underflow, nextSampleStartIndex, nextSampleETASecs] = PsychPortAudio('FillBuffer', pahandle, bufferdata [, streamingrefill=0][, startIndex=Append]);";
    //                          1          2                     3                                                 1         2             3                    4
    static char synopsisString[] =
    "Fill audio data playback buffer of a PortAudio audio device. 'pahandle' is the handle of the device "
    "whose buffer is to be filled.\n"
    #if PSYCH_LANGUAGE == PSYCH_MATLAB
    "'bufferdata' is usually a matrix with audio data in double() or single() format. "
    "Each row of the matrix specifies one sound channel, each column one sample for each channel. "
    #else
    "'bufferdata' is usually a NumPy 2D matrix with audio data in (ideally) float32 format, or also float64 format. "
    "Each column of the matrix specifies one sound channel, each row one sample for each channel. "
    #endif
    "Only floating point values are supported. Samples need to be in range -1.0 to +1.0, with 0.0 for silence. This is "
    "intentionally a very restricted interface. For lowest latency and best timing we want you to provide audio "
    "data exactly at the optimal format and sample rate, so the driver can save computation time and latency for "
    "expensive sample rate conversion, sample format conversion, and bounds checking/clipping.\n"
    "Instead of a matrix, you can also pass in the bufferhandle of an audio buffer as 'bufferdata'. This buffer "
    "must have been created beforehand via PsychPortAudio('CreateBuffer', ...). Its content must satisfy the "
    "same constraints as in case of passing a matrix. The content will be copied from the given buffer "
    "to the standard audio buffer, so it is safe to delete that source buffer if you want.\n"
    "'streamingrefill' optional: If set to 1, ask the driver to refill the buffer immediately while playback "
    "is active. You can think of this as appending the audio data to the audio data already present in the buffer. "
    "This is useful for streaming playback or for creating live audio feedback loops. However, the current implementation "
    "doesn't really append the audio data. Instead it replaces already played audio data with your new data. This means "
    "that if you try to refill more than what has been actually played, this function will wait until enough storage space "
    "is available. A 'streamingrefill' flag of 2 will always refill immediately, ie., without waiting for sufficient buffer "
    "space to become available, even if this causes audible artifacts or some sound data to be overwritten. This is useful "
    "for a few very special audio feedback tricks, only use if you really know what you're doing!\n"
    "It will also fail if you try to refill more than the total buffer capacity. Default is to not do "
    "streaming refills, i.e., the buffer is filled in one batch while playback is stopped. Such a refill will also "
    "reset any playloop setting done via the 'SetLoop' subfunction to the full size of the refilled buffer.\n"
    "If the 'streamingrefill' flag is non-zero and the optional 'startIndex' argument is provided, then the refilling "
    "of the buffer will happen at the provided linear sample index 'startIndex'. If the argument is omitted, new data "
    "will be appended at the end of the current soundbuffers content. The 'startIndex' argument is ignored if no streaming "
    "refill is requested.\n"
    "\nOptionally the function returns the following values:\n"
    "'underflow' A flag: If 1 then the audio buffer underflowed because you didn't refill it in time, ie., some audible "
    "glitches were present in playback and your further playback timing is screwed.\n"
    "'nextSampleStartIndex' This is the absolute index in samples since start of playback of the sample that would "
    "follow after the last sample you added during this 'FillBuffer' call, ie., the first sample during a successive "
    "'FillBuffer' call.\n"
    "'nextSampleETASecs' This value is undefined (NaN) if playback isn't running. During a streaming refill, it contains "
    "the predicted audio onset time in seconds of the sample with index 'nextSampleStartIndex'. Please note that this "
    "prediction can accumulate a prediction error if your buffer is so large that it contains samples that will only "
    "playback far in the future.\n";

    static char seeAlsoString[] = "Open GetDeviceSettings ";

    PsychPABuffer* inbuffer;
    int inbufferhandle = 0;
    float*  indatafloat = NULL;
    psych_bool userfloat = FALSE;
    psych_int64 inchannels, insamples, p;
    size_t buffersize;
    psych_int64 totalplaycount;
    double*    indata = NULL;
    float*  outdata = NULL;
    int pahandle   = -1;
    int streamingrefill = 0;
    int underrun = 0;
    double currentTime, etaSecs;
    psych_int64 startIndex = 0;
    double tBehind = 0.0;
    psych_bool c_layout = PsychUseCMemoryLayoutIfOptimal(TRUE);

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(4));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(2)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(3));    // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgRequired, &pahandle);
    if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");
    if ((audiodevices[pahandle].opmode & kPortAudioPlayBack) == 0) PsychErrorExitMsg(PsychError_user, "Audio device has not been opened for audio playback, so this call doesn't make sense.");

    // Bufferhandle instead of input data matrix provided?
    if (PsychCopyInIntegerArg(2, kPsychArgAnything, &inbufferhandle) && (inbufferhandle > 0)) {
        // Seems so. Double check:
        inbuffer = PsychPAGetAudioBuffer(inbufferhandle);

        // Assign properties:
        inchannels = inbuffer->outchannels;
        insamples  = inbuffer->outputbuffersize / sizeof(float) / inchannels;
        indatafloat = inbuffer->outputbuffer;
    }
    else {
        // Regular double matrix with sound data from runtime?
        if (!PsychAllocInDoubleMatArg64(2, kPsychArgAnything, &inchannels, &insamples, &p, &indata)) {
            // Or regular float matrix instead?
            PsychAllocInFloatMatArg64(2, kPsychArgRequired, &inchannels, &insamples, &p, &indatafloat);
            userfloat = TRUE;
        }

        if (p != 1)
            PsychErrorExitMsg(PsychError_user, "Audio data matrix must be a 2D matrix, but this one is not a 2D matrix!");

        // Swap inchannels <-> insamples to take transposed 2D matrix of C vs. Fortran layout into account:
        if (c_layout) {
            p = inchannels;
            inchannels = insamples;
            insamples = p;
        }
    }

    if (inchannels != audiodevices[pahandle].outchannels) {
        printf("PTB-ERROR: Audio device %i has %i output channels, but provided matrix has non-matching number of %i %s.\n",
               pahandle, (int) audiodevices[pahandle].outchannels, (int) inchannels, (c_layout) ? "columns" : "rows");
        if (c_layout)
            PsychErrorExitMsg(PsychError_user, "Number of columns of audio data matrix doesn't match number of output channels of selected audio device.\n");
        else
            PsychErrorExitMsg(PsychError_user, "Number of rows of audio data matrix doesn't match number of output channels of selected audio device.\n");
    }

    if (insamples < 1) PsychErrorExitMsg(PsychError_user, "You must provide at least 1 sample in your audio buffer!");

    // Get optional streaming refill flag:
    PsychCopyInIntegerArg(3, kPsychArgOptional, &streamingrefill);

    // Full refill or streaming refill?
    if (streamingrefill <= 0) {
        // Standard refill with possible buffer reallocation. Engine needs to be
        // stopped, full reset of engine at refill:

        // Wait for playback on this stream to finish, before refilling it:
        PsychPALockDeviceMutex(&audiodevices[pahandle]);
        while (audiodevices[pahandle].state > 0) {
            // Wait for a state-change before reevaluating the .state:
            PsychPAWaitForChange(&audiodevices[pahandle]);
        }

        // Device is idle, we hold the lock. We can safely drop the lock here and still modify
        // device data, as none of this will get touched by the engine in idle state:
        PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);

        // Ok, everything sane, fill the buffer:
        buffersize = sizeof(float) * (size_t) (inchannels * insamples);
        if (audiodevices[pahandle].outputbuffer && (audiodevices[pahandle].outputbuffersize != buffersize)) {
            free(audiodevices[pahandle].outputbuffer);
            audiodevices[pahandle].outputbuffer = NULL;
            audiodevices[pahandle].outputbuffersize = 0;
        }

        if (audiodevices[pahandle].outputbuffer == NULL) {
            audiodevices[pahandle].outputbuffersize = buffersize;
            audiodevices[pahandle].outputbuffer = (float*) malloc(buffersize);
            if (audiodevices[pahandle].outputbuffer==NULL) PsychErrorExitMsg(PsychError_outofMemory, "Out of system memory when trying to allocate audio buffer.");
        }

        // Reset play position:
        audiodevices[pahandle].playposition = 0;

        outdata = audiodevices[pahandle].outputbuffer;
        if (indata || userfloat) {
            if (indata) {
                // Copy the data, convert it from double to float:
                while(buffersize) {
                    *(outdata++) = (float) (PA_ANTICLAMPGAIN *  *(indata++));
                    buffersize-=sizeof(float);
                }
            } else {
                // Copy the data:
                while(buffersize) {
                    *(outdata++) = (float) (PA_ANTICLAMPGAIN *  *(indatafloat++));
                    buffersize-=sizeof(float);
                }
            }
        }
        else {
            // Data copy from internal audio buffer (already in float format and premultiplied with anti-clamp gain):
            memcpy(outdata, indatafloat, buffersize);
        }

        // Reset write position to end of buffer:
        audiodevices[pahandle].writeposition = (psych_int64) inchannels * insamples;

        // Elapsed count of played out samples must be zero as engine is stopped and will restart sometime after this call:
        totalplaycount = 0;

        // Current playout time undefined when playback is stopped:
        currentTime = PsychGetNanValue();

        // Reset playback loop to full buffer:
        audiodevices[pahandle].loopStartFrame = 0;
        audiodevices[pahandle].loopEndFrame = (audiodevices[pahandle].outputbuffersize / sizeof(float) / audiodevices[pahandle].outchannels) - 1;
    }
    else {
        // Streaming refill while playback is running:

        // Get optional startIndex for new writePosition, if any:
        if (PsychCopyInIntegerArg64(4, kPsychArgOptional, &startIndex)) {
            // New writePosition provided:
            if (startIndex < 0) PsychErrorExitMsg(PsychError_user, "Invalid 'startIndex' provided. Must be greater or equal to zero.");

            // Set writePosition to given startIndex: We can do this without lock held, because writeposition is only
            // touched by us, ie., the main thread, but never by one of the PortAudio threads:
            audiodevices[pahandle].writeposition = (psych_int64) inchannels * (psych_int64) startIndex;
        }

        // Engine stopped? [No need to mutex-lock, as engine can't change from state 0 to other state without our intervention]
        if (audiodevices[pahandle].state == 0) PsychErrorExitMsg(PsychError_user, "Audiodevice not in playback mode! Can't do a streaming buffer refill while stopped.");

        // No buffer allocated? [No need to mutex-lock, see above]
        if (audiodevices[pahandle].outputbuffer == NULL) PsychErrorExitMsg(PsychError_user, "No audio buffer allocated! You must call this method once before start of playback to initially allocate a buffer of sufficient size.");

        // Buffer of sufficient size for a streaming refill of this amount?
        buffersize = sizeof(float) * (size_t) ((psych_int64) inchannels * (psych_int64) insamples);
        if (audiodevices[pahandle].outputbuffersize < (psych_int64) buffersize) PsychErrorExitMsg(PsychError_user, "Total capacity of audio buffer is too small for a refill of this size! Allocate an initial buffer of at least the size of the biggest refill.");

        // Need to lock b'cause of 'playposition':
        PsychPALockDeviceMutex(&audiodevices[pahandle]);

        // Check for buffer underrun:
        if (audiodevices[pahandle].writeposition < audiodevices[pahandle].playposition) {
            underrun = 1;
            tBehind = (double) audiodevices[pahandle].playposition - (double) audiodevices[pahandle].writeposition;
        }

        // Boundary conditions met. Can we refill immediately or do we need to wait for playback
        // position to progress far enough? We skip this test if the streamingrefill flag is > 1:
        while ((streamingrefill < 2) && (audiodevices[pahandle].state > 0) && (!underrun) && (((audiodevices[pahandle].outputbuffersize / (psych_int64) sizeof(float)) - (audiodevices[pahandle].writeposition - audiodevices[pahandle].playposition) - (psych_int64) inchannels) <= (inchannels * insamples))) {
            // Sleep a bit, drop the lock throughout sleep:
            PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);
            // TODO: We could do better here by predicting how long it will take at least until we're ready to refill,
            // but a perfect solution would require quite a bit of effort... ...Something for a really boring afternoon.
            PsychYieldIntervalSeconds(yieldInterval);
            PsychPALockDeviceMutex(&audiodevices[pahandle]);

            // Recheck for buffer underrun:
            if (audiodevices[pahandle].writeposition < audiodevices[pahandle].playposition) {
                underrun = 1;
                tBehind = (double) audiodevices[pahandle].playposition - (double) audiodevices[pahandle].writeposition;
            }
        }

        // Exit with lock held...

        // Have we left the while-loop because the engine stopped? In that case we won't
        // be able to ever get the needed headroom and need to error-out:
        if (audiodevices[pahandle].state == 0) {
            // Ohoh...
            PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);
            PsychErrorExitMsg(PsychError_user, "Audiodevice no longer in playback mode (Auto stopped?!?)! Can't continue a streaming buffer refill while stopped. Check your code!");
        }

        // Ok, device locked and enough headroom for batch streaming refill:

        // Copy the data, convert it from double to float, take ringbuffer wraparound into account:
        if (indata || userfloat) {
            if (indata) {
                while(buffersize > 0) {
                    // Fetch next sample and copy it to matrix:
                    audiodevices[pahandle].outputbuffer[(audiodevices[pahandle].writeposition % (audiodevices[pahandle].outputbuffersize / sizeof(float)))] = (float) (PA_ANTICLAMPGAIN *  *(indata++));

                    // Update sample write counter:
                    audiodevices[pahandle].writeposition++;

                    // Decrement copy counter:
                    buffersize-=sizeof(float);
                }
            }
            else {
                while(buffersize > 0) {
                    // Fetch next sample and copy it to matrix:
                    audiodevices[pahandle].outputbuffer[(audiodevices[pahandle].writeposition % (audiodevices[pahandle].outputbuffersize / sizeof(float)))] = (float) (PA_ANTICLAMPGAIN *  *(indatafloat++));

                    // Update sample write counter:
                    audiodevices[pahandle].writeposition++;

                    // Decrement copy counter:
                    buffersize-=sizeof(float);
                }
            }
        }
        else {
            // Data copy from internal audio buffer (already in float format and premultiplied with anti-clamp gain):
            while(buffersize > 0) {
                // Fetch next sample and copy it to matrix:
                audiodevices[pahandle].outputbuffer[(audiodevices[pahandle].writeposition % (audiodevices[pahandle].outputbuffersize / sizeof(float)))] = *(indatafloat++);

                // Update sample write counter:
                audiodevices[pahandle].writeposition++;

                // Decrement copy counter:
                buffersize-=sizeof(float);
            }
        }

        // Retrieve total count of played out samples from engine:
        totalplaycount = audiodevices[pahandle].totalplaycount;

        // Retrieve corresponding timestamp of last playout:
        currentTime = audiodevices[pahandle].currentTime;

        // Check for buffer underrun:
        if (audiodevices[pahandle].writeposition < audiodevices[pahandle].playposition) {
            underrun = 1;
            tBehind = (double) audiodevices[pahandle].playposition - (double) audiodevices[pahandle].writeposition;
        }

        // Drop lock here, no longer needed:
        PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);

        if ((underrun > 0) && (verbosity > 1)) {
            printf("PsychPortAudio-WARNING: Underrun of audio playback buffer detected during streaming refill at approximate play position %f secs [%f msecs behind]. Sound will be skipped, timing may be wrong and audible glitches may occur!\n",
                   ((double) audiodevices[pahandle].playposition / ((double) audiodevices[pahandle].outchannels * (double) audiodevices[pahandle].streaminfo->sampleRate)) , tBehind / ((double) audiodevices[pahandle].outchannels * (double) audiodevices[pahandle].streaminfo->sampleRate) * 1000.0);
        }
    }

    // Copy out underrun flag:
    PsychCopyOutDoubleArg(1, FALSE, (double) underrun);

    // Copy out number of submitted sample frames, aka the absolute sample frame startindex for next 'Fillbuffer' call.
    PsychCopyOutDoubleArg(2, FALSE, (double) audiodevices[pahandle].writeposition / (double) audiodevices[pahandle].outchannels);

    // Compute and return predicted estimated time of arrival for playout of first sample for next 'FillBuffer' call, aka almost of the last sample in this 'FillBuffer' call:
    etaSecs = currentTime + ((((double) audiodevices[pahandle].writeposition - (double) totalplaycount) / (double) audiodevices[pahandle].outchannels) / (double) audiodevices[pahandle].streaminfo->sampleRate);
    PsychCopyOutDoubleArg(3, FALSE, etaSecs);

    // Buffer ready.
    return(PsychError_none);
}

/* PsychPortAudio('RefillBuffer') - Refill existing audio outputbuffer of a device with data.
 */
PsychError PSYCHPORTAUDIORefillBuffer(void)
{
    static char useString[] = "PsychPortAudio('RefillBuffer', pahandle [, bufferhandle=0], bufferdata [, startIndex=0]);";
    static char synopsisString[] =
    "Refill part of an audio data playback buffer of a PortAudio audio device. 'pahandle' is the handle of the device "
    "whose buffer is to be filled.\n"
    "'bufferhandle' is the handle of the buffer: Use a handle of zero for the standard "
    "buffer created and accessed via 'FillBuffer'.\n"
    #if PSYCH_LANGUAGE == PSYCH_MATLAB
    "'bufferdata' is a matrix with audio data in double() or single() format. "
    "Each row of the matrix specifies one sound channel, each column one sample for each channel. "
    #else
    "'bufferdata' is usually a NumPy 2D matrix with audio data in (ideally) float32 format, or also float64 format. "
    "Each column of the matrix specifies one sound channel, each row one sample for each channel. "
    #endif
    "Only floating point values are supported. Samples need to be in range -1.0 to +1.0, with 0.0 for silence. This is "
    "intentionally a very restricted interface. For lowest latency and best timing we want you to provide audio "
    "data exactly at the optimal format and sample rate, so the driver can save computation time and latency for "
    "expensive sample rate conversion, sample format conversion, and bounds checking/clipping.\n"
    "Instead of a matrix, you can also pass in the bufferhandle of an audio buffer as 'bufferdata'. This buffer "
    "must have been created beforehand via PsychPortAudio('CreateBuffer', ...). Its content must satisfy the "
    "same constraints as in case of passing a matrix. The content will be copied from the given buffer "
    "to the standard audio buffer, so it is safe to delete that source buffer if you want.\n"
    "'startIndex' optional: Defines the first sample frame within the buffer where refill should start. "
    "By default, refilling starts at the beginning of the buffer - at sample frame 0. 'startIndex' allows to "
    "start refilling at some offset.\n"
    "Please note that 'RefillBuffer' can't resize an existing buffer - you can't fill in more data than the "
    "current buffer capacity permits. If you want to add more sound, you'll need to use 'FillBuffer' or "
    "create a new buffer of proper capacity.\n"
    "'RefillBuffer' can be used any time on a buffer, even if the buffer is currently playing, allowing for "
    "on-the-fly replacement of content. However be careful to avoid the currently played section, or you'll "
    "hear audio artifacts. For streaming out audio content in a glitch-free way, you may want to use the "
    "'streamingrefill' option of the 'FillBuffer' subfunction instead.\n";

    static char seeAlsoString[] = "Open FillBuffer GetStatus ";

    PsychPABuffer* buffer = NULL;
    PsychPABuffer* inbuffer = NULL;

    psych_int64 inchannels, insamples, p;
    size_t buffersize, outbuffersize;
    double*    indata = NULL;
    int inbufferhandle = 0;
    float*  indatafloat = NULL;
    psych_bool userfloat = FALSE;
    float*  outdata = NULL;
    int pahandle   = -1;
    int bufferhandle = 0;
    psych_int64 startIndex = 0;
    psych_bool c_layout = PsychUseCMemoryLayoutIfOptimal(TRUE);

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(4));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(3)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(0));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgRequired, &pahandle);
    if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");
    if ((audiodevices[pahandle].opmode & kPortAudioPlayBack) == 0) PsychErrorExitMsg(PsychError_user, "Audio device has not been opened for audio playback, so this call doesn't make sense.");

    PsychCopyInIntegerArg(2, kPsychArgOptional, &bufferhandle);
    // Check for valid bufferhandle:
    if (bufferhandle < 0) PsychErrorExitMsg(PsychError_user, "Invalid audio 'bufferhandle' provided.");

    // If it is a non-zero handle, try to dereference from dynamic buffer:
    if (bufferhandle > 0) {
        // Deref bufferHandle: Issue error if no buffer with such a handle exists:
        buffer = PsychPAGetAudioBuffer(bufferhandle);

        // Validate matching output channel count:
        if (buffer->outchannels != audiodevices[pahandle].outchannels) {
            printf("PsychPortAudio-ERROR: Audio channel count %i of audiobuffer with handle %i doesn't match channel count %i of audio device!\n",
                   (int) buffer->outchannels, bufferhandle, (int) audiodevices[pahandle].outchannels);
            PsychErrorExitMsg(PsychError_user, "Target audio buffer 'bufferHandle' has an audio channel count that doesn't match channels of audio device!");
        }
    }

    // Bufferhandle instead of input data matrix provided?
    if (PsychCopyInIntegerArg(3, kPsychArgAnything, &inbufferhandle) && (inbufferhandle > 0)) {
        // Seems so. Double check:
        inbuffer = PsychPAGetAudioBuffer(inbufferhandle);

        // Assign properties:
        inchannels = inbuffer->outchannels;
        insamples = inbuffer->outputbuffersize / sizeof(float) / inchannels;
        indatafloat = inbuffer->outputbuffer;
    }
    else {
        // Regular double matrix with sound data from runtime:
        if (!PsychAllocInDoubleMatArg64(3, kPsychArgAnything, &inchannels, &insamples, &p, &indata)) {
            // Or regular float matrix instead:
            PsychAllocInFloatMatArg64(3, kPsychArgRequired, &inchannels, &insamples, &p, &indatafloat);
            userfloat = TRUE;
        }

        if (p != 1)
            PsychErrorExitMsg(PsychError_user, "Audio data matrix must be a 2D matrix, but this one is not a 2D matrix!");

        // Swap inchannels <-> insamples to take transposed 2D matrix of C vs. Fortran layout into account:
        if (c_layout) {
            p = inchannels;
            inchannels = insamples;
            insamples = p;
        }
    }

    if (inchannels != audiodevices[pahandle].outchannels) {
        printf("PTB-ERROR: Audio device %i has %i output channels, but provided matrix has non-matching number of %i %s.\n",
               pahandle, (int) audiodevices[pahandle].outchannels, (int) inchannels, (c_layout) ? "columns" : "rows");
        if (c_layout)
            PsychErrorExitMsg(PsychError_user, "Number of columns of audio data matrix doesn't match number of output channels of selected audio device.\n");
        else
            PsychErrorExitMsg(PsychError_user, "Number of rows of audio data matrix doesn't match number of output channels of selected audio device.\n");
    }

    if (insamples < 1) PsychErrorExitMsg(PsychError_user, "You must provide at least 1 sample for refill of your audio buffer!");

    // Get optional startIndex:
    PsychCopyInIntegerArg64(4, kPsychArgOptional, &startIndex);
    if (startIndex < 0) PsychErrorExitMsg(PsychError_user, "Invalid 'startIndex' provided. Must be greater or equal to zero.");

    // Assign bufferpointer based on bufferhandle:
    if (bufferhandle > 0) {
        // Generic buffer:
        outdata = buffer->outputbuffer;
        outbuffersize = (size_t) buffer->outputbuffersize;
    }
    else {
        // Standard playout buffer:
        outdata = audiodevices[pahandle].outputbuffer;
        outbuffersize = (size_t) audiodevices[pahandle].outputbuffersize;
    }

    // Buffer exists?
    if (outdata == NULL) PsychErrorExitMsg(PsychError_user, "No such buffer with given 'bufferhandle', or buffer not yet created!");

    // Compute required buffersize for copying all data from given startIndex:
    buffersize = sizeof(float) * (size_t) inchannels * ((size_t) insamples + (size_t) startIndex);

    // Buffer of sufficient size?
    if (buffersize > outbuffersize) {
        // Nope, too small: Adapt 'buffersize' to allowable maximum amount:
        if (verbosity > 1) printf("PsychPortAudio: WARNING: In 'RefillBuffer' for bufferhandle %i at startindex %i: Insufficient\nbuffersize %i for %i new audioframes starting at given startindex.\nWill truncate to maximum possible.\n", bufferhandle, (int) startIndex, (int) (outbuffersize / (sizeof(float) * inchannels)), (int) insamples);
        buffersize = outbuffersize;
        buffersize -= sizeof(float) * (size_t) inchannels * (size_t) startIndex;
    }
    else {
        // Big enough:
        buffersize = sizeof(float) * (size_t) inchannels * (size_t) insamples;
    }

    // Map startIndex to offset in buffer:
    outdata += (size_t) inchannels * (size_t) startIndex;

    // Ok, everything sane, fill the buffer: 'buffersize' iterations into 'outdata':
    //fprintf(stderr, "buffersize = %i\n", buffersize);

    if (indata || userfloat) {
        if (indata) {
            // Copy the data, convert it from double to float:
            while(buffersize > 0) {
                *(outdata++) = (float) (PA_ANTICLAMPGAIN *  *(indata++));
                buffersize-=sizeof(float);
            }
        }
        else {
            // Copy the data:
            while(buffersize > 0) {
                *(outdata++) = (float) (PA_ANTICLAMPGAIN *  *(indatafloat++));
                buffersize-=sizeof(float);
            }
        }
    }
    else {
        // Data copy from internal audio buffer (already in float format and premultiplied with anti-clamp gain):
        memcpy(outdata, indatafloat, buffersize);
    }

    // Done.
    return(PsychError_none);
}

/* PsychPortAudio('DeleteBuffer') - Delete a dynamic audio outputbuffer.
 */
PsychError PSYCHPORTAUDIODeleteBuffer(void)
{
    static char useString[] = "result = PsychPortAudio('DeleteBuffer'[, bufferhandle] [, waitmode]);";
    static char synopsisString[] =
    "Delete an existing dynamic audio data playback buffer.\n"
    "'bufferhandle' is the handle for the buffer to delete. If it is omitted, all "
    "buffers will be deleted. 'waitmode' defines what happens if a buffer shall be "
    "deleted that is currently in use, i.e., part of the audio playback schedule "
    "of an active audio device. The default of zero will simply return without deleting "
    "the buffer. A setting of 1 will wait until the buffer can be safely deleted.\n";

    static char seeAlsoString[] = "Open FillBuffer GetStatus ";

    int bufferhandle = 0;
    int waitmode = 0;
    int rc = 0;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(2));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(0)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    // Get optional waitmode with default zero:
    PsychCopyInIntegerArg(2, FALSE, &waitmode);

    if (PsychCopyInIntegerArg(1, FALSE, &bufferhandle)) {
        // Specific bufferhandle provided for deletion:

        // Check if handle valid and error out if so:
        (void) PsychPAGetAudioBuffer(bufferhandle);

        // Valid: Try to delete buffer:
        rc = PsychPADeleteAudioBuffer(bufferhandle, waitmode);
    }
    else {
        // No specific handle: Try to delete all buffers:
        if (PsychPAUpdateBufferReferences()) {
            // At least one buffer locked. What do do?
            if (waitmode == 0) {
                // Just fail -> No op.
                rc = 0;
            }
            else {
                // Retry until it works:
                while (PsychPAUpdateBufferReferences()) PsychYieldIntervalSeconds(yieldInterval);
                rc = 1;
            }
        }
        else {
            rc = 1;
        }

        // Really delete all buffers if rc == 1:
        if (rc == 1) PsychPADeleteAllAudioBuffers();
    }

    // Return status:
    PsychCopyOutDoubleArg(1, FALSE, (double) rc);

    // Done.
    return(PsychError_none);
}

/* PsychPortAudio('CreateBuffer') - Create and fill dynamic audio outputbuffer of a device with data.
 */
PsychError PSYCHPORTAUDIOCreateBuffer(void)
{
    static char useString[] = "bufferhandle = PsychPortAudio('CreateBuffer' [, pahandle], bufferdata);";
    static char synopsisString[] =
    "Create a new dynamic audio data playback buffer for a PortAudio audio device and fill it with initial data.\n"
    "Return a 'bufferhandle' to the new buffer. 'pahandle' is the optional handle of the device "
    "whose buffer is to be filled.\n"
    #if PSYCH_LANGUAGE == PSYCH_MATLAB
    "'bufferdata' is a matrix with audio data in double() or single() format. "
    "Each row of the matrix specifies one sound channel, each column one sample for each channel. "
    #else
    "'bufferdata' is usually a NumPy 2D matrix with audio data in (ideally) float32 format, or also float64 format. "
    "Each column of the matrix specifies one sound channel, each row one sample for each channel. "
    #endif
    "Only floating point values are supported. Samples need to be in range -1.0 to +1.0, with 0.0 for silence. This is "
    "intentionally a very restricted interface. For lowest latency and best timing we want you to provide audio "
    "data exactly at the optimal format and sample rate, so the driver can save computation time and latency for "
    "expensive sample rate conversion, sample format conversion, and bounds checking/clipping.\n\n"
    "You can refill the buffer anytime via the PsychPortAudio('RefillBuffer') call.\n"
    "You can delete the buffer via the PsychPortAudio('DeleteBuffer') call, once it is not used anymore. \n"
    "You can attach the buffer to an audio playback schedule for actual audio playback via the "
    "PsychPortAudio('AddToSchedule') call.\n"
    "The same buffer can be attached to and used by multiple audio devices simultaneously, or multiple "
    "times within one or more playback schedules. ";

    static char seeAlsoString[] = "Open FillBuffer GetStatus ";

    PsychPABuffer* buffer;
    psych_int64 inchannels, insamples, p;
    size_t buffersize;
    double*    indata = NULL;
    float* indatafloat = NULL;
    float*  outdata = NULL;
    int pahandle   = -1;
    int bufferhandle = 0;
    psych_bool c_layout = PsychUseCMemoryLayoutIfOptimal(TRUE);

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(2));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(0)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    // Get data matrix with initial buffer content:
    if (!PsychAllocInDoubleMatArg64(2, kPsychArgAnything, &inchannels, &insamples, &p, &indata)) {
        // Or regular float matrix instead:
        PsychAllocInFloatMatArg64(2, kPsychArgRequired, &inchannels, &insamples, &p, &indatafloat);
    }

    if (p != 1)
        PsychErrorExitMsg(PsychError_user, "Audio data matrix must be a 2D matrix, but this one is not a 2D matrix!");

    // Swap inchannels <-> insamples to take transposed 2D matrix of C vs. Fortran layout into account:
    if (c_layout) {
        p = inchannels;
        inchannels = insamples;
        insamples = p;
    }

    // If the optional pahandle is provided...
    if (PsychCopyInIntegerArg(1, kPsychArgOptional, &pahandle)) {
        // ...then we use it to validate the configuration of the datamatrix for the new
        // buffer against the requirements of that audiodevice pahandle:
        if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");
        if ((audiodevices[pahandle].opmode & kPortAudioPlayBack) == 0) PsychErrorExitMsg(PsychError_user, "Audio device has not been opened for audio playback, so this call doesn't make sense.");

        if (inchannels != audiodevices[pahandle].outchannels) {
            printf("PTB-ERROR: Audio device %i has %i output channels, but provided matrix has non-matching number of %i %s.\n",
                pahandle, (int) audiodevices[pahandle].outchannels, (int) inchannels, (c_layout) ? "columns" : "rows");
            if (c_layout)
                PsychErrorExitMsg(PsychError_user, "Number of columns of audio data matrix doesn't match number of output channels of selected audio device.\n");
            else
                PsychErrorExitMsg(PsychError_user, "Number of rows of audio data matrix doesn't match number of output channels of selected audio device.\n");
        }
    }

    if (inchannels < 1) PsychErrorExitMsg(PsychError_user, "You must provide at least a vector for creation of at least one audio channel in your audio buffer!");
    if (insamples < 1) PsychErrorExitMsg(PsychError_user, "You must provide at least 1 sample for creation of your audio buffer!");

    // Create buffer and assign bufferhandle:
    bufferhandle = PsychPACreateAudioBuffer(inchannels, insamples);

    // Deref bufferHandle:
    buffer = PsychPAGetAudioBuffer(bufferhandle);
    outdata = buffer->outputbuffer;
    buffersize = sizeof(float) * (size_t) inchannels * (size_t) insamples;

    if (indata) {
        // Copy the data, convert it from double to float:
        while(buffersize > 0) {
            *(outdata++) = (float) (PA_ANTICLAMPGAIN *  *(indata++));
            buffersize-=sizeof(float);
        }
    } else {
        // Copy the float data:
        while(buffersize > 0) {
            *(outdata++) = (float) (PA_ANTICLAMPGAIN *  *(indatafloat++));
            buffersize-=sizeof(float);
        }
    }

    // Return bufferhandle:
    PsychCopyOutDoubleArg(1, FALSE, (double) bufferhandle);

    // Done.
    return(PsychError_none);
}

/* PsychPortAudio('GetAudioData') - Retrieve captured audio data.
 */
PsychError PSYCHPORTAUDIOGetAudioData(void)
{
    #if PSYCH_LANGUAGE == PSYCH_MATLAB
    static char useString[] = "[audiodata absrecposition overflow cstarttime] = PsychPortAudio('GetAudioData', pahandle [, amountToAllocateSecs][, minimumAmountToReturnSecs][, maximumAmountToReturnSecs][, singleType=0]);";
    #else
    static char useString[] = "[audiodata absrecposition overflow cstarttime] = PsychPortAudio('GetAudioData', pahandle [, amountToAllocateSecs][, minimumAmountToReturnSecs][, maximumAmountToReturnSecs][, singleType=1]);";
    #endif
    static char synopsisString[] =
    "Retrieve captured audio data from a audio device. 'pahandle' is the handle of the device "
    "whose data is to be retrieved.\n"
    #if PSYCH_LANGUAGE == PSYCH_MATLAB
    "'audiodata' is a matrix with audio data in floating point format. "
    "Each row of the matrix returns one sound channel, each column one sample for each channel. "
    #else
    "'audiodata' is a NumPy 2D matrix with audio data in float32 floating point format by default. "
    "Each column of the matrix returns one sound channel, each row one sample for each channel. "
    #endif
    "Returned samples are in range -1.0 to +1.0, with 0.0 for silence. This is "
    "intentionally a very restricted interface. For lowest latency and best timing we want you to accept audio "
    "data exactly at the optimal format and sample rate, so the driver can save computation time and latency for "
    "expensive sample rate conversion, sample format conversion, and bounds checking/clipping.\n"
    "You must call this function once before start of capture operations to allocate an internal buffer "
    "that stores captured audio data inbetween your periodic calls. Provide 'amountToAllocateSecs' as "
    "requested buffersize in seconds. After start of capture you must call this function periodically "
    "at least every 'amountToAllocateSecs' seconds to drain the internal buffer into your "
    "matrix 'audiodata'. If you fail to call the function frequently enough, sound data will get lost!\n"
    "'minimumAmountToReturnSecs' optional minimum amount of recorded data to return at each call. The "
    "driver will only return control to your script when it was able to collect at least that amount "
    "of seconds of sound data - or if the capture engine was stopped. If you don't set this parameter, "
    "the driver will return immediately, giving you whatever amount of sound data was available - including "
    "an empty matrix if nothing was available.\n"
    "'maximumAmountToReturnSecs' allows you to optionally restrict the amount of returned sound data to "
    "a specific duration in seconds. By default, you'll get whatever is available.\n"
    "If you provide both, 'minimumAmountToReturnSecs' and 'maximumAmountToReturnSecs' and set them to equal "
    "values (but significantly lower than the 'amountToAllocateSecs' buffersize!!) then you'll always "
    "get an 'audiodata' matrix back that is of a fixed size. This may be convenient for postprocessing "
    "in the scripting language. It may also reduce or avoid memory fragmentation...\n"
    #if PSYCH_LANGUAGE == PSYCH_MATLAB
    "'singleType' if set to 1 will return a sound data matrix of single() type instead of double() type. "
    "By default, double() type is returned. single() type matrices only consume half as much memory as "
    "double() type matrices, without loss of audio precision for up to 24-Bit ADC hardware.\n"
    #else
    "'singleType' if set to 1 will return a sound data matrix of float32 type instead of float64 type. "
    "By default, float32 type is returned, as float32 matrices only consume half as much memory as "
    "float64 matrices, usually without loss of audio precision for up to 24-Bit ADC hardware.\n"
    #endif
    "\n"
    "\nOptional return arguments other than 'audiodata':\n\n"
    "'absrecposition' is the absolute position (in samples) of the first column in the returned data matrix, "
    "assuming that sample zero was the very first recorded sample in this session. The count is reset each time "
    "you start a new capture session via call to PsychPortAudio('Start').\n"
    "Each call to this function will return a new chunk of recorded sound data. The 'absrecposition' provides "
    "you with absolute matrix column indices to stitch together the results of all calls into one seamless "
    "recording if you want. 'overflow' if this flag is zero then everything went fine. If it is one then you "
    "didn't manage to call this function frequent enough, the capacity of the internal recording buffer was "
    "exceeded and therefore you lost captured sound data, i.e., there is a gap in your recording. When "
    "initially allocating the internal buffer, make sure to allocate it big enough so it is able to easily "
    "store all recorded data inbetween your calls to 'GetAudioData'. Example: You expect to call this routine "
    "once every second in your trial loop, then allocate a sound buffer of at least 2 seconds for some security "
    "headroom. If you know that the recording time of each recording has an upper bound then you can allocate "
    "an internal buffer of sufficient size and fetch the buffer all at once at the end of a recording.\n"
    "'ctsstarttime' this is an estimate of the system time (in seconds) when the very first sample of this "
    "recording was captured by the sound input of your hardware. This is only a rough estimate, not to be "
    "trusted down to the millisecond level, at least not without former careful calibration of your setup!\n";

    static char seeAlsoString[] = "Open GetDeviceSettings ";

    //int inchannels, insamples, p, maxSamples;
    psych_int64 insamples, maxSamples;
    size_t buffersize;
    double*    indata = NULL;
    float*  indatafloat = NULL;
    int pahandle   = -1;
    double allocsize;
    double minSecs, maxSecs, minSamples;
    int overrun = 0;
    int singleType = (PSYCH_LANGUAGE == PSYCH_MATLAB) ? 0 : 1;
    psych_bool c_layout = PsychUseCMemoryLayoutIfOptimal(TRUE);

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(5));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(4));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgRequired, &pahandle);
    if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");
    if ((audiodevices[pahandle].opmode & kPortAudioCapture) == 0) PsychErrorExitMsg(PsychError_user, "Audio device has not been opened for audio capture, so this call doesn't make sense.");

    buffersize = (size_t) audiodevices[pahandle].inputbuffersize;

    // Copy in optional amount of buffer memory to allocate for internal recording ringbuffer:
    allocsize = 0;
    PsychCopyInDoubleArg(2, kPsychArgOptional, &allocsize);

    // Internal capture ringbuffer already allocated?
    if (buffersize == 0) {
        // Nope. This call needs to allocate it. We know the engine is idle in this case, because it
        // can't be started or brought out of idle if buffersize == 0 due to checks in 'Start' and 'RescheduleStart'.
        if (allocsize <= 0) PsychErrorExitMsg(PsychError_user, "You must first call this function with a positive 'amountToAllocateSecs' argument to allocate internal bufferspace first!");

        // Ready to allocate ringbuffer memory outside this if-clause...
    }
    else {
        // Buffer already allocated. Need to realloc?
        if (allocsize > 0) {
            // Calling script wants to reallocate the buffer. Check if this is possible here:

            // Test 1: Engine running? [No need to mutex lock: If running, we fail safely. If state = 0 it can't transition to > 0 behind our back]
            if (audiodevices[pahandle].state > 0) PsychErrorExitMsg(PsychError_user, "Tried to resize internal buffer while recording engine is running! You must stop recording before resizing the buffer!");

            // Test 2: Pending samples to read from current ringbuffer? Engine is idle, so we can safely access device data lock-free...
            if (audiodevices[pahandle].readposition < audiodevices[pahandle].recposition) PsychErrorExitMsg(PsychError_user, "Tried to resize internal buffer without emptying it beforehand. You must drain the buffer before resizing it!");

            // Ok, reallocation allowed, as engine is idle. Delete old buffer:
            audiodevices[pahandle].inputbuffersize = 0;
            free(audiodevices[pahandle].inputbuffer);
            audiodevices[pahandle].inputbuffer = NULL;

            // At this point we are ready to re-allocate ringbuffer outside this if-clause...
        }
    }

    // Still (re-)allocation wanted?
    if (allocsize > 0) {
        // We know the engine is idle if we reach this point, so no need to acquire locks or check state...

        // Calculate needed buffersize in samples: Convert allocsize in seconds to size in bytes:
        audiodevices[pahandle].inputbuffersize = sizeof(float) * ((psych_int64) (allocsize * audiodevices[pahandle].streaminfo->sampleRate)) * audiodevices[pahandle].inchannels;
        audiodevices[pahandle].inputbuffer = (float*) calloc(1, (size_t) audiodevices[pahandle].inputbuffersize);
        if (audiodevices[pahandle].inputbuffer == NULL) PsychErrorExitMsg(PsychError_outofMemory, "Free system memory exhausted when trying to allocate audio recording buffer!");

        // This was an (re-)allocation call, so no data is pending in the buffer.
        // Therefore we don't return any data, just reset the counters:
        audiodevices[pahandle].recposition = 0;
        audiodevices[pahandle].readposition = 0;
        return(PsychError_none);
    }

    // This is not an allocation call, but a real data fetch call:

    // Get optional "minimum amount to return" argument:
    // We default to "whatever we can get" ie. zero seconds.
    minSecs = 0;
    PsychCopyInDoubleArg(3, kPsychArgOptional, &minSecs);

    // Get optional "maximum amount to return" argument:
    // We default to "whatever we can get" ie. infinite seconds.
    maxSecs = 0;
    PsychCopyInDoubleArg(4, kPsychArgOptional, &maxSecs);

    // Get optional singleType flag:
    PsychCopyInIntegerArg(5, kPsychArgOptional, &singleType);
    if (singleType < 0 || singleType > 1) PsychErrorExitMsg(PsychError_user, "'singleType' flag must be zero or one!");

    // The engine is potentially running, so we need to mutex-lock our accesses...
    PsychPALockDeviceMutex(&audiodevices[pahandle]);

    // How much samples are available in ringbuffer to fetch?
    insamples = (psych_int64) (audiodevices[pahandle].recposition - audiodevices[pahandle].readposition);

    // Convert amount of available data into seconds and check if our minimum
    // requirements are fulfilled:
    if (minSecs > 0) {
        // Convert seconds to samples:
        minSamples = minSecs * ((double) audiodevices[pahandle].streaminfo->sampleRate) * ((double) audiodevices[pahandle].inchannels) + ((double) audiodevices[pahandle].inchannels);

        // Bigger than buffersize? That would be a no no...
        if (((psych_int64) (minSamples * sizeof(float))) > audiodevices[pahandle].inputbuffersize) {
            PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);
            PsychErrorExitMsg(PsychError_user, "Invalid 'minimumAmountToReturnSecs' parameter: The requested minimum is bigger than the whole capture buffer size!'");
        }

        // Loop until either request is fullfillable or the device gets stopped - in which
        // case we'll never be able to fullfill the request...
        while (((double) insamples < minSamples) && (audiodevices[pahandle].state > 0)) {
            // Compute amount of time to elapse before request could be fullfilled:
            minSecs = (minSamples - (double) insamples) / ((double) audiodevices[pahandle].inchannels) / ((double) audiodevices[pahandle].streaminfo->sampleRate);
            // Ok, required data will be available earliest in 'minSecs' seconds. Sleep until then with lock dropped:
            PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);
            PsychWaitIntervalSeconds(minSecs);
            PsychPALockDeviceMutex(&audiodevices[pahandle]);

            // We've slept at least the estimated amount of required time. Recalculate amount
            // of available sound data and check again...
            insamples = (psych_int64) (audiodevices[pahandle].recposition - audiodevices[pahandle].readposition);
        }
    }

    // Lock held here...

    // Never ever fetch the samples for the last sampleframe. We do not want to fetch
    // a possibly not yet updated or incomplete sample frame. Leave this to next call
    // of this function. Well, unless state is zero == engine stopped. In that case we
    // know that the playhead won't move anymore and we can safely fetch all remaining
    // data.
    if (audiodevices[pahandle].state > 0) {
        insamples = insamples - (insamples % audiodevices[pahandle].inchannels);
        insamples-= audiodevices[pahandle].inchannels;
    }

    // Can unlock here: The remainder of the routine doesn't touch any critical device variables anymore,
    // only variables that aren't modified by the engine, or not used/touched by engine.
    // Well, theoretically the engine could overwrite the portion of the buffer we're going to
    // read out if we stall massively and the capturebuffer is way too "undersized", but in that
    // case the user code is screwed anyway and it (or the system) needs to be fixed...
    PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);

    insamples = (insamples < 0) ? 0 : insamples;
    buffersize = (size_t) insamples * sizeof(float);

    // Buffer "overflow" detected?
    if ((psych_int64) buffersize > audiodevices[pahandle].inputbuffersize) {
        // Ok, the buffer did overrun and captured data was lost. Limit returned data
        // to buffersize and set the overrun flag, optionally output a warning:
        buffersize = (size_t) audiodevices[pahandle].inputbuffersize;
        insamples = buffersize / sizeof(float);

        // Set overrun flag:
        overrun = 1;

        if (verbosity > 1) printf("PsychPortAudio-WARNING: Overflow of audio capture buffer detected. Some sound data will be lost!\n");
    }

    // Limitation of returned amount of data wanted?
    if (maxSecs > 0) {
        // Yes. Convert maximum amount in seconds to maximum amount in samples:
        maxSamples = (psych_int64) (ceil(maxSecs * ((double) audiodevices[pahandle].streaminfo->sampleRate)) * ((double) audiodevices[pahandle].inchannels));
        // Clamp insamples to that value, if neccessary:
        if (insamples > maxSamples) {
            insamples = maxSamples;
            buffersize = (size_t) insamples * sizeof(float);
        }
    }

    if (singleType & 1) {
        // Allocate output float matrix with matching number of channels and samples:
        if (c_layout)
            PsychAllocOutFloatMatArg(1, FALSE, insamples / audiodevices[pahandle].inchannels, audiodevices[pahandle].inchannels, 1, &indatafloat);
        else
            PsychAllocOutFloatMatArg(1, FALSE, audiodevices[pahandle].inchannels, insamples / audiodevices[pahandle].inchannels, 1, &indatafloat);
    }
    else {
        // Allocate output double matrix with matching number of channels and samples:
        if (c_layout)
            PsychAllocOutDoubleMatArg(1, FALSE, insamples / audiodevices[pahandle].inchannels, audiodevices[pahandle].inchannels, 1, &indata);
        else
            PsychAllocOutDoubleMatArg(1, FALSE, audiodevices[pahandle].inchannels, insamples / audiodevices[pahandle].inchannels, 1, &indata);
    }

    // Copy out absolute sample read position of first sample in buffer:
    PsychCopyOutDoubleArg(2, FALSE, (double) (audiodevices[pahandle].readposition / audiodevices[pahandle].inchannels));

    // Copy the data, convert it from float to double: Take ringbuffer wraparound into account:
    if (indatafloat) {
        // Copy to float/single matrix:
        while(buffersize > 0) {
            // Fetch next sample and copy it to matrix:
            *(indatafloat++) = (float) audiodevices[pahandle].inputbuffer[(audiodevices[pahandle].readposition % (audiodevices[pahandle].inputbuffersize / sizeof(float)))];

            // Update sample read counter:
            audiodevices[pahandle].readposition++;

            // Decrement copy counter:
            buffersize-=sizeof(float);
        }
    }
    else {
        // Copy to double matrix:
        while(buffersize > 0) {
            // Fetch next sample and copy it to matrix:
            *(indata++) = (double) audiodevices[pahandle].inputbuffer[(audiodevices[pahandle].readposition % (audiodevices[pahandle].inputbuffersize / sizeof(float)))];

            // Update sample read counter:
            audiodevices[pahandle].readposition++;

            // Decrement copy counter:
            buffersize-=sizeof(float);
        }
    }

    // Copy out overrun flag:
    PsychCopyOutDoubleArg(3, FALSE, (double) overrun);

    // Return capture timestamp in system time of first captured sample in this session. This is a bit problematic:
    // In full-duplex mode, at least OS/X doesn't return separate timestamps, so we'll provide the playback onset time
    // instead - the best we can do. In pure capture mode we get a capture timestamp...
    PsychCopyOutDoubleArg(4, FALSE, (audiodevices[pahandle].captureStartTime > 0) ? audiodevices[pahandle].captureStartTime : audiodevices[pahandle].startTime);

    // Buffer ready.
    return(PsychError_none);
}

/* PsychPortAudio('RescheduleStart') - Set new start time for an already running audio device via PortAudio.
 */
PsychError PSYCHPORTAUDIORescheduleStart(void)
{
    static char useString[] = "startTime = PsychPortAudio('RescheduleStart', pahandle, when [, waitForStart=0] [, repetitions] [, stopTime]);";
    static char synopsisString[] =
    "Modify requested start time 'when' of an already started PortAudio audio device.\n"
    "After you've started an audio device via the 'Start' subfunction, but *before* the "
    "device has really started playback (because the 'when' time provided to the 'Start' "
    "method is still far in the future), you can use this function to reschedule the start for "
    "a different 'when' time - including a value of zero for immediate start.\n"
    "\n"
    "The 'pahandle' is the handle of the device to start. Starting a "
    "device means: Start playback of output devices, start recording on capture device, do both on "
    "full duplex devices. 'waitForStart' if set to 1 will wait until device has really started, default "
    "is to continue immediately, ie. only schedule start of device. 'when' Requested time, when device "
    "should start. Defaults to zero, i.e. start immediately. If set to a non-zero system time, PTB will "
    "do its best to start the device at the requested time, but the accuracy of start depends on the "
    "operating system, audio hardware and system load. If 'waitForStart' is set to non-zero value, ie "
    "if PTB should wait for sound onset, then the optional return argument 'startTime' will contain an "
    "estimate of when the first audio sample hit the speakers, i.e., the real start time.\n"
    "Please note that the 'when' value always refers to playback, so it defines the starttime of "
    "playback. The start time of capture is related to the start time of playback in duplex mode, "
    "but it isn't the same. In pure capture mode (without playback), 'when' will be ignored and "
    "capture always starts immediately. See the help for subfunction 'GetStatus' for more info on "
    "the meaning of the different timestamps.\n"
    "The 'repetitions' parameter will change the number of playback repetitions if provided. The "
    "value for 'repetitions' from the PsychPortAudio('Start') function will be used if the parameter "
    "is omitted. See explanation in the 'Start' function for its meaning.\n"
    "'stopTime' is an optional override for the 'stopTime' parameter from the 'Start' function, see "
    "explanations there.\n";

    static char seeAlsoString[] = "Open";

    int pahandle= -1;
    int waitForStart = 0;
    double when = 0.0;
    double repetitions = -1;
    double stopTime = -1;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(5));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(2)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgRequired, &pahandle);
    if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");
    if ((audiodevices[pahandle].opmode & kPortAudioMonitoring) == 0) {
        // Not in monitoring mode: We must have in/outbuffers allocated:
        if ((audiodevices[pahandle].opmode & kPortAudioPlayBack) && (audiodevices[pahandle].outputbuffer == NULL) && (audiodevices[pahandle].schedule == NULL)) PsychErrorExitMsg(PsychError_user, "Sound outputbuffer doesn't contain any sound to play?!?");
        if ((audiodevices[pahandle].opmode & kPortAudioCapture) && (audiodevices[pahandle].inputbuffer == NULL)) PsychErrorExitMsg(PsychError_user, "Sound inputbuffer not prepared/allocated for capture?!?");
    }

    // Get new required 'when' start time:
    PsychCopyInDoubleArg(2, kPsychArgRequired, &when);
    if (when < 0) PsychErrorExitMsg(PsychError_user, "Invalid setting for 'when'. Valid values are zero or greater.");

    PsychCopyInIntegerArg(3, kPsychArgOptional, &waitForStart);
    if (waitForStart < 0) PsychErrorExitMsg(PsychError_user, "Invalid setting for 'waitForStart'. Valid values are zero or greater.");

    // Get new optional 'repetitions' count:
    if (PsychCopyInDoubleArg(4, kPsychArgOptional, &repetitions)) {
        // Argument provided: Range-Check and assign it:
        if (repetitions < 0) PsychErrorExitMsg(PsychError_user, "Invalid setting for 'repetitions'. Valid values are zero or greater.");
    }
    else {
        repetitions = -1;
    }

    // Get new optional 'stopTime':
    if (PsychCopyInDoubleArg(5, kPsychArgOptional, &stopTime)) {
        // Argument provided: Range-Check and assign it:
        if (stopTime <= when && (stopTime < DBL_MAX)) PsychErrorExitMsg(PsychError_user, "Invalid setting for 'stopTime'. Valid values are greater than 'when' starttime.");
    }
    else {
        stopTime = -1;
    }

    // Audio engine running? That is the minimum requirement for this function to work:
    if (!Pa_IsStreamActive(audiodevices[pahandle].stream)) PsychErrorExitMsg(PsychError_user, "Audio device not started. You need to call the 'Start' function first!");

    // Lock the device:
    PsychPALockDeviceMutex(&audiodevices[pahandle]);

    // Whatever the current scheduled starttime is, override it to be infinity:
    audiodevices[pahandle].reqStartTime = DBL_MAX;

    // Reset any pending requests:
    audiodevices[pahandle].reqstate = 255;

    // Engine is running. Is it in a state ready for rescheduling a start?

    // In runMode zero it must be in hot-standby as immediately after a 'Start' in order to be reschedulable:
    if ((audiodevices[pahandle].runMode == 0) && (audiodevices[pahandle].state != 1)) {
        PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);
        PsychErrorExitMsg(PsychError_user, "Audio device not started and waiting. You need to call the 'Start' function first with an infinite 'when' time or a 'when' time in the far future!");
    }

    // In runMode 1 the device itself is always running and has to be in a logically stopped/idle (=0) state or hotstandby for rescheduling.
    if ((audiodevices[pahandle].runMode == 1) && (audiodevices[pahandle].state > 1)) {
        PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);
        PsychErrorExitMsg(PsychError_user, "Audio device not idle. Make sure it is idle first, e.g., by proper use of the 'Stop' function or by checking its 'Active' state via the 'GetStatus' function!");
    }

    // Audio engine is in a proper state for rescheduling now:

    // New repetitions provided?
    if (repetitions >=0) {
        // Set number of requested repetitions: 0 means loop forever, default is 1 time.
        audiodevices[pahandle].repeatCount = (repetitions == 0) ? -1 : repetitions;
    }

    // New stopTime provided?
    if (stopTime >= 0) audiodevices[pahandle].reqStopTime = stopTime;

    // Reset statistics:
    audiodevices[pahandle].xruns = 0;
    audiodevices[pahandle].noTime = 0;
    audiodevices[pahandle].captureStartTime = 0;
    audiodevices[pahandle].startTime = 0.0;
    audiodevices[pahandle].estStopTime = 0;
    audiodevices[pahandle].currentTime = 0;
    audiodevices[pahandle].schedule_pos = 0;

    // Reset recorded samples counter:
    audiodevices[pahandle].recposition = 0;

    // Reset read samples counter: This will discard possibly not yet fetched data.
    audiodevices[pahandle].readposition = 0;

    // Reset play position:
    audiodevices[pahandle].playposition = 0;

    // Reset total count of played out samples:
    audiodevices[pahandle].totalplaycount = 0;

    // Setup new rescheduled target start time:
    audiodevices[pahandle].reqStartTime = when;

    if (audiodevices[pahandle].runMode == 1) {
        // Set the state to hot-standby to actually make this scheduling request active:
        audiodevices[pahandle].state = 1;
    }

    // Safety check for deadlock avoidance with waiting slaves:
    if ((waitForStart > 0) && (audiodevices[pahandle].opmode & kPortAudioIsSlave) &&
        (!Pa_IsStreamActive(audiodevices[pahandle].stream) || Pa_IsStreamStopped(audiodevices[pahandle].stream) ||
        audiodevices[audiodevices[pahandle].pamaster].state < 1)) {
        // We are a slave that shall wait for start, but the master audio device hasn't even
        // started its engine. This looks like a deadlock to avoid:
        PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);
        printf("PTB-ERROR: Failed to reschedule start of audio device %i. You asked me to wait for start to happen, but the\n", pahandle);
        printf("PTB-ERROR: master audio device %i isn't active yet, so i would hang forever! Aborting!!\n", audiodevices[pahandle].pamaster);
        PsychErrorExitMsg(PsychError_user, "Asked to 'waitForStart' of a slave device, but associated master device not even started! Deadlock avoided!");
    }

    if (waitForStart>0) {
        // Wait for real start of device: We enter the first while() loop iteration with
        // the device lock still held from above, so the while() loop will iterate at
        // least once...
        while (audiodevices[pahandle].state == 1 && Pa_IsStreamActive(audiodevices[pahandle].stream)) {
            // Wait for a state-change before reevaluating the .state:
            PsychPAWaitForChange(&audiodevices[pahandle]);
        }

        // Device has started (potentially even already finished for very short sounds!)
        // In any case we have a valid .startTime to wait for in the device struct,
        // which won't change by itself anymore, so it is safe to access unlocked.
        // Unlock device:
        PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);

        // Ok, relevant audio buffer with real sound onset submitted to engine.
        // We now have an estimate of real sound onset in startTime, wait until
        // then:
        PsychWaitUntilSeconds(audiodevices[pahandle].startTime);

        // Engine should run now. Return real onset time:
        PsychCopyOutDoubleArg(1, kPsychArgOptional, audiodevices[pahandle].startTime);
    }
    else {
        // Unlock device: This will trigger actual start at next paCallback() invocation by engine:
        PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);

        // Return empty zero timestamp to signal that this info is not available:
        PsychCopyOutDoubleArg(1, kPsychArgOptional, 0.0);
    }

    return(PsychError_none);
}


/* PsychPortAudio('Start') - Start an audio device via PortAudio.
 */
PsychError PSYCHPORTAUDIOStartAudioDevice(void)
{
    static char useString[] = "startTime = PsychPortAudio('Start', pahandle [, repetitions=1] [, when=0] [, waitForStart=0] [, stopTime=inf] [, resume=0]);";
    static char synopsisString[] =
    "Start a PortAudio audio device. The 'pahandle' is the handle of the device to start. Starting a "
    "device means: Start playback of output devices, start recording on capture device, do both on "
    "full duplex devices. 'waitForStart' if set to 1 will wait until device has really started, default "
    "is to continue immediately, ie. only schedule start of device. 'when' Requested time, when device "
    "should start. Defaults to zero, i.e. start immediately. If set to a non-zero system time, PTB will "
    "do its best to start the device at the requested time, but the accuracy of start depends on the "
    "operating system, audio hardware and system load. If 'waitForStart' is set to non-zero value, ie "
    "if PTB should wait for sound onset, then the optional return argument 'startTime' will contain an "
    "estimate of when the first audio sample hit the speakers, i.e., the real start time.\n"
    "Please note that the 'when' value always refers to playback, so it defines the starttime of "
    "playback. The start time of capture is related to the start time of playback in duplex mode, "
    "but it isn't the same. In pure capture mode (without playback), 'when' will be ignored and "
    "capture always starts immediately. See the help for subfunction 'GetStatus' for more info on "
    "the meaning of the different timestamps.\n"
    "The 'repetitions' parameter defines how often the playback of the sound data should be repeated. "
    "A setting of zero will cause infinite repetitions, ie., until manually stopped via the 'Stop' "
    "subfunction. A positive setting will cause the provided number of repetitions to happen. The "
    "default setting is 1, ie., play exactly once, then stop. Fractional values are allowed, e.g, "
    "1.5 for one and a half repetition.\n"
    "The optional parameter 'stopTime' allows to set a specific system time when sound playback "
    "should stop by itself at latest, regardless if the requested number of 'repetitions' has "
    "completed. PTB will do its best to stop sound at exactly that time, see comments about the "
    "'when' parameter - The same mechanism is used, with the same restrictions.\n"
    "The optional parameter 'resume' if set to 1, allows to resume playback at the position it "
    "was last stopped, instead of starting at the beginning again. By default, playback starts "
    "at the beginning.\n";

    static char seeAlsoString[] = "Open";

    PaError err;
    int pahandle= -1;
    int waitForStart = 0;
    int resume = 0;
    double repetitions = 1;
    double when = 0.0;
    double stopTime = DBL_MAX;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(6));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgRequired, &pahandle);
    if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");

    // If this is a master device, we refault repetitions to zero, ie., infinite playback/capture.
    // Makes much more sense for a master device:
    if (audiodevices[pahandle].opmode & kPortAudioIsMaster) repetitions = 0;

    PsychCopyInDoubleArg(2, kPsychArgOptional, &repetitions);
    if (repetitions < 0) PsychErrorExitMsg(PsychError_user, "Invalid setting for 'repetitions'. Valid values are zero or greater.");

    PsychCopyInDoubleArg(3, kPsychArgOptional, &when);
    if (when < 0) PsychErrorExitMsg(PsychError_user, "Invalid setting for 'when'. Valid values are zero or greater.");

    PsychCopyInIntegerArg(4, kPsychArgOptional, &waitForStart);
    if (waitForStart < 0) PsychErrorExitMsg(PsychError_user, "Invalid setting for 'waitForStart'. Valid values are zero or greater.");

    PsychCopyInDoubleArg(5, kPsychArgOptional, &stopTime);
    if (stopTime <= when && (stopTime < DBL_MAX)) PsychErrorExitMsg(PsychError_user, "Invalid setting for 'stopTime'. Valid values are greater than 'when' starttime.");

    PsychCopyInIntegerArg(6, kPsychArgOptional, &resume);
    if (resume < 0) PsychErrorExitMsg(PsychError_user, "Invalid setting for 'resume'. Valid values are zero or greater.");

    if ((audiodevices[pahandle].opmode & kPortAudioMonitoring) == 0) {
        // Not in monitoring mode: We must have in/outbuffers allocated:
        if ((audiodevices[pahandle].opmode & kPortAudioPlayBack) && (audiodevices[pahandle].outputbuffer == NULL) && (audiodevices[pahandle].schedule == NULL)) PsychErrorExitMsg(PsychError_user, "Sound outputbuffer doesn't contain any sound to play?!?");
        if ((audiodevices[pahandle].opmode & kPortAudioCapture) && (audiodevices[pahandle].inputbuffer == NULL)) PsychErrorExitMsg(PsychError_user, "Sound inputbuffer not prepared/allocated for capture?!?");
    }

    // Make sure current state is zero, aka fully stopped and engine is really stopped: Output a warning if this looks like an
    // unintended "too early" restart: [No need to mutex-lock here, as iff these .state setting is not met,
    // then we are good and they can't change by themselves behind our back -- paCallback() can't change .state to > 0]
    if ((audiodevices[pahandle].state > 0) && Pa_IsStreamActive(audiodevices[pahandle].stream)) {
        if (verbosity > 1) {
            printf("PsychPortAudio-WARNING: 'Start' method on audiodevice %i called, although playback on device not yet completely stopped.\nWill forcefully restart with possible audible artifacts or timing glitches.\nCheck your playback timing or use the 'Stop' function properly!\n", pahandle);
        }
    }

    // Safeguard: If the stream is not stopped in runMode 0, do it now:
    if (!Pa_IsStreamStopped(audiodevices[pahandle].stream)) {
        if (audiodevices[pahandle].runMode == 0) Pa_StopStream(audiodevices[pahandle].stream);
    }

    // Mutex-lock here: Needed if engine already/still running in runMode1, doesn't hurt if engine is stopped
    PsychPALockDeviceMutex(&audiodevices[pahandle]);

    // Reset statistics values:
    audiodevices[pahandle].batchsize = 0;
    audiodevices[pahandle].xruns = 0;
    audiodevices[pahandle].paCalls = 0;
    audiodevices[pahandle].noTime = 0;
    audiodevices[pahandle].captureStartTime = 0;
    audiodevices[pahandle].startTime = 0.0;
    audiodevices[pahandle].reqStopTime = stopTime;
    audiodevices[pahandle].estStopTime = 0;
    audiodevices[pahandle].currentTime = 0;
    if (!resume) audiodevices[pahandle].schedule_pos = 0;

    // Reset recorded samples counter:
    audiodevices[pahandle].recposition = 0;

    // Reset read samples counter: This will discard possibly not yet fetched data.
    audiodevices[pahandle].readposition = 0;

    // Reset play position:
    if (!resume) audiodevices[pahandle].playposition = 0;

    // Reset total count of played out samples:
    if (!resume) audiodevices[pahandle].totalplaycount = 0;

    // Set number of requested repetitions: 0 means loop forever, default is 1 time.
    audiodevices[pahandle].repeatCount = (repetitions == 0) ? -1 : repetitions;

    // Reset any pending requests:
    audiodevices[pahandle].reqstate = 255;

    // Setup target start time:
    audiodevices[pahandle].reqStartTime = when;

    // Mark state as "hot-started":
    audiodevices[pahandle].state = 1;

    if (!(audiodevices[pahandle].opmode & kPortAudioIsSlave)) {
        // Engine running?
        if (!Pa_IsStreamActive(audiodevices[pahandle].stream) || Pa_IsStreamStopped(audiodevices[pahandle].stream)) {
            // Try to start stream if the engine isn't running, either because it is the very
            // first call to 'Start' in any runMode, or because the engine got stopped in
            // preparation for a restart in runMode zero. Need to drop the lock during
            // call to Pa_StartStream to avoid deadlock:
            PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);

            // Safeguard: If the stream is not stopped, do it now:
            if (!Pa_IsStreamStopped(audiodevices[pahandle].stream)) Pa_StopStream(audiodevices[pahandle].stream);

            // Start engine:
            if ((err=Pa_StartStream(audiodevices[pahandle].stream))!=paNoError) {
                printf("PTB-ERROR: Failed to start audio device %i. PortAudio reports this error: %s \n", pahandle, Pa_GetErrorText(err));
                PsychErrorExitMsg(PsychError_system, "Failed to start PortAudio audio device.");
            }

            // Reacquire lock:
            PsychPALockDeviceMutex(&audiodevices[pahandle]);
        }
    }
    // From here on, the engine is running, and we have the mutex-lock. Unless we're a slave,
    // then we have the lock but the engine may not be running yet.

    // Safety check for deadlock avoidance with waiting slaves:
    if ((waitForStart > 0) && (audiodevices[pahandle].opmode & kPortAudioIsSlave) &&
        (!Pa_IsStreamActive(audiodevices[pahandle].stream) || Pa_IsStreamStopped(audiodevices[pahandle].stream) ||
        audiodevices[audiodevices[pahandle].pamaster].state < 1)) {
        // We are a slave that shall wait for start, but the master audio device hasn't even
        // started its engine. This looks like a deadlock to avoid:
        PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);
        printf("PTB-ERROR: Failed to start audio device %i. You asked me to wait for start to happen, but the\n", pahandle);
        printf("PTB-ERROR: master audio device %i isn't active yet, so i would hang forever! Aborting!!\n", audiodevices[pahandle].pamaster);
        PsychErrorExitMsg(PsychError_user, "Asked to 'waitForStart' of a slave device, but associated master device not even started! Deadlock avoided!");
    }

    // Wait for real start of playback/capture?
    if (waitForStart > 0) {
        // Device will be in state == 1 until playback really starts:
        // We need to enter the first while() loop iteration with
        // the device lock held from above, so the while() loop will iterate at
        // least once...
        while (audiodevices[pahandle].state == 1 && Pa_IsStreamActive(audiodevices[pahandle].stream)) {
            // Wait for a state-change before reevaluating the .state:
            PsychPAWaitForChange(&audiodevices[pahandle]);
        }

        // Device has started (potentially even already finished for very short sounds!)
        // In any case we have a valid .startTime to wait for in the device struct,
        // which won't change by itself anymore, so it is safe to access unlocked.
        // Unlock device:
        PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);

        // Ok, relevant audio buffer with real sound onset submit to engine.
        // We now have an estimate of real sound onset in startTime, wait until
        // then:
        PsychWaitUntilSeconds(audiodevices[pahandle].startTime);

        // Engine should run now. Return real onset time:
        PsychCopyOutDoubleArg(1, kPsychArgOptional, audiodevices[pahandle].startTime);
    }
    else {
        // Unlock device:
        PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);

        // Return empty zero timestamp to signal that this info is not available:
        PsychCopyOutDoubleArg(1, kPsychArgOptional, 0.0);
    }

    return(PsychError_none);
}

/* PsychPortAudio('Stop') - Stop an audio device via PortAudio.
 */
PsychError PSYCHPORTAUDIOStopAudioDevice(void)
{
    static char useString[] = "[startTime endPositionSecs xruns estStopTime] = PsychPortAudio('Stop', pahandle [, waitForEndOfPlayback=0] [, blockUntilStopped=1] [, repetitions] [, stopTime]);";
    static char synopsisString[] =
    "Stop a PortAudio audio device. The 'pahandle' is the handle of the device to stop.\n"
    "'waitForEndOfPlayback' - If set to 1, this method will wait until playback of the "
    "audio stream finishes by itself. This only makes sense if you perform playback with "
    "a defined playback duration. The flag will be ignored when infinite repetition is "
    "requested (as playback would never stop by itself, resulting in a hang) and if no "
    "scheduled 'stopTime' has been set, either in this call or in the 'Start' function. "
    "The setting will be also ignored if this is a pure recording session.\n"
    "A setting of 0 (which is the default) requests stop of playback without waiting for all "
    "the sound to be finished playing. Sound may continue to be played back for multiple "
    "milliseconds after this call, as this is a polite request to the hardware to finish up.\n"
    "A setting of 2 requests abortion of playback and/or capture as soon as possible with your "
    "sound hardware, even if this creates audible artifacts etc. Abortion may or may not be faster "
    "than a normal stop, this depends on your specific hardware, but our driver tries as hard as "
    "possible to get the hardware to shut up. In a worst-case setting, the hardware would continue "
    "to playback the sound data that is stored in its internal buffers, so the latency for stopping "
    "sound would be the the same as the latency for starting sound as quickly as possible. E.g., "
    "a 2nd generation Intel MacBook Pro seems to have a stop-delay of roughly 5-6 msecs under "
    "optimal conditions (e.g., buffersize = 64, frequency=96000, OS/X 10.4.10, no other sound apps running).\n"
    "A setting of 3 will not try to stop the playback now: You can use this setting if you just "
    "want to use this function to wait until playback stops by itself, e.g. because the set "
    "number of 'repetitions' or the set 'stopTime' has been reached. You can also use this to "
    "just change the 'repetitions' or 'stopTime' settings without waiting for anything.\n"
    "The optional parameter 'blockUntilStopped' defines if the subfunction should wait until "
    "sound processing has really stopped (at a setting of 1, which is the default), or if the "
    "function should return with minimal delay after only scheduling a stop at a zero setting. "
    "If 'waitForEndOfPlayback' is set to 1, then 'blockUntilStopped' is meaningless and the function "
    "will always block until the stop is completed.\n"
    "The optional parameter 'stopTime' allows to set a defined system time when playback should "
    "stop by itself. Similar, the optional 'repetitions' setting allows to change the number of "
    "repetitions after which playback should stop by itself. These settings allow you to override "
    "the same settings made during a call to the 'Start' or 'RescheduleStart' function.\n"
    "The optional return argument 'startTime' returns an estimate of when the stopped "
    "stream actually started its playback and/or recording. Its the same timestamp as the one "
    "returned by the start command when executed in waiting mode. 'endPositionSecs' is the final "
    "playback/recording position in seconds. 'xruns' is the number of buffer over- or underruns. "
    "This should be zero if the playback operation was glitch-free, however a zero value doesn't "
    "imply glitch free operation, as the glitch detection algorithm can miss some types of glitches. "
    "The optional return argument 'estStopTime' returns an estimate of when playback likely stopped.\n"
    "The return arguments are undefined if you set the 'blockUntilStopped' flag to zero.\n";

    static char seeAlsoString[] = "Open GetDeviceSettings ";

    PaError err;
    int pahandle= -1;
    int waitforend = 0;
    int blockUntilStopped = 1;
    double stopTime = -1;
    double repetitions = -1;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(5));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(4));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgRequired, &pahandle);
    if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");

    // Get optional wait-flag:
    PsychCopyInIntegerArg(2, kPsychArgOptional, &waitforend);

    // Get optional blockUntilStopped-flag:
    PsychCopyInIntegerArg(3, kPsychArgOptional, &blockUntilStopped);

    // Get new optional 'repetitions' count:
    if (PsychCopyInDoubleArg(4, kPsychArgOptional, &repetitions)) {
        // Argument provided: Range-Check and assign it:
        if (repetitions < 0) PsychErrorExitMsg(PsychError_user, "Invalid setting for 'repetitions'. Valid values are zero or greater.");
    }
    else {
        repetitions = -1;
    }

    // Get optional stopTime:
    if (PsychCopyInDoubleArg(5, kPsychArgOptional, &stopTime)) {
        if ((stopTime <= audiodevices[pahandle].reqStartTime) && (stopTime < DBL_MAX)) PsychErrorExitMsg(PsychError_user, "Invalid setting for 'stopTime'. Valid values are greater than previously set 'when' starttime.");
    }
    else {
        stopTime = -1;
    }

    // Lock device:
    PsychPALockDeviceMutex(&audiodevices[pahandle]);

    // New repetitions provided?
    if (repetitions >=0) {
        // Set number of requested repetitions: 0 means loop forever, default is 1 time.
        audiodevices[pahandle].repeatCount = (repetitions == 0) ? -1 : repetitions;
    }

    // New stopTime provided?
    if (stopTime > 0) {
        // Yes. Quickly assign it:
        audiodevices[pahandle].reqStopTime = stopTime;
    }

    // Wait for automatic stop of playback if requested: This only makes sense if we
    // are in playback mode and not in infinite playback mode! Would not make sense in
    // looping mode (infinite repetitions of playback) or pure recording mode. It is also
    // allowed if we have infinite repetitions set, but a finite stopTime is defined, so
    // the engine will eventually stop by itself. Same goes for an operative schedule which
    // will run empty if not regularly updated:
    if ((waitforend == 1) && Pa_IsStreamActive(audiodevices[pahandle].stream) && (audiodevices[pahandle].state > 0) &&
        (audiodevices[pahandle].opmode & kPortAudioPlayBack) && ((audiodevices[pahandle].repeatCount != -1) || (audiodevices[pahandle].schedule) || (audiodevices[pahandle].reqStopTime < DBL_MAX))) {
        while ( ((audiodevices[pahandle].runMode == 0) && Pa_IsStreamActive(audiodevices[pahandle].stream) && (audiodevices[pahandle].state > 0)) ||
            ((audiodevices[pahandle].runMode == 1) && (audiodevices[pahandle].state > 0))) {

            // Wait for a state-change before reevaluating:
            PsychPAWaitForChange(&audiodevices[pahandle]);
        }
    }

        // Lock held here in any case...

    if (waitforend == 3) {
        // No immediate stop request: This was only either a query for end of playback,
        // or a call to simply set new 'stopTime' or 'repetitions' parameters on the fly.
        // Unlock the device and skip stop requests...
        PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);
    }
    else {
        // Some real immediate stop request wanted:
        // Soft stop requested (as opposed to fast stop)?
        if (waitforend!=2) {
            // Softstop: Try to stop stream:
            if (audiodevices[pahandle].state > 0) {
                // Stream running. Request a stop of stream, to be honored by playback thread:
                audiodevices[pahandle].reqstate = 0;
            }

            // Drop lock, so request can get through...
            PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);

            // If blockUntilStopped is non-zero, then explicitely stop as well:
            if ((blockUntilStopped > 0) && (audiodevices[pahandle].runMode == 0) && (!Pa_IsStreamStopped(audiodevices[pahandle].stream)) && (err=Pa_StopStream(audiodevices[pahandle].stream))!=paNoError) {
                printf("PTB-ERROR: Failed to stop audio device %i. PortAudio reports this error: %s \n", pahandle, Pa_GetErrorText(err));
                PsychErrorExitMsg(PsychError_system, "Failed to stop PortAudio audio device.");
            }
        }
        else {
            // Faststop: Try to abort stream. Skip if already stopped/not yet started:

            // Stream active?
            if (audiodevices[pahandle].state > 0) {
                // Yes. Set the 'state' flag to signal our IO-Thread not to push any audio
                // data anymore, but only zeros for silence and to paAbort asap:
                audiodevices[pahandle].reqstate = 3;
            }

            // Drop lock, so request can get through...
            PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);

            // If blockUntilStopped is non-zero, then send abort request to hardware:
            if ((blockUntilStopped > 0) && (audiodevices[pahandle].runMode == 0) && (!Pa_IsStreamStopped(audiodevices[pahandle].stream)) && ((err=Pa_AbortStream(audiodevices[pahandle].stream))!=paNoError)) {
                printf("PTB-ERROR: Failed to abort audio device %i. PortAudio reports this error: %s \n", pahandle, Pa_GetErrorText(err));
                PsychErrorExitMsg(PsychError_system, "Failed to fast stop (abort) PortAudio audio device.");
            }
        }
    }

    // No lock held here...

    // Wait for real stop:
    if (blockUntilStopped > 0) {
        // Lock device:
        PsychPALockDeviceMutex(&audiodevices[pahandle]);

        // Wait for stop / idle:
        if (Pa_IsStreamActive(audiodevices[pahandle].stream)) {
            while ( ((audiodevices[pahandle].runMode == 0) && Pa_IsStreamActive(audiodevices[pahandle].stream) && (audiodevices[pahandle].state > 0)) ||
                ((audiodevices[pahandle].runMode == 1) && (audiodevices[pahandle].state > 0))) {

                // Wait for a state-change before reevaluating:
                PsychPAWaitForChange(&audiodevices[pahandle]);
            }
        }

        // We are stopped/idle, with lock held:

        // Need to update stream state and reqstate manually here, as the Pa_Stop/AbortStream()
        // requests may have stopped the paCallback() thread before it could update/honor state/reqstate by himself.

        // Mark state as stopped:
        audiodevices[pahandle].state = 0;

        // Reset request to none:
        audiodevices[pahandle].reqstate = 255;

        // Can unlock here, as the fields we're interested in will remain static with an idle/stopped engine:
        PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);

        // Copy out our estimate of when playback really started for the just stopped stream:
        PsychCopyOutDoubleArg(1, kPsychArgOptional, audiodevices[pahandle].startTime);

        // Copy out final playback position (secs) since start:
        PsychCopyOutDoubleArg(2, kPsychArgOptional, ((double)(audiodevices[pahandle].playposition / audiodevices[pahandle].outchannels)) / (double) audiodevices[pahandle].streaminfo->sampleRate);

        // Copy out number of buffer over-/underflows since start:
        PsychCopyOutDoubleArg(3, kPsychArgOptional, audiodevices[pahandle].xruns);

        // Copy out estimated stopTime:
        PsychCopyOutDoubleArg(4, kPsychArgOptional, audiodevices[pahandle].estStopTime);
    }
    else {
        // No block until stopped. That means we won't have meaningful return arguments available.
        // Just return dummy args to signal this:
        PsychCopyOutDoubleArg(1, kPsychArgOptional, -1);

        // Copy out final playback position (secs) since start:
        PsychCopyOutDoubleArg(2, kPsychArgOptional, -1);

        // Copy out number of buffer over-/underflows since start:
        PsychCopyOutDoubleArg(3, kPsychArgOptional, -1);

        // Copy out estimated stopTime:
        PsychCopyOutDoubleArg(4, kPsychArgOptional, -1);
    }

    return(PsychError_none);
}

/* PsychPortAudio('GetStatus') - Return current status of stream.
 */
PsychError PSYCHPORTAUDIOGetStatus(void)
{
    static char useString[] = "status = PsychPortAudio('GetStatus', pahandle);";
    static char synopsisString[] =
    "Returns 'status', a struct with status information about the current state of device 'pahandle'.\n"
    "The struct contains the following fields:\n"
    "Active: Can be 1 if playback or recording is active, or 0 if playback/recording is stopped or not yet started.\n"
    "RequestedStartTime: Is the requested start time of the audio stream after start of playback/recording. "
    "StartTime: Is the real start time of the audio stream after start of playback/recording. If both, playback and "
    "recording are active (full-duplex mode), it is the start time of sound playback, ie an estimate of when the first "
    "sample hit the speakers. Same goes for pure playback. \n"
    "CaptureStartTime: Start time of audio capture (if any is active) - an estimate of when the first sound sample was "
    "captured. In pure capture mode, this is nearly identical to StartTime, but whenever playback is active, StartTime and "
    "CaptureStartTime will differ. CaptureStartTime doesn't take the user provided 'LatencyBias' into account.\n"
    "RequestedStopTime: The requested stop / sound offset time for playback of sounds, as selected via the 'stopTime' "
    "paramter in the 'Start', 'RescheduleStart' or 'Stop' function. This will show a very large (~ infinite) value "
    "if no stop time has been sprecified.\n"
    "EstimatedStopTime: Estimated time when sound playback has stopped. This is an estimate of when exactly the last "
    "audio sample will leave the speaker. The value is zero as long as the estimate isn't available. Due to the latency "
    "involved in sound playback, the value may become available a few msecs before or after actual sound offset.\n"
    "CurrentStreamTime: Estimate of when the most recently submitted sample will hit the speaker. This corresponds "
    "roughly to 'PositionSecs' below, but in absolute realtime.\n"
    "ElapsedOutSamples: Total number of samples played out since start of playback. This count increments monotonically "
    "from start of playback to stop of playback. This denotes the absolute sample position that will hit the speaker "
    "at time 'CurrentStreamTime'. \n"
    "PositionSecs is an estimate of the current stream playback position in seconds within the current playback loop "
    "of the current buffer. it's not totally accurate, because "
    "it measures how much sound has been submitted to the sound system, not how much sound has left the "
    "speakers, i.e., it doesn't take driver and hardware latency into account.\n"
    "SchedulePosition: Current position in a running schedule, if any.\n"
    "XRuns: Number of dropouts due to buffer overrun or underrun conditions. This is not perfectly reliable, "
    "as the algorithm can miss some dropouts. Iow.: A non-zero or increasing value means that audio glitches "
    "during playback or capture happened, but a zero or constant value doesn't mean everything was glitch-free, "
    "because some glitches can't get reliably detected on some operating systems or audio hardware.\n"
    "TotalCalls, TimeFailed and BufferSize are only for debugging of PsychPortAudio itself.\n"
    "CPULoad: How much load does the playback engine impose on the CPU? Values can range from 0.0 = 0% "
    "to 1.0 for 100%. Values close to 1.0 indicate that your system can't handle the load and timing glitches "
    "or sound glitches are likely. In such a case, try to reduce the load on your system.\n"
    "PredictedLatency: Is the latency in seconds of your driver+hardware combo. It tells you, "
    "how far ahead of time a sound device must be started ahead of the requested onset time via "
    "PsychPortAudio('Start'...) to make sure it actually starts playing in time. High quality systems like "
    "Linux or MacOS/X may allow values as low as 5 msecs or less on standard hardware. Other operating "
    "systems may require dozens or hundreds of milliseconds of headstart. Caution: In full-duplex mode, "
    "this value only refers to the latency on the sound output, not in the sound input! Also, this is just "
    "an estimate, not 100% reliable.\n"
    "LatencyBias: Is an additional bias setting you can impose via PsychPortAudio('LatencyBias', pahandle, bias); "
    "in case our drivers estimate is a bit off. Allows fine-tuning.\n"
    "SampleRate: Is the sampling rate for playback/recording in samples per second (Hz).\n"
    "OutDeviceIndex: Is the deviceindex of the playback device, or -1 if not opened for playback. "
    "You can pass OutDeviceIndex to PsychPortAudio('GetDevices', [], OutDeviceIndex); to query information "
    "about the device.\n"
    "InDeviceIndex: Is the deviceindex of the capture device, or -1 if not opened for capture.\n"
    "RecordedSecs: Is the total amount of recorded sound data (in seconds) since start of capture.\n"
    "ReadSecs: Is the total amount of sound data (in seconds) that has been fetched from the internal buffer. "
    "The difference between RecordedSecs and ReadSecs is the amount of recorded sound data pending for retrieval. ";

    static char seeAlsoString[] = "Open GetDeviceSettings ";
    PsychGenericScriptType     *status;
    double currentTime;
    psych_int64 playposition, totalplaycount, recposition;

    const char *FieldNames[]={    "Active", "State", "RequestedStartTime", "StartTime", "CaptureStartTime", "RequestedStopTime", "EstimatedStopTime", "CurrentStreamTime", "ElapsedOutSamples", "PositionSecs", "RecordedSecs", "ReadSecs", "SchedulePosition",
        "XRuns", "TotalCalls", "TimeFailed", "BufferSize", "CPULoad", "PredictedLatency", "LatencyBias", "SampleRate",
        "OutDeviceIndex", "InDeviceIndex" };
    int pahandle = -1;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(1));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgRequired, &pahandle);
    if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");

    PsychAllocOutStructArray(1, kPsychArgOptional, -1, 23, FieldNames, &status);

    // Ok, in a perfect world we should hold the device mutex while querying all the device state.
    // However, we don't: This reduces lock contention at the price of a small chance that the
    // fetched information is not 100% up to date / that this is not an atomic snapshot of state.
    //
    // Instead we only hold the lock to get the most crucial values atomically, then release and get the rest
    // while not holding the lock:
    PsychPALockDeviceMutex(&audiodevices[pahandle]);
    currentTime = audiodevices[pahandle].currentTime;
    totalplaycount = audiodevices[pahandle].totalplaycount;
    playposition = audiodevices[pahandle].playposition;
    recposition = audiodevices[pahandle].recposition;
    PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);

    // Atomic snapshot for remaining fields would only be needed for low-level debugging, so who cares?
    PsychSetStructArrayDoubleElement("Active", 0, (audiodevices[pahandle].state >= 2) ? 1 : 0, status);
    PsychSetStructArrayDoubleElement("State", 0, audiodevices[pahandle].state, status);
    PsychSetStructArrayDoubleElement("RequestedStartTime", 0, audiodevices[pahandle].reqStartTime, status);
    PsychSetStructArrayDoubleElement("StartTime", 0, audiodevices[pahandle].startTime, status);
    PsychSetStructArrayDoubleElement("CaptureStartTime", 0, audiodevices[pahandle].captureStartTime, status);
    PsychSetStructArrayDoubleElement("RequestedStopTime", 0, audiodevices[pahandle].reqStopTime, status);
    PsychSetStructArrayDoubleElement("EstimatedStopTime", 0, audiodevices[pahandle].estStopTime, status);
    PsychSetStructArrayDoubleElement("CurrentStreamTime", 0, currentTime, status);
    PsychSetStructArrayDoubleElement("ElapsedOutSamples", 0, ((double)(totalplaycount / audiodevices[pahandle].outchannels)), status);
    PsychSetStructArrayDoubleElement("PositionSecs", 0, ((double)(playposition / audiodevices[pahandle].outchannels)) / (double) audiodevices[pahandle].streaminfo->sampleRate, status);
    PsychSetStructArrayDoubleElement("RecordedSecs", 0, ((double)(recposition / audiodevices[pahandle].inchannels)) / (double) audiodevices[pahandle].streaminfo->sampleRate, status);
    PsychSetStructArrayDoubleElement("ReadSecs", 0, ((double)(audiodevices[pahandle].readposition / audiodevices[pahandle].inchannels)) / (double) audiodevices[pahandle].streaminfo->sampleRate, status);
    PsychSetStructArrayDoubleElement("SchedulePosition", 0, audiodevices[pahandle].schedule_pos, status);
    PsychSetStructArrayDoubleElement("XRuns", 0, audiodevices[pahandle].xruns, status);
    PsychSetStructArrayDoubleElement("TotalCalls", 0, audiodevices[pahandle].paCalls, status);
    PsychSetStructArrayDoubleElement("TimeFailed", 0, audiodevices[pahandle].noTime, status);
    PsychSetStructArrayDoubleElement("BufferSize", 0, (double) audiodevices[pahandle].batchsize, status);
    PsychSetStructArrayDoubleElement("CPULoad", 0, (Pa_IsStreamActive(audiodevices[pahandle].stream)) ? Pa_GetStreamCpuLoad(audiodevices[pahandle].stream) : 0.0, status);
    PsychSetStructArrayDoubleElement("PredictedLatency", 0, audiodevices[pahandle].predictedLatency, status);
    PsychSetStructArrayDoubleElement("LatencyBias", 0, audiodevices[pahandle].latencyBias, status);
    PsychSetStructArrayDoubleElement("SampleRate", 0, audiodevices[pahandle].streaminfo->sampleRate, status);
    PsychSetStructArrayDoubleElement("OutDeviceIndex", 0, audiodevices[pahandle].outdeviceidx, status);
    PsychSetStructArrayDoubleElement("InDeviceIndex", 0, audiodevices[pahandle].indeviceidx, status);
    return(PsychError_none);
}

/* PsychPortAudio('Verbosity') - Set level of verbosity.
 */
PsychError PSYCHPORTAUDIOVerbosity(void)
{
    static char useString[] = "oldlevel = PsychPortAudio('Verbosity' [,level]);";
    static char synopsisString[] =
    "Set level of verbosity for error/warning/status messages. 'level' optional, new level "
    "of verbosity. 'oldlevel' is the old level of verbosity. The following levels are "
    "supported: 0 = Shut up. 1 = Print errors, 2 = Print also warnings, 3 = Print also some info, "
    "4 = Print more useful info (default), >5 = Be very verbose (mostly for debugging the driver itself). ";
    static char seeAlsoString[] = "Open GetDeviceSettings ";

    int level= -1;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(1));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(0)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    PsychCopyInIntegerArg(1, kPsychArgOptional, &level);
    if (level < -1) PsychErrorExitMsg(PsychError_user, "Invalid level of verbosity provided. Valid are levels of zero and greater.");

    // Return current/old level:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) verbosity);

    // Set new level, if one was provided:
    if (level > -1) {
        verbosity = level;

        #if PSYCH_SYSTEM == PSYCH_LINUX
            // Set an error handler for ALSA debug output/errors to stop the spewage of utterly
            // pointless ALSA warning messages to stderr. At verbosity <= 5 we sent ALSA chatter
            // to a dummy error handler. At levels > 5 we disable our error handler, so ALSA
            // chatter goes to stderr...
            if (verbosity <= 5)
                snd_lib_error_set_handler(ALSAErrorHandler);
            else
                snd_lib_error_set_handler(NULL);
        #endif
    }

    return(PsychError_none);
}

/* PsychPortAudio('GetOpenDeviceCount') - Get number of open audio devices.
 */
PsychError PSYCHPORTAUDIOGetOpenDeviceCount(void)
{
    static char useString[] = "count = PsychPortAudio('GetOpenDeviceCount');";
    static char synopsisString[] = "Return the number of currently open audio devices.\n";
    static char seeAlsoString[] = "Open ";

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(0));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(0)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    // Return count:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) audiodevicecount);

    return(PsychError_none);
}

/* PsychPortAudio('LatencyBias') - Set a manual bias for the latencies we operate on.
 */
PsychError PSYCHPORTAUDIOLatencyBias(void)
{
    static char useString[] = "oldbias = PsychPortAudio('LatencyBias', pahandle [,biasSecs]);";
    static char synopsisString[] =
    "Set audio output latency bias in seconds to 'biasSecs' and/or return old bias for a device "
    "'pahandle'. The device must be open for this setting to take effect. It is reset to zero at "
    "each reopening of the device. PsychPortAudio computes a latency value for the expected latency "
    "of an audio output device to get its timing right. If this latency value is slightly off for "
    "some reason, you can provide a bias value with this function to correct the computed value.\n\n"
    "Please note that this 'biasSecs' setting is applied to either audio playback or audio capture "
    "if your device is configured in simplex mode, either for playback *or* capture. If your device "
    "is opened in full-duplex mode for simultaneous capture and playback, then the bias value is only "
    "applied to the playback timing and timestamps, but not to the capture timing and timestamps!\n\n"
    "See the online help for PsychPortAudio for more in depth explanation of latencies. ";

    static char seeAlsoString[] = "Open GetDeviceSettings ";

    double bias= DBL_MAX;
    int pahandle = -1;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(2));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgRequired, &pahandle);
    if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");

    // Copy in optional new bias value:
    PsychCopyInDoubleArg(2, kPsychArgOptional, &bias);

    // Return current/old bias:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, audiodevices[pahandle].latencyBias);

    // Set new bias, if one was provided:
    if (bias!=DBL_MAX) {
        if (audiodevices[pahandle].opmode & kPortAudioIsSlave) PsychErrorExitMsg(PsychError_user, "Change of latency bias is not allowed on slave devices! Set it on associated master device.");
        if (Pa_IsStreamActive(audiodevices[pahandle].stream) && (audiodevices[pahandle].state > 0)) PsychErrorExitMsg(PsychError_user, "Tried to change 'biasSecs' while device is active! Forbidden!");
        audiodevices[pahandle].latencyBias = bias;
    }

    return(PsychError_none);
}

/* PsychPortAudio('Volume') - Set volume per device.
 */
PsychError PSYCHPORTAUDIOVolume(void)
{
    static char useString[] = "[oldMasterVolume, oldChannelVolumes] = PsychPortAudio('Volume', pahandle [, masterVolume][, channelVolumes]);";
    static char synopsisString[] =
    "Set audio output volume(s) and/or return old volumes for a device 'pahandle'.\n"
    "The device must be open for this setting to take effect. It is initially set to 1.0. "
    "On a regular audio device or master device, you can only set the global 'masterVolume' for the "
    "whole device. On a slave device, you can optionally define a per-channel volume vector 'channelVolumes'. "
    "E.g., if your slave device has two output channels, a 'channelVolumes' "
    "vector of [0.9 ; 0.25] would set the first channel to 90% intensity, the 2nd channel to "
    "25% intensity. The elements of 'channelVolumes' and the 'masterVolume' can be any value. The value for a device or a "
    "channel is simply multiplied with each output sample before output to the hardware, so "
    "that a value of 1.0 (the default setting for each device and channel) will pass samples "
    "unmodified, a value of 0.5 would reduce intensity to 50%, a value of 1.25 would amplify "
    "to 125% of the original value, or a value of -1 would invert the signal.\n"
    "Be careful with amplification (values > 1). If the final sample leaves the valid output "
    "range of -1.0 to 1.0, you may hear ugly auditory clipping artifacts and possibly damage "
    "your ears or speakers! In 'normal' mode, PsychPortAudio clamps output samples to the valid "
    "range of -1 to 1, but in high timing precision/low latency mode, clamping is not performed "
    "in order to save computation time.\n\n"
    "Caution: All volume settings are silently ignored on regular audio devices or master devices "
    "if audio monitoring mode is active. On a slave device, per-channel volumes would be still applied "
    "in monitoring mode, but not the masterVolume. This is for efficiency reasons. If you need volume "
    "settings in monitoring mode, then open a master device, attach a slave device to it in monitoring "
    "mode, then set the 'channelVolumes' of that slave device to control volume of monitoring.\n\n"
    "Caution: The volume settings only affect what is sent to your computers sound system. The "
    "operating system or audio hardware itself may apply additional volume settings, gains, filters "
    "etc., depending on your specific setup. These are only controllable via external system specific tools, "
    "PsychPortAudio doesn't know about such additional settings and can't control them in any way!\n\n"
    "If you need to modulate output volume over time, you can also attach an additional slave "
    "device whose opMode is set to act as an AMmodulator to a master or slave device, see help for 'OpenSlave'. "
    "Such a slave device will not output sound itself, but use the data stored in its playback "
    "buffers to modulate the amplitude of its parent slave, or of the summed signal of all previously "
    "attached slaves for a master over time. Example: You create a master device, then 'OpenSlave' three regular slave "
    "devices, then 'OpenSlave' a AMmodulator slave device. During playback, the sound signals "
    "of the three regular slaves will be mixed together. The combined signal will be multiplied "
    "by the per-sample volume values provided by the AMmodulator slave device, thereby modulating "
    "the amplitude or ''acoustic envelope'' of the mixed signals. The resulting signal will be "
    "played by the master device. If you'd attach more regular slave devices after the AMmodulator, "
    "their signal would not get modulated, but simply added to the modulated signal of the first "
    "three devices.\n\n"
    "You can also modulate only the signals of a specific slave device, by attaching the modulator "
    "to that slave device in the 'OpenSlave' call. You'd simply pass the handle of the slave that "
    "should be modulated to 'OpenSlave', instead of the handle of a master device.\n\n";

    static char seeAlsoString[] = "Open ";

    double masterVolume;
    double *channelVolumes;
    int m, n, p, i;
    int pahandle = -1;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(3));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(2));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgRequired, &pahandle);
    if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");
    if ((audiodevices[pahandle].opmode & kPortAudioPlayBack) == 0) PsychErrorExitMsg(PsychError_user, "Audio device has not been opened for audio playback, so this call doesn't make sense.");

    // Return old masterVolume:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) audiodevices[pahandle].masterVolume);

    // Copy in optional new masterVolume: We can assign this without locking, as a single scalar
    // 32-bit float assignment is basically atomic:
    if (PsychCopyInDoubleArg(2, kPsychArgOptional, &masterVolume)) audiodevices[pahandle].masterVolume = (float) masterVolume;

    // Slaves can have per-channel volumes:
    if (audiodevices[pahandle].opmode & kPortAudioIsSlave) {
        // Slave device: Need whole vector:

        // Copy out old settings:
        PsychAllocOutDoubleMatArg(2, kPsychArgOptional, 1, audiodevices[pahandle].outchannels, 1, &channelVolumes);
        for (i = 0; i < audiodevices[pahandle].outchannels; i++) channelVolumes[i] = (double) audiodevices[pahandle].outChannelVolumes[i];

        // Get optional new settings:
        if (PsychAllocInDoubleMatArg(3, kPsychArgOptional, &m, &n, &p, &channelVolumes)) {
            // Valid?
            if (m * n != audiodevices[pahandle].outchannels || p != 1) PsychErrorExitMsg(PsychError_user, "Invalid channelVolumes vector for audio slave device provided. Number of elements doesn't match number of audio output channels!");

            // Assign, but with device mutex of master device held, so we don't update in
            // the middle of a mix cycle for our slave device:
            PsychPALockDeviceMutex(&audiodevices[audiodevices[pahandle].pamaster]);
            for (i = 0; i < audiodevices[pahandle].outchannels; i++) audiodevices[pahandle].outChannelVolumes[i]  = (float) channelVolumes[i];
            PsychPAUnlockDeviceMutex(&audiodevices[audiodevices[pahandle].pamaster]);
        }
    }
    else {
        // Superfluous channelVolumes for a non-slave provided?
        if (PsychAllocInDoubleMatArg(3, kPsychArgOptional, &m, &n, &p, &channelVolumes)) {
            PsychErrorExitMsg(PsychError_user, "Invalid channelVolumes vector for a non-slave device provided. Only slave devices accept this vector!");
        }
    }

    return(PsychError_none);
}

/* PsychPortAudio('GetDevices') - Enumerate all available sound devices.
 */
PsychError PSYCHPORTAUDIOGetDevices(void)
{
    static char useString[] = "devices = PsychPortAudio('GetDevices' [, devicetype] [, deviceIndex]);";
    static char synopsisString[] =
    "Returns 'devices', an array of structs, one struct for each available PortAudio device.\n\n"
    "If the optional parameter 'deviceIndex' is provided and the optional parameter 'devicetype' "
    "is set to [], then only returns a single struct with information about the device with index "
    "'deviceIndex'.\n\n"
    "Each struct contains information about its associated PortAudio device. The optional "
    "parameter 'devicetype' can be used to enumerate only devices of a specific class: \n"
    #ifdef PTB_USE_ASIO
    "1=Windows/DirectSound, 2=Windows/MME, 3=Windows/ASIO, 11=Windows/WDMKS, 13=Windows/WASAPI, "
    #else
    "1=Windows/DirectSound, 2=Windows/MME, 11=Windows/WDMKS, 13=Windows/WASAPI, "
    #endif
    "8=Linux/ALSA, 7=Linux/OSS, 12=Linux/JACK, 5=MacOSX/CoreAudio.\n\n"
    "On OS/X you'll usually only see devices for the CoreAudio API, a first-class audio subsystem. "
    "On Linux you may have the choice between ALSA, JACK and OSS. ALSA or JACK provide very low "
    "latencies and very good timing, OSS is an older system which is less capable but not very "
    "widespread in use anymore. On MS-Windows you'll have the choice between up to 5 different "
    "audio subsystems:\n"
    #ifdef PTB_USE_ASIO
    "- If you buy a sound card with ASIO drivers, then you can pick that API for low latency. It "
    "should give you comparable performance to OS/X or Linux.\n"
    #endif
    "WASAPI (on Windows-Vista and later), or WDMKS (on Windows-2000/XP) should provide ok latency.\n"
    "DirectSound is the next worst choice if you have hardware with DirectSound support.\n"
    "If everything else fails, you'll be left with MME, a completely unusable API for precise or "
    "low latency timing.\n"
    "\n";

    static char seeAlsoString[] = "Open GetDeviceSettings ";
    PsychGenericScriptType     *devices;
    const char *FieldNames[]= { "DeviceIndex", "HostAudioAPIId", "HostAudioAPIName", "DeviceName", "NrInputChannels", "NrOutputChannels",
                                "LowInputLatency", "HighInputLatency", "LowOutputLatency", "HighOutputLatency",  "DefaultSampleRate", "xxx" };
    int devicetype = -1;
    int deviceindex = -1;
    int count = 0;
    int i, ic, filteredcount;
    const PaDeviceInfo* padev = NULL;
    const PaHostApiInfo* hainfo = NULL;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(2));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(0)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgOptional, &devicetype);
    if (devicetype < -1) PsychErrorExitMsg(PsychError_user, "Invalid 'devicetype' provided. Valid are values of zero and greater.");

    PsychCopyInIntegerArg(2, kPsychArgOptional, &deviceindex);
    if (deviceindex < -1) PsychErrorExitMsg(PsychError_user, "Invalid 'deviceindex' provided. Valid are values of zero and greater.");

    // Provided deviceIndex overrides potentially provided deviceType, if any:
    if (deviceindex >= 0 && devicetype >=0) PsychErrorExitMsg(PsychError_user, "Provided 'deviceindex' and 'devicetype'! This is forbidden. Provide one or the other.");

    // Query number of devices and allocate out struct array:
    count = (int) Pa_GetDeviceCount();
    if (count > 0) {
        filteredcount = count;
        if (devicetype!=-1) {
            filteredcount = 0;
            // Filtering by host API requested: Do it.
            for (i=0; i<count; i++) {
                padev = Pa_GetDeviceInfo((PaDeviceIndex) i);
                hainfo = Pa_GetHostApiInfo(padev->hostApi);
                if (hainfo->type == devicetype) filteredcount++;
            }
        }

        // Only return one single struct if specific deviceindex given:
        if (deviceindex >= 0) filteredcount = 1;

        // Alloc output struct array:
        PsychAllocOutStructArray(1, kPsychArgOptional, filteredcount, 11, FieldNames, &devices);
    }
    else {
        PsychErrorExitMsg(PsychError_user, "PTB-ERROR: PortAudio can't detect any supported sound device on this system.");
    }

    // Iterate through device list:
    ic = 0;
    for (i=0; i<count; i++) {
        // Return info about deviceindex i if it matches the selected deviceindex or if
        // all devices should be returned:
        if ((deviceindex == -1) || (deviceindex == i)) {
            // Get info about deviceindex i:
            padev = Pa_GetDeviceInfo((PaDeviceIndex) i);
            hainfo = Pa_GetHostApiInfo(padev->hostApi);

            // Return info if devicetype doesn't matter or if it matches the required one:
            if ((devicetype==-1) || (hainfo->type == devicetype)) {
                // Fill slot ic of struct array with info of deviceindex i:
                PsychSetStructArrayDoubleElement("DeviceIndex", ic, i, devices);
                PsychSetStructArrayDoubleElement("HostAudioAPIId", ic, hainfo->type, devices);
                PsychSetStructArrayStringElement("HostAudioAPIName", ic, (char*) (hainfo->name), devices);
                PsychSetStructArrayStringElement("DeviceName", ic, (char*) (padev->name), devices);
                PsychSetStructArrayDoubleElement("NrInputChannels", ic, padev->maxInputChannels, devices);
                PsychSetStructArrayDoubleElement("NrOutputChannels", ic, padev->maxOutputChannels, devices);
                PsychSetStructArrayDoubleElement("LowInputLatency", ic, padev->defaultLowInputLatency, devices);
                PsychSetStructArrayDoubleElement("HighInputLatency", ic, padev->defaultHighInputLatency, devices);
                PsychSetStructArrayDoubleElement("LowOutputLatency", ic, padev->defaultLowOutputLatency, devices);
                PsychSetStructArrayDoubleElement("HighOutputLatency", ic, padev->defaultHighOutputLatency, devices);
                PsychSetStructArrayDoubleElement("DefaultSampleRate", ic, padev->defaultSampleRate, devices);
                // PsychSetStructArrayDoubleElement("xxx", ic, 0, devices);
                ic++;
            }
        }
    }

    return(PsychError_none);
}

/* PsychPortAudio('RunMode') - Select a different mode of operation.
 */
PsychError PSYCHPORTAUDIORunMode(void)
{
    static char useString[] = "oldRunMode = PsychPortAudio('RunMode', pahandle [,runMode]);";
    static char synopsisString[] =
    "Set general run mode to 'runMode' and/or return old runMode for a device 'pahandle'.\n"
    "The device must be open for this setting to take effect and playback must be stopped if "
    "one wants to change the setting. If playback isn't stopped, it will be forcefully stopped. "
    "Change of runMode is not supported on slave devices. "
    "The current/old runMode is returned in the optional return value 'oldRunMode'. "
    "'runMode' is the optional new runmode: At device open time, the runMode defaults to one. "
    "In mode zero, the audio hardware and all internal processing are completely stopped at end "
    "of audio playback. This reduces system ressource usage (both hardware and computation time), "
    "but may cause slightly longer latencies when re'Start'ing the device for playback of a different "
    "sound, e.g., a sequence of sounds. In 'runMode' 1, the audio hardware and processing don't shut "
    "down at the end of audio playback. Instead, everything remains active in a ''hot standby'' state. "
    "This allows to very quickly (with low latency) restart sound playback via the 'RescheduleStart' "
    "function. The downside is a permanent use of system ressources even if no sound is playing. "
    "Future runMode settings may provide more interesting options, stay tuned...\n";

    static char seeAlsoString[] = "Start Stop RescheduleStart ";

    int runMode = -1;
    int pahandle = -1;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(2));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgRequired, &pahandle);
    if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");

    // Copy in optional runMode value:
    PsychCopyInIntegerArg(2, kPsychArgOptional, &runMode);

    // Return current/old runMode:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, audiodevices[pahandle].runMode);

    // Set new runMode, if one was provided:
    if (runMode != -1) {
        // Doing this is a no-no on slaves:
        if (audiodevices[pahandle].opmode & kPortAudioIsSlave) PsychErrorExitMsg(PsychError_user, "Change of runmode is not allowed on slave devices!");

        // Stop engine if it is running:
        if (!Pa_IsStreamStopped(audiodevices[pahandle].stream)) Pa_StopStream(audiodevices[pahandle].stream);

        // Reset state:
        audiodevices[pahandle].state = 0;
        audiodevices[pahandle].reqstate = 255;

        if (runMode < 0 || runMode > 1) PsychErrorExitMsg(PsychError_user, "Invalid 'runMode' provided. Must be 0 or 1!");

        // Assign new runMode:
        audiodevices[pahandle].runMode = runMode;
    }

    return(PsychError_none);
}

/* PsychPortAudio('SetLoop') - Define a subrange of samples in the audio playback buffer to play.
 */
PsychError PSYCHPORTAUDIOSetLoop(void)
{
    static char useString[] = "PsychPortAudio('SetLoop', pahandle[, startSample=0][, endSample=max][, UnitIsSeconds=0]);";
    static char synopsisString[] =
    "Restrict audio playback to a subrange of sound samples in the current audio playback buffer for "
    "audio device 'pahandle'. The device must be open and a soundbuffer created via 'FillBuffer' "
    "for this setting to take effect. The setting is reset at each non-streaming 'FillBuffer' call.\n"
    "'startSample' defines the first sample to be played, it defaults to zero, ie. the very first in "
    "the buffer. 'endSample' the last sample to be played, it defaults to the end of the buffer.\n"
    "Calling this function without any parameter except 'pahandle' will therefore reset the playloop "
    "to the full buffer size. It is possible but not advisable to change the loop range while "
    "sound playback is active, as this may cause audible glitches or artifacts in the playback.\n"
    "By default, 'startSample' and 'endSample' are specified in units of samples, but if the "
    "optional flag 'UnitIsSeconds' is set to a non-zero value, then the given range is "
    "interpreted in units of seconds. ";

    static char seeAlsoString[] = "FillBuffer Start Stop RescheduleStart ";

    double startSample, endSample, sMultiplier;
    psych_int64 maxSample;
    int unitIsSecs;
    int pahandle = -1;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(4));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(0));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgRequired, &pahandle);
    if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");
    if ((audiodevices[pahandle].opmode & kPortAudioPlayBack) == 0) PsychErrorExitMsg(PsychError_user, "Audio device has not been opened for audio playback, so this call doesn't make sense.");

    unitIsSecs = 0;
    PsychCopyInIntegerArg(4, kPsychArgOptional, &unitIsSecs);
    sMultiplier = (unitIsSecs > 0) ? (double) audiodevices[pahandle].streaminfo->sampleRate : 1.0;

    // Compute maxSample the maximum possible sampleframe index for given soundbuffer:
    maxSample = (audiodevices[pahandle].outputbuffersize / sizeof(float) / audiodevices[pahandle].outchannels) - 1;

    // Copy in optional startSample:
    startSample = 0;
    PsychCopyInDoubleArg(2, kPsychArgOptional, &startSample);
    if (startSample < 0) PsychErrorExitMsg(PsychError_user, "Invalid 'startSample' provided. Must be greater or equal to zero!");
    startSample *= sMultiplier;

    // Copy in optional endSample:
    if (PsychCopyInDoubleArg(3, kPsychArgOptional, &endSample)) {
        endSample *= sMultiplier;
        if (endSample > maxSample) PsychErrorExitMsg(PsychError_user, "Invalid 'endSample' provided. Must be no greater than total buffersize!");
    }
    else {
        endSample = (double) maxSample;
    }

    if (endSample < startSample) PsychErrorExitMsg(PsychError_user, "Invalid 'endSample' provided. Must be greater or equal than 'startSample'!");

    // Ok, range is valid. Assign it:
    PsychPALockDeviceMutex(&audiodevices[pahandle]);
    audiodevices[pahandle].loopStartFrame = (psych_int64) startSample;
    audiodevices[pahandle].loopEndFrame = (psych_int64) endSample;
    PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);

    return(PsychError_none);
}

/* PsychPortAudio('EngineTunables') - Set tunable low-level engine parameters
 */
PsychError PSYCHPORTAUDIOEngineTunables(void)
{
    static char useString[] = "[oldyieldInterval, oldMutexEnable, lockToCore1, audioserver_autosuspend] = PsychPortAudio('EngineTunables' [, yieldInterval] [, MutexEnable] [, lockToCore1] [, audioserver_autosuspend]);";
    static char synopsisString[] =
    "Return, and optionally set low-level tuneable driver parameters.\n"
    "The driver must be idle, ie., no audio device must be open, if you want to change tuneables! "
    "These tuneable parameters usually have reasonably chosen defaults and you should only "
    "need to change them to work around bugs or flaws in your operating system, sound hardware or drivers, "
    "or if you have very unusual needs or setups. Only touch these if you know what you're doing, probably "
    "after consultation with the Psychtoolbox forum or Wiki. Some of these have potential to cause serious "
    "system malfunctions if not selected properly!\n\n"
    "'yieldInterval' - If the driver has to perform polling operations, it will release the cpu for "
    "yieldInterval seconds inbetween unsuccessful polling iterations. Valid range is 0.0 to 0.1 secs, with "
    "a reasonable default of 0.001 secs ie. 1 msec.\n"
    "'MutexEnable' - Enable (1) or Disable (0) internal mutex locking of driver data structures to prevent "
    "potential race-conditions between internal processing threads. Locking is enabled by default. Only "
    "disable locking to work around seriously broken audio device drivers or system setups and be aware "
    "that this may have unpleasant side effects and can cause all kinds of malfunctions by itself!\n"
    "'lockToCore1' - Deprecated: Enable (1) or Disable (0) locking of all audio engine processing threads to cpu core 1 "
    "on Microsoft Windows systems. By default threads are locked to cpu core 1 to avoid problems with "
    "timestamping due to bugs in some microprocessors clocks and in Microsoft Windows itself. If you're "
    "confident/certain that your system is bugfree wrt. to its clocks and want to get a bit more "
    "performance out of multi-core machines, you can disable this. You must perform this setting before "
    "you open the first audio device the first time, otherwise the setting might be ignored. In the current "
    "driver this setting is silently ignored, as a new method of handling this has been implemented.\n"
    "'audioserver_autosuspend' - Enable (1) or Disable (0) automatic suspending of running desktop "
    "audio servers, e.g., PulseAudio, while PsychPortAudio is active. Default is (1) - suspend while "
    "PsychPortAudio does its thing. Desktop sound servers like the commonly used PulseAudio server "
    "can interfere with low level audio device access and low-latency / high-precision audio timing. "
    "For this reason it is a good idea to switch them to standby (suspend) while a PsychPortAudio "
    "session is active. Sometimes this isn't needed or not even desireable. Therefore this option "
    "allows to inhibit this automatic suspending of audio servers.\n";

    static char seeAlsoString[] = "Open ";

    int mutexenable, mylockToCore1, mysuspend;
    double myyieldInterval;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(4));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(0)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(4));    // The maximum number of outputs

    // Make sure no settings are changed while an audio device is open:
    if ((PsychGetNumInputArgs() > 0) && (audiodevicecount > 0)) PsychErrorExitMsg(PsychError_user, "Tried to change low-level engine parameter while at least one audio device is open! Forbidden!");

    // Return current/old audioserver_suspend:
    PsychCopyOutDoubleArg(4, kPsychArgOptional, (double) ((pulseaudio_autosuspend) ? 1 : 0));

    // Get optional new audioserver_suspend:
    if (PsychCopyInIntegerArg(4, kPsychArgOptional, &mysuspend)) {
        if (mysuspend < 0 || mysuspend > 1) PsychErrorExitMsg(PsychError_user, "Invalid setting for 'audioserver_autosuspend' provided. Valid are 0 and 1.");
        pulseaudio_autosuspend = (mysuspend > 0) ? TRUE : FALSE;
        if (verbosity > 3) printf("PsychPortAudio: INFO: Automatic suspending of desktop audio servers %s.\n", (pulseaudio_autosuspend) ? "enabled" : "disabled");
    }

    // Make sure PortAudio is online: Must be done after setup of audioserver_suspend!
    PsychPortAudioInitialize();

    // Return old yieldInterval:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, yieldInterval);

    // Get optional new yieldInterval:
    if (PsychCopyInDoubleArg(1, kPsychArgOptional, &myyieldInterval)) {
        if (myyieldInterval < 0 || myyieldInterval > 0.1) PsychErrorExitMsg(PsychError_user, "Invalid setting for 'yieldInterval' provided. Valid are between 0.0 and 0.1 seconds.");
        yieldInterval = myyieldInterval;
        if (verbosity > 3) printf("PsychPortAudio: INFO: Engine yieldInterval changed to %lf seconds.\n", yieldInterval);
    }

    // Return current/old mutexenable:
    PsychCopyOutDoubleArg(2, kPsychArgOptional, (double) ((uselocking) ? 1 : 0));

    // Get optional new mutexenable:
    if (PsychCopyInIntegerArg(2, kPsychArgOptional, &mutexenable)) {
        if (mutexenable < 0 || mutexenable > 1) PsychErrorExitMsg(PsychError_user, "Invalid setting for 'MutexEnable' provided. Valid are 0 and 1.");
        uselocking = (mutexenable > 0) ? TRUE : FALSE;
        if (verbosity > 3) printf("PsychPortAudio: INFO: Engine Mutex locking %s.\n", (uselocking) ? "enabled" : "disabled");
    }

    // Return current/old lockToCore1:
    PsychCopyOutDoubleArg(3, kPsychArgOptional, (double) ((lockToCore1) ? 1 : 0));

    // Get optional new lockToCore1:
    if (PsychCopyInIntegerArg(3, kPsychArgOptional, &mylockToCore1)) {
        if (mylockToCore1 < 0 || mylockToCore1 > 1) PsychErrorExitMsg(PsychError_user, "Invalid setting for 'lockToCore1' provided. Valid are 0 and 1.");
        lockToCore1 = (mylockToCore1 > 0) ? TRUE : FALSE;
        if (verbosity > 3) printf("PsychPortAudio: INFO: Locking of all engine threads to cpu core 1 %s.\n", (lockToCore1) ? "enabled" : "disabled");
    }

    return(PsychError_none);
}

/* PsychPortAudio('UseSchedule') - Enable & Create, or disable and destroy a playback schedule.
 */
PsychError PSYCHPORTAUDIOUseSchedule(void)
{
    static char useString[] = "PsychPortAudio('UseSchedule', pahandle, enableSchedule [, maxSize = 128]);";
    static char synopsisString[] =
    "Enable or disable use of a preprogrammed schedule for audio playback on audio device 'pahandle'.\n"
    "Schedules are similar to playlists on your favorite audio player. A schedule allows to define a sequence "
    "of distinct sounds to play in succession. When PsychPortAudio('Start') is called, processing of "
    "the schedule begins and the first programmed sound snippet is played, followed by the 2nd, 3rd, ... "
    "until the whole schedule has been played once and playback stops. You can add new sound snippets to the "
    "schedule while playback is running for uninterrupted playback of long sequences of sound.\n"
    "Each sound snippet or slot in the schedule defines a soundbuffer to play back via a 'bufferhandle', "
    "within the buffer a subsegment (a so called playloop) defined by start- and endpoint, and a number "
    "of repetitions for that playloop.\n"
    "This subfunction allows to either enable use of schedules by setting the 'enableSchedule' flag to 1, "
    "in which case a schedule with a maximum of 'maxSize' distinct slots is created ('maxSize' defaults to 128 slots), "
    "or to disable use of schedules by setting 'enableSchedule' to 0, in which case PsychPortAudio reverts "
    "back to its normal playback behaviour and an existing schedule is deleted.\n"
    "A 'enableSchedule' setting of 2 will reset an existing schedule, ie. clear it of all its entries, "
    "so it is ready to be rewritten with new entries. You should reset and rewrite a schedule each "
    "time after playback/processing of a schedule has finished or has been stopped.\n"
    "A 'enableSchedule' setting of 3 will reactivate an existing schedule, ie. prepare it for a replay.\n"
    "See the subfunction 'AddToSchedule' on how to populate the schedule with actual entries.\n";

    static char seeAlsoString[] = "FillBuffer Start Stop RescheduleStart AddToSchedule";

    int pahandle = -1;
    int enableSchedule;
    int maxSize = 128;
    unsigned int j;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(3));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(2)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(0));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgRequired, &pahandle);
    if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");
    if ((audiodevices[pahandle].opmode & kPortAudioPlayBack) == 0) PsychErrorExitMsg(PsychError_user, "Audio device has not been opened for audio playback, so this call doesn't make sense.");

    // Make sure the device is fully idle: We can check without mutex held, as a device which is
    // already idle (state == 0) can't switch by itself out of idle state (state > 0), neither
    // can an inactive stream start itself.
    if ((audiodevices[pahandle].state > 0) && Pa_IsStreamActive(audiodevices[pahandle].stream)) PsychErrorExitMsg(PsychError_user, "Tried to enable/disable audio schedule while audio device is active. Forbidden! Call 'Stop' first.");

    // At this point the deivce is idle and will remain so during this routines execution,
    // so it won't touch any of the schedule related variables and we can manipulate them
    // without a need to lock the device:

    // Get required enable flag:
    PsychCopyInIntegerArg(2, kPsychArgRequired, &enableSchedule);
    if (enableSchedule < 0 || enableSchedule > 3)  PsychErrorExitMsg(PsychError_user, "Invalid 'enableSchedule' provided. Must be 0, 1, 2 or 3!");

    // Get the optional maxSize parameter:
    PsychCopyInIntegerArg(3, kPsychArgOptional, &maxSize);
    if (maxSize < 1) PsychErrorExitMsg(PsychError_user, "Invalid 'maxSize' provided. Must be greater than zero!");

    // Revival of existing schedule requested?
    if (enableSchedule == 3) {
        if (NULL == audiodevices[pahandle].schedule) {
            PsychErrorExitMsg(PsychError_user, "'enableSchedule' == 3 requested to revive current schedule, but no such schedule exists! You must create it first.");
        }

        // Reset current position in schedule to start:
        audiodevices[pahandle].schedule_pos = 0;

        for (j = 0; j < audiodevices[pahandle].schedule_size; j++) {
            // Slot occupied?
            if (audiodevices[pahandle].schedule[j].mode & 1) {
                // Reactivate this slot to pending:
                audiodevices[pahandle].schedule[j].mode |= 2;
            }
        }

        // Done.
        return(PsychError_none);
    }


    // Reset of existing schedule requested?
    if ((enableSchedule == 2) && (audiodevices[pahandle].schedule)) {
        // Yes: Simply set requested size to current size, this will trigger
        // a memset() to zero below, instead of a realloc or alloc, which is
        // exactly what we want:
        maxSize = audiodevices[pahandle].schedule_size;
    }

    // Release an already existing schedule: This will take care of both,
    // disabling use of schedules if this is a disable call, and resetting
    // of an existing schedule if this is an enable call following another
    // enable call:
    if (audiodevices[pahandle].schedule) {
        // Schedule already exists: Is this by any chance an enable call and
        // the requested size of the new schedule matches the size of the current
        // one?
        if (enableSchedule && (audiodevices[pahandle].schedule_size == maxSize)) {
            // Yes! Have a schedule of exactly wanted size, no need to free and
            // realloc - We simply clear it out by zero filling:
            memset((void*) audiodevices[pahandle].schedule, 0, (size_t) (maxSize * sizeof(PsychPASchedule)));
        }
        else {
            // No. Release old schedule...
            free(audiodevices[pahandle].schedule);
            audiodevices[pahandle].schedule = NULL;
            audiodevices[pahandle].schedule_size = 0;
        }
    }

    // Reset current position in schedule to start and size to zero in any case:
    audiodevices[pahandle].schedule_pos = 0;
    audiodevices[pahandle].schedule_writepos = 0;

    // Enable/Reset request?
    if (enableSchedule && (NULL == audiodevices[pahandle].schedule)) {
        // Enable request - Allocate proper schedule:
        audiodevices[pahandle].schedule_size = 0;
        audiodevices[pahandle].schedule = (PsychPASchedule*) calloc(maxSize, sizeof(PsychPASchedule));
        if (audiodevices[pahandle].schedule == NULL) PsychErrorExitMsg(PsychError_outofMemory, "Insufficient free system memory when trying to create a schedule!");

        // Assign new size:
        audiodevices[pahandle].schedule_size = maxSize;
    }

    // Done.
    return(PsychError_none);
}

/* PsychPortAudio('AddToSchedule') - Add command slots to a playback schedule.
 */
PsychError PSYCHPORTAUDIOAddToSchedule(void)
{
    static char useString[] = "[success, freeslots] = PsychPortAudio('AddToSchedule', pahandle [, bufferHandle=0][, repetitions=1][, startSample=0][, endSample=max][, UnitIsSeconds=0][, specialFlags=0]);";
    //                                                                                  1           2                 3                4                5                6                  7
    static char synopsisString[] =
    "Add a new item to an existing schedule for audio playback on audio device 'pahandle'.\n"
    "The schedule must have been created and enabled already via a previous call to 'UseSchedule'. "
    "The function returns if the addition of a new item was successfull via the return argument "
    "'success' (1=Success, 0=Failed), and the number of remaining free slots in 'freeslots'. "
    "Failure to add an item can happen if the schedule is full. If playback is running, you can "
    "simply retry after some time, because eventually the playback will consume and thereby free "
    "at least one slot in the schedule. If playback is stopped and you get this failure, you should "
    "reallocate the schedule with a bigger size via a proper call to 'UseSchedule'.\n"
    "Please note that after playback/processing of a schedule has finished by itself, or due to "
    "'Stop'ping the playback via the stop function, you should clear or reactivate the schedule and rewrite "
    "it, otherwise results at next call to 'Start' may be undefined. You can clear/reactivate a schedule "
    "efficiently without resizing it by calling 'UseSchedule' with an enableFlag of 2 or 3.\n\n"
    "The following optional paramters can be used to define the new slot in the schedule:\n"
    "'bufferHandle' Handle of the audio buffer which should be used for playback of this slot. "
    "The default value zero will play back the standard audio buffer created by a call to 'FillBuffer'.\n"
    "'repetitions' How often should playback of this slot be repeated. Fractional positive values are "
    "allowed, the value zero (ie. infinite repetition) is not allowed in this driver release.\n"
    "'startSample' and 'endSample' define a playback loop - a subsegment of the audio buffer to which "
    "playback should be restricted, and 'UnitIsSeconds' tells if the given loop boundaries are expressed "
    "in audio sample frames, or in seconds realtime. See the help for the 'SetLoop' function for more "
    "explanation.\n\n"
    "'specialFlags' is an optional parameter to pass in some special mode flags to alter processing "
    "of a slot: Defaults to zero. If the value '1' is used (or added), then a slot is not automatically "
    "disabled after it has been used up, but it will be reused on repeated execution of the schedule. "
    "You'll need to set this flag on all slots in a schedule if you want the schedule to auto-repeat "
    "without the need for manual reset commands."
    "\n\n"
    "This function can also be used to sneak special command slots into the schedule:\n"
    "If you specify a negative number for the 'bufferHandle' argument, then this actually "
    "defines a command slot instead of a regular playback slot, and the number is a command code "
    "that defines the action to take. For timing related actions, the 'repetitions' parameter is "
    "abused as a 'tWhen' time specification in seconds.\n"
    "Following command codes are currently defined: (Add numbers to combine options)\n"
    "1   = Pause audio playback (and capture) immediately, resume it at a given 'tWhen' target time.\n"
    "2   = (Re-)schedule the end of playback for given 'tWhen' target time.\n"
    "You must specify the type of 'tWhen' if you specify 1 or 2 as commandCode:"
    "+4  = 'tWhen' is an absolute GetSecs() style system time.\n"
    "+8  = 'tWhen' is a time delta to the last requested start time of playback.\n"
    "+16 = 'tWhen' is a time delta to the last actual start time of playback.\n"
    "+32 = 'tWhen' is a time delta to the last requested end time of playback. Broken & Defective!!\n"
    "+64 = 'tWhen' is a time delta to the last actual end time of playback. Broken & Defective!!\n"
    "\n"
    "E.g., you want to (re)start playback at a certain time, then you'd set 'bufferHandle' to -5, because "
    "command code would be 1 + 4 == 5, so negated it is -5. Then you'd specify the requested time in the "
    "'repetitions' parameter as an absolute time in seconds.\n\n";

    static char seeAlsoString[] = "FillBuffer Start Stop RescheduleStart UseSchedule";

    PsychPASchedule* slot;
    PsychPABuffer* buffer;
    int    slotid;
    double startSample, endSample, sMultiplier;
    psych_int64 maxSample;
    int unitIsSecs;
    int pahandle = -1;
    int bufferHandle = 0;
    unsigned int commandCode = 0;
    int specialFlags = 0;
    double repetitions = 1;
    int success = 0;
    int freeslots = 0;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(7));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(2));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgRequired, &pahandle);
    if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");
    if ((audiodevices[pahandle].opmode & kPortAudioPlayBack) == 0) PsychErrorExitMsg(PsychError_user, "Audio device has not been opened for audio playback, so this call doesn't make sense.");

    // Make sure there is a schedule available:
    if (audiodevices[pahandle].schedule == NULL) PsychErrorExitMsg(PsychError_user, "You tried to AddToSchedule, but use of schedules is disabled! Call 'UseSchedule' first to enable them.");

    // Get optional bufferhandle:
    PsychCopyInIntegerArg(2, kPsychArgOptional, &bufferHandle);

    // if (bufferHandle < 0) PsychErrorExitMsg(PsychError_user, "Invalid 'bufferHandle' provided. Must be greater or equal to zero, and a handle to an existing buffer!");
    if (bufferHandle < 0) {
        // This isn't a bufferHandle but a command code, ergo this ain't a
        // audio playout buffer, but a command buffer:
        commandCode = -bufferHandle;
        bufferHandle = 0;

        // Child protection: Must give a timespec type specifier if some time related action is requested:
        if ((commandCode & (1 | 2)) && !(commandCode & (4 | 8 | 16 | 32 | 64))) PsychErrorExitMsg(PsychError_user, "Invalid commandCode provided: You requested scheduled (re)start or end of operation, but didn't provide any of the required timespec-type specifiers!");
    }

    // If it is a non-zero handle, try to dereference from dynamic buffer:
    if (bufferHandle > 0) {
        // Deref bufferHandle: Issue error if no buffer with such a handle exists:
        buffer = PsychPAGetAudioBuffer(bufferHandle);

        // Validate matching output channel count:
        if (buffer->outchannels != audiodevices[pahandle].outchannels) {
            printf("PsychPortAudio-ERROR: Audio channel count %i of audiobuffer with handle %i doesn't match channel count %i of audio device!\n", (int) buffer->outchannels, bufferHandle, (int) audiodevices[pahandle].outchannels);
            PsychErrorExitMsg(PsychError_user, "Referenced audio buffer 'bufferHandle' has an audio channel count that doesn't match channels of audio device!");
        }
    }

    // Get optional repetitions. We abuse it also for the tWhen parameter if this
    // is a commandbuffer instead of an audio buffer:
    PsychCopyInDoubleArg(3, kPsychArgOptional, &repetitions);
    if ((repetitions < 0) && (commandCode == 0)) PsychErrorExitMsg(PsychError_user, "Invalid 'repetitions' provided. Must be a positive or zero number!");

    // Get loop parameters, if any:
    unitIsSecs = 0;
    PsychCopyInIntegerArg(6, kPsychArgOptional, &unitIsSecs);
    sMultiplier = (unitIsSecs > 0) ? (double) audiodevices[pahandle].streaminfo->sampleRate : 1.0;

    // Set maxSample to maximum integer: The scheduler (aka PsychPAProcessSchedule()) will test at runtime if the playloop extends
    // beyond valid playbuffer boundaries and clamp to end-of-buffer if needed, so this is safe:
    // Ok, not quite the maximum 64 bit signed integer, but 2^32 counts less. Why? Because we assign
    // maxSample to a double variable below, then that back to a int64. Due to limited precision of
    // the double data type, roundoff errors would cause a INT64_MAX to overflow the representable
    // max value of 2^63 for psych_int64, thereby causing storage of wrapped around negative value,
    // and the whole logic would blow up! Keep some (overly large) security margin to prevent this.
    // Oh, and we need to divide by max number of device channels, as otherwise some multiplication
    // in the audio schedule processing may overflow and wreak havoc:
    maxSample = ((INT64_MAX - (psych_int64) INT32_MAX)) / ((psych_int64) MAX_PSYCH_AUDIO_CHANNELS_PER_DEVICE);

    // Copy in optional startSample:
    startSample = 0;
    PsychCopyInDoubleArg(4, kPsychArgOptional, &startSample);
    if (startSample < 0) PsychErrorExitMsg(PsychError_user, "Invalid 'startSample' provided. Must be greater or equal to zero!");
    startSample *= sMultiplier;

    // Copy in optional endSample:
    if (PsychCopyInDoubleArg(5, kPsychArgOptional, &endSample)) {
        endSample *= sMultiplier;
        if (endSample > maxSample) PsychErrorExitMsg(PsychError_user, "Invalid 'endSample' provided. Must be no greater than total buffersize!");
    }
    else {
        endSample = (double) maxSample;
    }

    if (endSample < startSample) PsychErrorExitMsg(PsychError_user, "Invalid 'endSample' provided. Must be greater or equal than 'startSample'!");

    // Copy in optional specialFlags:
    PsychCopyInIntegerArg(7, kPsychArgOptional, &specialFlags);

    // All settings validated and ready to initialize a slot in the schedule:

    // Lock device:
    PsychPALockDeviceMutex(&audiodevices[pahandle]);

    // Map writepos to slotindex:
    slotid = audiodevices[pahandle].schedule_writepos % audiodevices[pahandle].schedule_size;

    // Enough unoccupied space in schedule? Ie., is this slot free (either never used, or already consumed and ready for recycling)?
    if ((audiodevices[pahandle].schedule[slotid].mode & 2) == 0) {
        // Fill slot:
        slot = (PsychPASchedule*) &(audiodevices[pahandle].schedule[slotid]);
        slot->mode = 1 | 2 | ((specialFlags & 1) ? 4 : 0);
        slot->bufferhandle   = bufferHandle;
        slot->repetitions    = (commandCode == 0) ? ((repetitions == 0) ? -1 : repetitions) : 0.0;;
        slot->loopStartFrame = (psych_int64) startSample;
        slot->loopEndFrame   = (psych_int64) endSample;
        slot->command         = commandCode;
        slot->tWhen             = (commandCode > 0) ? repetitions : 0.0;

        // Advance write position for next update iteration:
        audiodevices[pahandle].schedule_writepos++;

        // Recompute number of free slots:
        if (audiodevices[pahandle].schedule_size >= (audiodevices[pahandle].schedule_writepos - audiodevices[pahandle].schedule_pos)) {
            freeslots = audiodevices[pahandle].schedule_size - (audiodevices[pahandle].schedule_writepos - audiodevices[pahandle].schedule_pos);
        }
        else {
            freeslots = 0;
        }

        success = 1;
    }
    else {
        // Nope. No free slot:
        success = 0;
        freeslots = 0;
    }

    // Unlock device:
    PsychPAUnlockDeviceMutex(&audiodevices[pahandle]);

    // Return optional result code:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) success);

    // Return optional remaining number of free slots:
    PsychCopyOutDoubleArg(2, kPsychArgOptional, (double) freeslots);

    return(PsychError_none);
}

/* PsychPortAudio('SetOpMode') - Change opmode of an already opened device.
 */
PsychError PSYCHPORTAUDIOSetOpMode(void)
{
    static char useString[] = "oldOpMode = PsychPortAudio('SetOpMode', pahandle [, opModeOverride]);";
    static char synopsisString[] =
    "Override basic mode of operation of an open audio device 'pahandle' and/or return old/current mode.\n"
    "The device must be open for this setting to take effect and operations must be stopped. "
    "If the device isn't stopped, it will be forcefully stopped. "
    "The current/old opMode is returned in the optional return value 'oldOpMode'. "
    "'opModeOverride' is the optional new opMode: At device open time, the initial opMode is assigned "
    "via the 'mode' parameter of PsychPortAudio('Open', ...); and defaults to audio playback if omitted.\n\n"
    "The mode can only be changed within certain constraints. It is not possible to switch a device "
    "from playback to capture mode, or from simplex to full-duplex mode after it has been opened. It "
    "is possible to change other special opMode flags though, e.g., one can enable/disable the live (software based) "
    "low-latency monitoring mode of full-duplex devices by adding or not adding the value '4' to the opMode value.\n\n"
    "The function will ignore settings that can't be changed while the device is open.\n";

    static char seeAlsoString[] = "Start Stop RescheduleStart Open Close";

    // Mode bits defined in ignoreMask can not be changed/overriden by this function:
    // Currently the untouchables are device playback/capture/full-/halfduplex configuration, as
    // these require a full close/reopen cycle in PortAudio, as well as any master/slave setting:
    const int ignoreMask = (kPortAudioPlayBack | kPortAudioCapture | kPortAudioFullDuplex | kPortAudioIsMaster | kPortAudioIsSlave);
    int opMode = -1;
    int pahandle = -1;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(2));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    PsychCopyInIntegerArg(1, kPsychArgRequired, &pahandle);
    if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided.");

    // Copy in optional opMode value:
    PsychCopyInIntegerArg(2, kPsychArgOptional, &opMode);

    // Return current/old opMode:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, audiodevices[pahandle].opmode);

    // Set new opMode, if one was provided:
    if (opMode != -1) {
        // Stop engine if it is running:
        if (!Pa_IsStreamStopped(audiodevices[pahandle].stream)) Pa_StopStream(audiodevices[pahandle].stream);

        // Reset state:
        audiodevices[pahandle].state = 0;
        audiodevices[pahandle].reqstate = 255;

        if (opMode < 0) PsychErrorExitMsg(PsychError_user, "Invalid 'opModeOverride' provided. Check the 'mode' parameter in the help for PsychPortAudio('Open', ...)!");

        // Make sure that number of capture and playback channels is the same and device is in full-duplex mode for fast monitoring/feedback mode:
        if (opMode & kPortAudioMonitoring) {
            if (((audiodevices[pahandle].opmode & kPortAudioFullDuplex) != kPortAudioFullDuplex) || (audiodevices[pahandle].outchannels != audiodevices[pahandle].inchannels)) {
                PsychErrorExitMsg(PsychError_user, "Fast monitoring/feedback mode selected, but device is not in full-duplex mode or number of capture and playback channels differs! They must be the same for this mode!");
            }
        }

        // Combine old current mode with override mode, mask out untouchable bits:
        opMode = (audiodevices[pahandle].opmode & ignoreMask) | (opMode & (~ignoreMask));

        // Assign new opMode:
        audiodevices[pahandle].opmode = opMode;
    }

    return(PsychError_none);
}

/* PsychPortAudio('DirectInputMonitoring') - Enable/Disable or reconfigure direct input monitoring.
 */
PsychError PSYCHPORTAUDIODirectInputMonitoring(void)
{
    static char useString[] = "result = PsychPortAudio('DirectInputMonitoring', pahandle, enable [, inputChannel = -1][, outputChannel = 0][, gainLevel = 0.0][, stereoPan = 0.5]);";
    //                                                                          1         2         3                    4                    5                  6
    static char synopsisString[] =
    "Change the current settings for the \"direct input monitoring\" feature on device 'pahandle'.\n\n"
    #if (PSYCH_SYSTEM == PSYCH_OSX) && !defined(paMacCoreChangeDeviceParameters)
    "NOTE: The current PsychPortAudio driver supports direct input monitoring only on MacOSX with some pro-sound hardware, "
    "but not with the builtin sound chip. Also this feature is completely untested by the developer so far.\n\n"
    #else
    "NOTE: This function is not actually functional in this PsychPortAudio driver build. It is only here for backwards compatibility.\n\n"
    #endif
    "The device must be open for this setting to take effect. Changed settings may or may not "
    "persist across closing and opening the device, this is hardware dependent and not to be relied on.\n"
    "So-called \"Zero latency direct input monitoring\" is a hardware feature of some modern "
    "(and usually higher end) soundcards. It allows to directly feed audio signals that are received "
    "at the audio input connectors of the soundcard back to the output connectors, without any extended "
    "intermediate processing of the audio signals by either the sound hardware or the host computer and "
    "its software. Due to this direct signal path, which only applies selectable amplification and "
    "some stereo panning and rerouting, the feedback latency from input to output (e.g, microphone to "
    "speakers) is as minimal as technically possible. On many high-end cards it is instantaneous!\n\n"
    "The 'enable' flag is mandatory: If set to zero, monitoring will be disabled for the given 'inputChannel'. "
    "A setting of one will enable input live-monitoring of the given 'inputChannel' to the given 'outputChannel' with "
    "the selected other settings.\n"
    "All following settings are optional and have reasonable defaults. Depending on your hardware, some or all of them "
    "may be silently ignored by your sound hardware.\n"
    "The optional 'inputChannel' argument specifies which input audio channels monitoring settings should be modified. "
    "If omitted or set to -1, all input channels settings will be modified, or at least tried to be modified.\n"
    "The optional 'outputChannel' specifies the index of the base-channel of a channel stereo-pair to which the 'inputChannel' "
    "should be routed. It must be an even number like 0, 2, 4, .... If omitted, channel 0, i.e., the first output channel "
    "stereo pair will be used.\n"
    "The optional 'gainLevel' defines the desired amplifier gain for the routed signal. The value should be negative for "
    "signal attenuation (i.e., negative gain) and positive for amplification (i.e., positive gain). "
    "As specification of gains is not standardized across operating systems, the numbers you'll have to pass in for a desired "
    "effect will vary across operating systems and audio hardware. However, the default setting of zero tries to set a neutral "
    "gain of zero decibel - the signal is passed through without change in intensity. The setting may get completely ignored or "
    "only approximately implemented by given hardware. Double-check your results!\n"
    "The optional 'stereoPan' parameter allows to select panning between the two output channels of a selected stereo output "
    "channel pair if the hardware allows that. Range 0.0 - 1.0 selects between left-channel and right channel, with the default "
    "of 0.5 selecting a centered output with equal distribution to both channels.\n\n"
    "In the optional return argument 'result' the function returns a status code to report if the requested change could be carried "
    "out successfully. A value of zero means success. A value of 1 means some error, e.g., invalid parameters specified. A value of "
    "2 means that your combination of operating system, sound system, soundcard device driver and soundcard hardware does not support "
    "direct input monitoring, at least not for the given configuration. A setting of 3 means that your PortAudio driver plugin does "
    "not support the feature."
    "\n";

    static char seeAlsoString[] = "Open GetDeviceSettings ";

    int pahandle = -1;
    int enable, inputChannel, outputChannel, rc;
    double gain, stereoPan;
    const PaDeviceInfo* padev = NULL;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(6));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(2)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));    // The maximum number of outputs

    // Make sure PortAudio is online:
    PsychPortAudioInitialize();

    // Get mandatory device handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &pahandle);
    if (pahandle < 0 || pahandle>=MAX_PSYCH_AUDIO_DEVS || audiodevices[pahandle].stream == NULL) PsychErrorExitMsg(PsychError_user, "Invalid audio device handle provided. No such device with that handle open!");

    // Get mandatory enable flag:
    PsychCopyInIntegerArg(2, kPsychArgRequired, &enable);
    if (enable < 0 || enable > 1) PsychErrorExitMsg(PsychError_user, "Invalid enable flag provided. Must be zero or one for on or off!");

    // Copy in optional inputChannel id:
    if (PsychCopyInIntegerArg(3, kPsychArgOptional, &inputChannel)) {
        // Find out how many real input channels the device has and check provided index against them:
        padev = Pa_GetDeviceInfo((PaDeviceIndex) audiodevices[pahandle].indeviceidx);
        if (inputChannel < -1 || inputChannel >= (int) padev->maxInputChannels) PsychErrorExitMsg(PsychError_user, "Invalid inputChannel provided. No such input channel available on device!");
    }
    else {
        inputChannel = -1;
    }

    // Copy in optional outputChannel id:
    if (PsychCopyInIntegerArg(4, kPsychArgOptional, &outputChannel)) {
        // Find out how many real output channels the device has and check provided index against them:
        padev = Pa_GetDeviceInfo((PaDeviceIndex) audiodevices[pahandle].outdeviceidx);
        if (outputChannel < 0 || outputChannel >= (int) padev->maxOutputChannels) PsychErrorExitMsg(PsychError_user, "Invalid outputChannel provided. No such outputChannel channel available on device!");
    }
    else {
        outputChannel = 0;
    }

    // Copy in optional gain:
    gain = 0.0;
    PsychCopyInDoubleArg(5, kPsychArgOptional, &gain);

    // Copy in optional stereoPan:
    stereoPan = 0.5;
    PsychCopyInDoubleArg(6, kPsychArgOptional, &stereoPan);
    if (stereoPan < 0 || stereoPan > 1) PsychErrorExitMsg(PsychError_user, "Invalid stereoPan provided. Not in valid range 0.0 to 1.0!");

    // Default result code is "totally unsupported by our driver":
    rc = 3;

    // Feature currently only supported on OS/X with our custom build / extended Portaudio variant, not with Portaudio upstream:
    #if (PSYCH_SYSTEM == PSYCH_OSX) && !defined(paMacCoreChangeDeviceParameters)
    if (audiodevices[pahandle].hostAPI == paCoreAudio) {
        // CoreAudio device: Is the Portaudio plugin recent enough to support our custom directmonitoring interface?
        if (strstr(Pa_GetVersionText(), "WITH-DIM")) {
            PaError rcp;
            // Plugin supports the API, so at least we can safely call it without crashing.
            // Lower the fail level to rc = 2, can't fail because of our deficiencies anymore:
            if (verbosity > 4)
                printf("PsychPortAudio('DirectInputMonitoring'): Calling with padev=%i (%p), enable = %i, in=%i, out=%i, gain=%f, pan=%f.\n", pahandle, audiodevices[pahandle].stream, enable, inputChannel, outputChannel, gain, stereoPan);

            rcp = Pa_DirectInputMonitoring(audiodevices[pahandle].stream, enable, inputChannel, outputChannel, gain, stereoPan);
            switch (rcp) {
                case paNoError:
                    // Success!
                    rc = 0;
                    break;

                case paInvalidFlag:
                    // Some invalid request:
                    rc = 1;
                    break;

                case paBadIODeviceCombination:
                    // Unsupported by device:
                    rc = 2;
                    break;

                default:
                    // Default to unknown failure:
                    rc = 1;
            }
            if ((verbosity > 1) && (rc > 0)) {
                switch(rc) {
                    case 1:
                        printf("PsychPortAudio('DirectInputMonitoring'): Failed to change monitoring settings for call with padev=%i (%p), enable = %i, in=%i, out=%i, gain=%f, pan=%f.\n", pahandle, audiodevices[pahandle].stream, enable, inputChannel, outputChannel, gain, stereoPan);
                        printf("PsychPortAudio('DirectInputMonitoring'): Could be due to invalid parameters? Or some other unknown error. Rerun with PsychPortAudio('Verbosity', 10); for possible clues.\n");
                        break;

                    case 2:
                        printf("PsychPortAudio('DirectInputMonitoring'): Failed to change monitoring settings for call with padev=%i (%p), enable = %i, in=%i, out=%i, gain=%f, pan=%f.\n", pahandle, audiodevices[pahandle].stream, enable, inputChannel, outputChannel, gain, stereoPan);
                        printf("PsychPortAudio('DirectInputMonitoring'): DirectInputMonitoring is not supported by the target audio hardware.\n");
                        break;

                    default:
                        printf("PsychPortAudio('DirectInputMonitoring'): Failed to change monitoring settings for call with padev=%i (%p), enable = %i, in=%i, out=%i, gain=%f, pan=%f.\n", pahandle, audiodevices[pahandle].stream, enable, inputChannel, outputChannel, gain, stereoPan);
                        printf("PsychPortAudio('DirectInputMonitoring'): Reason unknown.\n");
                }
            }
        }
        else {
            if (verbosity > 1 && (audiodevices[pahandle].hostAPI == paCoreAudio))
                printf("PsychPortAudio('DirectInputMonitoring'): Your Portaudio plugin is too old to support this feature!\n");
        }
    }
    else {
        if (verbosity > 3) printf("PsychPortAudio('DirectInputMonitoring'): Tried to call, but feature not supported on this system.\n");
    }
    #else
    // Standard upstream Portaudio:
    if (verbosity > 3) printf("PsychPortAudio('DirectInputMonitoring'): Tried to call, but feature not supported on your operating system.\n");
    #endif

    // Return success status:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, rc);

    return(PsychError_none);
}
