/*
 *    IOPort.c
 *
 *    PROJECTS:
 *
 *        IOPort only.
 *
 *    AUTHORS:
 *
 *        Allen.Ingling@nyu.edu                   awi
 *        mario.kleiner.de@gmail.com              mk
 *
 *    PLATFORMS:
 *
 *        All.
 *
 *    HISTORY:
 *
 *        12/20/2004    awi       Wrote it.
 *        04/10/2008    mk        Started to extend/rewrite it to become a full-fledged generic I/O driver (serial port, parallel port, etc...).
 *
 *    DESCRIPTION:
 *
 */

//begin include once
#ifndef PSYCH_IS_INCLUDED_IOPort
#define PSYCH_IS_INCLUDED_IOPort

//project includes
#include "Psych.h"
#include "PsychTimeGlue.h"
#if PSYCH_SYSTEM != PSYCH_WINDOWS
#include "PsychSerialUnixGlue.h"
#else
#include "PsychSerialWindowsGlue.h"
#endif

// Special filterflags for input processing:
#define kPsychIOPortCMUPSTFiltering             1            // Filtering for CMU/PST button boxes.
#define kPsychIOPortCRLFFiltering               2            // Filtering for USB/32 Bitwhacker with StickOS.
#define kPsychIOPortAsyncLineBufferFiltering    4            // Filtering for emulation of line-buffering, like in "cooked" Unixish canonical input processing.

// Types of Input/Output port we support:
#define KPsychIOPortNone        0                // No port: This indicates a free slot.
#define kPsychIOPortSerial      1                // Serial port.

typedef struct PsychPortIORecord {
    unsigned int        portType;       // Type of I/O port, see defines above.
    void*               device;         // Opaque pointer to struct with device specific data - Different types need different structs...
} PsychPortIORecord;

// Operating system specific glue functions:
PsychSerialDeviceRecord* PsychIOOSOpenSerialPort(const char* portSpec, const char* configString, char* errmsg);
void PsychIOOSCloseSerialPort(PsychSerialDeviceRecord* device);
PsychError PsychIOOSConfigureSerialPort(PsychSerialDeviceRecord* device, const char* configString);
int PsychIOOSWriteSerialPort(PsychSerialDeviceRecord* device, void* writedata, unsigned int amount, int blocking, char* errmsg, double* timestamp);
int PsychIOOSReadSerialPort(PsychSerialDeviceRecord* device, void** readdata, unsigned int amount, int blocking, char* errmsg, double* timestamp);
int PsychIOOSBytesAvailableSerialPort(PsychSerialDeviceRecord* device);
void PsychIOOSFlushSerialPort(PsychSerialDeviceRecord* device);
void PsychIOOSPurgeSerialPort(PsychSerialDeviceRecord* device);
void PsychIOOSShutdownSerialReaderThread(PsychSerialDeviceRecord* device);

// Public subfunction prototypes
PsychError MODULEVersion(void);
#if PSYCH_SYSTEM == PSYCH_OSX
PsychError IOPORTSetPort(void);
PsychError IOPORTInit(void);
#endif
PsychError IOPORTDisplaySynopsis(void);
PsychError IOPORTVerbosity(void);

// Generic functions for port access -- The ones supposed to be useful for all types of ports:
PsychError IOPORTClose(void);
PsychError IOPORTCloseAll(void);
PsychError IOPORTRead(void);
PsychError IOPORTWrite(void);
PsychError IOPORTBytesAvailable(void);
PsychError IOPORTPurge(void);
PsychError IOPORTFlush(void);

// Serial port specific functions:
PsychError IOPORTOpenSerialPort(void);
PsychError IOPORTConfigureSerialPort(void);

// Initialize usage info -- function overview:
const char** InitializeSynopsis(void);

// Exit function: Called at shutdown time:
PsychError PsychExitIOPort(void);
// Init function: Called at startup time:
PsychError PsychInitIOPort(void);
// Close function: Closes port referenced by 'handle':
PsychError PsychCloseIOPort(int handle);
// Write function:
int PsychWriteIOPort(int handle, void* writedata, unsigned int amount, int blocking, char* errmsg, double* timestamp);
int    PsychReadIOPort(int handle, void** readbuffer, unsigned int amount, int blocking, char* errmsg, double* timestamp);
int PsychBytesAvailableIOPort(int handle);
void PsychPurgeIOPort(int handle);
void PsychFlushIOPort(int handle);

//end include once
#endif
