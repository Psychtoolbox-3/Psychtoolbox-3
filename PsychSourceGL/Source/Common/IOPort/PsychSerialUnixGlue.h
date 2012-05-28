/*
	Psychtoolbox3/PsychSourceGL/Source/Common/IOPort/PsychSerialGlue.h		
    
	PROJECTS: 
	
		IOPort for now.
  
	AUTHORS:
	
		mario.kleiner at tuebingen.mpg.de	mk
  
	PLATFORMS:	
	
		All Unix (aka Posix) systems: OS/X and Linux.
    
	HISTORY:

		04/10/2008	mk		Initial implementation.
 
	DESCRIPTION:
	
		This is the operating system dependent "glue code" header file for access to serial ports for the
		Unices, ie., for GNU/Linux and Apple MacOS/X. It is used by the higher-level serial port
		routines to abstract out operating system dependencies.

		The code is shared and #ifdef'ed if needed, because most of it is identical for Linux and OS/X.
*/

#ifndef PSYCH_IS_INCLUDED_SerialUnixGlue
#define PSYCH_IS_INCLUDED_SerialUnixGlue

#include <unistd.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <errno.h>
#include <paths.h>
#include <termios.h>
#include <sysexits.h>
#include <sys/param.h>
#include <sys/select.h>
#include <sys/time.h>
#include <time.h>
#include <pthread.h>

// OS/X specific includes and structures:
#if PSYCH_SYSTEM == PSYCH_OSX

#include <AvailabilityMacros.h>
#include <CoreFoundation/CoreFoundation.h>
#include <IOKit/IOKitLib.h>
#include <IOKit/serial/IOSerialKeys.h>
#include <IOKit/serial/ioss.h>
#include <IOKit/IOBSD.h>

// End OS/X specific includes...
#endif

// Linux specific includes and structures:
#if PSYCH_SYSTEM == PSYCH_LINUX
// None yet.
#endif

typedef struct PsychSerialDeviceRecord {
	char				portSpec[1000];			// Name string of the device file.
	int					fileDescriptor;			// Device handle.
	struct termios		OriginalTTYAttrs;		// Stores original settings of device to allow restore on close.
	unsigned char*		readBuffer;				// Pointer to memory buffer for reading data.
	unsigned int		readBufferSize;			// Size of readbuffer.
	double				readTimeout;			// Backup copy of current read timeout value.
	double				pollLatency;			// Seconds to sleep between spin-wait polls in 'Read'.
	pthread_t			readerThread;			// Thread handle for background reading thread.
	pthread_mutex_t		readerLock;				// Primary lock.
	int					readerThreadWritePos;	// Position of next data write for readerThread.
	int					clientThreadReadPos;	// Position of next data read from main thread.
	int					readGranularity;		// Amount of bytes to request per blocking read call in readerThread.
	int					isBlockingBackgroundRead;  // 1 = Blocking background read, 0 = Polling operation.
	double*				timeStamps;				// Buffer for async-read timestamps. Size = readBufferSize / readGranularity Bytes.
	int					bounceBufferSize;		// Size of bounceBuffer in Bytes.
	unsigned char*		bounceBuffer;			// Bouncebuffer.
	unsigned int		readFilterFlags;		// Special flags to enable certain postprocessing operations on read data.
	int					asyncReadBytesCount;	// Counter of total bytes read via async thread so far. [Updates not mutex protected!]
	unsigned char		lineTerminator;			// Line terminator byte, if any.
	unsigned char		cookedMode;				// Cooked input processing mode active? Set to 1 if so.
	int					dontFlushOnWrite;		// If set to 1, don't tcdrain() after blocking writes, otherwise do.
	double				triggerWhen;			// Target time for trigger byte emission.
} PsychSerialDeviceRecord;

#endif
