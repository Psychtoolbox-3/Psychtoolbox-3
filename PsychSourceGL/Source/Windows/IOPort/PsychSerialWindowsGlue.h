/*
	Psychtoolbox3/PsychSourceGL/Source/Windows/IOPort/PsychSerialWindowsGlue.h		
    
	PROJECTS: 
	
		IOPort for now.
  
	AUTHORS:
	
		mario.kleiner at tuebingen.mpg.de	mk
  
	PLATFORMS:	
	
		All Windows (aka Win32) systems: Windows 2000, Windows-XP, Windows-Vista, ...
    
	HISTORY:

		05/10/2008	mk		Initial implementation.
 
	DESCRIPTION:
	
		This is the operating system dependent "glue code" layer for access to serial ports for the
		Microsoft Windows platforms. It is used by the higher-level serial port
		routines to abstract out operating system dependencies.
*/

#ifndef PSYCH_IS_INCLUDED_SerialWindowsGlue
#define PSYCH_IS_INCLUDED_SerialWindowsGlue

#include <time.h>

typedef struct PsychSerialDeviceRecord {
	char				portSpec[1000];			// Name string of the device file - typically "COM0:", but we alloc large space for virtual device files etc.
	HANDLE				fileDescriptor;			// Device handle.
	DCB					OriginalTTYAttrs;		// Stores original settings of device to allow restore on close.
	COMMTIMEOUTS		Originaltimeouts;		// dto.
	COMMTIMEOUTS		timeouts;				// Structure with Win-specific timeout settings.
	unsigned char*		readBuffer;				// Pointer to memory buffer for reading data.
	unsigned int		readBufferSize;			// Size of readbuffer.
	double				readTimeout;			// Backup copy of current read timeout value.
	double				pollLatency;			// Seconds to sleep between spin-wait polls in 'Read'.	
	unsigned int        num_output_pending;		// Number of pending bytes for writeout.
	psych_thread		readerThread;			// Thread handle for background reading thread.
	psych_mutex			readerLock;				// Primary lock.
	unsigned int		abortThreadReq;			// Set to 1 to request a thread to abort.
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
	unsigned char		eventCharEnabled;		// 0 = Waiting for "event character received" event is disabled. Otherwise it is enabled.
	int					dontFlushOnWrite;		// If set to 1, don't tcdrain() after blocking writes, otherwise do.	
} PsychSerialDeviceRecord;

#endif
