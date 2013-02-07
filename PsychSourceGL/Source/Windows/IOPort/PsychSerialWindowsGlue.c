/*
	Psychtoolbox3/PsychSourceGL/Source/Windows/IOPort/PsychSerialWindowsGlue.c		
    
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

#include "IOPort.h"

// Internal function:
int PsychIOOSCheckError(PsychSerialDeviceRecord* device, char* inerrmsg);

// Externally defined level of verbosity:
extern int verbosity;

int PsychSerialWindowsGlueAsyncReadbufferBytesAvailable(PsychSerialDeviceRecord* device)
{	
	int navail = 0;
	
	// Lock data structure:
	PsychLockMutex(&(device->readerLock));
	
	// Compute amount of pending data for readout:
	navail = (device->readerThreadWritePos - device->clientThreadReadPos);

	// Unlock data structure:
	PsychUnlockMutex(&(device->readerLock));

	// Return it.
	return(navail);
}

void* PsychSerialWindowsGlueReaderThreadMain(void* deviceToCast)
{
	char errmsg[1024];
	int rc, nread;
	int navail;
	int tmpcurpos, naccumread;
	unsigned char lastcharacter, lineterminator;	
	double dt, oldt;
	double t;
	
	COMSTAT dstatus;
	DWORD   dummy;

	// Get a handle to our device struct: These pointers must not be NULL!!!
	PsychSerialDeviceRecord* device = (PsychSerialDeviceRecord*) deviceToCast;

	// Try to raise our priority: We ask to switch ourselves (NULL) to priority class 1 aka
	// round robin realtime scheduling, with a tweakPriority of +1, ie., raise the relative
	// priority level by +1 wrt. to the current level:
	if ((rc = PsychSetThreadPriority(NULL, 1, 1)) > 0) {
		if (verbosity > 0) fprintf(stderr, "PTB-ERROR: In IOPort:PsychSerialWindowsGlueReaderThreadMain(): Failed to switch to realtime priority [%i]!\n", rc);
	}

	// Perform initial check for device errors, clear error state:
	errmsg[0] = 0;
	if ((rc=PsychIOOSCheckError(device, &errmsg[0]))>0) {
		if (verbosity > 0) printf("PTB-ERROR: In startup of IOPort:PsychSerialWindowsGlueReaderThreadMain(): Errorcheck on device %s returned errorcode (%d) [%s]! Trying to recover.\n", device->portSpec, rc, errmsg);
	}

	// Init reference timestamp of last scan:
	PsychGetAdjustedPrecisionTimerSeconds(&oldt);
	dt = 0;

	// Main loop: Runs until external thread cancellation signal device->abortThreadReq set to non-zero:
	while (0 == device->abortThreadReq) {
		if (verbosity > 9) printf("IOPort-DEBUG: In PsychSerialWindowsGlueReaderThreadMain(): New loop iteration, abortflag = %i, checking abort signal.\n", device->abortThreadReq); 
		PsychTestCancelThread(&(device->readerThread));
		if (verbosity > 9) printf("IOPort-DEBUG: In PsychSerialWindowsGlueReaderThreadMain(): New loop iteration, abortflag = %i, abort signal negative.\n", device->abortThreadReq); 

		errmsg[0] = 0;
        if ((rc=PsychIOOSCheckError(device, &errmsg[0]))>0) {
			if (verbosity > 0) printf("PTB-ERROR: In IOPort:PsychSerialWindowsGlueReaderThreadMain(): Errorcheck on device %s returned errorcode (%d) [%s]! Trying to recover.\n", device->portSpec, rc, errmsg);
        }
        
		// Polling read requested?
		if (device->isBlockingBackgroundRead == 0) {
			// Polling operation:
			
			// Enough data available for read of requested granularity?
			if (ClearCommError(device->fileDescriptor, &dummy, &dstatus)==0) {
				if (verbosity > 0) fprintf(stderr, "PTB-ERROR: In IOPort:PsychSerialWindowsGlueReaderThreadMain(): ClearCommError() bytes available query on device %s returned error code (%d)!\n", device->portSpec, GetLastError());
			}
			navail = (int) dstatus.cbInQue;
			
			// If not, we sleep for some time, then retry:
			while(navail < device->readGranularity) {
				// Cancellation point: Test for cancel request and terminate, if so:
				if (device->abortThreadReq > 0) return(0);
				PsychTestCancelThread(&(device->readerThread));

				// Yield...
				PsychYieldIntervalSeconds(device->pollLatency);

				// ... and retest:
				if (ClearCommError(device->fileDescriptor, &dummy, &dstatus)==0) {
					if (verbosity > 0) fprintf(stderr, "PTB-ERROR: In IOPort:PsychSerialWindowsGlueReaderThreadMain(): ClearCommError() bytes available query on device %s returned error code (%d)!\n", device->portSpec, GetLastError());
				}
				navail = (int) dstatus.cbInQue;
			}
		}
		else {
			// Non-polling operation. We perform a blocking read on the device.
			// Assuming the masterthread doesn't perform any non-blocking operations,
			// which would screw us or it, this is safe. The read() call is a thread cancellation
			// point, so we can be safely aborted by the masterthread even if blocked on
			// input.

			// Set filedescriptor to blocking mode:
			device->timeouts.ReadIntervalTimeout = (int) (device->readTimeout * 1000 + 0.5);
			device->timeouts.ReadTotalTimeoutMultiplier = device->timeouts.ReadIntervalTimeout;
			device->timeouts.ReadTotalTimeoutConstant = 0;
			
			if (SetCommTimeouts(device->fileDescriptor, &(device->timeouts)) == 0) {
				if (verbosity > 0) fprintf(stderr, "PTB-ERROR: In IOPort:PsychSerialWindowsGlueReaderThreadMain(): Error calling SetCommTimeouts on device %s for blocking read - (%d).\n", device->portSpec, GetLastError());
			}
		}

		// Zerofill the space we're gonna read in next read() request, so we have a nice
		// clean buffersegment in case of a short-read, e.g., in cooked mode on end-of-line:
		memset(&(device->readBuffer[(device->readerThreadWritePos) % (device->readBufferSize)]), 0, device->readGranularity);
		
		// Async linebuffered read op?
		if (device->readFilterFlags & kPsychIOPortAsyncLineBufferFiltering) {
			// Emulation of linebuffered readop, similar to Unix cooked, canonical input processing mode:
			tmpcurpos = (device->readerThreadWritePos) % (device->readBufferSize);
			lineterminator = device->lineTerminator;
			lastcharacter  = lineterminator + 1;
			naccumread = 0;
			t = DBL_MIN;

			// Repeat reading until a maximum of readGranularity bytes has been read or until
			// lineterminator character is encountered, whatever comes first:
			while ( (naccumread < device->readGranularity) && (lastcharacter != lineterminator) ) {
                // Check for termination request:
                if (device->abortThreadReq > 0) return(0);
				PsychTestCancelThread(&(device->readerThread));

				// Read 1 Byte:
				if (ReadFile(device->fileDescriptor, &lastcharacter, 1, &nread, NULL) == 0) {
					// Some error:
					nread = -1;
				}
				
				if (nread > 0) {
					// Put it into readBuffer:
					device->readBuffer[tmpcurpos++] = lastcharacter;

					// Get representative timestamp of "start of line terminated new line" for first received byte:
					if (0 == naccumread) PsychGetAdjustedPrecisionTimerSeconds(&t);

					// Increment count of received bytes:
					naccumread++;
				}
				else {
					// 1 Byte read failed!

					// Diagnostic:
					if (verbosity > 5) fprintf(stderr, "PTB-ERROR: In IOPort:PsychSerialWindowsGlueReaderThreadMain(): Linebuffered read: Failed to read 1 byte of data for reason [%d] at relative position %i! Padding...\n", GetLastError(), naccumread);
					
					// Break out of while() loop:
					break;
				}
			}
			
			// Done with this read quantum, either due to error, line-terminator reached, or readGranularity bytes stored.

			// Increment serial bytes received counter:
			device->asyncReadBytesCount += (naccumread > 0) ? naccumread : 0;			
		}
		else {
			// Standard non-linebuffered readop:
			
			// How much data to read? We use the last 8-Bytes of a readGranularity quantum for our
			// serial bytes counter and our dt to last scan time if kPsychIOPortCMUPSTFiltering is active:
			naccumread = (device->readFilterFlags & kPsychIOPortCMUPSTFiltering) ? (device->readGranularity - 8) : device->readGranularity;
			if (naccumread < 0) naccumread = 0;
			
			// Enough data available. Read it!
			if (ReadFile(device->fileDescriptor, &(device->readBuffer[(device->readerThreadWritePos) % (device->readBufferSize)]), naccumread, &nread, NULL) == 0) {
				// Some error:
				if (verbosity > 0) fprintf(stderr, "PTB-ERROR: In IOPort:PsychSerialWindowsGlueReaderThreadMain(): Error during blocking read from device %s - (%d).\n", device->portSpec, GetLastError());
			}
			
			// Sufficient amount read?
			if (nread != naccumread) {
                // Should not happen, unless device is in cooked (canonical) input processing mode, where any read() will
                // at most return the content of a single line of buffered input, regardless of requested amount. In this
                // case it is not only a perfectly valid and expected result, but also safe, because our memset() above will
                // zero-fill the remainder of the buffer with a defined value. For this reason we only output an error
                // if high verbosity level is selected for debug output.

                // One special case: Zero bytes read due to read timeout exceeded. In this case we just skip and retry,
                // instead of storing a complete zero-filled block of readGranularity bytes. Padding is nice, but returning
                // completely empty data, which is indistinguishable from a sequence of zeros is not good. This caused serious
                // pain when trying to receive one-byte scanner triggers, and it timed out due to scan not yet started and the
                // valid scanner trigger bytes are value==zero bytes, indistinguishable by usercode from timeout. This injected
                // "ghost triggers" from the usercodes perspective, very bad! Note that the actual timeout case is handled by
                // the high level read code by independent checking for timeout. The high level code will hit timeout if we
                // hit this "continue" too often, and will return a empty [] variable, reliably signalling to calling code that
                // timeout occured - something that can be cleanly handled by caller.
                if (nread == 0) continue;

                if (verbosity > 5) fprintf(stderr, "PTB-ERROR: In IOPort:PsychSerialWindowsGlueReaderThreadMain(): Failed to read %i bytes of data for unknown reason (Got only %i bytes)! Padding...\n", naccumread, nread);
            }
			
			// Take read completion timestamp:
			PsychGetAdjustedPrecisionTimerSeconds(&t);
			
			// Compute timedelta to last read() scan:
			dt = t - oldt;
			oldt = t;
			
			// Increment serial bytes received counter:
			device->asyncReadBytesCount += (nread > 0) ? nread : 0;
			
			// Filtermode for filtering out CR and LF characters active (e.g., for UBW32-Bitwhacker with StickOS)?
			if ((device->readFilterFlags & kPsychIOPortCRLFFiltering) &&
				((device->readBuffer[(device->readerThreadWritePos) % (device->readBufferSize)] == 10) ||
				 (device->readBuffer[(device->readerThreadWritePos) % (device->readBufferSize)] == 13))) {
				// Current read byte is code 10 or 13 aka CR or LF. Discard & Skip:
				continue;
			}
			
			// Filtermode for CMU button box or PST button box enabled?
			if (device->readFilterFlags & kPsychIOPortCMUPSTFiltering) {
				// Special input data filter for the CMU button box and the PST button box.
				// Both boxes are hillarious masterpieces of totally braindamaged protocol design.
				// They send a continous stream of status bytes, at a rate of 1000 Hz (!?!), regardless
				// if the status of the box has changed or not, instead of just sending a status update
				// when actually something has changed. This creates a lot of load on the host computer
				// and a s***load of redundant data. As these shoddy beasts are still sold to customers,
				// and quite widespread, we implement special filtering. We check each received byte if
				// it matches its predecessor. If so, we discard it, as it is redundant.
				if ((device->readerThreadWritePos > 0) && (device->readBuffer[(device->readerThreadWritePos) % (device->readBufferSize)] == lastcharacter)) {
					// Current read byte value is identical to last stored value.
					// --> No status change, therefore no reason to store this redundant value.
					// We skip processing and wait for the next byte:
					continue;
				}
				else {
					// Store current character as "lastcharacter" reference for next iteration:
					lastcharacter = device->readBuffer[(device->readerThreadWritePos) % (device->readBufferSize)];
					
					// Store new counter as a 32-bit unsigned int, which may possibly be not 32-bit boundary aligned
					// on the target architecture!
					*((unsigned int*) &(device->readBuffer[(device->readerThreadWritePos+1) % (device->readBufferSize)])) = (unsigned int) device->asyncReadBytesCount;
					// Store dt as a 32 bit unsigned int: It contains dt in microseconds - That resolution should be more than sufficient!
					*((unsigned int*) &(device->readBuffer[(device->readerThreadWritePos+5) % (device->readBufferSize)])) = (unsigned int) (dt * 1e6);
				}
			}
		}	// End of regular non-linebuffered readop.

		// Store timestamp for this read chunk of data:
		device->timeStamps[(device->readerThreadWritePos / device->readGranularity) % (device->readBufferSize / device->readGranularity)] = t;
		
		// Try to lock, block until available if not available:
		if ((rc=PsychLockMutex(&(device->readerLock)))) {
			// This could potentially kill Matlab, as we're printing from outside the main interpreter thread.
			// Use fprintf() instead of the overloaded printf() (aka mexPrintf()) in the hope that we don't
			// wreak havoc -- maybe it goes to the system log, which should be safer...
			fprintf(stderr, "PTB-ERROR: In IOPort:PsychSerialWindowsGlueReaderThreadMain(): mutex_lock failed.\n");
			
			// Commit suicide:
			return(0);
		}

		// Update linear write pointer:
		device->readerThreadWritePos += device->readGranularity;
		
		// Need to unlock the mutex:
		if ((rc=PsychUnlockMutex(&(device->readerLock)))) {
			// This could potentially kill Matlab, as we're printing from outside the main interpreter thread.
			// Use fprintf() instead of the overloaded printf() (aka mexPrintf()) in the hope that we don't
			// wreak havoc -- maybe it goes to the system log, which should be safer...
			fprintf(stderr, "PTB-ERROR: In IOPort:PsychSerialWindowsGlueReaderThreadMain(): Last mutex_unlock in termination failed.\n");
			
			// Commit suicide:
			return(0);
		}

		// Next iteration...
	}
	
	// Go and die peacefully...
	return(0);
}

void PsychIOOSShutdownSerialReaderThread( PsychSerialDeviceRecord* device)
{
	if (device->readerThread) {
		// Signal the thread that it should cancel:
		device->abortThreadReq = 1;
		
		if (verbosity > 6) printf("IOPort-DEBUG: In PsychIOOSShutdownSerialReaderThread(): Calling PsychAbortThread()...\n"); 
		PsychAbortThread(&(device->readerThread));
		
		// Wait for it to die:
		if (verbosity > 6) printf("IOPort-DEBUG: In PsychIOOSShutdownSerialReaderThread(): Calling PsychDeleteThread()...\n"); 
		PsychDeleteThread(&(device->readerThread));

		// Mark it as dead:
		device->readerThread = NULL;
		
		// Release the mutex:
		if (verbosity > 6) printf("IOPort-DEBUG: In PsychIOOSShutdownSerialReaderThread(): Calling PsychDestroyMutex()...\n"); 
		PsychDestroyMutex(&(device->readerLock));

		// Reset cancel signal:
		device->abortThreadReq = 0;

		// Release timestamp buffer:
		free(device->timeStamps);
		device->timeStamps = NULL;
	}
	
	return;
}

/* PsychIOOSOpenSerialPort()
 *
 * Open a serial port device and configure it.
 *
 * portSpec - String with the port name / device name of the serial port device.
 * configString - String with port configuration parameters.
 * errmsg - Pointer to char[] buffer in which error messages should be returned, if any.
 *
 * On success, allocate a PsychSerialDeviceRecord with all relevant settings,
 * return a pointer to it.
 *
 * Otherwise abort with error message.
 */
PsychSerialDeviceRecord* PsychIOOSOpenSerialPort(const char* portSpec, const char* configString, char* errmsg)
{
    HANDLE			fileDescriptor = INVALID_HANDLE_VALUE;
    DCB				options;
	PsychSerialDeviceRecord* device = NULL;
	psych_bool			usererr = FALSE;
	unsigned int	myerrno;
	
	// Init errmsg error message to empty == no error:
	errmsg[0] = 0;

    // Open the serial port read/write, for exclusive access with normal attributes and default other settings:
    fileDescriptor = CreateFile(portSpec, GENERIC_READ | GENERIC_WRITE, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (fileDescriptor == INVALID_HANDLE_VALUE)
    {
		myerrno = GetLastError();
		if (myerrno == ERROR_FILE_NOT_FOUND) {
			sprintf(errmsg, "Error opening serial port device %s - No such serial port device exists! (%d) [ENOENT].\n", portSpec, myerrno);
			usererr = TRUE;
		}
		else if (myerrno == ERROR_SHARING_VIOLATION) {
			sprintf(errmsg, "Error opening serial port device %s - The serial port is already open in this or some other application, close it first! (Code = %d) [ERROR_SHARING_VIOLATION].\n", portSpec, myerrno);
			usererr = TRUE;
		}
		else if (myerrno == ERROR_ACCESS_DENIED) {
			sprintf(errmsg, "Error opening serial port device %s - The system didn't grant access (Code = %d) [ERROR_ACCESS_DENIED]. Maybe a permission or configuration problem?\n", portSpec, myerrno);
			usererr = TRUE;
		}
		else {
			sprintf(errmsg, "Error opening serial port device %s - Errorcode (%d).\n", portSpec, myerrno);
		}
        goto error;
    }
    
    // Create the device struct and init it:
	device = calloc(1, sizeof(PsychSerialDeviceRecord));
	device->fileDescriptor = INVALID_HANDLE_VALUE;
	device->readBuffer = NULL;
	device->readBufferSize = 0;
	
    // Get the current options and save them so we can restore the default settings later.
	device->OriginalTTYAttrs.DCBlength = sizeof(DCB);
    if (GetCommState(fileDescriptor, &(device->OriginalTTYAttrs)) == 0) {
        sprintf(errmsg, "Error getting original device settings for device %s - (%d).\n", portSpec, GetLastError());
        goto error;
    }

	// Get and save the current timeout settings as well:
    if (GetCommTimeouts(fileDescriptor, &(device->timeouts)) == 0) {
        sprintf(errmsg, "Error getting current device comm timeouts for device %s - (%d).\n", portSpec, GetLastError());
		goto error;
    }
	
    if (GetCommTimeouts(fileDescriptor, &(device->Originaltimeouts)) == 0) {
        sprintf(errmsg, "Error getting original device comm timeouts for device %s - (%d).\n", portSpec, GetLastError());
		goto error;
    }

	strncpy(device->portSpec, portSpec, sizeof(device->portSpec));
	device->fileDescriptor = fileDescriptor;
    options = device->OriginalTTYAttrs;

	// Ok, the device is basically open. Call the reconfiguration routine with the setting
	// string. It will do all further setup work:
	if (PsychError_none != PsychIOOSConfigureSerialPort(device, configString)) {
		// Something failed. Error handling...
        sprintf(errmsg, "Error changing device settings for device %s - (%d).\n", portSpec, GetLastError());
		usererr = TRUE;
        goto error;
	}

	if (device->readBuffer == NULL) {
        sprintf(errmsg, "Error for device %s - No InputBuffer allocated! You must specify the 'InputBuffer' size in the configuration.\n", portSpec);
		usererr = TRUE;
        goto error;
	}

	// Setup for blocking serial port write operations: Tell that a transmit-buffer-empty event shall be signalled:
	// This will set the EV_TXEMPTY whenever a write operation completes (and resets it at start of each new write
	// operation!). The WaitCommEvent() routine in PsychIOOSWriteSerialPort() will wait for such events when a
	// blocking write is requested by usercode.
	//
	// N.B.: This used to be part of the blocking write routine in PsychIOOSWriteSerialPort(), but we
	// learned that this call is awfully slow for unknown reasons, and it seems to sync us to the USB
	// duty-cycle on Serial-over-USB ports, which is not what we want! As switching this on and off on
	// demand would be also awfully slow, and as it doesn't do any harm for non-blocking writes AFAIK,
	// we now always unconditionally enable it already here during device setup, where extra delays of
	// 1-2 msecs don't matter.
	//
	// If eventCharEnabled is set, we also enable EV_RXFLAG events...
	if (SetCommMask(device->fileDescriptor, ((device->eventCharEnabled > 0) ? (EV_TXEMPTY | EV_RXFLAG) : EV_TXEMPTY))==0) {
		sprintf(errmsg, "Error during setup for blocking write operations to device %s during call to SetCommMask(EV_TXEMPTY) - (%d).\n", portSpec, GetLastError());
		goto error;
	}

    // Success! Return Pointer to new device structure:
    return(device);
    
    // Failure path: Called on error with error message in errmsg:
error:

		if (fileDescriptor != INVALID_HANDLE_VALUE)
		{
			// File open.
			if (device) {
				// Settings retrieved or made?
				if (device->fileDescriptor != INVALID_HANDLE_VALUE) {
					// Have original settings available: Try to restore them:
					if (SetCommState(fileDescriptor, &options) == 0) {
						if (verbosity > 1) printf("WARNING: In error handling: Error restoring tty attributes %s - (%d).\n", portSpec, GetLastError());
					}
					
					// Restore original timeout settings as well:
					if (SetCommTimeouts(fileDescriptor, &(device->Originaltimeouts)) == 0) {
						if (verbosity > 1) printf("WARNING: In error handling: Error restoring comm timeouts for device %s - (%d).\n", portSpec, GetLastError());
					}				
				}
				
				// Release read buffer, if any:
				if (device->readBuffer) free(device->readBuffer);
				
				// Release device struct:
				free(device);
			}
			
			// Close file:
			CloseHandle(fileDescriptor);
		}
    
	// Return with error message:
	if (verbosity > 0) {
		PsychErrorExitMsg(((usererr) ? PsychError_user : PsychError_system), errmsg);
	}
	
	return(NULL);
}

/* PsychIOOSCloseSerialPort()
 *
 * Close serial port connection/device referenced by given 'device' record.
 * Free the device record afterwards.
 */
void PsychIOOSCloseSerialPort(PsychSerialDeviceRecord* device)
{
	if (device == NULL) PsychErrorExitMsg(PsychError_internal, "NULL-Ptr instead of valid device pointer!");
	
	PsychIOOSShutdownSerialReaderThread(device);
	
	// Drain all send-buffers:
    // Block until all written output has been sent from the device.
    // Note that this call is simply passed on to the serial device driver. 
	// See tcsendbreak(3) ("man 3 tcsendbreak") for details.
    if ((!device->dontFlushOnWrite) && (FlushFileBuffers(device->fileDescriptor) == 0))
    {
		if (verbosity > 1) printf("IOPort: WARNING: While trying to close serial port %s : Error waiting for drain - (%d).\n", device->portSpec, GetLastError());
    }
    
    // Traditionally it is good practice to reset a serial port back to
    // the state in which you found it. This is why the original settings
    // were saved.
    if (SetCommState(device->fileDescriptor, &(device->OriginalTTYAttrs)) == 0)
    {
		if (verbosity > 1) printf("IOPort: WARNING: While trying to close serial port %s : Could not restore original port settings - (%d).\n", device->portSpec, GetLastError());
    }

	// Restore original timeout settings as well:
	if (SetCommTimeouts(device->fileDescriptor, &(device->Originaltimeouts)) == 0) {
		if (verbosity > 1) printf("IOPort: WARNING: While trying to close serial port %s : Could not restore original timeout settings - (%d).\n", device->portSpec, GetLastError());
	}
				
	// Close device:
    CloseHandle(device->fileDescriptor);

	// Release read buffer, if any:
	if (device->readBuffer) free(device->readBuffer);
	if (device->bounceBuffer) free(device->bounceBuffer);

	// Release memory for device struct:
	free(device);
	
	// Done.
	return;
}

/* PsychIOOSConfigureSerialPort()
 *
 * (Re-)configure serial port connection/device referenced by given 'device' record.
 * The relevant parameters are stored in 'configString'.
 *
 * Returns success- or failure status.
 */
PsychError PsychIOOSConfigureSerialPort( PsychSerialDeviceRecord* device, const char* configString)
{
    DCB				options;
	char*			p;
	float			infloat;
	int				inint, inint2, rc;
	psych_bool			updatetermios = FALSE;

	// Init DCBlength field to size of structure! Fixed as of 31.12.2008 -- May have caused
	// trouble in some setups with earlier releases!
	memset(&options, 0, sizeof(DCB));
	options.DCBlength = sizeof(DCB);
	
	// Retrieve current settings:
    if (GetCommState(device->fileDescriptor, &options) == 0)
    {
        if (verbosity > 0) printf("Error getting current serial port device settings for device %s - (%d).\n", device->portSpec, GetLastError());
		return(PsychError_system);
    }
	
	// Set opmode to binary mode, give a warning if 'Cooked' mode is requested, as only binary mode is supported on Windows:
	options.fBinary = 1;
	if (strstr(configString, "ProcessingMode=Cooked") && (verbosity > 1)) printf("IOPort: Warning: Requested 'ProcessingMode=Cooked' not supported on Windows. Will use binary mode aka raw mode.\n");

    // Print the current input and output baud rates.
    // See tcsetattr(4) ("man 4 tcsetattr") for details.
    if (verbosity > 3) {
		printf("IOPort-Info: Configuration for device %s:\n", device->portSpec);
		printf("IOPort-Info: Current baud rate is %d\n", (int) options.BaudRate);
    }
	
	if ((p = strstr(configString, "ReceiveTimeout="))) {
		// Set timeouts for read operations:
		if ((1!=sscanf(p, "ReceiveTimeout=%f", &infloat)) || (infloat < 0)) {
			if (verbosity > 0) printf("Invalid parameter for ReceiveTimeout set! Typo, or negative value provided.\n");
			return(PsychError_user);
		}
		else {
			// Set timeout: It is in granularity of milliseconds.
			// Clamp to a minimum of 1 millisecond if a non-zero timeout is requested:
			if ((infloat < 0.001) && (infloat > 0)) {
				if (verbosity > 1) printf("IOPort: Warning: Requested per-byte 'ReceiveTimeout' value %f secs is positive but smaller than supported minimum of 1 millisecond on Windows. Changed value to 0.001 secs.\n", infloat);
				infloat = 0.001f;
			}

            // Make sure we don't change non-mutex-protected variables behind the
            // back of our readerThread by only allowing this function to be called
            // with inactive thread:
            if (device->readerThread) {
                if (verbosity > 0) printf("Assigned ReceiveTimeout= while background read operations are already enabled! Disable first via 'StopBackgroundRead'!\n");
                return(PsychError_user);			
            }
			
            device->readTimeout = infloat;
		}
	}
	
	if ((p = strstr(configString, "SendTimeout="))) {
		// Set timeouts for write operations:
		if ((1!=sscanf(p, "SendTimeout=%f", &infloat)) || (infloat < 0)) {
			if (verbosity > 0) printf("Invalid parameter for SendTimeout set! Typo, or negative value provided.\n");
			return(PsychError_user);
		}
		else {
			// Set timeout: It is in granularity of milliseconds.
			// Clamp to a minimum of 1 millisecond if a non-zero timeout is requested:
			if ((infloat < 0.001) && (infloat > 0)) {
				if (verbosity > 1) printf("IOPort: Warning: Requested per-byte 'SendTimeout' value %f secs is positive but smaller than supported minimum of 1 millisecond on Windows. Changed value to 0.001 secs.\n", infloat);
				infloat = 0.001f;
			}

			device->timeouts.WriteTotalTimeoutMultiplier = (int) (infloat * 1000 + 0.5);
			device->timeouts.WriteTotalTimeoutConstant = 0;
			
			if (SetCommTimeouts(device->fileDescriptor, &(device->timeouts)) == 0) {
				if (verbosity > 0) printf("IOPort: Error: Failed to set requested write timeout %f seconds on serial port %s : SetCommTimeouts() failed with error code %d.\n", infloat, device->portSpec, GetLastError());
				return(PsychError_system);
			}
		}
	}

	// Assign line terminator if any:
	if ((p = strstr(configString, "Terminator="))) {
		if (1!=sscanf(p, "Terminator=%i", &inint)) {
			if (verbosity > 0) printf("Invalid parameter for Terminator= set!\n");
			return(PsychError_invalidIntegerArg);
		}

        // Make sure we don't change non-mutex-protected variables behind the
        // back of our readerThread by only allowing this function to be called
        // with inactive thread:
        if (device->readerThread) {
            if (verbosity > 0) printf("Assigned Terminator= while background read operations are already enabled! Disable first via 'StopBackgroundRead'!\n");
            return(PsychError_user);			
        }

		device->lineTerminator = (unsigned char) inint;
		
		// We also assign the line terminator as Event Character: This allows
		// some line-delimited input buffering with the lineTerminator as
		// delimiter in synchronous reads.
		// It may also set the lineTerminator as event character and enable it
		// in Serial-USB converters from FTDI, thereby allowing to drastically
		// lower receive latency on reception of lineTerminator delimited data
		// packets:
		if (inint >= 0) {
			// Enable event character:
			options.EvtChar = (char) inint;
			device->eventCharEnabled = 1;
		}
		else {
			// Disable use of event characters (maybe):
			options.EvtChar = (char) 255;
			device->eventCharEnabled = 0;
		}
		updatetermios = TRUE;

		// Need to update the event mask as well, depending on the new enable setting
		// for event characters. See comments in PsychIOOSOpenSerialPort(); for SetCommMask() call as well!
		if (SetCommMask(device->fileDescriptor, ((device->eventCharEnabled > 0) ? (EV_TXEMPTY | EV_RXFLAG) : EV_TXEMPTY))==0) {
			if (verbosity > 0) printf("IOPort: Error: SetCommMask() call for 'Terminator' setup failed on device %s - (%d).\n", device->portSpec, GetLastError());
			return(PsychError_system);
		}
	}

	// Assign receiver enable state, if any:
	if ((p = strstr(configString, "ReceiverEnable="))) {
		if (1!=sscanf(p, "ReceiverEnable=%i", &inint)) {
			if (verbosity > 0) printf("Invalid parameter for ReceiverEnable= set!\n");
			return(PsychError_invalidIntegerArg);
		}

		// 'ReceiverEnable' is so far unsupported on Windows. Don't know if the OS has any
		// way to allow to change the setting.
		/*
		if (inint > 0) {
			// Enable receiver:
			
		}
		else {
			// Disable receiver:
		}

		updatetermios = TRUE;
		*/
	}

	// Set common baud rate for send and receive:
	if ((p = strstr(configString, "BaudRate="))) {
		if ((1!=sscanf(p, "BaudRate=%i", &inint)) || (inint < 1)) {
			if (verbosity > 0) printf("Invalid parameter %i for BaudRate set!\n", inint);
			return(PsychError_invalidIntegerArg);
		}
		else {
			// Set common speed for input- and output queues:
			options.BaudRate = inint;
			updatetermios = TRUE;
		}
	}

	// Generation and handling of parity bits:
	if ((p = strstr(configString, "Parity="))) {
		// Clear all parity settings:
		options.fParity = 0;
		options.Parity  = NOPARITY;
		
		if (strstr(p, "Parity=Even")) {
			options.fParity = 1;
			options.Parity  = EVENPARITY;
		}
		else 
		if (strstr(p, "Parity=Odd")) {
			options.fParity = 1;
			options.Parity  = ODDPARITY;
		}
		else
		if (strstr(p, "Parity=None")) {
			// Nothing to do.
			options.fParity = 0;
			options.Parity  = NOPARITY;
		}
		else {
			// Invalid parity spec:
			if (verbosity > 0) printf("Invalid parity setting %s not accepted! (Valid are None, Even and Odd", p);
			return(PsychError_user);
		}
		updatetermios = TRUE;		
	}
	
	// Handling of Break conditions: Not reproducible in a OS/X & Linux compatible way on MS-Windows to my knowledge,
	// so we simply ignore it for now. Imho even on OS/X and Linux this parameter is not very useful...
	//	if ((p = strstr(configString, "BreakBehaviour="))) {
	//		// Clear all break settings:
	//		options.c_iflag &= ~(IGNBRK | BRKINT);
	//		
	//		if (strstr(configString, "BreakBehaviour=Ignore")) {
	//			// Break signals will be ignored:
	//			options.c_iflag |= IGNBRK;
	//		}
	//		else 
	//		if (strstr(configString, "BreakBehaviour=Flush")) {
	//			// A break signal will flush the input- and output-queues:
	//			options.c_iflag |= BRKINT;
	//		}
	//		else
	//		if (strstr(configString, "BreakBehaviour=Zero")) {
	//			// Don't change anything, ie. IGNBRK and BRKINT not set:
	//			// A break signal will be output as single zero byte.
	//		}
	//		else {
	//			// Invalid spec:
	//			printf("Invalid break behaviour %s not accepted (Valid: Ignore, Flush, or Zero)!", p);
	//			return(PsychError_user);
	//		}
	//		updatetermios = TRUE;
	//	}

	// Handling of data bits:
	if ((p = strstr(configString, "DataBits="))) {
		// Clear all databit settings:
		options.ByteSize = 0;
		
		if (strstr(p, "DataBits=5")) {
			options.ByteSize = 5;
		}
		else 
		if (strstr(p, "DataBits=6")) {
			options.ByteSize = 6;
		}
		else
		if (strstr(p, "DataBits=7")) {
			options.ByteSize = 7;
		}
		else
		if (strstr(p, "DataBits=8")) {
			options.ByteSize = 8;
		}
		else
		if (strstr(p, "DataBits=16")) {
			options.ByteSize = 16;
		}
		else {
			// Invalid spec:
			if (verbosity > 0) printf("Invalid setting for data bits %s not accepted! (Valid: 5, 6, 7, 8 or 16 databits)", p);
			return(PsychError_user);
		}
		updatetermios = TRUE;	
	}
	
	// Handling of stop bits:
	if ((p = strstr(configString, "StopBits="))) {
		// Clear all stopbit settings:
		options.StopBits = 0;
		
		if (strstr(p, "StopBits=1")) {
			options.StopBits = ONESTOPBIT;
		}
		else 
		if (strstr(p, "StopBits=2")) {
			options.StopBits = TWOSTOPBITS;
		}
		else {
			// Invalid spec:
			if (verbosity > 0) printf("Invalid setting for stop bits %s not accepted! (Valid: 1 or 2 stopbits)", p);
			return(PsychError_user);
		}
		updatetermios = TRUE;	
	}
	
	// Handling of flow control:
	if ((p = strstr(configString, "FlowControl="))) {
		// Clear all flow control settings:
		options.fOutX = 0;
		options.fInX = 0;
		options.fRtsControl = 0;
		options.fDtrControl = 0;
		options.fOutxCtsFlow = 0;
		options.fOutxDsrFlow = 0;

		if (strstr(p, "FlowControl=None")) {
			// Set above already...
		}
		else 
		if (strstr(p, "FlowControl=Software")) {
			options.fOutX = 1;
			options.fInX = 1;
		}
		else
		if (strstr(p, "FlowControl=Hardware")) {
			options.fRtsControl = RTS_CONTROL_HANDSHAKE;
			options.fOutxCtsFlow = 1;
		}
		else {
			// Invalid spec:
			if (verbosity > 0) printf("Invalid setting for flow control %s not accepted! (Valid: None, Software, Hardware)", p);
			return(PsychError_user);
		}
		updatetermios = TRUE;		
	}

    // Print the new input and output baud rates.
    if (verbosity > 3) {
		// Output new baud rate:
		printf("IOPort-Info: Baud rate changed to %d\n", (int) options.BaudRate);
    }
    
	// Handling of DTR :
	if ((p = strstr(configString, "DTR="))) {
		if (strstr(p, "DTR=0")) {
			options.fDtrControl = DTR_CONTROL_DISABLE;
		}
		else
		if (strstr(p, "DTR=1")) {
			options.fDtrControl = DTR_CONTROL_ENABLE;
		}
		else {
			// Invalid spec:
			if (verbosity > 0) printf("Invalid setting for DTR %s not accepted! (Valid: 1 or 0)", p);
			return(PsychError_user);
		}		
		updatetermios = TRUE;		
	}
	
	// Handling of RTS :
	if ((p = strstr(configString, "RTS="))) {
		if (strstr(p, "RTS=0")) {
			options.fRtsControl = RTS_CONTROL_DISABLE;
		}
		else
		if (strstr(p, "RTS=1")) {
			options.fRtsControl = RTS_CONTROL_ENABLE;
		}
		else {
			// Invalid spec:
			if (verbosity > 0) printf("Invalid setting for RTS %s not accepted! (Valid: 1 or 0)", p);
			return(PsychError_user);
		}		
		updatetermios = TRUE;		
	}

    // Cause the new options to take effect immediately.
    if (updatetermios && (SetCommState(device->fileDescriptor, &options) == 0))
    {
        if (verbosity > 0) printf("Error setting new serial port configuration attributes for device %s - (%d).\n", device->portSpec, GetLastError());
        return(PsychError_system);
    }

	// Set input buffer size for receive ops:
	if ((p = strstr(configString, "InputBufferSize="))) {
        // Make sure we don't change non-mutex-protected variables behind the
        // back of our readerThread by only allowing this function to be called
        // with inactive thread:
        if (device->readerThread) {
            if (verbosity > 0) printf("Assigned InputBufferSize= while background read operations are already enabled! Disable first via 'StopBackgroundRead'!\n");
            return(PsychError_user);			
        }
        
		if ((1!=sscanf(p, "InputBufferSize=%i", &inint)) || (inint < 1)) {
			if (verbosity > 0) printf("Invalid parameter for InputBufferSize set! Typo or requested buffer size smaller than 1 byte.\n");
			return(PsychError_invalidIntegerArg);
		}
		else {
			// Set InputBufferSize:
			if (device->readBuffer) {
				// Buffer already exists. Try to realloc:
				p = (char*) realloc(device->readBuffer, inint);
				if (p == NULL) {
					// Realloc failed:
					if (verbosity > 0) printf("Reallocation of Inputbuffer of device %s for new size %i failed!", device->portSpec, inint);
					return(PsychError_outofMemory);
				}
				// Worked. Assign new values:
				device->readBuffer = (unsigned char*) p;
				device->readBufferSize = (unsigned int) inint;
			}
			else {
				// No Buffer yet. Allocate
				p = malloc(inint);
				if (p == NULL) {
					// Realloc failed:
					if (verbosity > 0) printf("Allocation of Inputbuffer of device %s of size %i failed!", device->portSpec, inint);
					return(PsychError_outofMemory);
				}
				
				// Worked. Assign new values:
				device->readBuffer = (unsigned char*) p;
				device->readBufferSize = (unsigned int) inint;
			}
		}
	}

	// Set buffer sizes of the hardware itself, ie., the driver-internal buffers:
	if ((p = strstr(configString, "HardwareBufferSizes="))) {
		if ((2!=sscanf(p, "HardwareBufferSizes=%i,%i", &inint, &inint2)) || (inint < 1) || (inint2 < 1)) {
			if (verbosity > 0) printf("Invalid parameter for HardwareBufferSizes set! Typo or at least one of the requested buffer sizes smaller than 1 byte.\n");
			return(PsychError_invalidIntegerArg);
		}
		else {
			// Set HardwareBufferSizes: This functionality is only available on the Windows platform.
			// It is only a hint or polite request to the underlying driver. The driver is free to
			// ignore the request and choose whatever size or buffering strategy it finds appropriate:
			if (!SetupComm(device->fileDescriptor, inint, inint2)) {
				// Failed:
				if (verbosity > 0) printf("Setup of explicit HardwareBufferSizes for device %s of inputsize %i and outputsize %i failed with error code %d!\n", device->portSpec, inint, inint2, GetLastError());
				return(PsychError_system);
			}
		}
	}

	if ((p = strstr(configString, "PollLatency="))) {
		// Set polling latency for busy-waiting read operations:
		if ((1!=sscanf(p, "PollLatency=%f", &infloat)) || (infloat < 0)) {
			if (verbosity > 0) printf("Invalid parameter for PollLatency set! Typo, or negative value provided.\n");
			return(PsychError_user);
		}
		else {
            // Make sure we don't change non-mutex-protected variables behind the
            // back of our readerThread by only allowing this function to be called
            // with inactive thread:
            if (device->readerThread) {
                if (verbosity > 0) printf("Assigned PollLatency= while background read operations are already enabled! Disable first via 'StopBackgroundRead'!\n");
                return(PsychError_user);			
            }

			device->pollLatency = infloat;
		}
	}

	if ((p = strstr(configString, "BlockingBackgroundRead="))) {
		if (1!=sscanf(p, "BlockingBackgroundRead=%i", &inint)) {
			if (verbosity > 0) printf("Invalid parameter for BlockingBackgroundRead= set!\n");
			return(PsychError_invalidIntegerArg);
		}

        // Make sure we don't change non-mutex-protected variables behind the
        // back of our readerThread by only allowing this function to be called
        // with inactive thread:
        if (device->readerThread) {
            if (verbosity > 0) printf("Assigned BlockingBackgroundRead= while background read operations are already enabled! Disable first via 'StopBackgroundRead'!\n");
            return(PsychError_user);			
        }

		device->isBlockingBackgroundRead = inint;
	}

	if ((p = strstr(configString, "ReadFilterFlags="))) {
		if (1!=sscanf(p, "ReadFilterFlags=%i", &inint)) {
			if (verbosity > 0) printf("Invalid parameter for ReadFilterFlags= set!\n");
			return(PsychError_invalidIntegerArg);
		}

        // Make sure we don't change non-mutex-protected variables behind the
        // back of our readerThread by only allowing this function to be called
        // with inactive thread:
        if (device->readerThread) {
            if (verbosity > 0) printf("Assigned ReadFilterFlags= while background read operations are already enabled! Disable first via 'StopBackgroundRead'!\n");
            return(PsychError_user);			
        }
		
        device->readFilterFlags = (unsigned int) inint;
	}

	// Stop a background reader?
	if ((p = strstr(configString, "StopBackgroundRead"))) {
		PsychIOOSShutdownSerialReaderThread(device);
	}

	// Async background read via parallel thread requested?
	if ((p = strstr(configString, "StartBackgroundRead="))) {
		if (1!=sscanf(p, "StartBackgroundRead=%i", &inint)) {
			if (verbosity > 0) printf("Invalid parameter for StartBackgroundRead set!\n");
			return(PsychError_invalidIntegerArg);
		}
		else {
			// Set data-fetch granularity for single background read requests:
			if (inint < 1) {
				if (verbosity > 0) printf("Invalid StartBackgroundRead fetch granularity of %i bytes provided. Must be at least 1 byte!\n", inint);
				return(PsychError_invalidIntegerArg);
			}
			
			if (device->readBufferSize < (unsigned int) inint) {
				if (verbosity > 0) printf("Invalid StartBackgroundRead fetch granularity of %i bytes provided. Bigger than current Inputbuffer size of %i bytes!\n", inint, device->readBufferSize);
				return(PsychError_invalidIntegerArg);
			}

			if ((device->readBufferSize % inint) != 0) {
				if (verbosity > 0) printf("Invalid StartBackgroundRead fetch granularity of %i bytes provided. Inputbuffer size of %i bytes is not an integral multiple of fetch granularity as required!\n", inint, device->readBufferSize);
				return(PsychError_invalidIntegerArg);
			}

			if (device->readerThread) {
				if (verbosity > 0) printf("Called StartBackgroundRead, but background read operations are already enabled! Disable first via 'StopBackgroundRead'!\n");
				return(PsychError_user);			
			}
			
			// Setup data structures:
			device->asyncReadBytesCount = 0;			
			device->readerThreadWritePos = 0;
			device->clientThreadReadPos  = 0;
			device->readGranularity = inint;
			
			// Allocate sufficiently large timestamp buffer:
			device->timeStamps = (double*) calloc(sizeof(double), device->readBufferSize / device->readGranularity);
			
			// Create & Init the mutex:
			if ((rc=PsychInitMutex(&(device->readerLock)))) {
				printf("PTB-ERROR: In StartBackgroundReadCould(): Could not create readerLock mutex lock.\n");
				return(PsychError_system);
			}
            
            // Perform lock->unlock mutex sequence to inject some memory ordering barriers here, so all our
            // settings are picked up by the newborn thread:
            if ((rc=PsychLockMutex(&(device->readerLock))) || (rc=PsychUnlockMutex(&(device->readerLock)))) {
				printf("PTB-ERROR: In StartBackgroundRead(): Could not lock + unlock readerLock mutex lock [%s].\n", strerror(rc));
				return(PsychError_system);
            }
			
			// Create and startup thread:
			if ((rc=PsychCreateThread(&(device->readerThread), NULL, PsychSerialWindowsGlueReaderThreadMain, (void*) device))) {
				printf("PTB-ERROR: In StartBackgroundReadCould(): Could not create background reader thread.\n");
				return(PsychError_system);
			}
		}
	}

	if ((p = strstr(configString, "DontFlushOnWrite="))) {
		if (1!=sscanf(p, "DontFlushOnWrite=%i", &inint)) {
			if (verbosity > 0) printf("Invalid parameter for DontFlushOnWrite= set!\n");
			return(PsychError_invalidIntegerArg);
		}
		device->dontFlushOnWrite = inint;
	}

	// Done.
	return(PsychError_none);
}

/* PsychIOOSWriteSerialPort()
 *
 * Write data to serial port:
 * writedata = Ptr. to data buffer of uint8 values to write.
 * amount = Buffersize aka amount of bytes to write.
 * blocking = 0 --> Async, 1 --> Flush buffer and wait for write completion.
 * errmsg = Pointer to char array where error messages should be put to.
 * timestamp = Pointer to a 4 element vectorof double values where write-completion timestamps should be put to.
 * [0] = Final write completion timestamp.
 * [1] = Timestamp immediately before write() submission.
 * [2] = Timestamp immediately after write() submission.
 * [3] = Timestamp of last check for write completion.
 *
 * Returns number of bytes written, or -1 on error.
 */
int PsychIOOSWriteSerialPort(PsychSerialDeviceRecord* device, void* writedata, unsigned int amount, int blocking, char* errmsg, double* timestamp)
{
	int nwritten;
	DWORD EvtMask;
	
	// Nonblocking mode?
	if (blocking <= 0) {
		// Non-blocking write:

		// Ok, this is problematic: There isn't a simple way to make the WriteFile non-blocking, so for now, our
		// "non-blocking write" is a (probably blocking under certain conditions) WriteFile() op, just with the
		// FlushFileBuffer call missing.

		// Write the data: Take pre- and postwrite timestamps.
		PsychGetAdjustedPrecisionTimerSeconds(&timestamp[1]);
		if (WriteFile(device->fileDescriptor, writedata, amount, &nwritten, NULL) == 0)
		{
			sprintf(errmsg, "Error during non-blocking write to device %s - (%d).\n", device->portSpec, GetLastError());
			return(-1);
		}
		PsychGetAdjustedPrecisionTimerSeconds(&timestamp[2]);
	}
	else {
		// Blocking write:

		// Write the data: Take pre- and postwrite timestamps.
		PsychGetAdjustedPrecisionTimerSeconds(&timestamp[1]);
		if (WriteFile(device->fileDescriptor, writedata, amount, &nwritten, NULL) == 0)
		{
			sprintf(errmsg, "Error during blocking write to device %s - (%d).\n", device->portSpec, GetLastError());
			return(-1);
		}
		PsychGetAdjustedPrecisionTimerSeconds(&timestamp[2]);
		
		// Flush the write buffer and wait for write completion on physical hardware:
		if ((!device->dontFlushOnWrite) && (FlushFileBuffers(device->fileDescriptor) == 0)) {
			sprintf(errmsg, "Error during blocking write to device %s while draining the write buffers - (%d).\n", device->portSpec, GetLastError());
			return(-1);
		}
		
		// Special polling mode to wait for transmit completion instead of WaitCommEvent()?
		if (blocking == 2) {
			// Yes. Use a tight polling loop which spin-waits until the driver output queue is empty (contains zero pending bytes):
			device->num_output_pending = nwritten;
			while (device->num_output_pending > 0) {
				// Take timestamp for this iteration:
				PsychGetAdjustedPrecisionTimerSeconds(&timestamp[3]);

				// Poll: This updates num_output_pending as side effect:
				if (PsychIOOSCheckError(device, errmsg) > 0) {
					return(-1);
				}
			}
		}
		else {
			// Take timestamp for completeness although it doesn't make much sense in the blocking case:
			PsychGetAdjustedPrecisionTimerSeconds(&timestamp[3]);

			// Wait until a event mask with EV_TXEMPTY event signalled is received:
			EvtMask = 0;
			while (0 == (EvtMask & EV_TXEMPTY)) {
				// Really wait for write completion, ie., until the last byte has left the transmitter:
				// (N.B.: No need to clear or init EvtMask, it is a pure return parameter)
				if (WaitCommEvent(device->fileDescriptor, (LPDWORD) &EvtMask, NULL)==0) {
					sprintf(errmsg, "Error during blocking write to device %s during call to WaitCommEvent(EV_TXEMPTY) - (%d).\n", device->portSpec, GetLastError());
					return(-1);
				}
			}
		}
	}
	
	// Take write completion timestamp:
	PsychGetAdjustedPrecisionTimerSeconds(timestamp);

	// Check for communication errors, acknowledge them, clear the error state and return error message, if any:
	if (PsychIOOSCheckError(device, errmsg)) return(-1); 
	
	// Write successfully completed if we reach this point. Return:
	errmsg[0] = 0;

	return(nwritten);
}

int PsychIOOSReadSerialPort( PsychSerialDeviceRecord* device, void** readdata, unsigned int amount, int blocking, char* errmsg, double* timestamp)
{
	double timeout;
	int raPos, i;	
	int nread = 0;	
	*readdata = NULL;

	// Clamp 'amount' of data to be read to receive buffer size:
	if (amount > device->readBufferSize) {
		// Too much. Is amount unspecified aka INT_MAX? In that case,
		// we'll clamp it to readbuffers size. Otherwise we abort with
		// error.
		if (amount == INT_MAX) {
			// Clamp to buffersize:
			amount = device->readBufferSize;
		}
		else {
			// amount specified and impossible to satisfy request at
			// current readbuffer size. Abort with error and tell
			// user to resize readbuffer:
			sprintf(errmsg, "Amount of requested data %i is more than device %s can satisfy, as its input buffer is too small (%i bytes).\nSet a bigger readbuffer size please.\n", amount, device->portSpec, device->readBufferSize);
			return(-1);			
		}		
	}
	
	// Nonblocking mode?
	if (blocking <= 0) {
		// Yes.
		
		// Background read active?
		if (device->readerThread) {
			// Just return what is available in the internal buffer:
			nread = PsychSerialWindowsGlueAsyncReadbufferBytesAvailable(device);
		}
		else {	
			// No, regular sync read: Set COMMTIMEOUT for non-blocking read mode:		
			device->timeouts.ReadIntervalTimeout = MAXDWORD;
			device->timeouts.ReadTotalTimeoutConstant = 0;
			device->timeouts.ReadTotalTimeoutMultiplier = 0;
			
			if (SetCommTimeouts(device->fileDescriptor, &(device->timeouts)) == 0)
			{
				sprintf(errmsg, "Error calling SetCommTimeouts on device %s for non-blocking read - (%d).\n", device->portSpec, GetLastError());
				return(-1);
			}
			
			// Read the data, at most 'amount' bytes, nonblocking:
			if (ReadFile(device->fileDescriptor, device->readBuffer, amount, &nread, NULL) == 0)
			{			
				// Some error:
				sprintf(errmsg, "Error during non-blocking read from device %s - (%d).\n", device->portSpec, GetLastError());
				return(-1);
			}
		}
	}
	else {
		// Nope. This is a blocking read request:

		// Background read active?
		if (device->readerThread) {
			// Poll until requested amount of data is available:
			// TODO FIXME: Could do this in a more friendly way by using condition
			// variables, so we sleep until async reader thread signals availability
			// of sufficient data...
			PsychGetAdjustedPrecisionTimerSeconds(&timeout);
			*timestamp = timeout;
			timeout+=device->readTimeout;
			
			while((*timestamp < timeout) && (PsychSerialWindowsGlueAsyncReadbufferBytesAvailable(device) < (int) amount)) {
				PsychGetAdjustedPrecisionTimerSeconds(timestamp);
				PsychYieldIntervalSeconds(device->pollLatency);
			}
			
			// Return amount of available data:
			nread = PsychSerialWindowsGlueAsyncReadbufferBytesAvailable(device);			
		}
		else {
			// Nope. Setup timeouts for blocking read with user code spec'd timeout parameters:
			// We set all fields to the same value, so reception of a single byte should not take longer than device->readTimeout
			// seconds, neither absolute, nor on average. ReadIntervalTimeout is equivalent to VTIME on Unix, its an inter-byte
			// timer. As this only gets triggered after receipt of first byte, we need to use ReadTotalTimeoutXXX as well, so
			// even if nothing is received, we time out after (readTimeout * amount) time units.
			// All these timeouts are per byte, not absolute for the whole read op, so the real maximum duration for the whole
			// read is the given readTimeout multiplied by the total number of requested bytes.
			// This is conformant to the Unix behaviour where the timeout is also per byte, so it scales with the total
			// number of requested bytes worst case.
			//
			// Granularity of setting is milliseconds, so we need to convert our seconds value to milliseconds, with proper
			// rounding applied:
			device->timeouts.ReadIntervalTimeout = (int) (device->readTimeout * 1000 + 0.5);
			device->timeouts.ReadTotalTimeoutMultiplier = device->timeouts.ReadIntervalTimeout;
			device->timeouts.ReadTotalTimeoutConstant = 0;
			
			if (SetCommTimeouts(device->fileDescriptor, &(device->timeouts)) == 0)
			{
				sprintf(errmsg, "Error calling SetCommTimeouts on device %s for blocking read - (%d).\n", device->portSpec, GetLastError());
				return(-1);
			}

			// Skip this polling for first byte if BlockingBackgroundRead is non-zero:
			if (device->isBlockingBackgroundRead == 0) {			
				if (amount > 1) {
					// Need to poll for arrival of first byte, as the interbyte ReadIntervalTimeout timer is only armed after
					// reception of first byte:
					PsychGetAdjustedPrecisionTimerSeconds(&timeout);
					*timestamp = timeout;
					timeout+=device->readTimeout;
					
					while((*timestamp < timeout) && (PsychIOOSBytesAvailableSerialPort(device) < 1)) {
						PsychGetAdjustedPrecisionTimerSeconds(timestamp);
						// Due to the shoddy Windows scheduler, we can't use a nice sleep duration of 500 microseconds like
						// on Unices, but need to use a rather large 5 msecs so the system actually puts us to sleep
						// and we don't overload the machine in realtime scheduling:
						PsychYieldIntervalSeconds(device->pollLatency);
					}
					
					if (PsychIOOSBytesAvailableSerialPort(device) < 1) {
						// Ok, first byte didn't arrive within one timeout period:
						return(0);
					}
				}
			}

			// Read the data, exactly 'amount' bytes, blocking, unless timeout occurs:
			if (ReadFile(device->fileDescriptor, device->readBuffer, amount, &nread, NULL) == 0)
			{			
				// Some error:
				sprintf(errmsg, "Error during blocking read from device %s - (%d).\n", device->portSpec, GetLastError());
				return(-1);
			}
		} // End Sync blocking read.
	} // End blocking read.
	
	// Read successfully completed if we reach this point. Clear error message:
	errmsg[0] = 0;
	
	// Was this a fetch-op from an active background read?
	if (device->readerThread) {
		// Yes.

		// Check for buffer overflow:
		if (nread > (int) device->readBufferSize) {
			sprintf(errmsg, "Error: Readbuffer overflow for background read operation on device %s. Flushing buffer to recover. At least %i bytes of input data have been lost, expect data corruption!\n", device->portSpec, nread);

			// Flush readBuffer - Try to get a fresh start...
			PsychLockMutex(&(device->readerLock));
			
			// Set read pointer to current write pointer, effectively emptying the buffer:
			device->clientThreadReadPos = device->readerThreadWritePos;
			
			// Unlock data structure:
			PsychUnlockMutex(&(device->readerLock));

			// Return error code:
			return(-1);
		}
		
		// Clamp available amount to requested (or maximum allowable) amount:
		nread = (nread > (int) amount) ? (int) amount : nread;
		
		// Compute startindex in buffer for readout operation:
		raPos = (device->clientThreadReadPos) % (device->readBufferSize);
		
		// Exceeding the end of the buffer?
		if (raPos + nread - 1 >= (int) device->readBufferSize) {
			// Yes: Need an intermediate bounce-buffer to make this safe:
			if (device->bounceBufferSize < nread) {
				free(device->bounceBuffer);
				device->bounceBuffer = NULL;
				device->bounceBufferSize = 0;
			}
			
			if (NULL == device->bounceBuffer) {
				// (Re-)allocate bounceBuffer. Allocate at least 4096 Bytes, ie., one memory page.
				// Anything smaller would be just plain inefficient speed-wise:
				device->bounceBufferSize = (nread < 4096) ? 4096 : nread;
				device->bounceBuffer = (unsigned char*) calloc(sizeof(unsigned char), device->bounceBufferSize);
			}
			
			// Copy data to bounce-buffer, take wraparound in readBuffer into account:
			for(i = 0; i < nread; i++, raPos++) device->bounceBuffer[i] = device->readBuffer[raPos % (device->readBufferSize)];
			*readdata = (void*) device->bounceBuffer;
		}
		else {
			// No: Can directly copy from readBuffer, starting at raPos:
			*readdata = (void*) &(device->readBuffer[raPos]);
		}

		// Retrieve timestamp for this read chunk of data:
		*timestamp = device->timeStamps[(device->clientThreadReadPos / device->readGranularity) % (device->readBufferSize / device->readGranularity)];

		// Update of read-pointer:
		device->clientThreadReadPos += nread;
	}
	else {
		// No. Standard read:
		
		// Take timestamp of read completion:
		PsychGetAdjustedPrecisionTimerSeconds(timestamp);

		// Check for communication errors, acknowledge them, clear the error state and return error message:
		if (PsychIOOSCheckError(device, errmsg)) return(-1);

		// Assign returned data:
		*readdata = (void*) device->readBuffer;
	}

	return(nread);
}

int PsychIOOSCheckError(PsychSerialDeviceRecord* device, char* inerrmsg)
{
	COMSTAT dstatus;
	DWORD   estatus;
	char errmsg[100];
	errmsg[0]=0;
	
	if (ClearCommError(device->fileDescriptor, &estatus, &dstatus)==0) {
		if (verbosity > 0) printf("Error during ClearCommError() on device %s returned (%d)\n", device->portSpec, GetLastError());
		if (inerrmsg) sprintf(inerrmsg, "Error during ClearCommError() on device %s returned (%d)\n", device->portSpec, GetLastError());
		return(0xdeadbeef);
	}
	
	if (estatus > 0 && inerrmsg) {
		if (estatus & CE_BREAK) { sprintf(errmsg, "IOPort: Break condition on receive line detected.\n"); strcat(inerrmsg, errmsg); }
		if (estatus & CE_FRAME) { sprintf(errmsg, "IOPort: Data packet framing error detected.\n"); strcat(inerrmsg, errmsg); }
		if (estatus & CE_OVERRUN) { sprintf(errmsg, "IOPort: Character buffer overrun detected.\n"); strcat(inerrmsg, errmsg); }
		if (estatus & CE_RXOVER) { sprintf(errmsg, "IOPort: Receive buffer overflow detected.\n"); strcat(inerrmsg, errmsg); }
		if (estatus & CE_RXPARITY) { sprintf(errmsg, "IOPort: Parity error detected.\n"); strcat(inerrmsg, errmsg); }
	}
	
	// Store number of pending characters to write out:
	device->num_output_pending = (unsigned int) dstatus.cbOutQue;
	
	return((int) estatus);
}

int PsychIOOSBytesAvailableSerialPort(PsychSerialDeviceRecord* device)
{
	COMSTAT dstatus;
	DWORD   dummy;

	if (device->readerThread) {
		// Async reader poll: Calculate what is available in the internal buffer:
		return(PsychSerialWindowsGlueAsyncReadbufferBytesAvailable(device));
	}
	else {
		if (ClearCommError(device->fileDescriptor, &dummy, &dstatus)==0) {
			if (verbosity > 0) printf("Error during 'BytesAvailable': ClearCommError() bytes available query on device %s returned (%d)\n", device->portSpec, GetLastError());
			return(-1);
		}
	}
	return((int) dstatus.cbInQue);
}

void PsychIOOSFlushSerialPort(PsychSerialDeviceRecord* device)
{	
	if (FlushFileBuffers(device->fileDescriptor)==0) {
		if (verbosity > 0) printf("Error during 'Flush': FlushFileBuffers() on device %s returned (%d)\n", device->portSpec, GetLastError());
	}
	
	return;
}

void PsychIOOSPurgeSerialPort( PsychSerialDeviceRecord* device)
{
	if (PurgeComm(device->fileDescriptor, PURGE_RXABORT | PURGE_RXCLEAR | PURGE_TXABORT | PURGE_TXCLEAR)==0) {
		if (verbosity > 0) printf("Error during 'Purge': PurgeComm() on device %s returned (%d)\n", device->portSpec, GetLastError());
	}
	
	if (device->readerThread) {
		// Purge the input buffer of async reader thread as well:

		// Lock data structure:
		PsychLockMutex(&(device->readerLock));
		
		// Set read pointer to current write pointer, effectively emptying the buffer:
		// It is important to not modify the write pointer, only the read pointer, otherwise
		// we would have an ugly race-condition within the writer thread, as the writer thread
		// does not mutex-lock during reading the writeposition pointer, only during updating
		// the writeposition pointer!
		device->clientThreadReadPos = device->readerThreadWritePos;
		
		// Unlock data structure:
		PsychUnlockMutex(&(device->readerLock));
	}

	return;
}
