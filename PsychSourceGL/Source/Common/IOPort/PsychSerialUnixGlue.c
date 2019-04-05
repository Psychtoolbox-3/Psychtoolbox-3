/*
    Psychtoolbox3/PsychSourceGL/Source/Common/IOPort/PsychSerialUnixGlue.c

    PROJECTS:

        IOPort for now.

    AUTHORS:

        mario.kleiner.de@gmail.com      mk

    PLATFORMS:

        All Unix (aka Posix) systems: OS/X and Linux.

    HISTORY:

        04/10/2008      mk      Initial implementation.
        04/03/2011      mk      Audited (but not tested) to be 64-bit safe. Doesn't take advantage of > 2 GB
                                buffers yet. This will require some minor changes.

    DESCRIPTION:

        This is the operating system dependent "glue code" layer for access to serial ports for the
        Unices, ie., for GNU/Linux and Apple MacOS/X. It is used by the higher-level serial port
        routines to abstract out operating system dependencies.

        The code is shared and #ifdef'ed if needed, because most of it is identical for Linux and OS/X.
*/

#include "Psych.h"

// Need to exclude this files source-code from compile on Win32:
#if PSYCH_SYSTEM != PSYCH_WINDOWS

#include "IOPort.h"

#if PSYCH_SYSTEM == PSYCH_LINUX
// Defines serial_struct and ASYNC_LOW_LATENCY
#include <linux/serial.h>
#endif

// Externally defined level of verbosity:
extern int verbosity;

// Map numeric baud rate to Posix constant:
static int BaudToConstant(int inint)
{
    // Need to map from numeric values to predefined constants. On OS/X the
    // constants are identical to the input values, but on Linux they aren't, so better safe
    // than sorry, although it looks as if Linux accepts the "raw" numbers properly as well...
    switch(inint) {
        case 50: inint = B50; break;
        case 75: inint = B75; break;
        case 110: inint = B110; break;
        case 134: inint = B134; break;
        case 150: inint = B150; break;
        case 200: inint = B200; break;
        case 300: inint = B300; break;
        case 600: inint = B600; break;
        case 1200: inint = B1200; break;
        case 1800: inint = B1800; break;
        case 2400: inint = B2400; break;
        case 4800: inint = B4800; break;
        case 9600: inint = B9600; break;
        case 19200: inint = B19200; break;
        case 38400: inint = B38400; break;
        case 57600: inint = B57600; break;
        case 115200: inint = B115200; break;
        case 230400: inint = B230400; break;
        #if PSYCH_SYSTEM == PSYCH_OSX
        // Only defined on OS/X:
        case 7200: inint = B7200; break;
        case 14400: inint = B14400; break;
        case 28800: inint = B28800; break;
        case 76800: inint = B76800; break;
        #endif

        default:
            // Pass inint unchanged in the hope that the OS will accept it:
            if (verbosity > 1) printf("IOPort: Non-standard BaudRate %i specified. We'll see if the OS likes this. Double-check the settings!\n", inint);
    }

    return(inint);
}

// Map Posix constant to numeric baud rate:
static int ConstantToBaud(int inint)
{
    // Need to map predefined constants to numeric values On OS/X the
    // constants are identical to the input values, but on Linux they aren't, so better safe
    // than sorry, although it looks as if Linux accepts the "raw" numbers properly as well...
    switch(inint) {
        case B50: inint = 50; break;
        case B75: inint = 75; break;
        case B110: inint = 110; break;
        case B134: inint = 134; break;
        case B150: inint = 150; break;
        case B200: inint = 200; break;
        case B300: inint = 300; break;
        case B600: inint = 600; break;
        case B1200: inint = 1200; break;
        case B1800: inint = 1800; break;
        case B2400: inint = 2400; break;
        case B4800: inint = 4800; break;
        case B9600: inint = 9600; break;
        case B19200: inint = 19200; break;
        case B38400: inint = 38400; break;
        case B57600: inint = 57600; break;
        case B115200: inint = 115200; break;
        case B230400: inint = 230400; break;
        #if PSYCH_SYSTEM == PSYCH_OSX
        // Only defined on OS/X:
        case B7200: inint = 7200; break;
        case B14400: inint = 14400; break;
        case B28800: inint = 28800; break;
        case B76800: inint = 76800; break;
        #endif
        default:
            // Pass inint unchanged in the hope that the OS will accept it:
            if (verbosity > 1) printf("IOPort: Non-standard BaudRate %i detected. Let's see if this makes sense...\n", inint);
    }

    return(inint);
}

int PsychSerialUnixGlueFcntl(PsychSerialDeviceRecord* device, int newVal)
{
    int rc;

    // Query current setting:
    rc = fcntl(device->fileDescriptor, F_GETFL, 0);

    // Only perform fcntl() call if current setting isn't already equal to requested new setting newVal,
    // because fcntl() calls can be very time-consuming on some devices, e.g., USB-Serial ports:
    if (rc != newVal) { rc = fcntl(device->fileDescriptor, F_SETFL, newVal); }

    // Return new setting or error code:
    return(rc);
}

int PsychSerialUnixGlueSetBlockingMinBytes(PsychSerialDeviceRecord* device, int minBytes)
{
    struct termios    options;
    int rc;

    // Retrieve current termios settings:
    if (tcgetattr(device->fileDescriptor, &options) == -1)
    {
        if (verbosity > 0) fprintf(stderr, "IOPort: Error getting current serial port device settings for device %s - %s(%d).\n", device->portSpec, strerror(errno), errno);
        return(-1);
    }

    rc = minBytes;

    // Clamp amount smaller or equal one to zero: Clamping negative values to zero is obviously needed
    // to prevent passing illegal values (valid range 0-255). The amount == 1 --> amount = 0 case is special.
    // Normally the interbyte timer VTIME will only get armed after reception of 1st byte, so asking for
    // VMIN=1 for amount == 1 would effectively disable the VTIME interbyte timeout semantic that we want.
    // However the special combination VTIME = someinterbytetimeout and VTIME = 0 will (see man termios)
    // ask the system to wait until exactly 1 Byte is available or VTIME time elapses (absolute time, no
    // interbyte timeout anymore!), whatever comes first. The net result is exactly the interbyte VTIME timeout
    // semantics we want for a single byte receive.
    if (minBytes <= 1) {
        rc = (minBytes < 1) ? 0 : 1;
        minBytes = 0;
    }

    // Clamp at max. 255 Bytes, the maximum for blocking reads: This is a bit inconvenient. Unix limits
    // the maximum amount of Bytes one can wait for in a blocking receive to 255, because the VMIN
    // parameter is only 1 Byte in size. We have to clamp to 255 if 255 is exceeded and just hope that
    // this doesn't lead to undefined behaviour in usercode. Luckily, must data transfers in typical
    // cognitive science applications of serial ports are smaller than 255 Bytes, e.g., 1 Byte for
    // shoddy response boxes or for TTL trigger receivers, 1 - 8 Bytes for good response boxes, and
    // maybe anywhere between 4 - 20 Bytes for eyetracker stuff and similar equipment. So i guess usercode
    // would seldomly hit this special case. We output a warning anyway if a high level of verbosity is
    // selected.
    if (minBytes > 255) {
        minBytes = 255;
        rc = 255;
    }

    // Change the number of bytes that must be received for blocking read()'s before read() returns.
    // read() will return if it can fetch at least that many bytes (we want 'minBytes' bytes at least),
    // or if the user specified timeout occurs (as set in Open or Configure routine):
    // N.B.: We do this lazily, only executing the call if the current and requested setting differ,
    // as the tcsetattr() call may be expensive on some system/driver combinations...
    if (options.c_cc[VMIN] != (unsigned char) minBytes) {
        options.c_cc[VMIN] = (unsigned char) minBytes;

        // Cause the new options to take effect immediately.
        if (tcsetattr(device->fileDescriptor, TCSANOW, &options) == -1)
        {
            if (verbosity > 0) fprintf(stderr, "IOPort: Error setting new serial port configuration attributes for device %s - %s(%d).\n", device->portSpec, strerror(errno), errno);
            return(-1);
        }
    }

    // Return really set value:
    return(rc);
}

int PsychSerialUnixGlueAsyncReadbufferBytesAvailable(PsychSerialDeviceRecord* device)
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

void* PsychSerialUnixGlueReaderThreadMain( void* deviceToCast)
{
    int rc, nread, oldstate;
    int navail;
    int tmpcurpos, naccumread;
    unsigned char lastcharacter, lineterminator;
    double dt, oldt;
    double t;
    int doBlockingRead = 0;

    // Get a handle to our device struct: These pointers must not be NULL!!!
    PsychSerialDeviceRecord* device = (PsychSerialDeviceRecord*) deviceToCast;

    // Assign a name to ourselves, for debugging:
    PsychSetThreadName("IOPortSerialRd");

    // Try to raise our priority: We ask to switch ourselves (NULL) to priority class 2 aka
    // realtime scheduling, with a tweakPriority of +1, ie., raise the relative
    // priority level by +1 wrt. to the current level:
    if ((rc = PsychSetThreadPriority(NULL, 2, 1)) > 0) {
        if (verbosity > 0) printf("PTB-ERROR: In IOPort:PsychSerialUnixGlueReaderThreadMain(): Failed to switch to realtime priority [%s]!\n", strerror(rc));
    }

    // Init reference timestamp of last scan:
    PsychGetAdjustedPrecisionTimerSeconds(&oldt);
    dt = 0;

    // Main loop: Runs until external thread cancellation:
    while (1) {
        // Test for explicit cancellation by mother-thread:
        PsychTestCancelThread(&(device->readerThread));

        // Polling read requested?
        doBlockingRead = device->isBlockingBackgroundRead;

        if (doBlockingRead == 0) {
            // Polling operation:

            // Enough data available for read of requested granularity?
            // If not, we sleep for some time, then retry:
            ioctl(device->fileDescriptor, FIONREAD, &navail);

            while(navail < device->readGranularity) {
                PsychTestCancelThread(&(device->readerThread));
                PsychWaitIntervalSeconds(device->pollLatency);
                ioctl(device->fileDescriptor, FIONREAD, &navail);
            }
        }
        else {
            // Non-polling operation. We perform a blocking read on the device.
            // Assuming the masterthread doesn't perform any non-blocking operations,
            // which would screw us or it, this is safe. The read() call is a thread cancellation
            // point, so we can be safely aborted by the masterthread even if blocked on
            // input.

            // Set filedescriptor to blocking mode:
            // Clear the O_NONBLOCK flag so subsequent I/O will block.
            // See fcntl(2) ("man 2 fcntl") for details.
            // We do this in a cached, lazy way b'cause it can be expensive:
            PsychSerialUnixGlueFcntl(device, 0);
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

            // Setup minimum byte counter for 1 Byte blocking reads:
            if (doBlockingRead > 0) PsychSerialUnixGlueSetBlockingMinBytes(device, 1);

            // Repeat reading until a maximum of readGranularity bytes has been read or until
            // lineterminator character is encountered, whatever comes first:
            while ( (naccumread < device->readGranularity) && (lastcharacter != lineterminator) ) {
                // Read 1 Byte:
                if ((nread = read(device->fileDescriptor, &lastcharacter, 1)) > 0) {
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
                    if (nread == -1 && errno == EAGAIN) {
                        // Timeout.
                        if (verbosity > 5) fprintf(stderr, "PTB-ERROR: In IOPort:PsychSerialUnixGlueReaderThreadMain(): Linebuffered read: Failed to read 1 byte of data due to read-timeout at relative position %i! Padding...\n", naccumread);
                    }
                    else {
                        if (verbosity > 5) fprintf(stderr, "PTB-ERROR: In IOPort:PsychSerialUnixGlueReaderThreadMain(): Linebuffered read: Failed to read 1 byte of data for reason [%s] at relative position %i! Padding...\n", strerror(errno), naccumread);
                    }

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

            // Setup minimum byte counter for 'naccumread' Bytes blocking reads:
            if (doBlockingRead > 0) {
                PsychSerialUnixGlueSetBlockingMinBytes(device, naccumread);

                // More than 255 Bytes requested for one read() call?
                if (naccumread > 255) {
                    // The Unix blocking VMIN/VTIME mechanism can only handle blocking/waiting
                    // for at most 255 bytes, ie., it can't handle this request. We need to resort
                    // to a polling wait for naccumread bytes :-(

                    // Enough data available for read of requested granularity?
                    // If not, we sleep for some time, then retry:
                    ioctl(device->fileDescriptor, FIONREAD, &navail);

                    while(navail < naccumread) {
                        PsychTestCancelThread(&(device->readerThread));
                        PsychWaitIntervalSeconds(device->pollLatency);
                        ioctl(device->fileDescriptor, FIONREAD, &navail);
                    }

                    // Ok, we've got our share of bytes...
                }
            }

            // Enough data available. Read it!
            if ((nread = read(device->fileDescriptor, &(device->readBuffer[(device->readerThreadWritePos) % (device->readBufferSize)]), naccumread)) != naccumread) {
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

                if (verbosity > 5) fprintf(stderr, "PTB-ERROR: In IOPort:PsychSerialUnixGlueReaderThreadMain(): Failed to read %i bytes of data for unknown reason (Got only %i bytes)! Padding...\n", naccumread, nread);
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
        }    // End of regular non-linebuffered readop.

        // Prevent our cancellation:
        pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &oldstate);

        // Store timestamp for this read chunk of data:
        device->timeStamps[(device->readerThreadWritePos / device->readGranularity) % (device->readBufferSize / device->readGranularity)] = t;

        // Try to lock, block until available if not available:
        if ((rc=PsychLockMutex(&(device->readerLock)))) {
            // This could potentially kill Matlab, as we're printing from outside the main interpreter thread.
            // Use fprintf() instead of the overloaded printf() (aka mexPrintf()) in the hope that we don't
            // wreak havoc -- maybe it goes to the system log, which should be safer...
            fprintf(stderr, "PTB-ERROR: In IOPort:PsychSerialUnixGlueReaderThreadMain(): mutex_lock failed  [%s].\n", strerror(rc));

            // Commit suicide:
            return(NULL);
        }

        // Update linear write pointer:
        device->readerThreadWritePos += device->readGranularity;

        // Need to unlock the mutex:
        if ((rc=PsychUnlockMutex(&(device->readerLock)))) {
            // This could potentially kill Matlab, as we're printing from outside the main interpreter thread.
            // Use fprintf() instead of the overloaded printf() (aka mexPrintf()) in the hope that we don't
            // wreak havoc -- maybe it goes to the system log, which should be safer...
            fprintf(stderr, "PTB-ERROR: In IOPort:PsychSerialUnixGlueReaderThreadMain(): Last mutex_unlock in termination failed  [%s].\n", strerror(rc));

            // Commit suicide:
            return(NULL);
        }

        // Reenable cancellation:
        pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, &oldstate);

        // Next iteration...
    }

    // Go and die peacefully...
    return(NULL);
}

void PsychIOOSShutdownSerialReaderThread(PsychSerialDeviceRecord* device)
{
    if (device->readerThread) {
        // Cancel the thread:
        PsychAbortThread(&(device->readerThread));

        // Wait for it to die:
        PsychDeleteThread(&(device->readerThread));

        // Mark it as dead:
        device->readerThread = (psych_thread) NULL;

        // Release the mutex:
        PsychDestroyMutex(&(device->readerLock));

        // Release timestamp buffer:
        free(device->timeStamps);
        device->timeStamps = NULL;
    }

    return;
}

/* PsychSerialUnixGlueJLTriggerThreadMain() -- Timed trigger byte emission.
 * This is a proof-of-concept prototype for the async-task API of IOPort,
 * not for public use by regular users!
 *
 * Emits a 0xff triggerbyte over serial port at time device->triggerWhen.
 *
 */
void* PsychSerialUnixGlueJLTriggerThreadMain( void* deviceToCast)
{
    int rc;
    unsigned char writedata = 0xff;
    double timestamp[4];
    char errmsg[256];
    errmsg[0] = 0;

    // Get a handle to our device struct: These pointers must not be NULL!!!
    PsychSerialDeviceRecord* device = (PsychSerialDeviceRecord*) deviceToCast;

    // Try to raise our priority: We ask to switch ourselves (NULL) to priority class 2 aka
    // realtime scheduling, with a tweakPriority of +2, ie., raise the relative
    // priority level by +2 wrt. to the current level:
    if ((rc = PsychSetThreadPriority(NULL, 2, 2)) > 0) {
        if (verbosity > 0) printf("PTB-ERROR: In IOPort:PsychSerialUnixGlueJLTriggerThreadMain(): Failed to switch to realtime priority [%s]!\n", strerror(rc));
    }

    // Detach ourselves, so parent thread doesn't need to pthread_join() us:
    if (pthread_detach(pthread_self())) printf("PTB-ERROR: In IOPort:PsychSerialUnixGlueJLTriggerThreadMain(): Failed to detach!\n");

    // Wait for our target time to come...
    PsychWaitUntilSeconds(device->triggerWhen);

    if (1 != PsychIOOSWriteSerialPort(device, &writedata, 1, 1, errmsg, timestamp)) {
        printf("PTB-ERROR: In IOPort:PsychSerialUnixGlueJLTriggerThreadMain(): Failed to write triggerbyte!\n");
    }

    // Good enough?
    if ((verbosity > 3) && (timestamp[0] - device->triggerWhen > 0.003)) printf("PTB-WARNING: In IOPort:PsychSerialUnixGlueJLTriggerThreadMain(): Trigger emission delayed by up to %f msecs wrt. to deadline!\n", (float) 1000.0 * (timestamp[0] - device->triggerWhen));

    // Store write completion timestamp:
    device->triggerWhen = timestamp[0];

    // Go and die peacefully...
    return(NULL);
}

/* PsychIOOSOpenSerialPort()
 *
 * Open a serial port device and configure it.
 *
 * portSpec - String with the port name / device name of the serial port device.
 * configString - String with port configuration parameters.
 * errmsg - Pointer to char[] buffer in which error messages should be returned, if any.
 * On success, allocate a PsychSerialDeviceRecord with all relevant settings,
 * return a pointer to it.
 *
 * Otherwise abort with error message.
 */
PsychSerialDeviceRecord* PsychIOOSOpenSerialPort(const char* portSpec, const char* configString, char* errmsg)
{
    int fileDescriptor = -1;
    struct termios options;
    PsychSerialDeviceRecord* device = NULL;
    psych_bool usererr = FALSE;

    // Init errmsg error message to empty == no error:
    errmsg[0] = 0;

    // Open the serial port read/write, with no controlling terminal, and don't wait for a connection.
    // The O_NONBLOCK flag also causes subsequent I/O on the device to be non-blocking.
    // See open(2) ("man 2 open") for details.
    fileDescriptor = open(portSpec, O_RDWR | O_NOCTTY | O_NONBLOCK);
    if (fileDescriptor == -1) {
        if (errno == ENOENT) {
            sprintf(errmsg, "Error opening serial port device %s - No such serial port device exists! (%d) [ENOENT].\n", portSpec, errno);
            usererr = TRUE;
        }
        else if (errno == EBUSY || errno == EPERM) {
            sprintf(errmsg, "Error opening serial port device %s - The serial port is already open, close it first! %s(%d) [EBUSY EPERM]. Could be a permission problem as well.\n", portSpec, strerror(errno), errno);
            usererr = TRUE;
        }
        else {
            sprintf(errmsg, "Error opening serial port device %s - %s(%d).\n", portSpec, strerror(errno), errno);
        }
        goto error;
    }

    // Note that open() follows POSIX semantics: multiple open() calls to the same file will succeed
    // unless the TIOCEXCL ioctl is issued. This will prevent additional opens except by root-owned
    // processes.
    // See tty(4) ("man 4 tty") and ioctl(2) ("man 2 ioctl") for details.
    if (ioctl(fileDescriptor, TIOCEXCL) == -1) {
        if (verbosity > 1) printf("IOPort-Warning: Error setting exclusive access (TIOCEXCL) on device %s - %s(%d).\n", portSpec, strerror(errno), errno);
        // goto error;
    }

    // Create the device struct and init it:
    device = calloc(1, sizeof(PsychSerialDeviceRecord));
    device->fileDescriptor = -1;
    device->readBuffer = NULL;
    device->readBufferSize = 0;
    device->readerThread = (psych_thread) NULL;
    device->lineTerminator = _POSIX_VDISABLE;

    // Get the current options and save them so we can restore the default settings later.
    if (tcgetattr(fileDescriptor, &(device->OriginalTTYAttrs)) == -1) {
        sprintf(errmsg, "Error getting original device settings for device %s - %s(%d).\n", portSpec, strerror(errno), errno);
        goto error;
    }

    strncpy(device->portSpec, portSpec, sizeof(device->portSpec));
    device->fileDescriptor = fileDescriptor;
    options = device->OriginalTTYAttrs;

    // Ok, the device is basically open. Call the reconfiguration routine with the setting
    // string. It will do all further setup work:
    if (PsychError_none != PsychIOOSConfigureSerialPort(device, configString)) {
        // Something failed. Error handling...
        sprintf(errmsg, "Error changing device settings for device %s - %s(%d).\n", portSpec, strerror(errno), errno);
        usererr = TRUE;
        goto error;
    }

    if (device->readBuffer == NULL) {
        sprintf(errmsg, "Error for device %s - No InputBuffer allocated! You must specify the 'InputBuffer' size in the configuration.\n", portSpec);
        usererr = TRUE;
        goto error;
    }

    // printf("MAX_INPUT = %i / %i \n", (int) fpathconf(device->fileDescriptor, _PC_MAX_INPUT), (int) fpathconf(device->fileDescriptor, _PC_MAX_CANON));

    // Success! Return Pointer to new device structure:
    return(device);

    // Failure path: Called on error with error message in errmsg:
error:

    if (fileDescriptor != -1) {
        // File open.
        if (device) {
            // Settings retrieved or made?
            if (device->fileDescriptor != -1) {
                // Have original settings available: Try to restore them:
                if (tcsetattr(fileDescriptor, TCSANOW, &options) == -1) {
                    if (verbosity > 1) printf("WARNING: In error handling: Error restoring tty attributes %s - %s(%d).\n", portSpec, strerror(errno), errno);
                }
            }

            // Release read buffer, if any:
            if (device->readBuffer) free(device->readBuffer);

            // Release device struct:
            free(device);
        }

        // Close file:
        close(fileDescriptor);
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
    if ((!device->dontFlushOnWrite) && (tcdrain(device->fileDescriptor) == -1)) {
        if (verbosity > 1) printf("IOPort: WARNING: While trying to close serial port: Error waiting for drain - %s(%d).\n", strerror(errno), errno);
    }

    // Traditionally it is good practice to reset a serial port back to
    // the state in which you found it. This is why the original termios struct
    // was saved.
    if (tcsetattr(device->fileDescriptor, TCSANOW, &(device->OriginalTTYAttrs)) == -1) {
        if (verbosity > 1) printf("IOPort: WARNING: While trying to close serial port: Could not restore original port settings - %s(%d).\n", strerror(errno), errno);
    }

    // Close device:
    close(device->fileDescriptor);

    // Release read buffer and/or bounceBuffer, if any:
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
PsychError PsychIOOSConfigureSerialPort(PsychSerialDeviceRecord* device, const char* configString)
{
#if PSYCH_SYSTEM == PSYCH_LINUX
    struct serial_struct serialstruct;
#endif
    int rc;
    struct termios options;
    int handshake;
    char* p;
    float infloat;
    int inint;
    unsigned long mics = 0UL;
    psych_bool updatetermios = FALSE;

    // The serial port attributes such as timeouts and baud rate are set by modifying the termios
    // structure and then calling tcsetattr() to cause the changes to take effect. Note that the
    // changes will not become effective without the tcsetattr() call.
    // See tcsetattr(4) ("man 4 tcsetattr") for details.

    // Retrieve current settings:
    if (tcgetattr(device->fileDescriptor, &options) == -1) {
        if (verbosity > 0) printf("Error getting current serial port device settings for device %s - %s(%d).\n", device->portSpec, strerror(errno), errno);
        return(PsychError_system);
    }

    // Print the current input and output baud rates.
    // See tcsetattr(4) ("man 4 tcsetattr") for details.
    if (verbosity > 4) {
        printf("IOPort-Info: Configuration for device %s:\n", device->portSpec);
        printf("IOPort-Info: Current input baud rate is %d\n", ConstantToBaud((int) cfgetispeed(&options)));
        printf("IOPort-Info: Current output baud rate is %d\n", ConstantToBaud((int) cfgetospeed(&options)));
    }

    // Set raw input (non-canonical) mode: No line-based processing or editing of
    // received data via special control characters, unless a processing mode of
    // cooked was explicitely requested:
    if ((p=strstr(configString, "ProcessingMode="))!=NULL) {
        if (strstr(p, "ProcessingMode=Cooked") == NULL) {
            // Disable cooked mode, switch port to raw ops:
            cfmakeraw(&options);
            updatetermios = TRUE;
            device->cookedMode = 0;
        }
        else {
            // Enable cooked mode:
            device->cookedMode = 1;
            updatetermios = TRUE;

            // Disable any kind of output processing in cooked mode.
            // We only use the input processing facilities of cooked mode!
            options.c_oflag &= ~OPOST;

            // If no lineTerminator is set, assign it to be the Unix default of ASCII 10 aka LF aka NL:
            if (device->lineTerminator == _POSIX_VDISABLE) device->lineTerminator = (unsigned char) 10;
        }
    }

    // No minimum number of bytes to receive by default. This will get overriden
    // in the synchronous and asynchronous (background thread driven) functions
    // prior to read() system calls if blocking reads are requested:
    options.c_cc[VMIN] = 0;

    if ((p = strstr(configString, "ReceiveTimeout="))) {
        // Set timeouts for read operations:
        // See tcsetattr(4) ("man 4 tcsetattr") and termios(4) ("man 4 termios") for details.
        if ((1!=sscanf(p, "ReceiveTimeout=%f", &infloat)) || (infloat < 0)) {
            if (verbosity > 0) printf("Invalid parameter for ReceiveTimeout set! Typo, or negative value provided.\n");
            return(PsychError_user);
        }
        else {
            // Set timeout: It is in granularity of 1/10th seconds, so we need to quantize to 10th of seconds:
            inint = (int) (infloat * 10 + 0.5);
            // Clamp to a minimum of 0.1 secs if a non-zero timeout is requested:
            if ((inint <= 0) && (infloat > 0)) {
                if (verbosity > 1) printf("IOPort: Warning: Requested per-byte 'ReceiveTimeout' value %f secs is positive but smaller than supported minimum of 0.1 secs on Unix. Changed value to 0.1 secs.\n", infloat);
                infloat = 0.1;
                inint = 1;
            }

            // Clamp to maximum of 25.5 seconds aka 255 10th of seconds. As Unix defines
            // VTIME as 1 Byte, we can't go higher than that.
            if (inint > 255) {
                if (verbosity > 1) printf("IOPort: Warning: Requested per-byte 'ReceiveTimeout' value %f secs is bigger than supported maximum of 25.5 seconds on Unix. Changed value to 25.5 secs.\n", infloat);
                infloat = 25.5;
                inint = 255;
            }

            // Make sure we don't change non-mutex-protected variables behind the
            // back of our readerThread by only allowing this function to be called
            // with inactive thread:
            if (device->readerThread) {
                if (verbosity > 0) printf("Assigned ReceiveTimeout= while background read operations are already enabled! Disable first via 'StopBackgroundRead'!\n");
                return(PsychError_user);
            }

            options.c_cc[VTIME] = inint;
            device->readTimeout = infloat;
            updatetermios = TRUE;
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

        // We also assign the line terminator as VEOL, which will be effective in
        // canonical input processing mode, and - maybe - as event character on
        // OS/X et al. for FTDI serial converter chips: A negative value for
        // LineTerminator disables the feature:
        options.c_cc[VEOL] = (inint >= 0) ? (unsigned char) inint : _POSIX_VDISABLE;
        // printf("ISNOW %x \n", options.c_cc[VEOL]);
        updatetermios = TRUE;
    }

    // Assign receiver enable state, if any:
    if ((p = strstr(configString, "ReceiverEnable="))) {
        if (1!=sscanf(p, "ReceiverEnable=%i", &inint)) {
            if (verbosity > 0) printf("Invalid parameter for ReceiverEnable= set!\n");
            return(PsychError_invalidIntegerArg);
        }

        if (inint > 0) {
            // Enable receiver:
            options.c_cflag |= CREAD;
        }
        else {
            // Disable receiver:
            options.c_cflag &= ~CREAD;
        }

        updatetermios = TRUE;
    }

    // Set common baud rate for send and receive:
    if ((p = strstr(configString, "BaudRate="))) {
        if (1!=sscanf(p, "BaudRate=%i", &inint)) {
            if (verbosity > 0) printf("Invalid parameter for BaudRate set!\n");
            return(PsychError_invalidIntegerArg);
        }
        else {
            // Set common speed for input- and output queues:
            inint = BaudToConstant(inint);
            if (-1==cfsetspeed(&options, inint)) {
                if (verbosity > 0) printf("Invalid BaudRate %i not accepted! (%s)", inint, strerror(errno));
                return(PsychError_invalidIntegerArg);
            }
            updatetermios = TRUE;
        }
    }

    // Generation and handling of parity bits:
    if ((p = strstr(configString, "Parity="))) {
        // Clear all parity settings:
        options.c_cflag &= ~(PARENB | PARODD);
        options.c_iflag &= ~INPCK;

        if (strstr(p, "Parity=Even")) {
            options.c_cflag |= PARENB;
            options.c_iflag |= INPCK;
        }
        else
        if (strstr(p, "Parity=Odd")) {
            options.c_cflag |= (PARENB | PARODD);
            options.c_iflag |= INPCK;
        }
        else
        if (strstr(p, "Parity=None")) {
            // Nothing to do.
        }
        else {
            // Invalid parity spec:
            if (verbosity > 0) printf("Invalid parity setting %s not accepted! (Valid are None, Even and Odd", p);
            return(PsychError_user);
        }
        updatetermios = TRUE;
    }

    // Handling of Break conditions:
    if ((p = strstr(configString, "BreakBehaviour="))) {
        // Clear all break settings:
        options.c_iflag &= ~(IGNBRK | BRKINT);

        if (strstr(p, "BreakBehaviour=Ignore")) {
            // Break signals will be ignored:
            options.c_iflag |= IGNBRK;
        }
        else
        if (strstr(p, "BreakBehaviour=Flush")) {
            // A break signal will flush the input- and output-queues:
            options.c_iflag |= BRKINT;
        }
        else
        if (strstr(p, "BreakBehaviour=Zero")) {
            // Don't change anything, ie. IGNBRK and BRKINT not set:
            // A break signal will be output as single zero byte.
        }
        else {
            // Invalid spec:
            if (verbosity > 0) printf("Invalid break behaviour %s not accepted (Valid: Ignore, Flush, or Zero)!", p);
            return(PsychError_user);
        }
        updatetermios = TRUE;
    }

    // Handling of data bits:
    if ((p = strstr(configString, "DataBits="))) {
        // Clear all databit settings:
        options.c_cflag &= ~(CS5 | CS6 | CS7 | CS8);

        if (strstr(p, "DataBits=5")) {
            options.c_cflag |= CS5;
        }
        else
        if (strstr(p, "DataBits=6")) {
            options.c_cflag |= CS6;
        }
        else
        if (strstr(p, "DataBits=7")) {
            options.c_cflag |= CS7;
        }
        else
        if (strstr(p, "DataBits=8")) {
            options.c_cflag |= CS8;
        }
        else {
            // Invalid spec:
            if (verbosity > 0) printf("Invalid setting for data bits %s not accepted! (Valid: 5, 6, 7 or 8 databits)", p);
            return(PsychError_user);
        }
        updatetermios = TRUE;
    }

    // Handling of stop bits:
    if ((p = strstr(configString, "StopBits="))) {
        // Clear all stopbit settings:
        options.c_cflag &= ~CSTOPB;

        if (strstr(p, "StopBits=1")) {
            // Nothing to do, this is already set above.
        }
        else
        if (strstr(p, "StopBits=2")) {
            options.c_cflag |= CSTOPB;
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
        options.c_cflag &= ~CRTSCTS;
        options.c_iflag &= ~(IXON | IXOFF);

        if (strstr(p, "FlowControl=None")) {
            // Set above already...
        }
        else
        if (strstr(p, "FlowControl=Software")) {
            options.c_iflag |= (IXON | IXOFF);
        }
        else
        if (strstr(p, "FlowControl=Hardware")) {
            options.c_cflag |= CRTSCTS;
        }
        else {
            // Invalid spec:
            if (verbosity > 0) printf("Invalid setting for flow control %s not accepted! (Valid: None, Software, Hardware)", p);
            return(PsychError_user);
        }
        updatetermios = TRUE;
    }

    // Print the new input and output baud rates. Note that the IOSSIOSPEED ioctl interacts with the serial driver
    // directly bypassing the termios struct. This means that the following two calls will not be able to read
    // the current baud rate if the IOSSIOSPEED ioctl was used but will instead return the speed set by the last call
    // to cfsetspeed.

    // Retrieve baudrate and check for equal rate on input- and output queue:
    if (verbosity > 1) {
        inint = cfgetispeed(&options);
        if (inint != (int) cfgetospeed(&options)) printf("IOPort: Warning: Hmm, new input- and output baudrates %i vs. %i don't match!? May or may not be a problem...\n", ConstantToBaud((int) inint), ConstantToBaud((int) cfgetospeed(&options)));
    }

    // Output new baud rate:
    if (verbosity > 4) printf("IOPort-Info: Baud rate changed to %d\n", (int) ConstantToBaud(inint));

    // Cause the new options to take effect immediately.
    if (updatetermios && (tcsetattr(device->fileDescriptor, TCSANOW, &options) == -1)) {
        if (verbosity > 0) printf("Error setting new serial port configuration attributes for device %s - %s(%d).\n", device->portSpec, strerror(errno), errno);
        return(PsychError_system);
    }

    // To set the modem handshake lines, use the following ioctls.
    // See tty(4) ("man 4 tty") and ioctl(2) ("man 2 ioctl") for details.

    // Retrieve current handshake state aka DTR, RTS, CTS, DSR et al.:
    handshake = 0;
    if (ioctl(device->fileDescriptor, TIOCMGET, &handshake) == -1) {
        if (verbosity > 0) printf("Error getting lines status for device %s - %s(%d).\n", device->portSpec, strerror(errno), errno);
        if (strstr(configString, "Lenient") == NULL) return(PsychError_system);
    }
    else if (verbosity > 4) {
        printf("IOPort-Info: Handshake lines currently set to %d : ", handshake);
        printf("DTR=%i : DSR=%i : RTS=%i : CTS=%i\n", (handshake & TIOCM_DTR) ? 1:0, (handshake & TIOCM_DSR) ? 1:0, (handshake & TIOCM_RTS) ? 1:0, (handshake & TIOCM_CTS) ? 1:0);
    }

    // Handling of DTR :
    if ((p = strstr(configString, "DTR="))) {
        if (strstr(p, "DTR=0")) {
            handshake&= ~TIOCM_DTR;
        }
        else
        if (strstr(p, "DTR=1")) {
            handshake|= TIOCM_DTR;
        }
        else {
            // Invalid spec:
            if (verbosity > 0) printf("Invalid setting for DTR %s not accepted! (Valid: 1 or 0)", p);
            return(PsychError_user);
        }
    }

    // Handling of DSR :
    if ((p = strstr(configString, "DSR="))) {
        if (strstr(p, "DSR=0")) {
            handshake&= ~TIOCM_DSR;
        }
        else
        if (strstr(p, "DSR=1")) {
            handshake|= TIOCM_DSR;
        }
        else {
            // Invalid spec:
            if (verbosity > 0) printf("Invalid setting for DSR %s not accepted! (Valid: 1 or 0)", p);
            return(PsychError_user);
        }
    }

    // Handling of RTS :
    if ((p = strstr(configString, "RTS="))) {
        if (strstr(p, "RTS=0")) {
            handshake&= ~TIOCM_RTS;
        }
        else
        if (strstr(p, "RTS=1")) {
            handshake|= TIOCM_RTS;
        }
        else {
            // Invalid spec:
            if (verbosity > 0) printf("Invalid setting for RTS %s not accepted! (Valid: 1 or 0)", p);
            return(PsychError_user);
        }
    }

    // Handling of CTS :
    if ((p = strstr(configString, "CTS="))) {
        if (strstr(p, "CTS=0")) {
            handshake&= ~TIOCM_CTS;
        }
        else
        if (strstr(p, "CTS=1")) {
            handshake|= TIOCM_CTS;
        }
        else {
            // Invalid spec:
            if (verbosity > 0) printf("Invalid setting for CTS %s not accepted! (Valid: 1 or 0)", p);
            return(PsychError_user);
        }
    }

    if (ioctl(device->fileDescriptor, TIOCMSET, &handshake) == -1) {
        if (verbosity > 0) printf("Error setting handshake lines status for device %s - %s(%d).\n", device->portSpec, strerror(errno), errno);
        if (strstr(configString, "Lenient") == NULL) return(PsychError_system);
    }

#if PSYCH_SYSTEM == PSYCH_OSX
    // OS/X specific implementation:
    if ((p = strstr(configString, "ReceiveLatency="))) {
        // Set receive latency for low-level driver: Minimum latency between dequeue and notify
        // operation on reception of new input data. The granularity is microseconds...
        if ((1!=sscanf(p, "ReceiveLatency=%f", &infloat)) || (infloat < 0)) {
            if (verbosity > 0) printf("Invalid parameter for ReceiveLatency set! Typo, or negative value provided.\n");
            return(PsychError_user);
        }
        else {
            // Set timeout: It is in granularity microseconds, so we need to quantize to microseconds:
            mics = (unsigned long) (infloat * 1000000 + 0.5);

            // A value of zero would revert to operating system default, so we need to make sure
            // the mics value is at least 1 microsecond:
            if (mics == 0) mics = 1;

            // Set the receive latency in microseconds. Serial drivers use this value to determine how often to
            // dequeue characters received by the hardware. Most applications don't need to set this value: if an
            // app reads lines of characters, the app can't do anything until the line termination character has been
            // received anyway. The most common applications which are sensitive to read latency are MIDI and IrDA
            // applications.
            if (ioctl(device->fileDescriptor, IOSSDATALAT, &mics) == -1) {
                if (verbosity > 0) printf("Error setting receive latency for device %s to %lu microseconds - %s(%d).\n", device->portSpec, mics, strerror(errno), errno);
                if (strstr(configString, "Lenient") == NULL) return(PsychError_system);
            }
        }
    }
#endif

#if PSYCH_SYSTEM == PSYCH_LINUX
    // Linux specific implementation: Controls low-latency mode. If a value <= 0.001 seconds aka <= 1msec
    // is requested, we set the low-latency flag, otherwise we clear it:
    if ((p = strstr(configString, "ReceiveLatency="))) {
        // Set receive latency for low-level driver: Minimum latency between dequeue and notify
        // operation on reception of new input data. The granularity is microseconds...
        if ((1!=sscanf(p, "ReceiveLatency=%f", &infloat)) || (infloat < 0)) {
            if (verbosity > 0) printf("Invalid parameter for ReceiveLatency set! Typo, or negative value provided.\n");
            return(PsychError_user);
        }
        else {
            // Set timeout: It is in granularity microseconds, so we need to quantize to microseconds:
            mics = (infloat <= 0.001) ? 1 : 0;

            // Retrieve current settings:
            if (ioctl(device->fileDescriptor, TIOCGSERIAL, &serialstruct) == -1) {
                if (verbosity > 0) printf("Error setting receive latency for device %s to %s latency mode - %s(%d).\n", device->portSpec, ((mics) ? "low" : "normal"), strerror(errno), errno);
                if (strstr(configString, "Lenient") == NULL) return(PsychError_system);
            }

            // Clear or set low-latency flag:
            // Setting it enables whatever leads to low receive latency on
            // a given driver. In the FTDI serial-converter driver (ftdi_sio.c)
            // it forces the chips receive latency timer to its minimum value of
            // 1 msec, overriding whatever sysfs value is set. Clearing it will
            // enable use of the latency_timer value fro sysfs. Other drivers or
            // the higher level subsystem may make other/additional optimizations.
            if (mics) {
                serialstruct.flags |= ASYNC_LOW_LATENCY;
            } else {
                serialstruct.flags &= ~ASYNC_LOW_LATENCY;
            }

            // Set updated settings:
            if (ioctl(device->fileDescriptor, TIOCSSERIAL, &serialstruct) == -1) {
                if (verbosity > 0) printf("Error setting receive latency for device %s to %s latency mode - %s(%d).\n", device->portSpec, ((mics) ? "low" : "normal"), strerror(errno), errno);
                if (strstr(configString, "Lenient") == NULL) return(PsychError_system);
            }
        }
    }
#endif

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
                    if (verbosity > 0) printf("Reallocation of Inputbuffer of device %s for new size %i failed! (%s)", device->portSpec, inint, strerror(errno));
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
                    if (verbosity > 0) printf("Allocation of Inputbuffer of device %s of size %i failed! (%s)", device->portSpec, inint, strerror(errno));
                    return(PsychError_outofMemory);
                }

                // Worked. Assign new values:
                device->readBuffer = (unsigned char*) p;
                device->readBufferSize = (unsigned int) inint;
            }

            // Zerofill, so each page of memory gets touched once and we fault in
            // all memory pages to improve realtime behaviour of rt reader thread by
            // minimizing / avoiding page-faults:
            memset(device->readBuffer, 0, device->readBufferSize);
        }
    }

    if ((p = strstr(configString, "PollLatency="))) {
        // Set polling latency for busy-waiting read operations:

        // Make sure we don't change non-mutex-protected variables behind the
        // back of our readerThread by only allowing this function to be called
        // with inactive thread:
        if (device->readerThread) {
            if (verbosity > 0) printf("Assigned PollLatency= while background read operations are already enabled! Disable first via 'StopBackgroundRead'!\n");
            return(PsychError_user);
        }

        if ((1!=sscanf(p, "PollLatency=%f", &infloat)) || (infloat < 0)) {
            if (verbosity > 0) printf("Invalid parameter for PollLatency set! Typo, or negative value provided.\n");
            return(PsychError_user);
        }
        else {
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

            // Warn user if readGranularity is possibly to high for system to handle properly without weird side-effects:
            if ((device->readGranularity > 255) && (verbosity > 1)) printf("IOPort: WARNING: In call to 'StartBackgroundRead', requested read granularity of %i bytes exceeds maximum safe size of 255 Bytes.\nThis can cause malfunctions or unexpected behaviour/data loss on some systems with some device drivers!\n", device->readGranularity);

            // Allocate sufficiently large timestamp buffer:
            device->timeStamps = (double*) calloc(sizeof(double), device->readBufferSize / device->readGranularity);

            // Create & Init the mutex:
            if ((rc=PsychInitMutex(&(device->readerLock)))) {
                printf("PTB-ERROR: In StartBackgroundRead(): Could not create readerLock mutex lock [%s].\n", strerror(rc));
                return(PsychError_system);
            }

            // Perform lock->unlock mutex sequence to inject some memory ordering barriers here, so all our
            // settings are picked up by the newborn thread:
            if ((rc=PsychLockMutex(&(device->readerLock))) || (rc=PsychUnlockMutex(&(device->readerLock)))) {
                printf("PTB-ERROR: In StartBackgroundRead(): Could not lock + unlock readerLock mutex lock [%s].\n", strerror(rc));
                return(PsychError_system);
            }

            // Create and startup thread:
            if ((rc=PsychCreateThread(&(device->readerThread), NULL, PsychSerialUnixGlueReaderThreadMain, (void*) device))) {
                printf("PTB-ERROR: In StartBackgroundRead(): Could not create background reader thread [%s].\n", strerror(rc));
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

    // Proof-of-concept test code: Not for public use!
    // Async triggerbyte emission via parallel thread requested?
    if ((p = strstr(configString, "JLFireTrigger="))) {
        if (1!=sscanf(p, "JLFireTrigger=%f", &infloat)) {
            if (verbosity > 0) printf("Invalid parameter for JLFireTrigger set!\n");
            return(PsychError_user);
        }
        else {
            // Store target time in device struct:
            device->triggerWhen = (double) infloat;

            // Create and startup trigger thread: It will detach itself from us, do
            // its job and then die lonely and forgotten without us caring:
            psych_thread threadhandle;
            if ((rc=PsychCreateThread(&threadhandle, NULL, PsychSerialUnixGlueJLTriggerThreadMain, (void*) device))) {
                printf("PTB-ERROR: In JLFireTrigger(): Could not create background trigger thread [%s].\n", strerror(rc));
                return(PsychError_system);
            }
        }
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
    unsigned int lsr = 0;   // Serial transmitter line status register.
    int outqueue_pending;    // Pending bytes in output queue.

    // Nonblocking mode?
    if (blocking <= 0) {
        // Yep. Set filedescriptor to non-blocking mode:
        // Set the O_NONBLOCK flag so subsequent I/O will not block.
        // See fcntl(2) ("man 2 fcntl") for details.
        if (PsychSerialUnixGlueFcntl(device, O_NONBLOCK) == -1) {
            sprintf(errmsg, "Error setting O_NONBLOCK on device %s for non-blocking writes - %s(%d).\n", device->portSpec, strerror(errno), errno);
            return(-1);
        }

        // Write the data: Take pre- and postwrite timestamps.
        PsychGetAdjustedPrecisionTimerSeconds(&timestamp[1]);
        if ((nwritten = write(device->fileDescriptor, writedata, amount)) == -1) {
            sprintf(errmsg, "Error during write to device %s - %s(%d).\n", device->portSpec, strerror(errno), errno);
            return(-1);
        }
        PsychGetAdjustedPrecisionTimerSeconds(&timestamp[2]);
    }
    else {
        // Nope. Set filedescriptor to blocking mode:
        // Clear the O_NONBLOCK flag so subsequent I/O will block.
        // See fcntl(2) ("man 2 fcntl") for details.
        if (PsychSerialUnixGlueFcntl(device, 0) == -1) {
            sprintf(errmsg, "Error clearing O_NONBLOCK on device %s for blocking writes - %s(%d).\n", device->portSpec, strerror(errno), errno);
            return(-1);
        }

        // Write the data: Take pre- and postwrite timestamps.
        PsychGetAdjustedPrecisionTimerSeconds(&timestamp[1]);
        if ((nwritten = write(device->fileDescriptor, writedata, amount)) == -1) {
            sprintf(errmsg, "Error during write to device %s - %s(%d).\n", device->portSpec, strerror(errno), errno);
            return(-1);
        }
        PsychGetAdjustedPrecisionTimerSeconds(&timestamp[2]);

        // Special polling mode to wait for transmit completion instead of tcdrain() + wait?
        if (blocking == 2) {
            // Yes. Use a tight polling loop which spin-waits until the driver output queue is empty (contains zero pending bytes):
            outqueue_pending = 1;
            while (outqueue_pending > 0) {
                // Take timestamp for this iteration:
                PsychGetAdjustedPrecisionTimerSeconds(&timestamp[3]);

                // Poll:
                ioctl(device->fileDescriptor, TIOCOUTQ, &outqueue_pending);
            }
        }
        else if ((PSYCH_SYSTEM == PSYCH_LINUX) && (blocking == 3)) {
            // Yes. Use a tight polling loop which spin-waits on the transmitter idle flag:
            // TODO: This may be not what we want, because the transmitter might be temporarily
            // idle although there is data pending in the higher-level write buffers. A more
            // proper combo might be to have blocking==2 polling for empty queue followed by
            // this call, so both conditions must be met.
            #if PSYCH_SYSTEM == PSYCH_LINUX
                // ioctl supported by device?
                if (ioctl(device->fileDescriptor, TIOCSERGETLSR, &lsr)==-1) {
                    // Error! Can't use this method!
                    sprintf(errmsg, "Error during write to device %s while trying to use blocking mode 3 to wait for write completion - %s(%d).\nMaybe function unsupported by device??\n", device->portSpec, strerror(errno), errno);
                    return(-1);
                }

                // Seems to work. Go into the spin-waiting loop:
                while (!(lsr & TIOCSER_TEMT)) {
                    // Take timestamp for this iteration:
                    PsychGetAdjustedPrecisionTimerSeconds(&timestamp[3]);

                    // Poll:
                    ioctl(device->fileDescriptor, TIOCSERGETLSR, &lsr);
                }
            #else
                (void) lsr;
            #endif
        }
        else {
            // Take timestamp for completeness although it doesn't make much sense in the blocking case:
            PsychGetAdjustedPrecisionTimerSeconds(&timestamp[3]);

            // Flush the write buffer and wait for write completion on physical hardware:
            if ((!device->dontFlushOnWrite) && (tcdrain(device->fileDescriptor) == -1)) {
                sprintf(errmsg, "Error during write to device %s while draining the write buffers - %s(%d).\n", device->portSpec, strerror(errno), errno);
                return(-1);
            }
        }
    }

    // Write successfully completed if we reach this point. Take timestamp, clear error message, return:
    PsychGetAdjustedPrecisionTimerSeconds(&timestamp[0]);
    errmsg[0] = 0;

    return(nwritten);
}

int PsychIOOSReadSerialPort(PsychSerialDeviceRecord* device, void** readdata, unsigned int amount, int blocking, char* errmsg, double* timestamp)
{
    double timeout;
    int raPos, i;
    int nread = 0;
    int gotamount, reqamount;
    unsigned char* tmpbuffer;
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
        // Background read active?
        if (device->readerThread) {
            // Just return what is available in the internal buffer:
            nread = PsychSerialUnixGlueAsyncReadbufferBytesAvailable(device);
        }
        else {
            // Set filedescriptor to non-blocking mode:
            // Set the O_NONBLOCK flag so subsequent I/O will not block.
            // See fcntl(2) ("man 2 fcntl") for details.
            if (PsychSerialUnixGlueFcntl(device, O_NONBLOCK) == -1) {
                sprintf(errmsg, "Error setting O_NONBLOCK on device %s for non-blocking read - %s(%d).\n", device->portSpec, strerror(errno), errno);
                return(-1);
            }

            // Read the data, at most 'amount' bytes, nonblocking:
            if ((nread = read(device->fileDescriptor, device->readBuffer, amount)) == -1) {
                // Nothing to read? In nonblocking mode we simply return zero bytes read:
                if (errno == EAGAIN) return(0);

                // Some other error:
                sprintf(errmsg, "Error during non-blocking read from device %s - %s(%d).\n", device->portSpec, strerror(errno), errno);
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

            while((*timestamp < timeout) && (PsychSerialUnixGlueAsyncReadbufferBytesAvailable(device) < (int) amount)) {
                PsychGetAdjustedPrecisionTimerSeconds(timestamp);
                PsychWaitIntervalSeconds(device->pollLatency);
            }

            // Return amount of available data:
            nread = PsychSerialUnixGlueAsyncReadbufferBytesAvailable(device);
        }
        else {
            // Set filedescriptor to blocking mode:
            // Clear the O_NONBLOCK flag so subsequent I/O will block.
            // See fcntl(2) ("man 2 fcntl") for details.
            if (PsychSerialUnixGlueFcntl(device, 0) == -1) {
                sprintf(errmsg, "Error clearing O_NONBLOCK on device %s for blocking read - %s(%d).\n", device->portSpec, strerror(errno), errno);
                return(-1);
            }

            // Init start position for read to store data to the start of readBuffer:
            tmpbuffer = device->readBuffer;

            // Perform repeated read() requests in small chunks until the full 'amount'
            // of requested data has been read() (or error, timeout etc.). We need to
            // perform polling/blocking reads() in relatively small chunks because the
            // size of the serial port receive queue is limited to a few kilobytes at
            // best - often less, so data can easily get lost or we can deadlock on
            // reads() > maxqueuesize if we don't do it in < queuesize chunks:
            while (amount > 0) {
                // Perform pull iteration:

                // Setup minimum byte counter for 'naccumread' Bytes blocking reads:
                gotamount = PsychSerialUnixGlueSetBlockingMinBytes(device, amount);
                if ((int) amount != gotamount) {
                    // Didn't get what we wanted: Error!
                    // Unless this is just the case where we clamp to the 255 bytes max on Unix:
                    if (!(gotamount == 255 && amount > 255)) {
                        sprintf(errmsg, "Error setting wanted minimum amount of bytes %i on device %s for blocking read - %s(%d). Got %i instead!\n", amount, device->portSpec, strerror(errno), errno, gotamount);
                        return(-1);
                    }
                }

                // Skip this polling for first byte if BlockingBackgroundRead is non-zero:
                if (device->isBlockingBackgroundRead == 0) {
                    // Define timeout deadline for poll:
                    PsychGetAdjustedPrecisionTimerSeconds(&timeout);
                    *timestamp = timeout;

                    // Timeout is one interbyte timeout:
                    timeout += device->readTimeout;

                    // Wait polling until timeout or 1 byte available:
                    while((*timestamp < timeout) && (PsychIOOSBytesAvailableSerialPort(device) < 1)) {
                        PsychGetAdjustedPrecisionTimerSeconds(timestamp);
                        PsychWaitIntervalSeconds(device->pollLatency);
                    }

                    if (PsychIOOSBytesAvailableSerialPort(device) < 1) {
                        // Ok, first byte for this read() iteration didn't arrive within timeout period.
                        // Simply break out of while() loop to terminate this read.
                        break;
                    }
                }

                // Read the data, at most (and at least, unless timeout occurs) 'gotamount' bytes, blocking:
                reqamount = gotamount;
                if ((gotamount = read(device->fileDescriptor, tmpbuffer, reqamount)) == -1) {
                    sprintf(errmsg, "Error during blocking read from device %s - %s(%d).\n", device->portSpec, strerror(errno), errno);
                    return(-1);
                }
                // Successfully read nread >= 0 bytes.

                // Advance read() storebuffer pointer and decrement count of pending data to read:
                tmpbuffer+= gotamount;
                amount   -= gotamount;
                nread    += gotamount;

                // Is a line-terminator set, the device in cooked input processing mode and the last read character is the lineterminator?
                if ((gotamount > 0) && (device->cookedMode > 0) && (device->lineTerminator == *(tmpbuffer - 1))) {
                    // Line terminator detected in cooked mode: Break out of read() loop!
                    break;
                }

                // read() prematurely completed due to insufficient amount of data available and
                // interbyte-timeout occured? If so, we break out of the pull-loop:
                if (gotamount < reqamount) break;

                // Next pull loop iteration...
            } // Try to read() next chunk of data.

            // Done with this blocking read.

            // Reset minbytes back to zero:
            PsychSerialUnixGlueSetBlockingMinBytes(device, 0);
        }
    }

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
        nread = (nread > (int) amount) ? amount : nread;

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

        // Assign returned data:
        *readdata = (void*) device->readBuffer;
    }

    return(nread);
}

int PsychIOOSBytesAvailableSerialPort(PsychSerialDeviceRecord* device)
{
    int navail = 0;

    if (device->readerThread) {
        // Async reader poll: Calculate what is available in the internal buffer:
        navail = PsychSerialUnixGlueAsyncReadbufferBytesAvailable(device);
    }
    else {
        // Regular poll: Perform query to system:
        if (ioctl(device->fileDescriptor, FIONREAD, &navail)!=0) {
            if (verbosity > 0) printf("Error during 'BytesAvailable': ioctl - FIONREAD byte available query on device %s returned %s(%d)\n", device->portSpec, strerror(errno), errno);
            return(-1);
        }
    }

    return(navail);
}

void PsychIOOSFlushSerialPort(PsychSerialDeviceRecord* device)
{
    if (tcdrain(device->fileDescriptor)!=0) {
        if (verbosity > 0) printf("Error during 'Flush': tcdrain() on device %s returned %s(%d)\n", device->portSpec, strerror(errno), errno);
    }

    return;
}

void PsychIOOSPurgeSerialPort(PsychSerialDeviceRecord* device)
{
    if (tcflush(device->fileDescriptor, TCIOFLUSH)!=0) {
        if (verbosity > 0) printf("Error during 'Purge': tcflush(TCIFLUSH) on device %s returned %s(%d)\n", device->portSpec, strerror(errno), errno);
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

#endif
