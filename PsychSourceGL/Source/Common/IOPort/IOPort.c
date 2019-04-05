/*
    IOPort.c

    PROJECTS:

        IOPort only.

    AUTHORS:

        Allen.Ingling@nyu.edu                   awi
        mario.kleiner.de@gmail.com              mk

    PLATFORMS:

        All.

    HISTORY:

        12/20/2004    awi       Wrote it.
        04/10/2008    mk        Started to extend/rewrite it to become a full-fledged generic I/O driver (serial port, parallel port, etc...).

    DESCRIPTION:

*/

#include "IOPort.h"

// Maximum number of slots aka allowed open ports:
#define PSYCH_MAX_IOPORTS 100

// Level of verbosity:
int verbosity = 4;

// Bank of all port records:
static PsychPortIORecord portRecordBank[PSYCH_MAX_IOPORTS];

// Total count of open ports:
static unsigned int portRecordCount = 0;

// Array with function syntax strings:
#define MAX_SYNOPSIS_STRINGS 50
static const char *synopsisSYNOPSIS[MAX_SYNOPSIS_STRINGS];

// Initialize usage info -- function overview:
const char** InitializeSynopsis(void)
{
    int i=0;
    const char **synopsis = synopsisSYNOPSIS;  //abbreviate the long name

    synopsis[i++] = "IOPort - A Psychtoolbox driver for general port I/O:\n";

    synopsis[i++] = "\nGeneral information:\n";
    synopsis[i++] = "version = IOPort('Version');";
    synopsis[i++] = "oldlevel = IOPort('Verbosity' [,level]);";

    synopsis[i++] = "\nGeneral commands for all types of input/output ports:\n";
    synopsis[i++] = "IOPort('Close', handle);";
    synopsis[i++] = "IOPort('CloseAll');";
    synopsis[i++] = "[nwritten, when, errmsg, prewritetime, postwritetime, lastchecktime] = IOPort('Write', handle, data [, blocking=1]);";
    synopsis[i++] = "IOPort('Flush', handle);";
    synopsis[i++] = "[data, when, errmsg] = IOPort('Read', handle [, blocking=0] [, amount]);";
    synopsis[i++] = "navailable = IOPort('BytesAvailable', handle);";
    synopsis[i++] = "IOPort('Purge', handle);";

    synopsis[i++] = "\nCommands specific to serial ports:\n";
    synopsis[i++] = "[handle, errmsg] = IOPort('OpenSerialPort', port [, configString]);";
    synopsis[i++] = "IOPort('ConfigureSerialPort', handle, configString);";

    synopsis[i++] = NULL;  //this tells IOPORTDisplaySynopsis where to stop
    if (i > MAX_SYNOPSIS_STRINGS) {
        PrintfExit("%s: increase dimension of synopsis[] from %ld to at least %ld and recompile.",__FILE__,(long)MAX_SYNOPSIS_STRINGS,(long)i);
    }

    return(synopsisSYNOPSIS);
}

PsychError IOPORTVerbosity(void)
{
    static char useString[] = "oldlevel = IOPort('Verbosity' [,level]);";
    static char synopsisString[] =
        "Set level of verbosity for error/warning/status messages. 'level' optional, new level "
        "of verbosity. 'oldlevel' is the old level of verbosity. The following levels are "
        "supported: 0 = Shut up. 1 = Print errors, 2 = Print also warnings, 3 = Print also some info, "
        "4 = Print more useful info (default), >5 = Be very verbose (mostly for debugging the driver itself). ";
    static char seeAlsoString[] = "";

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
    if (level > -1) verbosity = level;

    return(PsychError_none);
}

// Display usage info, aka synopsis if IOPort is called without any subfunction arguments:
PsychError IOPORTDisplaySynopsis(void)
{
    int i;

    for (i = 0; synopsisSYNOPSIS[i] != NULL; i++) printf("%s\n",synopsisSYNOPSIS[i]);

    return(PsychError_none);
}

// Exit function: Called at shutdown time:
PsychError PsychExitIOPort(void)
{
    int i;

    // Close all open ports:
    for (i=0; i < PSYCH_MAX_IOPORTS; i++) {
        if (portRecordBank[i].portType) PsychCloseIOPort(i);
    }

    return(PsychError_none);
}

// Init function: Called at startup time:
PsychError PsychInitIOPort(void)
{
    int i;

    // Init array of IOPort records to empty:
    for (i=0; i < PSYCH_MAX_IOPORTS; i++) {
        portRecordBank[i].portType = KPsychIOPortNone;
    }

    // No ports open yet:
    portRecordCount = 0;

    return(PsychError_none);
}

// Receive pointer to PsychPortIORecord for given userspace port handle:
PsychPortIORecord* PsychGetPortIORecord(int handle)
{
    char errmsg[1000];
    errmsg[0] = 0;

    if (handle < 0 || handle >= PSYCH_MAX_IOPORTS || portRecordBank[handle].portType == KPsychIOPortNone) {
        sprintf(errmsg, "Invalid port handle %i provided. No such port open. Maybe you closed it beforehand?", handle);
        PsychErrorExitMsg(PsychError_user, errmsg);
    }

    // Valid handle, referring to open port: Return its data structure...
    return(&portRecordBank[handle]);
}

// Close port referenced by 'handle':
PsychError PsychCloseIOPort(int handle)
{
    PsychPortIORecord* portRecord;

    // Retrieve and assign open port for handle: Will check for invalid handles and closed ports...
    portRecord = PsychGetPortIORecord(handle);

    switch(portRecord->portType) {
        case kPsychIOPortSerial:
            // Close serial port:
            PsychIOOSCloseSerialPort(portRecord->device);
        break;

        default:
            PsychErrorExitMsg(PsychError_internal, "Unknown portType - Unsupported.");
    }

    // Reset port to "no device" aka closed:
    portRecord->device = NULL;
    portRecord->portType = KPsychIOPortNone;

    // Decrement count of open ports:
    portRecordCount--;

    return(PsychError_none);
}

// Write data to port:
int PsychWriteIOPort(int handle, void* writedata, unsigned int amount, int blocking, char* errmsg, double* timestamp)
{
    PsychPortIORecord* portRecord = PsychGetPortIORecord(handle);

    switch(portRecord->portType) {
        case kPsychIOPortSerial:
            // Write to serial port:
            return(PsychIOOSWriteSerialPort(portRecord->device, writedata, amount, blocking, errmsg, timestamp));
        break;

        default:
            PsychErrorExitMsg(PsychError_internal, "Unknown portType - Unsupported.");
    }

    // Not reached, just to make compiler happy:
    return(0);
}

int PsychReadIOPort(int handle, void** readbuffer, unsigned int amount, int blocking, char* errmsg, double* timestamp)
{
    PsychPortIORecord* portRecord = PsychGetPortIORecord(handle);

    switch(portRecord->portType) {
        case kPsychIOPortSerial:
            // Read from serial port:
            return(PsychIOOSReadSerialPort(portRecord->device, readbuffer, amount, blocking, errmsg, timestamp));
        break;

        default:
            PsychErrorExitMsg(PsychError_internal, "Unknown portType - Unsupported.");
    }

    // Not reached, just to make compiler happy:
    return(0);
}

int PsychBytesAvailableIOPort(int handle)
{
    PsychPortIORecord* portRecord = PsychGetPortIORecord(handle);

    switch(portRecord->portType) {
        case kPsychIOPortSerial:
            // Read from serial port:
            return(PsychIOOSBytesAvailableSerialPort(portRecord->device));
        break;

        default:
            PsychErrorExitMsg(PsychError_internal, "Unknown portType - Unsupported.");
    }

    // Not reached, just to make compiler happy:
    return(0);
}

void PsychPurgeIOPort(int handle)
{
    PsychPortIORecord* portRecord = PsychGetPortIORecord(handle);

    switch(portRecord->portType) {
        case kPsychIOPortSerial:
            // Purge serial port:
            PsychIOOSPurgeSerialPort(portRecord->device);
        break;

        default:
            PsychErrorExitMsg(PsychError_internal, "Unknown portType - Unsupported.");
    }
}

void PsychFlushIOPort(int handle)
{
    PsychPortIORecord* portRecord = PsychGetPortIORecord(handle);

    switch(portRecord->portType) {
        case kPsychIOPortSerial:
            // Purge serial port:
            PsychIOOSFlushSerialPort(portRecord->device);
        break;

        default:
            PsychErrorExitMsg(PsychError_internal, "Unknown portType - Unsupported.");
    }
}


PsychError IOPORTIOPort(void)
{
    return(PsychError_none);
}

// Close given port.
PsychError IOPORTClose(void)
{
    static char useString[] = "IOPort('Close', handle);";
    static char synopsisString[] =
        "Close a I/O port device 'handle'. This works for all device types.\n";
    static char seeAlsoString[] = "'CloseAll'";

    int handle;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(1));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(0));     // The maximum number of outputs

    // Get handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);

    // Close associated port:
    PsychCloseIOPort(handle);

    return(PsychError_none);
}

// Close all ports of all devices -- Effectively a shutdown.
PsychError IOPORTCloseAll(void)
{
    static char useString[] = "IOPort('CloseAll');";
    static char synopsisString[] =
        "Close all open I/O port devices, regardless of type. Release all ressources.\n";
    static char seeAlsoString[] = "'Close'";

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(0));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(0)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(0));     // The maximum number of outputs

    // Close all ports by calling our shutdown function:
    PsychExitIOPort();

    return(PsychError_none);
}

// Open a serial port on a serial port device:
PsychError IOPORTOpenSerialPort(void)
{
    static char useString[] = "[handle, errmsg] = IOPort('OpenSerialPort', port [, configString]);";
    static char synopsisString[] =
        "Open a serial port device, return a 'handle' to it.\n"
        "If a port can't be opened, the function will abort with error, unless the "
        "level of verbosity is set to zero, in which case the function will silently "
        "fail, but return an invalid (negative) handle to signal the failure to the "
        "calling script. The optional return argument 'errmsg' contains a text string "
        "which is either empty on success, or contains a descriptive error message. \n"
        "'port' is usually a name string that defines the serial port device "
        "to open. On MS-Windows this could be, e.g., 'COM1' or 'COM2' etc. On "
        "Apple OS/X, it is the path to a BSD device file, e.g., '/dev/cu.usbserial-FT3Z95V5' "
        "for a serial-over-USB device with unique id FT3Z95V5. On GNU/Linux it could be "
        "'/dev/ttyS0' for the first real serial port, or '/dev/ttyUSB0' for the first "
        "serial-over-USB device.\n\n"
        "The optional string 'configString' is a string with pairs of paramName=paramValue "
        "tokens, separated by a delimiter, e.g., a space. It allows to specify specific "
        "values 'paramValue' to specific serial port parameters 'paramName'. Not all "
        "parameters are supported by all operating systems, and all settings have reasonable "
        "defaults. Settings unknown to a specific operating system are ignored.\n"
        "The following is a list of (possibly) supported parameters with their defaults:\n\n"
        "Lenient -- If this keyword is present, then the driver will carry on on certain error conditions "
        "instead of aborting. This is sometimes neccessary for some special cases like virtual com ports or "
        "other non-standard setups.\n\n"
        "BaudRate=9600 -- The baud transmission rate of the connection. Standard baud rates include "
        "110, 300, 600, 1200, 2400, 4800, 9600, 14400, 19200, 38400, 57600, 115200, 128000 and "
        "256000 bits per second. Not all values may be supported by all operating systems and drivers.\n\n"
        "Parity=None   -- Type of parity checking: None, Even, Odd.\n\n"
        "DataBits=8    -- Number of data bits per packet: 5,6,7 or 8, on Windows also 16.\n\n"
        "StopBits=1    -- Number of stop bits per packet: 1 or 2.\n\n"
        "FlowControl=None  -- Type of flow control: None, Hardware (RTS/CTS lines), Software (XON/XOFF characters).\n\n"
        "ReceiverEnable=1 -- A non-zero setting will enable the serial receiver, a zero setting will disable the receiver. "
        "This setting may not be supported by all operating systems and hardware. In such a case, the receiver will always "
        "be enabled, irrespective of this setting. On MS-Windows, this setting is not available at all.\n\n"
        "Terminator=os default  -- Type of terminator, given as ASCII character value, e.g., 13 for char(13) aka CR or 10 for LF. Currently "
        "only used in async read mode (see 'StartBackgroundRead' below) if the 'ReadFilterFlags' are set to include value 4, "
        "or on OS/X and Linux in 'Cooked' processing mode as line delimiter. A setting of -1 will try to disable the line terminator.\n\n"
        "DTR=os default    -- Setting for 'Data Terminal Ready' pin: 0 or 1.\n\n"
        "RTS=os default    -- Setting for 'Request To Send' pin: 0 or 1.\n\n"
        "BreakBehaviour=Ignore -- Behaviour if a 'Break Condition' is detected on the line: Ignore, Flush, Zero. On Windows, this setting is ignored.\n\n"
        "OutputBufferSize=4096 -- Size of output buffer in bytes.\n\n"
        "InputBufferSize=4096 -- Size of input buffer in bytes. You can't read more than that amount per read command.\n\n"
        "HardwareBufferSizes=input,output -- Set size of the hardware driver internal input and output buffers in bytes. "
        "E.g., HardwareBufferSizes=32768,8192 would set the input buffer to 32768 bytes and the output buffer to 8192 bytes. "
        "This function is currently only supported on Windows, ignored on other systems. It is only a polite hint to the driver, "
        "the serial port driver is free to ignore the request and choose any buffer sizes or buffering strategy it finds appropriate. "
        "By default, this parameter is not set by IOPort and the hardware driver uses some built-in reasonable setting.\n\n"
        "The following timeout values are inter-byte timeouts. You specify how much time reception or "
        "transmission of a single byte is allowed to take. Timeout occurs if more than that time elapses "
        "between send/reception of two consecutive bytes or if the total amount of time exceeds the "
        "number of bytes, times the interbyte timeout value. A value of zero means not to use any timeout, "
        "in which case a blocking read or write may take forever. If a timeout occurs, the read or write "
        "operation will be aborted. \n"
        "Granularity of timeout settings is 100 msecs on OS/X and Linux, 1 msec on Windows, all values "
        "are rounded to the closest value matching that granularity. The minimal timeout is 100 msecs "
        "on OS/X and Linux, about 6 msecs on Windows.\n\n"
        "SendTimeout=1.0 -- Interbyte send timeout in seconds. Only used on Windows.\n\n"
        "ReceiveTimeout=1.0 -- Interbyte receive timeout in seconds.\n\n"
        "ReceiveLatency -- Latency in seconds for processing of new input bytes. Only used on OS/X and Linux for some devices.\n\n"
        "PollLatency=0.0005 (0.001 on Windows) -- Latency between polls in seconds for polling in some 'Read' operations.\n\n"
        "ProcessingMode=Raw -- Mode of input/output processing: Raw or Cooked. On Windows, only Raw (binary) mode is supported.\n\n"
        "DontFlushOnWrite=0 -- Do not flush the serial port write buffer at device close time or during blocking writes. "
        "This can be set to 1 to work around broken serial port drivers, but it may disrupt any kind of timing sensitive "
        "algorithms that interact with the serial port! Only use if you really know what you're doing!\n\n"
        "StartBackgroundRead=readGranularity -- Enable asynchronous background read operations on the port. "
        "A parallel background thread is started which tries to fetch 'readGranularity' bytes of data, "
        "polling the port every 'PollLatency' seconds for at least 'readGranularity' bytes of data. 'InputBufferSize' must be an "
        "integral multiple of 'readGranularity' for this to work. Later IOPort('Read') commands will pull collected data from "
        "the InputBuffer in quanta of at most 'readGranularity' bytes per invocation. This function is useful for background "
        "data collection from devices that stream some data at a constant rate. You set up background read, let the parallel "
        "thread do all data collection in the background and collect the data at the end of a session with a sequence of "
        "IOPort('Read') calls. This way, data collection doesn't clutter your main experiment script.\n\n"
        "BlockingBackgroundRead=0 -- Perform blocking background reads instead of polling reads, if set to 1.\n\n"
        "StopBackgroundRead -- Stop running background read operation, discard all pending data.\n\n"
        "ReadFilterFlags=0 -- Special flags to specify certain post-processing operations on read input data.\n"
        "* A setting of 1 will enable special filtering for serial input data from the CMU or PST response button boxes. "
        "  Redundant data bytes received will be discarded - only bytes that are different from their predecessor are stored. "
        "  All read data has a 4-Byte 32 bit count of total bytes read and a 4-Byte count of sampling delta in microseconds attached. "
        "  You should set 'readGranularity' = 9 for best effect with the CMU or PST button boxes or compatible devices. \n"
        "* A setting of 2 will filter out CR and LF character codes 10 and 13 from the inputstream.\n"
        "* A setting of 4 will implement simple line-buffering for async reads: Read up to 'readGranularity' bytes per iteration, "
        "  or until 'Terminator' character encountered, whatever comes first. Zero-Pad to full 'readGranularity' bytes in any case. "
        "  Read timestamps in this line-buffered mode correspond to the reception of the first byte of a line, not the last one!\n"
        "\n\n";

    static char seeAlsoString[] = "'CloseAll'";

    #if PSYCH_SYSTEM == PSYCH_WINDOWS
    // Difference to Unices: PollLatency defaults to 1 msecs instead of 0.5 msecs due to shoddy windows scheduler:
    static char defaultConfig[] = "BaudRate=9600 Parity=None DataBits=8 StopBits=1 FlowControl=None PollLatency=0.001 SendTimeout=1.0 ReceiveTimeout=1.0 ProcessingMode=Raw BreakBehaviour=Ignore OutputBufferSize=4096 InputBufferSize=4096 DontFlushOnWrite=0";
    #endif

    #if PSYCH_SYSTEM == PSYCH_OSX
    static char defaultConfig[] = "BaudRate=9600 Parity=None DataBits=8 StopBits=1 FlowControl=None PollLatency=0.0005 SendTimeout=1.0 ReceiveTimeout=1.0 ProcessingMode=Raw BreakBehaviour=Ignore OutputBufferSize=4096 InputBufferSize=4096 DontFlushOnWrite=0";
    #endif

    #if PSYCH_SYSTEM == PSYCH_LINUX
    static char defaultConfig[] = "BaudRate=9600 Parity=None DataBits=8 StopBits=1 FlowControl=None PollLatency=0.0005 SendTimeout=1.0 ReceiveTimeout=1.0 ProcessingMode=Raw BreakBehaviour=Ignore OutputBufferSize=4096 InputBufferSize=4096 DontFlushOnWrite=0";
    #endif

    char finalConfig[2000];
    char errmsg[1024];
    char* portSpec = NULL;
    char* configString = NULL;
    PsychSerialDeviceRecord* device = NULL;
    int handle;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(2));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(2));     // The maximum number of outputs

    // Get required portSpec:
    PsychAllocInCharArg(1, kPsychArgRequired, &portSpec);

    // Get the optional configString:
    if (!PsychAllocInCharArg(2, kPsychArgOptional, &configString)) {
        // No config string specified: Assign default string:
        sprintf(finalConfig, "%s", defaultConfig);
    }
    else {
        // Config string provided: Prepend it to default string. That way, all settings
        // spec'd by usercode will override the defaultConfig, but non-specified settings
        // will be provided by defaultConfig:
        sprintf(finalConfig, "%s %s", configString, defaultConfig);
    }

    // Search for a free slot:
    if (portRecordCount >= PSYCH_MAX_IOPORTS) PsychErrorExitMsg(PsychError_user, "Maximum number of open Input/Output ports exceeded.");

    // Iterate until end or free slot:
    for (handle=0; (handle < PSYCH_MAX_IOPORTS) && (portRecordBank[handle].portType); handle++);
    if (portRecordBank[handle].portType) PsychErrorExitMsg(PsychError_user, "Maximum number of open Input/Output ports exceeded.");
    // handle is index into our port record...

    // Call OS specific open routine for serial port:
    device = PsychIOOSOpenSerialPort(portSpec, finalConfig, errmsg);

    // Copy out optional errmsg string:
    PsychCopyOutCharArg(2, kPsychArgOptional, errmsg);

    if (device == NULL) {
        // Special case: Could not open port, but verbosity level is zero, no conventional
        // error return possible as this would clutter the console with output. Simply cancel
        // open op and return a negative handle to signal failure to user code:
        PsychCopyOutDoubleArg(1, kPsychArgRequired, -1);
        return(PsychError_none);
    }

    // If we reach this point, then device is the pointer to the class specific device struct and the
    // open operation was successfull. Build port struct:
    portRecordBank[handle].portType = kPsychIOPortSerial;
    portRecordBank[handle].device = (void*) device;
    portRecordCount++;

    // Return handle to new serial port object:
    PsychCopyOutDoubleArg(1, kPsychArgRequired, (double) handle);

    return(PsychError_none);
}

// Open a serial port on a serial port device:
PsychError IOPORTConfigureSerialPort(void)
{
    static char useString[] = "IOPort('ConfigureSerialPort', handle, configString);";
    static char synopsisString[] =
        "(Re-)Configure a serial port device, specified by 'handle'.\n"
        "The string 'configString' is a string with pairs of paramName=paramValue "
        "tokens, separated by a delimiter, e.g., a space. It allows to specify specific "
        "values 'paramValue' to specific serial port parameters 'paramName'. Not all "
        "parameters are supported by all operating systems, and all settings have reasonable "
        "defaults. Settings unknown to a specific operating system are ignored.\n\n"
        "See the help of 'OpenSerialPort' for possible settings.";

    static char seeAlsoString[] = "'OpenSerialPort'";
    char* configString = NULL;
    int handle;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(2));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(2)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(0));     // The maximum number of outputs

    // Get required portSpec:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);

    // Get the configString:
    PsychAllocInCharArg(2, kPsychArgRequired, &configString);

    // Return return value of configuration call:
    return(PsychIOOSConfigureSerialPort(PsychGetPortIORecord(handle)->device, configString));
}

PsychError IOPORTRead(void)
{
    static char useString[] = "[data, when, errmsg] = IOPort('Read', handle [, blocking=0] [, amount]);";
    static char synopsisString[] =
        "Read data from device, specified by 'handle'.\n"
        "Returned 'data' will be a row vector of read data bytes. 'when' will be a receive "
        "timestamp of when the data read was complete. 'errmsg' will be a human readable "
        "char string with an error message if any error occured, otherwise an empty string. "
        "The optional flag 'blocking' if set to 0 will ask the read function to not block, "
        "but return immediately: The read function will return whatever amount of data is "
        "currently available in the internal input queue, but at most 'amount' bytes if 'amount' "
        "is specified. If no data is available, it will return an empty matrix. This is the default.\n"
        "If 'blocking' is set to 1, you must specify the 'amount' of bytes to receive and the "
        "function will wait until that exact amount of data is available, then return it."
        "Even in blocking mode, the function will return if no data becomes available within "
        "the time period specified by the configuration setting 'ReadTimeout' (see help for "
        "'OpenSerialPort' for description). In some situations, the read function may return "
        "less data than requested, e.g., in specific error cases, where a read could only "
        "get partially satisfied.";

    static char seeAlsoString[] = "'Write', 'OpenSerialPort', 'ConfigureSerialPort'";

    char errmsg[1024];
    int handle, blocking, nread, amount, i;
    psych_uint8* readbuffer;
    double* outbuffer;
    double timestamp;
    errmsg[0] = 0;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(3));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(3));     // The maximum number of outputs

    // Get required port handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);

    // Get optional blocking flag: Defaults to 0 -- non-blocking.
    blocking = 0;
    PsychCopyInIntegerArg(2, kPsychArgOptional, &blocking);

    // Get optional maximum or exact amount to read:
    amount = INT_MAX;
    if (!PsychCopyInIntegerArg(3, kPsychArgOptional, &amount)) {
        // Not spec'd:
        if (blocking > 0) PsychErrorExitMsg(PsychError_user, "When issuing a 'Read' in blocking mode, you must specify the exact 'amount' to read, but 'amount' was omitted!");
    }

    if (amount < 0) PsychErrorExitMsg(PsychError_user, "Invalid (negative) 'amount' of data to read!");

    // Read data:
    nread = PsychReadIOPort(handle, (void**) &readbuffer, amount, blocking, errmsg, &timestamp);

    // Allocate outbuffer of proper size:
    PsychAllocOutDoubleMatArg(1, kPsychArgOptional, 1, ((nread >=0) ? nread : 0), 1, &outbuffer);

    // Copy to outbuffer: We copy our uint8 values to a double matrix. This is arguably a bit of a waste of
    // memory and bandwidth, but it simplifies interfacing with Octave, old Matlab versions and possible
    // future runtime environments a lot:
    if (nread > 0) for (i=0; i<nread; i++) outbuffer[i] = readbuffer[i];

    if (nread < 0 && verbosity > 0) printf("IOPort: Error: %s\n", errmsg);

    // Return timestamp and errmsg, if any:
    PsychCopyOutDoubleArg(2, kPsychArgOptional, timestamp);
    PsychCopyOutCharArg(3, kPsychArgOptional, errmsg);

    return(PsychError_none);
}

PsychError IOPORTWrite(void)
{
    static char useString[] = "[nwritten, when, errmsg, prewritetime, postwritetime, lastchecktime] = IOPort('Write', handle, data [, blocking=1]);";
    static char synopsisString[] =
        "Write data to device, specified by 'handle'.\n"
        "'data' must be a vector of data items to write, or a matrix (in which case data "
        "in the matrix will be transmitted in column-major order, ie., first the first "
        "column, then the 2nd column etc...), either with data elements of uint8 class "
        "or a (1 Byte per char) character string. The optional flag 'blocking' if "
        "set to 0 will ask the write function to not block, but return immediately, ie. "
        "data is sent/written in the background while your code continues to execute - There "
        "may be an arbitrary delay until data transmission is really finished. The default "
        "setting is 1, ie. blocking writes - The function waits until data transmission is really "
        "finished. You can also use blocking == 2 to request a different mode "
        "for blocking writes, where IOPort is polling for write-completion instead of "
        "a more cpu friendly wait. This may decrease latency for certain applications. "
        "Another even more agressive polling method is implemented via blocking == 3 on "
        "Linux systems with some limited set of hardware, e.g., real native serial ports.\n"
        "On systems without any support for specific polling modes, the 2 or 3 settings are treated "
        "as a standard blocking write.\n\n"
        "Optionally, the function returns the following return arguments:\n"
        "'nwritten' Number of bytes written -- Should match amount of data provided on success.\n"
        "'when' A timestamp of write completion: This is only meaningful in blocking mode!\n"
        "'errmsg' A system defined error message if something wen't wrong.\n"
        "The following three timestamps are for low-level debugging and special purpose:\n"
        "'prewritetime' A timestamp taken immediately before submitting the write request. "
        "'postwritetime' A timestamp taken immediately after submitting the write request. "
        "'lastchecktime' A timestamp taken at the time of last check for write completion if applicable. ";

    static char seeAlsoString[] = "";

    char errmsg[1024];
    int handle, blocking, m, n, p, nwritten;
    psych_uint8* inData = NULL;
    char* inChars = NULL;
    void* writedata = NULL;
    double timestamp[4] = {0, 0, 0, 0};
    errmsg[0] = 0;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(3));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(2)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(6));     // The maximum number of outputs

    // Get required port handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);

    // Get the data:
    switch(PsychGetArgType(2)) {
        case PsychArgType_uint8:
            PsychAllocInUnsignedByteMatArg(2, kPsychArgRequired, &m, &n, &p, &inData);
            if (p!=1 || m * n == 0) PsychErrorExitMsg(PsychError_user, "'data' is not a vector or 2D matrix, but some higher dimensional matrix!");
            n = m * n;
            writedata = (void*) inData;
        break;

        case PsychArgType_char:
            PsychAllocInCharArg(2, kPsychArgRequired, &inChars);
            n = strlen(inChars);
            writedata = (void*) inChars;
        break;

        default:
            n = 0;
            PsychErrorExitMsg(PsychError_user, "Invalid type for 'data' vector: Must be an uint8 or char vector.");
            return(PsychError_invalidArg_type);
    }

    // Get optional blocking flag: Defaults to one -- blocking.
    blocking = 1;
    PsychCopyInIntegerArg(3, kPsychArgOptional, &blocking);

    // Write data:
    nwritten = PsychWriteIOPort(handle, writedata, n, blocking, errmsg, &timestamp[0]);
    if (nwritten < 0 && verbosity > 0) printf("IOPort: Error: %s\n", errmsg);

    PsychCopyOutDoubleArg(1, kPsychArgOptional, nwritten);
    PsychCopyOutDoubleArg(2, kPsychArgOptional, timestamp[0]);
    PsychCopyOutCharArg(3, kPsychArgOptional, errmsg);
    PsychCopyOutDoubleArg(4, kPsychArgOptional, timestamp[1]);
    PsychCopyOutDoubleArg(5, kPsychArgOptional, timestamp[2]);
    PsychCopyOutDoubleArg(6, kPsychArgOptional, timestamp[3]);

    return(PsychError_none);
}

PsychError IOPORTBytesAvailable(void)
{
    static char useString[] = "navailable = IOPort('BytesAvailable', handle);";
    static char synopsisString[] = "Return number 'navailable' of data bytes available for readout from device specified by 'handle'.";
    static char seeAlsoString[] = "'Read'";

    int handle;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(1));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    // Get required port handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);

    // Copy out number of bytes available for read on port:
    PsychCopyOutDoubleArg(1, kPsychArgOptional, (double) PsychBytesAvailableIOPort(handle));

    return(PsychError_none);
}

PsychError IOPORTPurge(void)
{
    static char useString[] = "IOPort('Purge', handle);";
    static char synopsisString[] = "Purge all data queued for reading or writing from/to device specified by 'handle'. All unread or unwritten data is discarded.";
    static char seeAlsoString[] = "'Flush'";

    int handle;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(1));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    // Get required port handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);

    PsychPurgeIOPort(handle);

    return(PsychError_none);
}

PsychError IOPORTFlush(void)
{
    static char useString[] = "IOPort('Flush', handle);";
    static char synopsisString[] = "Flush all data queued for writeout to device specified by 'handle', wait for full write completion.";
    static char seeAlsoString[] = "'Flush'";

    int handle;

    // Setup online help:
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()) {PsychGiveHelp(); return(PsychError_none); };

    PsychErrorExit(PsychCapNumInputArgs(1));     // The maximum number of inputs
    PsychErrorExit(PsychRequireNumInputArgs(1)); // The required number of inputs
    PsychErrorExit(PsychCapNumOutputArgs(1));     // The maximum number of outputs

    // Get required port handle:
    PsychCopyInIntegerArg(1, kPsychArgRequired, &handle);

    PsychFlushIOPort(handle);

    return(PsychError_none);
}
