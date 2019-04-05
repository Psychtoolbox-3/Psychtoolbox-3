/*

    Psychtoolbox3/Source/Common/Screen/SCREENComputer.c

    AUTHORS:

        Allen.Ingling@nyu.edu           awi
        mario.kleiner@tuebingen.mpg.de  mk

    PLATFORMS:

        All.

    HISTORY:

        4/01/03     awi     Created.
        10/12/04    awi     Changed "SCREEN" to "Screen" in useString.
        10/21/04    awi     Added comment showing returned struct on OS X.
        1/11/04     awi     Fix bug in reporting memory size, reported by David Jones and uncovered in source by Mario Kleiner.
                            Cleaned up ReportSysctlError.

    DESCRIPTION:

        Returns information about the computer.

*/

#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[]= "comp=Screen('Computer')";
static char synopsisString[] =
        "Get information about the computer. The result is a struct holding information about your computer. "
        "Top-level flags in the returned struct are available on all operating systems and identify the "
        "operating system: 'macintosh', 'windows', 'osx', 'linux'. 'system' contains info about the operating "
        "system version. 'gstreamer' tells if GStreamer support is compiled into Screen(). 'supported' tells if "
        "the given operating system was supported by the given Psychtoolbox at the time of the Psychtoolbox "
        "release, ie., if the given system was tested and shown to work with the given Psychtoolbox, so there "
        "can be some expectation of Psychtoolbox working with this system. This doesn't mean that you will get "
        "support in case of problems though, as support for your system may have ended if it isn't a recent "
        "enough system. 'MACAddress' may contain the MAC address of your ethernet adapter.\n"
        "All other fields or additional sub-structures are platform-dependent, so their presence can not be "
        "relied on.\n"
        "SCREEN 'Computer' not longer supports the  obsolete usage:\n"
        "[model,owner,system,processor,cache,fpu,hz,busHz,vm,pci,emulating]=SCREEN('Computer')\n";
static char seeAlsoString[] = "";

//special includes for sysctl calls
#if PSYCH_SYSTEM == PSYCH_OSX
#include <sys/types.h>
#include <sys/sysctl.h>

#include <pwd.h>
#include <uuid/uuid.h>

//special include for SCDynamicStoreCopySpecific* functions
#include <sys/cdefs.h>
#include <CoreFoundation/CoreFoundation.h>
#include <SystemConfiguration/SCDynamicStore.h>
#include <SystemConfiguration/SCDynamicStoreCopySpecific.h>

//for getting the ethernet MAC address
#include "GetEthernetAddress.h"

// For OSX version query:
#include "PsychCocoaGlue.h"

static void ReportSysctlError(int errorValue)
{
    psych_bool  foundError;
    int sysctlErrors[] = {EFAULT, EINVAL, ENOMEM, ENOTDIR, EISDIR, EOPNOTSUPP, EPERM};
    int i, errorIndex, numSysctlErrors=7;
    char *sysctlErrorStrings[]={"EFAULT", "EINVAL", "ENOMEM", "ENOTDIR", "EISDIR", "EOPNOTSUPP", "EPERM", "UNRECOGNIZED"};

    if(errorValue == 0) return;

    foundError=0;
    errorIndex=7;

    for(i=0; i<numSysctlErrors; i++){
        if(errno==sysctlErrors[i]){
            foundError=1;
            errorIndex=i;
            break;
        }
    }

    PsychErrorExitMsg(PsychError_internal, sysctlErrorStrings[errorIndex]);
}

PsychError SCREENComputer(void) 
{
    const char *majorStructFieldNames[]={"macintosh", "windows", "osx" ,"linux", "kern", "hw", "processUserLongName", 
                                         "processUserShortName", "consoleUserName", "machineName", "localHostName", "location", "MACAddress", "system", "gstreamer", "supported" };
    const char *kernStructFieldNames[]={"ostype", "osrelease", "osrevision", "version","hostname"};
    const char *hwStructFieldNames[]={"machine", "model", "ncpu", "physmem", "usermem", "busfreq", "cpufreq"};
    int numMajorStructDimensions=-1, numKernStructDimensions=-1, numHwStructDimensions=-1;
    int numMajorStructFieldNames=16, numKernStructFieldNames=5, numHwStructFieldNames=7;
    PsychGenericScriptType	*kernStruct, *hwStruct, *majorStruct;
    //char tempStr[CTL_MAXNAME];   //this seems like a bug in Darwin, CTL_MAXNAME is shorter than the longest name.  
    char            tempStr[256], *ethernetMACStr;
    size_t          tempIntSize,  tempStrSize, tempULongIntSize;
    int             mib[2];
    int             tempInt;
    psych_uint64    tempULongInt;
    char            *tempStrPtr;
    CFStringRef     tempCFStringRef;
    psych_bool      stringSuccess;
    int             stringLengthChars, ethernetMACStrSizeBytes;
    long            gestaltResult;
    OSErr           gestaltError;
    int             i,strIndex, bcdDigit, lengthSystemVersionString;
    int             osMajor, osMinor, osBugfix;
    char            systemVersionStr[256];

    //all subfunctions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(0));

    //fill the major struct 
    PsychAllocOutStructArray(1, FALSE, numMajorStructDimensions, numMajorStructFieldNames, majorStructFieldNames, &majorStruct);
    PsychSetStructArrayDoubleElement("macintosh", 0, 0, majorStruct);
    PsychSetStructArrayDoubleElement("windows", 0, 0, majorStruct);
    PsychSetStructArrayDoubleElement("linux", 0, 0, majorStruct);
    PsychSetStructArrayDoubleElement("osx", 0, 1, majorStruct);

    // Official support status:
    PsychSetStructArrayStringElement("supported", 0, (char*) PsychSupportStatus(), majorStruct);

    // GStreamer availability and rough version:
    #if defined(PTB_USE_GSTREAMER)
        #if GST_CHECK_VERSION(1,0,0)
        PsychSetStructArrayDoubleElement("gstreamer", 0, 1 * 10000 + 0 * 100 + 0, majorStruct);
        #else
        PsychSetStructArrayDoubleElement("gstreamer", 0, 0 * 10000 + 10 * 100 + 0, majorStruct);
        #endif
    #else
        PsychSetStructArrayDoubleElement("gstreamer", 0, 0, majorStruct);
    #endif

    //fill the kern struct and implant it within the major struct
    PsychAllocOutStructArray(-1, FALSE, numKernStructDimensions, numKernStructFieldNames, kernStructFieldNames, &kernStruct);
    mib[0]=CTL_KERN;

    mib[1]=KERN_OSTYPE;
    tempStrSize=sizeof(tempStr);
    ReportSysctlError(sysctl(mib, 2, tempStr, &tempStrSize, NULL, 0));
    PsychSetStructArrayStringElement("ostype", 0, tempStr, kernStruct);

    mib[1]=KERN_OSRELEASE;
    tempStrSize=sizeof(tempStr);
    ReportSysctlError(sysctl(mib, 2, tempStr, &tempStrSize, NULL, 0));
    PsychSetStructArrayStringElement("osrelease", 0, tempStr, kernStruct);

    mib[1]=KERN_OSREV;
    tempIntSize=sizeof(tempInt);
    ReportSysctlError(sysctl(mib, 2, &tempInt, &tempIntSize, NULL, 0));
    PsychSetStructArrayDoubleElement("osrevision", 0, (double)tempInt, kernStruct);

    mib[1]=KERN_VERSION;
    tempStrSize=sizeof(tempStr);
    ReportSysctlError(sysctl(mib, 2, tempStr, &tempStrSize, NULL, 0));
    PsychSetStructArrayStringElement("version", 0, tempStr, kernStruct);

    mib[1]=KERN_HOSTNAME;
    tempStrSize=sizeof(tempStr);
    ReportSysctlError(sysctl(mib, 2, tempStr, &tempStrSize, NULL, 0));
    PsychSetStructArrayStringElement("hostname", 0, tempStr, kernStruct);
    PsychSetStructArrayStructElement("kern",0, kernStruct, majorStruct);

    //fill the hw struct and implant it within the major struct
    PsychAllocOutStructArray(-1, FALSE, numHwStructDimensions, numHwStructFieldNames, hwStructFieldNames, &hwStruct);
    mib[0]=CTL_HW;

    mib[1]=HW_MACHINE;
    tempStrSize=sizeof(tempStr);
    ReportSysctlError(sysctl(mib, 2, tempStr, &tempStrSize, NULL, 0));
    PsychSetStructArrayStringElement("machine", 0, tempStr, hwStruct);

    mib[1]=HW_MODEL;
    tempStrSize=sizeof(tempStr);
    ReportSysctlError(sysctl(mib, 2, tempStr, &tempStrSize, NULL, 0));
    PsychSetStructArrayStringElement("model", 0, tempStr, hwStruct);

    mib[1]=HW_NCPU;
    tempIntSize=sizeof(tempInt);
    ReportSysctlError(sysctl(mib, 2, &tempInt, &tempIntSize, NULL, 0));
    PsychSetStructArrayDoubleElement("ncpu", 0, (double)tempInt, hwStruct);

    mib[1]=HW_MEMSIZE;
    tempULongIntSize=sizeof(tempULongInt);
    tempULongInt = 0;
    ReportSysctlError(sysctl(mib, 2, &tempULongInt, &tempULongIntSize, NULL, 0));
    PsychSetStructArrayDoubleElement("physmem", 0, (double)tempULongInt, hwStruct);

    mib[1]=HW_USERMEM;
    tempULongIntSize=sizeof(tempULongInt);
    tempULongInt = 0;
    ReportSysctlError(sysctlbyname("hw.usermem", &tempULongInt, &tempULongIntSize, NULL, 0));
    PsychSetStructArrayDoubleElement("usermem", 0, (double)tempULongInt, hwStruct);

    mib[1]=HW_BUS_FREQ;
    tempULongIntSize=sizeof(tempULongInt);
    tempULongInt = 0;
    ReportSysctlError(sysctlbyname("hw.busfrequency", &tempULongInt, &tempULongIntSize, NULL, 0));
    PsychSetStructArrayDoubleElement("busfreq", 0, (double)tempULongInt, hwStruct);

    mib[1]=HW_CPU_FREQ;
    tempULongIntSize=sizeof(tempULongInt);
    tempULongInt = 0;
    ReportSysctlError(sysctlbyname("hw.cpufrequency", &tempULongInt, &tempULongIntSize, NULL, 0));
    PsychSetStructArrayDoubleElement("cpufreq", 0, (double)tempULongInt, hwStruct);
    PsychSetStructArrayStructElement("hw",0, hwStruct, majorStruct);

    //fill in the process user, console user and machine name in the root struct.
    tempCFStringRef = SCDynamicStoreCopyComputerName(NULL, NULL);
    if (tempCFStringRef) {
        stringLengthChars=(int) CFStringGetMaximumSizeForEncoding(CFStringGetLength(tempCFStringRef), kCFStringEncodingASCII);
        tempStrPtr=malloc(sizeof(char) * (stringLengthChars+1));
        stringSuccess= CFStringGetCString(tempCFStringRef, tempStrPtr, stringLengthChars+1, kCFStringEncodingASCII);

        if(stringSuccess) {
            PsychSetStructArrayStringElement("machineName", 0, tempStrPtr, majorStruct);
        }
        else {
            PsychSetStructArrayStringElement("machineName", 0, "UNKNOWN! QUERY FAILED DUE TO EMPTY OR PROBLEMATIC NAME.", majorStruct);
        }

        free(tempStrPtr);
        CFRelease(tempCFStringRef);
    }
    else {
        PsychSetStructArrayStringElement("machineName", 0, "UNKNOWN! QUERY FAILED DUE TO EMPTY OR PROBLEMATIC NAME.", majorStruct);
    }

    struct passwd* thisUser = getpwuid(getuid());
    if (thisUser) {
        PsychSetStructArrayStringElement("processUserShortName", 0, thisUser->pw_name, majorStruct);
    }
    else {
        PsychSetStructArrayStringElement("processUserShortName", 0, "UNKNOWN! QUERY FAILED DUE TO EMPTY OR PROBLEMATIC NAME.", majorStruct);
    }

    PsychSetStructArrayStringElement("processUserLongName", 0, PsychCocoaGetFullUsername(), majorStruct);

    tempCFStringRef= SCDynamicStoreCopyConsoleUser(NULL, NULL, NULL);
    if (tempCFStringRef) {
        stringLengthChars=(int) CFStringGetMaximumSizeForEncoding(CFStringGetLength(tempCFStringRef), kCFStringEncodingASCII);
        tempStrPtr=malloc(sizeof(char) * (stringLengthChars+1));
        stringSuccess= CFStringGetCString(tempCFStringRef, tempStrPtr, stringLengthChars+1, kCFStringEncodingASCII);

        if(stringSuccess) {
            PsychSetStructArrayStringElement("consoleUserName", 0, tempStrPtr, majorStruct);
        }
        else {
            PsychSetStructArrayStringElement("consoleUserName", 0, "UNKNOWN! QUERY FAILED DUE TO EMPTY OR PROBLEMATIC NAME.", majorStruct);
        }

        free(tempStrPtr);
        CFRelease(tempCFStringRef);
    }
    else {
        PsychSetStructArrayStringElement("consoleUserName", 0, "UNKNOWN! QUERY FAILED DUE TO EMPTY OR PROBLEMATIC NAME.", majorStruct);
    }

    tempCFStringRef= SCDynamicStoreCopyLocalHostName(NULL);
    if (tempCFStringRef) {
        stringLengthChars=(int) CFStringGetMaximumSizeForEncoding(CFStringGetLength(tempCFStringRef), kCFStringEncodingASCII);
        tempStrPtr=malloc(sizeof(char) * (stringLengthChars+1));
        stringSuccess= CFStringGetCString(tempCFStringRef, tempStrPtr, stringLengthChars+1, kCFStringEncodingASCII);
        if(stringSuccess) {
            PsychSetStructArrayStringElement("localHostName", 0, tempStrPtr, majorStruct);
        }
        else {
            PsychSetStructArrayStringElement("localHostName", 0, "UNKNOWN! QUERY FAILED DUE TO EMPTY OR PROBLEMATIC NAME.", majorStruct);
        }

        free(tempStrPtr);
        CFRelease(tempCFStringRef);
    }
    else {
        PsychSetStructArrayStringElement("localHostName", 0, "UNKNOWN! QUERY FAILED DUE TO EMPTY OR PROBLEMATIC NAME.", majorStruct);
    }

    tempCFStringRef= SCDynamicStoreCopyLocation(NULL);
    if (tempCFStringRef) {
        stringLengthChars=(int) CFStringGetMaximumSizeForEncoding(CFStringGetLength(tempCFStringRef), kCFStringEncodingASCII);
        tempStrPtr=malloc(sizeof(char) * (stringLengthChars+1));
        stringSuccess= CFStringGetCString(tempCFStringRef, tempStrPtr, stringLengthChars+1, kCFStringEncodingASCII);
        if(stringSuccess) {
            PsychSetStructArrayStringElement("location", 0, tempStrPtr, majorStruct);
        }
        else {
            PsychSetStructArrayStringElement("location", 0, "UNKNOWN! QUERY FAILED DUE TO EMPTY OR PROBLEMATIC NAME.", majorStruct);
        }

        free(tempStrPtr);
        CFRelease(tempCFStringRef);
    }
    else {
        PsychSetStructArrayStringElement("location", 0, "UNKNOWN! QUERY FAILED DUE TO EMPTY OR PROBLEMATIC NAME.", majorStruct);
    }

    //Add the ethernet MAC address of the primary ethernet interface to the stuct.  This can serve as a unique identifier for the computer.
    ethernetMACStrSizeBytes=GetPrimaryEthernetAddressStringLengthBytes(TRUE)+1;
    ethernetMACStr=(char*) malloc(sizeof(char) * ethernetMACStrSizeBytes);
    GetPrimaryEthernetAddressString(ethernetMACStr, TRUE, TRUE);
    PsychSetStructArrayStringElement("MACAddress", 0, ethernetMACStr, majorStruct);
    free(ethernetMACStr);

    //Add the system version string:
    PsychCocoaGetOSXVersion(&osMajor, &osMinor, &osBugfix);
    sprintf(systemVersionStr, "Mac OS %i.%i.%i", osMajor, osMinor, osBugfix);

    //embed it in the return struct
    PsychSetStructArrayStringElement("system", 0, systemVersionStr, majorStruct);

    return(PsychError_none);
}

#endif

#if PSYCH_SYSTEM == PSYCH_WINDOWS

// M$-Windows implementation of Screen('Computer'): This is very rudimentary for now.
// We only report the operating sytem type (="Windows") and MAC-Address, but don't report any more useful
// information. MAC query does not work yet - We do not have the neccessary libraries to compile the code :(
PsychError SCREENComputer(void)
{
    // Info struct for queries to OS:
    OSVERSIONINFO osvi;
    char versionString[256];

    // const char *majorStructFieldNames[]={"macintosh", "windows", "osx" ,"linux", "kern", "hw", "processUserLongName", 
    //      "processUserShortName", "consoleUserName", "machineName", "localHostName", "location", "MACAddress", "system" };
    const char *majorStructFieldNames[]={"macintosh", "windows", "osx" ,"linux", "system", "IsVistaClass", "gstreamer", "supported"};

    const char *kernStructFieldNames[]={"ostype", "osrelease", "osrevision", "version","hostname"};
    const char *hwStructFieldNames[]={"machine", "model", "ncpu", "physmem", "usermem", "busfreq", "cpufreq"};
    int numMajorStructDimensions=-1, numKernStructDimensions=-1, numHwStructDimensions=-1;
    int numMajorStructFieldNames=8, numKernStructFieldNames=5, numHwStructFieldNames=7;

    PsychGenericScriptType *majorStruct;
    //all subfunctions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(0));

    //fill the major struct 
    PsychAllocOutStructArray(1, FALSE, numMajorStructDimensions, numMajorStructFieldNames, majorStructFieldNames, &majorStruct);
    PsychSetStructArrayDoubleElement("macintosh", 0, 0, majorStruct);
    PsychSetStructArrayDoubleElement("windows", 0, 1, majorStruct);
    PsychSetStructArrayDoubleElement("linux", 0, 0, majorStruct);
    PsychSetStructArrayDoubleElement("osx", 0, 0, majorStruct);

    // Official support status:
    PsychSetStructArrayStringElement("supported", 0, (char*) PsychSupportStatus(), majorStruct);

    // GStreamer availability and rough version:
    #if defined(PTB_USE_GSTREAMER)
        #if GST_CHECK_VERSION(1,0,0)
        PsychSetStructArrayDoubleElement("gstreamer", 0, 1 * 10000 + 0 * 100 + 0, majorStruct);
        #else
        PsychSetStructArrayDoubleElement("gstreamer", 0, 0 * 10000 + 10 * 100 + 0, majorStruct);
        #endif
    #else
        PsychSetStructArrayDoubleElement("gstreamer", 0, 0, majorStruct);
    #endif

    // Query info about Windows version:
    versionString[0]=0;
    memset(&osvi, 0, sizeof(OSVERSIONINFO));
    osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
    GetVersionEx(&osvi);

    // Special case for Windows-10 and later, as GetVersionEx() doesn't report
    // version numbers faithfully beyond Windows 8, unless application manifest
    // would mark the app as Windows-8.1+ compatible. Fake a 10.0 version if this
    // is Windows 10 or later - should be good enough for our purposes.
    if (PsychOSIsMSWin10()) {
        osvi.dwMajorVersion = 10;
        osvi.dwMinorVersion = 0;
    }

    // Convert into string with major.minor.buildnumber - Name of service packs (max 128 chars) etc.:
    // Versions to products: 6.1 = Windows-7, 6.0  = Vista, 5.2 = Windows Server 2003, 5.1 = WindowsXP, 5.0 = Windows 2000, 4.x = NT
    sprintf(versionString, "NT-%i.%i.%i - %s", osvi.dwMajorVersion, osvi.dwMinorVersion, osvi.dwBuildNumber, (char*) osvi.szCSDVersion);
    PsychSetStructArrayStringElement("system", 0, versionString, majorStruct);

    PsychSetStructArrayDoubleElement("IsVistaClass", 0, (PsychIsMSVista() ? 1 : 0), majorStruct);

    return(PsychError_none);
}

#endif

#if PSYCH_SYSTEM == PSYCH_LINUX

// GNU/Linux implementation of Screen('Computer'):
#include <sys/ioctl.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <net/if.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <stdio.h>
// utsname for uname() so we can find out on which kernel we're running:
#include <sys/utsname.h>

PsychError SCREENComputer(void)
{
    //    const char *majorStructFieldNames[]={"macintosh", "windows", "osx" ,"linux", "kern", "hw", "processUserLongName", 
    //    "processUserShortName", "consoleUserName", "machineName", "localHostName", "location", "MACAddress", "system" };
    const char *majorStructFieldNames[]={"macintosh", "windows", "osx" ,"linux", "kern", "hw", "system", "machineName", "localHostName", "MACAddress", "gstreamer", "supported"};

    const char *kernStructFieldNames[]={"ostype", "osrelease", "osrevision", "version","hostname"};
    const char *hwStructFieldNames[]={"machine", "model", "ncpu", "physmem", "usermem", "busfreq", "cpufreq"};
    int numMajorStructDimensions=-1, numKernStructDimensions=-1, numHwStructDimensions=-1;
    int numMajorStructFieldNames=12, numKernStructFieldNames=5, numHwStructFieldNames=7;
    char ethernetMACStr[20];
    struct ifreq devea;
    int s, i;
    PsychGenericScriptType	*kernStruct, *hwStruct, *majorStruct;
    struct utsname unameresult;
    char tmpString[1024];

    //all subfunctions should have these two lines
    PsychPushHelp(useString, synopsisString, seeAlsoString);
    if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    PsychErrorExit(PsychCapNumOutputArgs(1));
    PsychErrorExit(PsychCapNumInputArgs(0));

    //fill the major struct 
    PsychAllocOutStructArray(1, FALSE, numMajorStructDimensions, numMajorStructFieldNames, majorStructFieldNames, &majorStruct);
    PsychSetStructArrayDoubleElement("macintosh", 0, 0, majorStruct);
    PsychSetStructArrayDoubleElement("windows", 0, 0, majorStruct);
    PsychSetStructArrayDoubleElement("linux", 0, 1, majorStruct);
    PsychSetStructArrayDoubleElement("osx", 0, 0, majorStruct);

    // Official support status:
    PsychSetStructArrayStringElement("supported", 0, (char*) PsychSupportStatus(), majorStruct);

    // GStreamer availability and rough version:
    #if defined(PTB_USE_GSTREAMER)
        #if GST_CHECK_VERSION(1,0,0)
        PsychSetStructArrayDoubleElement("gstreamer", 0, 1 * 10000 + 0 * 100 + 0, majorStruct);
        #else
        PsychSetStructArrayDoubleElement("gstreamer", 0, 0 * 10000 + 10 * 100 + 0, majorStruct);
        #endif
    #else
        PsychSetStructArrayDoubleElement("gstreamer", 0, 0, majorStruct);
    #endif

    // Query hardware MAC address of primary ethernet interface: This is a unique id of the computer,
    // good enough to disambiguate our statistics:
    sprintf(ethernetMACStr, "00:00:00:00:00:00");

    s = socket(PF_INET, SOCK_DGRAM, 0);
    if (s >= 0) {
        // Probe eth0 to eth99 for a valid network interface:
        for (i = 0; i < 100; i++) {
            sprintf(devea.ifr_name, "eth%d", i);
            if (ioctl(s, SIOCGIFHWADDR, &devea) >= 0) {
                // Success: Use its MAC address and be done:
                sprintf(ethernetMACStr, "%.2X:%.2X:%.2X:%.2X:%.2X:%.2X",
                        devea.ifr_ifru.ifru_hwaddr.sa_data[0]&0xff, devea.ifr_ifru.ifru_hwaddr.sa_data[1]&0xff,
                        devea.ifr_ifru.ifru_hwaddr.sa_data[2]&0xff, devea.ifr_ifru.ifru_hwaddr.sa_data[3]&0xff,
                        devea.ifr_ifru.ifru_hwaddr.sa_data[4]&0xff, devea.ifr_ifru.ifru_hwaddr.sa_data[5]&0xff);
                break;
            }
        }
        close(s);
    }

    PsychSetStructArrayStringElement("MACAddress", 0, ethernetMACStr, majorStruct);

    uname(&unameresult);
    sprintf(tmpString, "Linux %s", unameresult.release);
    PsychSetStructArrayStringElement("system", 0, tmpString, majorStruct);

    // Fill the kern struct and implant it within the major struct
    PsychAllocOutStructArray(-1, FALSE, numKernStructDimensions, numKernStructFieldNames, kernStructFieldNames, &kernStruct);
    PsychSetStructArrayStringElement("ostype", 0, unameresult.sysname, kernStruct);
    PsychSetStructArrayStringElement("osrelease", 0, unameresult.release, kernStruct);
    PsychSetStructArrayDoubleElement("osrevision", 0, 0, kernStruct);
    sprintf(tmpString, "%s Kernel Version %s: %s", unameresult.sysname, unameresult.release, unameresult.version);
    PsychSetStructArrayStringElement("version", 0, tmpString, kernStruct);
    PsychSetStructArrayStringElement("hostname", 0, unameresult.nodename, kernStruct);
    PsychSetStructArrayStructElement("kern",0, kernStruct, majorStruct);

    PsychSetStructArrayStringElement("localHostName", 0, unameresult.nodename, majorStruct);
    PsychSetStructArrayStringElement("machineName", 0, unameresult.nodename, majorStruct);

    // Fill the hw struct and implant it within the major struct:
    PsychAllocOutStructArray(-1, FALSE, numHwStructDimensions, numHwStructFieldNames, hwStructFieldNames, &hwStruct);
    PsychSetStructArrayStringElement("machine", 0, unameresult.machine, hwStruct);
    PsychSetStructArrayStringElement("model", 0, "", hwStruct);
    PsychSetStructArrayDoubleElement("ncpu", 0, (double) sysconf(_SC_NPROCESSORS_ONLN), hwStruct);
    PsychSetStructArrayDoubleElement("physmem", 0, ((double) sysconf(_SC_PHYS_PAGES)) * ((double) sysconf(_SC_PAGESIZE)), hwStruct);
    PsychSetStructArrayDoubleElement("usermem", 0, ((double) sysconf(_SC_AVPHYS_PAGES)) * ((double) sysconf(_SC_PAGESIZE)), hwStruct);
    PsychSetStructArrayDoubleElement("busfreq", 0, 0, hwStruct);
    PsychSetStructArrayDoubleElement("cpufreq", 0, 0, hwStruct);
    PsychSetStructArrayStructElement("hw",0, hwStruct, majorStruct);

    return(PsychError_none);
}

#endif
