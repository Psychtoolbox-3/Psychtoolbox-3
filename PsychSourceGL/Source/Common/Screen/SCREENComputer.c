/*

	Psychtoolbox3/Source/Common/Screen/SCREENComputer.c		

	AUTHORS:

		Allen.Ingling@nyu.edu			awi 
		mario.kleiner@tuebingen.mpg.de	mk
  
	PLATFORMS:	

		All.

	HISTORY:

		4/01/03		awi		Created. 
		10/12/04	awi		Changed "SCREEN" to "Screen" in useString.
		10/21/04	awi		Added comment showing returned struct on OS X.
		1/11/04		awi		Fix bug in reporting memory size, reported by David Jones and uncovered in source by Mario Kleiner.
							Cleaned up ReportSysctlError.

	DESCRIPTION:

		Returns information about the computer. 

	TO DO:

		Macro out the parts specfic to OS X or abstract them.
		Add OS X and Linux flags to the OS 9 and Windows version of Screen 'Computer'

		Migrate this out of Screen.  It makes no sense for it to be here. Rename it to 
		"ComputerInfo".  We already have "FontInfo".  Generally, Psychtoolbox functions
		which return a struct with a bunch of info should be named *Info.     

*/

/* This is the struct returned by the OS9 version of SCREEN 'Computer'.

	macintosh: 1
	windows: 0
	emulating: ''
	pci: 1
	vm: ''
	busHz: 133216625
	hz: 999999997
	fpu: ''
	cache: ''
	processor: 'G4'
	system: 'Mac OS 9.2.2'
	owner: ''
	model: 'Classic Mac OS Compatibility/999'

*/

/*

from sys/time.h

struct clockinfo {
	int	hz;		// clock frequency
	int	tick;		// micro-seconds per hz tick
	int	tickadj;	// clock skew rate for adjtime() 
	int	stathz;		// statistics clock frequency 
	int	profhz;		// profiling clock frequency 
};

*/

#include "Screen.h"

// If you change the useString then also change the corresponding synopsis string in ScreenSynopsis.c
static char useString[]= "comp=Screen('Computer')";
static char synopsisString[] =
        "Get information about the computer.  The result is a struct holding information about your computer. "
        "Top-level flags in the returned struct are available on all operating systems and identify the operating "
        "operating system: 'macintosh', 'windows', 'osx'.  All other fields in the returned struct are platform-dependent. \n"
        "\n"
        "OS X: results contains a hierarchial struct with major and minor fields names as with BSD's sysctl(3) MIB fields. \n"
        "\n"
        "SCREEN 'Computer' not longer supports the  obsolete usage: \n"
        "[model,owner,system,processor,cache,fpu,hz,busHz,vm,pci,emulating]=SCREEN('Computer')\n";
static char seeAlsoString[] = "";

//special includes for sysctl calls
#if PSYCH_SYSTEM == PSYCH_OSX
#include <sys/types.h>
#include <sys/sysctl.h>

//special include for SCDynamicStoreCopySpecific* functions
#include <sys/cdefs.h>
#include <CoreFoundation/CoreFoundation.h>
#include <SystemConfiguration/SCDynamicStore.h>
#include <SystemConfiguration/SCDynamicStoreCopySpecific.h>

//for getting the ethernet MAC address
#include "GetEthernetAddress.h"

static void ReportSysctlError(int errorValue)
{
	psych_bool	foundError;
    int sysctlErrors[]={EFAULT, EINVAL, ENOMEM, ENOTDIR, EISDIR, EOPNOTSUPP, EPERM};
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
	                                     "processUserShortName", "consoleUserName", "machineName", "localHostName", "location", "MACAddress", "system", "gstreamer" };
    const char *kernStructFieldNames[]={"ostype", "osrelease", "osrevision", "version","hostname"};
    const char *hwStructFieldNames[]={"machine", "model", "ncpu", "physmem", "usermem", "busfreq", "cpufreq"};
    int numMajorStructDimensions=1, numKernStructDimensions=1, numHwStructDimensions=1;
    int numMajorStructFieldNames=15, numKernStructFieldNames=5, numHwStructFieldNames=7;
    PsychGenericScriptType	*kernStruct, *hwStruct, *majorStruct;
    //char tempStr[CTL_MAXNAME];   //this seems like a bug in Darwin, CTL_MAXNAME is shorter than the longest name.  
    char						tempStr[256], *ethernetMACStr;
    size_t						tempIntSize,  tempStrSize, tempULongIntSize; 	
	int							mib[2];
	int							tempInt;
	psych_uint64                tempULongInt;
	char						*tempStrPtr;
	CFStringRef					tempCFStringRef;
	psych_bool						stringSuccess;
	int							stringLengthChars, ethernetMACStrSizeBytes;
	long						gestaltResult;
	OSErr						gestaltError;
	int							i,strIndex, bcdDigit, lengthSystemVersionString;
	SInt32						osMajor, osMinor, osBugfix;
	char						systemVersionStr[256];
	
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

    // GStreamer availability:
    #if defined(PTB_USE_GSTREAMER)
    PsychSetStructArrayDoubleElement("gstreamer", 0, 1, majorStruct);
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
    long long tempLongInt;
    tempULongIntSize=sizeof(tempLongInt);
    ReportSysctlError(sysctl(mib, 2, &tempLongInt, &tempULongIntSize, NULL, 0));
    PsychSetStructArrayDoubleElement("physmem", 0, (double)tempLongInt, hwStruct);

    mib[1]=HW_USERMEM;
    tempULongIntSize=sizeof(tempULongInt);
    ReportSysctlError(sysctl(mib, 2, &tempULongInt, &tempULongIntSize, NULL, 0));
    PsychSetStructArrayDoubleElement("usermem", 0, (double)tempULongInt, hwStruct);

    mib[1]=HW_BUS_FREQ;
    tempULongIntSize=sizeof(tempULongInt);
    ReportSysctlError(sysctl(mib, 2, &tempULongInt, &tempULongIntSize, NULL, 0));
    PsychSetStructArrayDoubleElement("busfreq", 0, (double)tempULongInt, hwStruct);

    mib[1]=HW_CPU_FREQ;
    tempULongIntSize=sizeof(tempULongInt);
    ReportSysctlError(sysctl(mib, 2, &tempULongInt, &tempULongIntSize, NULL, 0));
    PsychSetStructArrayDoubleElement("cpufreq", 0, (double)tempULongInt, hwStruct);
    PsychSetStructArrayStructElement("hw",0, hwStruct, majorStruct);

    //fill in the process user, console user and machine name in the root struct.
	tempCFStringRef= CSCopyMachineName();
	if (tempCFStringRef) {
		stringLengthChars=(int) CFStringGetMaximumSizeForEncoding(CFStringGetLength(tempCFStringRef), kCFStringEncodingUTF8);
		tempStrPtr=malloc(sizeof(char) * (stringLengthChars+1));
		stringSuccess= CFStringGetCString(tempCFStringRef, tempStrPtr, stringLengthChars+1, kCFStringEncodingUTF8);
		
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
	
	tempCFStringRef= CSCopyUserName(TRUE); //use short name
	if (tempCFStringRef) {		
		stringLengthChars=(int) CFStringGetMaximumSizeForEncoding(CFStringGetLength(tempCFStringRef), kCFStringEncodingUTF8);
		tempStrPtr=malloc(sizeof(char) * (stringLengthChars+1));
		stringSuccess= CFStringGetCString(tempCFStringRef, tempStrPtr, stringLengthChars+1, kCFStringEncodingUTF8);
		if(stringSuccess) {
            PsychSetStructArrayStringElement("processUserShortName", 0, tempStrPtr, majorStruct);
        }
        else {
            PsychSetStructArrayStringElement("processUserShortName", 0, "UNKNOWN! QUERY FAILED DUE TO EMPTY OR PROBLEMATIC NAME.", majorStruct);
        }
		
		free(tempStrPtr);
		CFRelease(tempCFStringRef);
	}
	else {
		PsychSetStructArrayStringElement("processUserShortName", 0, "UNKNOWN! QUERY FAILED DUE TO EMPTY OR PROBLEMATIC NAME.", majorStruct);
	}

	tempCFStringRef= CSCopyUserName(FALSE); //use long name
	if (tempCFStringRef) {		
		stringLengthChars=(int) CFStringGetMaximumSizeForEncoding(CFStringGetLength(tempCFStringRef), kCFStringEncodingUTF8);
		tempStrPtr=malloc(sizeof(char) * (stringLengthChars+1));
		stringSuccess= CFStringGetCString(tempCFStringRef, tempStrPtr, stringLengthChars+1, kCFStringEncodingUTF8);
		if(stringSuccess) {
            PsychSetStructArrayStringElement("processUserLongName", 0, tempStrPtr, majorStruct);
        }
        else {
            PsychSetStructArrayStringElement("processUserLongName", 0, "UNKNOWN! QUERY FAILED DUE TO EMPTY OR PROBLEMATIC NAME.", majorStruct);
        }

		free(tempStrPtr);
		CFRelease(tempCFStringRef);
	}
	else {
		PsychSetStructArrayStringElement("processUserLongName", 0, "UNKNOWN! QUERY FAILED DUE TO EMPTY OR PROBLEMATIC NAME.", majorStruct);
	}
	
	tempCFStringRef= SCDynamicStoreCopyConsoleUser(NULL, NULL, NULL);
	if (tempCFStringRef) {
		stringLengthChars=(int) CFStringGetMaximumSizeForEncoding(CFStringGetLength(tempCFStringRef), kCFStringEncodingUTF8);
		tempStrPtr=malloc(sizeof(char) * (stringLengthChars+1));
		stringSuccess= CFStringGetCString(tempCFStringRef, tempStrPtr, stringLengthChars+1, kCFStringEncodingUTF8);
		
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
		stringLengthChars=(int) CFStringGetMaximumSizeForEncoding(CFStringGetLength(tempCFStringRef), kCFStringEncodingUTF8);
		tempStrPtr=malloc(sizeof(char) * (stringLengthChars+1));
		stringSuccess= CFStringGetCString(tempCFStringRef, tempStrPtr, stringLengthChars+1, kCFStringEncodingUTF8);
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
		stringLengthChars=(int) CFStringGetMaximumSizeForEncoding(CFStringGetLength(tempCFStringRef), kCFStringEncodingUTF8);
		tempStrPtr=malloc(sizeof(char) * (stringLengthChars+1));
		stringSuccess= CFStringGetCString(tempCFStringRef, tempStrPtr, stringLengthChars+1, kCFStringEncodingUTF8);
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
	Gestalt(gestaltSystemVersionMajor, &osMajor);
	Gestalt(gestaltSystemVersionMinor, &osMinor);
	Gestalt(gestaltSystemVersionBugFix, &osBugfix);
	
	sprintf(systemVersionStr, "Mac OS %i.%i.%i", osMajor, osMinor, osBugfix);

	//embed it in the return struct
	PsychSetStructArrayStringElement("system", 0, systemVersionStr, majorStruct);

/*
	OLD DEAD Implementation, left for now as a reference...
	//Add the system version string:
	gestaltError=Gestalt(gestaltSystemVersion, &gestaltResult);

	//The result is a four-digit value stored in BCD in the lower 16-bits  of the result.  There are implicit decimal
	// points between the last three digis.  For example Mac OS 10.3.6 is:
	//
	//  0000 0000 0000 0000 0001 0000 0011 0110
	//                         1    0    3    6
	//                         1    0.   3.   6

	strIndex=0;
	//4th digit.
	bcdDigit=gestaltResult & 15;
	gestaltResult= gestaltResult>>4;
	strIndex=strIndex+sprintf(systemVersionStr+strIndex, "%i", bcdDigit);

	//decimal point
	strIndex=strIndex+sprintf(systemVersionStr+strIndex, "%s", ".");

	//3rd digit
	bcdDigit=gestaltResult & 15;
	gestaltResult= gestaltResult>>4;
	strIndex=strIndex+sprintf(systemVersionStr+strIndex, "%i", bcdDigit);

	//decimal point
	strIndex=strIndex+sprintf(systemVersionStr+strIndex, "%s", ".");

	//second digit
	//2nd digit.
	bcdDigit=gestaltResult & 15;
	gestaltResult= gestaltResult>>4;
	strIndex=strIndex+sprintf(systemVersionStr+strIndex, "%i", bcdDigit);

	//1st digit
	bcdDigit=gestaltResult & 15;
	gestaltResult= gestaltResult>>4;
	strIndex=strIndex+sprintf(systemVersionStr+strIndex, "%i", bcdDigit);

	//preface with "Mac OS "  
	strIndex=strIndex+sprintf(systemVersionStr+strIndex, "%s", " SO caM");

	//reverse to make it forward
	lengthSystemVersionString=strlen(systemVersionStr);
	for(i=0;i<lengthSystemVersionString;i++){
		systemVersionStrForward[lengthSystemVersionString-1-i]=systemVersionStr[i];
	}

	systemVersionStrForward[lengthSystemVersionString]='\0';
	//embed it in the return struct
	PsychSetStructArrayStringElement("system", 0, systemVersionStrForward, majorStruct);
*/

    return(PsychError_none);
}

#endif

#if PSYCH_SYSTEM == PSYCH_WINDOWS

//#include <Windows.h>
//#include <rpc.h>
//#include <rpcdce.h>
//#pragma comment(lib, "rpcrt4.lib")

// M$-Windows implementation of Screen('Computer'): This is very rudimentary for now.
// We only report the operating sytem type (="Windows") and MAC-Address, but don't report any more useful
// information. MAC query does not work yet - We do not have the neccessary libraries to compile the code :(
PsychError SCREENComputer(void)
{
	// Info struct for queries to OS:
	OSVERSIONINFO osvi;
	char	versionString[256];
	
    // const char *majorStructFieldNames[]={"macintosh", "windows", "osx" ,"linux", "kern", "hw", "processUserLongName", 
    //      "processUserShortName", "consoleUserName", "machineName", "localHostName", "location", "MACAddress", "system" };
    const char *majorStructFieldNames[]={"macintosh", "windows", "osx" ,"linux", "system", "MACAddress"};

    const char *kernStructFieldNames[]={"ostype", "osrelease", "osrevision", "version","hostname"};
    const char *hwStructFieldNames[]={"machine", "model", "ncpu", "physmem", "usermem", "busfreq", "cpufreq"};
    int numMajorStructDimensions=1, numKernStructDimensions=1, numHwStructDimensions=1;
    int numMajorStructFieldNames=5, numKernStructFieldNames=5, numHwStructFieldNames=7;
    // char ethernetMACStr[20];
    // unsigned char MACData[6];
    // UUID uuid;
    int i;

    PsychGenericScriptType	*kernStruct, *hwStruct, *majorStruct;
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

	// Query info about Windows version:
	versionString[0]=0;
	memset(&osvi, 0, sizeof(OSVERSIONINFO));
	osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
	GetVersionEx(&osvi);

	// Convert into string with major.minor.buildnumber - Name of service packs (max 128 chars) etc.:
	// Versions to products: 6.1 = Windows-7, 6.0  = Vista, 5.2 = Windows Server 2003, 5.1 = WindowsXP, 5.0 = Windows 2000, 4.x = NT
	sprintf(versionString, "%i.%i.%i - %s", osvi.dwMajorVersion, osvi.dwMinorVersion, osvi.dwBuildNumber, (char*) osvi.szCSDVersion);
	PsychSetStructArrayStringElement("system", 0, versionString, majorStruct);

    
    // Query hardware MAC address of primary ethernet interface: This is a unique id of the computer,
    // good enough to disambiguate our statistics:
    // sprintf(ethernetMACStr, "00:00:00:00:00:00");
   
    // Ask OS to create UUID. Windows uses the MAC address of primary interface
    // in bytes 2 to 7  to do this:
    //UuidCreateSequential( &uuid );    // Ask OS to create UUID
    //for (i=2; i<8; i++) MACData[i - 2] = uuid.Data4[i];
    //sprintf(ethernetMACStr, "%.2X:%.2X:%.2X:%.2X:%.2X:%.2X",
	 //   MACData[0] & 0xff, MACData[1] & 0xff,
	 //   MACData[2] & 0xff, MACData[3] & 0xff,
	 //   MACData[4] & 0xff, MACData[5] & 0xff);
	    	    	    
    // PsychSetStructArrayStringElement("MACAddress", 0, ethernetMACStr, majorStruct);

    return(PsychError_none);
}

#endif

#if PSYCH_SYSTEM == PSYCH_LINUX

// GNU/Linux implementation of Screen('Computer'): This is very rudimentary for now.
// We only report the operating sytem type (="Linux") and the MAC ethernet address of
// the primary interface, but don't report any more useful information.

#include <sys/ioctl.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <net/if.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <stdio.h>

PsychError SCREENComputer(void)
{
    //    const char *majorStructFieldNames[]={"macintosh", "windows", "osx" ,"linux", "kern", "hw", "processUserLongName", 
    //    "processUserShortName", "consoleUserName", "machineName", "localHostName", "location", "MACAddress", "system" };
    const char *majorStructFieldNames[]={"macintosh", "windows", "osx" ,"linux", "MACAddress"};

    const char *kernStructFieldNames[]={"ostype", "osrelease", "osrevision", "version","hostname"};
    const char *hwStructFieldNames[]={"machine", "model", "ncpu", "physmem", "usermem", "busfreq", "cpufreq"};
    int numMajorStructDimensions=1, numKernStructDimensions=1, numHwStructDimensions=1;
    int numMajorStructFieldNames=5, numKernStructFieldNames=5, numHwStructFieldNames=7;
    char ethernetMACStr[20];
    struct ifreq devea;
    int s;
    PsychGenericScriptType	*kernStruct, *hwStruct, *majorStruct;

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

    // Query hardware MAC address of primary ethernet interface: This is a unique id of the computer,
    // good enough to disambiguate our statistics:
    sprintf(ethernetMACStr, "00:00:00:00:00:00");

    s = socket(PF_INET, SOCK_DGRAM, 0);
    if (s>=0) {
      strcpy(devea.ifr_name, "eth0");
      if (ioctl(s, SIOCGIFHWADDR, &devea) >= 0) {
	sprintf(ethernetMACStr, "%.2X:%.2X:%.2X:%.2X:%.2X:%.2X",
		devea.ifr_ifru.ifru_hwaddr.sa_data[0]&0xff, devea.ifr_ifru.ifru_hwaddr.sa_data[1]&0xff,
		devea.ifr_ifru.ifru_hwaddr.sa_data[2]&0xff, devea.ifr_ifru.ifru_hwaddr.sa_data[3]&0xff,
		devea.ifr_ifru.ifru_hwaddr.sa_data[4]&0xff, devea.ifr_ifru.ifru_hwaddr.sa_data[5]&0xff);
      }

      close(s);
    }

    PsychSetStructArrayStringElement("MACAddress", 0, ethernetMACStr, majorStruct);
    
    return(PsychError_none);
}

#endif
