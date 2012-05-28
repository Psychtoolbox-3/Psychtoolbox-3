/*

  Psychtoolbox3/Source/Common/SCREENNull.c		

  AUTHORS:

  Allen.Ingling@nyu.edu		awi 

  PLATFORMS:	

  This file should build on any platform. 

  HISTORY:

  11/26/02  awi		Created. 

  DESCRIPTION:

  Arbitrary code to be executed for debugging/testing/hacking PTB. This subfunction is not
  documented or mentioned in official docs and only here for the purpose of internal testing
  and as a container for quick & dirty hacks and tests by the core developers.

  Therefore its syntax and semantics can change in any way, anytime without warning.

*/

#include "Screen.h"

#include "PsychGraphicsCardRegisterSpecs.h"

#if PSYCH_SYSTEM == PSYCH_OSX
#include <IOKit/IOKitLib.h>
#include <sys/time.h>
#endif

static char useString[] = "[[value1], [value2] ..]=SCREEN('Null',[value1],[value2],...);";
static char synopsisString[] = 
	"Special test function for Psychtoolbox testing and developments. Normal users: DO NOT USE!";
static char seeAlsoString[] = "";

PsychError SCREENNull(void) 
{
#define RADEON_D1CRTC_INTERRUPT_CONTROL 0x60DC

//#define RADEON_R500_GEN_INT_CNTL   0x100
//#define RADEON_R500_GEN_INT_STATUS 0x104
//#define RADEON_R500_GEN_INT_CNTL   0x040
//#define RADEON_R500_GEN_INT_STATUS 0x044
//#define RADEON_R500_GEN_INT_CNTL   0x200
//#define RADEON_R500_GEN_INT_STATUS 0x204

	const double defaultMatrix[] = {1.1, 1.2, 1.3, 1.4, 2.1, 2.2, 2.3, 2.4};
	const double defaultM=2, defaultN=4; 
	double tempValue; 
	double *array;
	int i, m,n, p, numInArgs, numOutArgs, numNamedOutArgs;
	char *str;
	PsychArgFormatType format;
	PsychWindowRecordType *windowRecord;
	const char defaultString[] = "I am the default string\n";

	//all sub functions should have these two lines
	PsychPushHelp(useString, synopsisString, seeAlsoString);
	if(PsychIsGiveHelp()){PsychGiveHelp();return(PsychError_none);};

    #if PSYCH_SYSTEM != PSYCH_WINDOWS
        // Test GPU low-level dithering control:
        int screenId, ditherEnable;
		PsychCopyInIntegerArg(1, TRUE, &screenId);
		PsychCopyInIntegerArg(2, TRUE, &ditherEnable);

		if (PsychPrefStateGet_Verbosity() > 2) printf("SetDithering: ScreenId %i : DitherSetting %i.\n", screenId, ditherEnable);
        PsychSetOutputDithering(NULL, screenId, (unsigned int) ditherEnable);


        unsigned int crtcoff[6] = { EVERGREEN_CRTC0_REGISTER_OFFSET, EVERGREEN_CRTC1_REGISTER_OFFSET, EVERGREEN_CRTC2_REGISTER_OFFSET, EVERGREEN_CRTC3_REGISTER_OFFSET, EVERGREEN_CRTC4_REGISTER_OFFSET, EVERGREEN_CRTC5_REGISTER_OFFSET };

//        printf("Pre:  AVIVO_DC_LUTA_BLACK_OFFSET_BLUE = %i\n", PsychOSKDReadRegister(screenId, AVIVO_DC_LUTA_BLACK_OFFSET_BLUE, NULL));
/*
        PsychOSKDWriteRegister(screenId, AVIVO_DC_LUTA_BLACK_OFFSET_BLUE + 0x0000, 0, NULL);
        PsychOSKDWriteRegister(screenId, AVIVO_DC_LUTA_BLACK_OFFSET_GREEN + 0x0000, 0, NULL);
        PsychOSKDWriteRegister(screenId, AVIVO_DC_LUTA_BLACK_OFFSET_RED + 0x0000, 0, NULL);

        PsychOSKDWriteRegister(screenId, AVIVO_DC_LUTA_WHITE_OFFSET_BLUE + 0x0000, 0xffff, NULL);
        PsychOSKDWriteRegister(screenId, AVIVO_DC_LUTA_WHITE_OFFSET_GREEN + 0x0000, 0xffff, NULL);
        PsychOSKDWriteRegister(screenId, AVIVO_DC_LUTA_WHITE_OFFSET_RED + 0x0000, 0xffff, NULL);

        PsychOSKDWriteRegister(screenId, AVIVO_DC_LUT_RW_INDEX, 0, NULL);        
        for (i=0; i < 256; i++) {
            m = i << 1;
            PsychOSKDWriteRegister(screenId, AVIVO_DC_LUT_30_COLOR, (m << 20) | (m << 10) | (m << 0), NULL);
        }
*/        
//        printf("Post: AVIVO_DC_LUTA_BLACK_OFFSET_BLUE = %i\n", PsychOSKDReadRegister(screenId, AVIVO_DC_LUTA_BLACK_OFFSET_BLUE, NULL));

        for (i=0; i < 6; i++) printf("EVERGREEN_DC_LUT_BLACK_OFFSET_BLUE[%i] = %i\n", i, PsychOSKDReadRegister(screenId, EVERGREEN_DC_LUT_BLACK_OFFSET_BLUE + crtcoff[i], NULL));

        
        return(PsychError_none);
    #endif

	#if PSYCH_SYSTEM == PSYCH_LINUX
/*		PsychAllocInWindowRecordArg(1, TRUE, &windowRecord);

		// Return current OpenML kPsychOpenMLDefective state:
		PsychCopyOutDoubleArg(1, FALSE, (int) (windowRecord->specialflags & kPsychOpenMLDefective));

		// Enable/Disable current OpenML workaround2:
		if (PsychCopyInIntegerArg(2, FALSE, &i)) {
			if (i > 0) windowRecord->specialflags |= kPsychOpenMLDefective;
			if (i == 0) windowRecord->specialflags &= ~kPsychOpenMLDefective;
		}

		if (PsychCopyInIntegerArg(3, FALSE, &i)) {
			windowRecord->specialflags = i;
		}
*/

		unsigned int regOffset, value, hi, lo;
		PsychCopyInIntegerArg(1, TRUE, &regOffset);
		value = PsychOSKDReadRegister(0, regOffset, NULL);

		hi = value >> 16;
		lo = value & 0xffff;

		if (PsychPrefStateGet_Verbosity() > 2) printf("%p :: hi = %i , lo = %i , val = %i\n", (void*) (size_t) regOffset, hi, lo, value);

		PsychCopyOutDoubleArg(1, FALSE, (double) value);
		PsychCopyOutDoubleArg(2, FALSE, (double) hi);
		PsychCopyOutDoubleArg(3, FALSE, (double) lo);

		return(PsychError_none);
	#endif
	

	#if PSYCH_SYSTEM == PSYCH_WINDOWS
		PsychAllocInWindowRecordArg(1, TRUE, &windowRecord);
		PsychOSFlipWindowBuffers(windowRecord);
		return(PsychError_none);
	#endif

	#if PSYCH_SYSTEM == PSYCH_OSX

// ----------- Linux DRI2/Radeon-KMS Scanline-Timestamping prototype for R500 core.

	#define RADEON_D1CRTC_STATUS_POSITION  0x60a0
	#define RADEON_D2CRTC_OFFSET           0x0800
	#define RADEON_D2CRTC_STATUS_POSITION  0x68a0
	#define RADEON_VBEAMPOSITION_BITMASK   0x1fff
	#define RADEON_HBEAMPOSITION_BITSHIFT  16
	#define AVIVO_D1CRTC_V_BLANK_START_END 0x6024
	#define AVIVO_D1CRTC_V_TOTAL           0x6020
	#define AVIVO_D1CRTC_H_TOTAL           0x6000
	
	int crtcid = 0, verbose = 0;
	PsychCopyInIntegerArg(1, FALSE, &crtcid);
	PsychCopyInIntegerArg(2, FALSE, &verbose);
	
	const int pedantic = 1;
	int vbl_start, vbl_end, vbl, htotal, vtotal, dotclock;
	psych_int64 linedur_ns, pixeldur_ns;
	psych_bool invbl = TRUE;
	int scanline, scanpixel, vblcount;
	unsigned int cardId;
	
	psych_uint32 crtco = (crtcid > 0) ? RADEON_D2CRTC_OFFSET : 0;
	struct timeval raw_time;
	struct timeval delta_time;	
	struct timeval vblank_time;

	// NVidia GPU? Check OpenGL vendor string and for recognized gpu-id:
	if ((strstr((char*) glGetString(GL_VENDOR), "NVIDIA")) && (PsychGetNVidiaGPUType(NULL) > 0)) {
		// Yes. Test beampos queries and vblank counter:
		if ((cardId = PsychGetNVidiaGPUType(NULL)) >= 0x50) {
			// NV50 or later, aka GeForce-8000 or later: Same registers for
			// everything shipping from NV50 up to at least NVC0 "Fermi":
			
			// Offset between crtc's is 0x800:
			crtco = (crtcid > 0) ? 0x800 : 0;
			
			// Lower 16 bits are horizontal scanout position, upper 16 bits are always zero:
			scanpixel = (PsychOSKDReadRegister(crtcid, 0x616344 + crtco, NULL)) & 0xFFFF;

			// Lower 16 bits are vertical scanout position (scanline), upper 16 bits are vblank counter:
			vblcount = (PsychOSKDReadRegister(crtcid, 0x616340 + crtco, NULL));
			scanline = vblcount & 0xFFFF;

			vblcount = (vblcount >> 16) & 0xFFFF;

			if (verbose > 1 || verbose < 0) {
//				printf("%i\n", (PsychOSKDReadRegister(crtcid, 0x616340 + crtco + 4 * verbose, NULL)));
				long int i;
				
				// crtc stride for CRTC control block (CRTC_VAL at offset 0xa00, stride 0x540)
				crtco = (crtcid > 0) ? 0x540 : 0;

//				for (i=0; i<=0xfffffc; i+=4) {
				for (i=0; i<0x540; i+=4) {
//					unsigned int base = 0x610000; // Display on NV50
//no					unsigned int base = 0x600000;
//					unsigned int base = 0x00008000;
					//crtco = 0x540;
					unsigned long int base = 0x610000 + 0xa00 + crtco;
					vblcount = PsychOSKDReadRegister(crtcid, base + i, NULL);
					unsigned int hi = vblcount >> 16;
					unsigned int lo = vblcount & 0xffff;
//					if ((vblcount > 0 && vblcount < 50) || (hi > 0 && hi < 50) || (hi> 750 && hi < 850) || (lo > 0 && lo < 50) || (lo> 750 && lo < 850)) {
//					if (abs((int)(hi + lo) - 38) < 5) {
if (1) {
						printf("%p :: hi = %i , lo = %i , val = %i\n", (void*) (base + i), hi, lo, vblcount);
					}
				}
//				printf("%i\n", (PsychOSKDReadRegister(crtcid, 0x0068080c + crtco, NULL)));
			}
			// SYNC_START_TO_BLANK_END high-word in CRTC_VAL block of NV50_PDISPLAY on NV-50 encodes
			// length of interval from vsync start line to vblank end line, the corrective offset we
			// need to subtract from scanline position. Scanline position measures positive distance from
			// vsync start line (== scanline 0).
			// To low-word likely encodes hsyncstart to hblank end length.
			int vbloffset = (PsychOSKDReadRegister(crtcid, 0x610000 + 0xa00 + 0xe8 + ((crtcid > 0) ? 0x540 : 0), NULL) >> 16) & 0xFFFF;
			// DISPLAY_TOTAL: Encodes VTOTAL in high-word, HTOTAL in low-word:
			int vtotal = (PsychOSKDReadRegister(crtcid, 0x610000 + 0xa00 + 0xf8 + ((crtcid > 0) ? 0x540 : 0), NULL) >> 16) & 0xFFFF;
			printf("head %i: vbloffset = %i : vtotal = %i\n", crtcid, vbloffset, vtotal);
			
		}
		else {
			// NV40 or earlier, aka GeForce-7000 or earlier: Same registers down to
			// earliest NVidia cards NV04 aka RivaTNT-1:

			// Offset between crtc's is 0x2000:
			crtco = (crtcid > 0) ? 0x2000 : 0;

			// Lower 12 bits are vertical scanout position, bit 16 is known to
			// indicate "in vblank" status. All other bits are always zero:
			vblcount = (PsychOSKDReadRegister(crtcid, 0x600808 + crtco, NULL));
			scanline = vblcount & 0xFFF;

			// Bit 4 after right-shift should indicate "in vblank", other bits
			// are always zero:
			vblcount = vblcount >> 12;

			// No support for readout of horizontal scanout position so far:
			scanpixel = 0;
                }

		if (verbose == 1) printf("NV-%x : CRTC %i : scanout x,y = %i , %i   -- vblCount %i\n", cardId, crtcid, scanpixel, scanline, vblcount);

		PsychCopyOutDoubleArg(1, FALSE, (double) scanpixel);
		PsychCopyOutDoubleArg(2, FALSE, (double) scanline);
		PsychCopyOutDoubleArg(3, FALSE, (double) vblcount);
	
		// Done with NVidia GPU test:
		return(PsychError_none);
	}

	// Must be a non-NVidia, hopefully an ATI/AMD piece, otherwise we die:

	// Get basis parameters:
	vbl = PsychOSKDReadRegister(crtcid, AVIVO_D1CRTC_V_BLANK_START_END + crtco, NULL);
	vbl_start = vbl & RADEON_VBEAMPOSITION_BITMASK;
	vbl_end = (vbl >> 16) & RADEON_VBEAMPOSITION_BITMASK;

	vtotal = PsychOSKDReadRegister(crtcid, AVIVO_D1CRTC_V_TOTAL + crtco, NULL);
	htotal = PsychOSKDReadRegister(crtcid, AVIVO_D1CRTC_H_TOTAL + crtco, NULL) + 1;

	// Dotclock shall be 160 Mhz:
	dotclock = 120 * 1e6;
	
	// Convert line length in pixel and video dot clock to line duration
	// in nanoseconds:
	linedur_ns = htotal * 1e9 / dotclock;
	pixeldur_ns = 1e9 / dotclock;
	
	if (verbose > 1) printf("CRTC %i : vbl_start %i : vbl_end %i : vtotal %i : htotal %i : linedur_ns = %i\n", crtcid, vbl_start, vbl_end, vtotal, htotal, linedur_ns);

	// Query current beampos:
	if (pedantic) {
		scanline  = PsychOSKDReadRegister(crtcid, RADEON_D1CRTC_STATUS_POSITION + crtco, NULL);
		scanpixel = (scanline >> RADEON_HBEAMPOSITION_BITSHIFT) & RADEON_VBEAMPOSITION_BITMASK;
		scanline  = scanline & RADEON_VBEAMPOSITION_BITMASK;
	}
	else {
		scanline = PsychOSKDReadRegister(crtcid, RADEON_D1CRTC_STATUS_POSITION + crtco, NULL) & RADEON_VBEAMPOSITION_BITMASK;
	}
	
	// Query current time:
	gettimeofday(&raw_time, NULL);

	// In vertical blank?
	if ((scanline < vbl_start) && (scanline >= vbl_end))
		invbl = FALSE;
	
    scanline = scanline - vbl_end;

    if (invbl && (scanline >= vbl_end))
        scanline = scanline - vtotal - 1;

	// Seconds field will always be zero for any video refresh
	// with a duration of less than 1 hour per frame:
	delta_time.tv_sec = 0;

	// Want vblank_time = raw_time - delta_time, but
	// split it into cases delta_time >= 0 aka scanline >=0
	// and delta_time < 0 aka scanline < 0.
	// This because timeradd/sub only operate on positive 2nd
	// operand:
	if (scanline >= 0) {
		// Convert scanline into microseconds elapsed:
		delta_time.tv_usec = scanline * linedur_ns / 1000;		
		if (pedantic) delta_time.tv_usec = (scanline * linedur_ns + scanpixel * pixeldur_ns) / 1000;
		timersub(&raw_time, &delta_time, &vblank_time);
	}
	else {
		// Convert scanline into microseconds to go:
		delta_time.tv_usec = -scanline * linedur_ns / 1000;
		if (pedantic) delta_time.tv_usec = (-scanline * linedur_ns + scanpixel * pixeldur_ns) / 1000;
		timeradd(&raw_time, &delta_time, &vblank_time);
	}

	// Map vblank_time to double value tUptime -- time in GetSecs() reference.
	double tUptime;
	PsychGetAdjustedPrecisionTimerSeconds(&tUptime);
	tUptime = ((double) raw_time.tv_sec + ((double) raw_time.tv_usec) / 1e6) - tUptime;
	tUptime = ((double) vblank_time.tv_sec + ((double) vblank_time.tv_usec) / 1e6) - tUptime;

	if (verbose > 0) printf("CRTC %i : At time %i.%i secs : Scanline %i : vblank = %i.%i secs.\n", crtcid, raw_time.tv_sec, raw_time.tv_usec, scanline, vblank_time.tv_sec, vblank_time.tv_usec);
	
	PsychCopyOutDoubleArg(1, FALSE, (double) tUptime);
	
	// Done.
	return(PsychError_none);

// ----------------------------------

	char pbuf[1024];
	int residualLength = 1024;

//	if (infoType == 2) {
//		// MMIO register Read for screenid "auxArg1", register offset "auxArg2":
//		PsychCopyInDoubleArg(3, TRUE, &auxArg1);
//		PsychCopyInDoubleArg(4, TRUE, &auxArg2);
//		PsychCopyOutDoubleArg(1, FALSE, (double) PsychOSKDReadRegister((int) auxArg1, (unsigned int) auxArg2, NULL));
//	}
//	
//	if (infoType == 3) {
//		// MMIO register Write for screenid "auxArg1", register offset "auxArg2", to value "auxArg3":
//		PsychCopyInDoubleArg(3, TRUE, &auxArg1);
//		PsychCopyInDoubleArg(4, TRUE, &auxArg2);
//		PsychCopyInDoubleArg(5, TRUE, &auxArg3);
//		PsychOSKDWriteRegister((int) auxArg1, (unsigned int) auxArg2, (unsigned int) auxArg3, NULL);
//	}


	//WORKS:	io_registry_entry_t myentry = IORegistryEntryFromPath(kIOMasterPortDefault, kIOServicePlane ":/AppleACPIPlatformExpert/PCI0@0/AppleACPIPCI/PEGP@1/IOPCI2PCIBridge/display@0/ATY,Wormy@0/ATY_Wormy/display0/AppleBacklightDisplay");
	//WORKS:	io_registry_entry_t myentry = IORegistryEntryFromPath(kIOMasterPortDefault, kIOServicePlane ":/AppleACPIPlatformExpert/PCI0@0/AppleACPIPCI/PEGP@1/IOPCI2PCIBridge/display@0/ATY,Wormy@0/ATY_Wormy/display0");
	io_registry_entry_t myentry = IORegistryEntryFromPath(kIOMasterPortDefault, kIOServicePlane ":/AppleACPIPlatformExpert/PCI0@0/AppleACPIPCI/PEGP@1/IOPCI2PCIBridge/display@0");
	if (myentry == MACH_PORT_NULL) PsychErrorExitMsg(PsychError_system, "PTB-DEBUG: FAILED TO FIND IORegistryEntry!\n");
	printf("GOT IORegistryEntry! >> ");

	if (1) {
		CFDataRef myprop = IORegistryEntryCreateCFProperty(myentry, CFSTR("AAPL00,Dither"), kCFAllocatorDefault, 0);
		if (myprop == NULL) PsychErrorExitMsg(PsychError_system, "PTB-DEBUG: FAILED TO FIND AAPL00,Dither key!!\n");
		m = CFDataGetLength(myprop);
		printf("Amount of data stored in key 'AAPL00,Dither' is %i bytes.\nData: ", m);
		psych_uint8* mypropdata = CFDataGetMutableBytePtr((CFMutableDataRef) myprop);
		if (mypropdata == NULL) PsychErrorExitMsg(PsychError_system, "PTB-DEBUG: FAILED TO GET MUTABLE DATA PTR OF AAPL00,Dither key!!\n");
		for (i=0; i < m; i++) printf("%x ", mypropdata[i]);
		
		mypropdata[0] = 0x1;
		//	FAILS with kIOReturnUnsupported error:	printf("\n COMMIT! RC is %x [kIOReturnUnsupported = %x] \n", (int) IORegistryEntrySetCFProperty(myentry, CFSTR("AAPL00,Dither"), myprop), kIOReturnUnsupported);
		io_connect_t myconnect;
		if (0) {
			// FAILS with kIOReturnUnsupported error:
			if (kIOReturnSuccess == IOServiceOpen(CGDisplayIOServicePort(CGMainDisplayID()), mach_task_self(), 0, &myconnect)) {
				printf("\n COMMIT-IOSERVICE2! RC is %x [kIOReturnUnsupported = %x] \n", (int) IOConnectSetCFProperty(myconnect, CFSTR("AAPL00,Dither"), myprop), kIOReturnUnsupported);
				IOServiceClose(myconnect);
			}
			else {
				printf("\nNOCON\n");
			}
		}

		if (1) {
			kern_return_t	kernResult; 
			io_service_t	service;
			io_connect_t	connect;
			io_iterator_t 	iterator;
			CFDictionaryRef	classToMatch;
			unsigned int				i;
			
			// Setup matching criterion to find our driver in the IORegistry device tree:
//WORKS but useless for our purpose:			
classToMatch = IOServiceNameMatching("ATIRadeonX1000");
//WORKS but fails to open:			classToMatch = IOServiceNameMatching("ATY,Wormy");
//WORKS, IS WHAT WE NEED but fails to open:			classToMatch = IOServiceNameMatching("display");
			if (classToMatch == NULL) {
				printf("PTB-DEBUG: IOServiceMatching() for Psychtoolbox kernel support driver returned a NULL dictionary. Kernel driver support disabled.\n");
				PsychErrorExitMsg(PsychError_system, "FAIL!");
			}
			
			// This creates an io_iterator_t of all instances of our driver that exist in the I/O Registry. Each installed graphics card
			// will get its own instance of a driver. The iterator allows to iterate over all instances:
			kernResult = IOServiceGetMatchingServices(kIOMasterPortDefault, classToMatch, &iterator);
			if (kernResult != KERN_SUCCESS) {
				printf("PTB-DEBUG: IOServiceGetMatchingServices for Psychtoolbox kernel support driver returned 0x%08x. Kernel driver support disabled.\n", kernResult);
				PsychErrorExitMsg(PsychError_system, "FAIL!");
			}
			
			// In a polished final version we would want to handle the case where more than one gfx-card is attached.
			// The iterator would return multiple instances of our driver and we need to decide which one to connect to.
			// For now, we do not handle this case but instead just get the first item from the iterator.
			service = IOIteratorNext(iterator);
			
			// Release the io_iterator_t now that we're done with it.
			IOObjectRelease(iterator);
			
			if (service == IO_OBJECT_NULL) {
				printf("PTB-INFO: Couldn't find requested IOService.\n");
				PsychErrorExitMsg(PsychError_system, "FAIL!");
			}
			
			// Try in a loop to open with different connect flags i = 0 to 0xffff:
			for(i=0; i< 0xffff; i++) {
				kernResult = IOServiceOpen(service, mach_task_self(), i, &myconnect);
				if (kIOReturnSuccess == kernResult) {
					printf("\nCOOLIO!\n");
									printf("\n COMMIT-IOSERVICE3! RC is %x [kIOReturnUnsupported = %x] \n", (int) IOConnectSetCFProperty(myconnect, CFSTR("AAPL00,Dither"), myprop), kIOReturnUnsupported);
					IOServiceClose(myconnect);
					break;
				}
				else {
					if (kernResult != kIOReturnUnsupported) printf("\nNOCON5(%i): %x\n", i, kernResult);
				}
			}
			if (kIOReturnSuccess != kernResult) printf("PTB-INFO: Couldn't open requested IOService.\n");	
		}
		
		CFRelease(myprop);
	}
	
	if (0) {
		CFStringRef myprop = IORegistryEntryCreateCFProperty(myentry, CFSTR("name"), kCFAllocatorDefault, 0);
		CFStringGetCString(myprop, pbuf, 1024, kCFStringEncodingASCII);
		printf("Value of 'name' is %s .\n", (char*) pbuf);
		CFRelease(myprop);
		//	IORegistryEntrySetCFProperty(myentry, CFSTR("name"), CFSTR("funkydisplay"));
		//	IORegistryEntrySetCFProperty(r, CFSTR("IORequestIdle"), kCFBooleanTrue);
	}

	// Release registry entry:
	IOObjectRelease(myentry);

	#endif

	return(PsychError_none);
	
	#if PSYCH_SYSTEM == PSYCH_OSX
		printf("PTB-DEBUG: Shutting down connection to PsychtoolboxKernelDriver, if any...\n");
		PsychOSShutdownPsychtoolboxKernelDriverInterface();
	#endif
	
	#if PSYCH_SYSTEM == PSYCH_LINUX
		printf("PTB-DEBUG: PreRADEON_R500_GEN_INT_CNTL: %x\n", PsychOSKDReadRegister(0, 0x040, NULL));
		PsychOSKDWriteRegister(0, 0x040, PsychOSKDReadRegister(0, 0x040, NULL) | 0x1, NULL);
		PsychWaitIntervalSeconds(0.1);
		printf("PTB-DEBUG: PostRADEON_R500_GEN_INT_CNTL: %x\n", PsychOSKDReadRegister(0, 0x040, NULL));

		printf("PTB-DEBUG: RADEON_R500_GEN_INT_STATUS: %x\n", PsychOSKDReadRegister(0, 0x044, NULL));
		printf("PTB-DEBUG: RADEON_R500_GEN_INT_CNTL: %x\n", PsychOSKDReadRegister(0, RADEON_R500_GEN_INT_CNTL, NULL));
		printf("PTB-DEBUG: RADEON_R500_GEN_INT_STATUS: %x\n", PsychOSKDReadRegister(0, RADEON_R500_GEN_INT_STATUS, NULL));
		printf("PTB-DEBUG: RADEON_D1CRTC_INTERRUPT_CONTROL: %x\n", PsychOSKDReadRegister(0, RADEON_D1CRTC_INTERRUPT_CONTROL, NULL));
	#endif

	return(PsychError_none);

	//demonstrate how we find the function and subfunction names
	//printf("Psychtoolbox function: %s, subfunction %s\n", PsychGetModuleName(), PsychGetFunctionName() );

	//copy all the input argument to their outputs if we have doubles, if not error.  
	numInArgs = PsychGetNumInputArgs();
	numOutArgs = PsychGetNumOutputArgs();
	numNamedOutArgs = PsychGetNumNamedOutputArgs();
	PsychErrorExit(PsychCapNumOutputArgs(numInArgs));

	/*
	printf("number of input arguments: %d\n", numInArgs);
	printf("number of output arguments: %d\n", numOutArgs);
	printf("number of named output arguments: %d\n", numNamedOutArgs);
	*/

	//iterate over each of the supplied inputs.  If the input is a two-dimensional array 
	//of doubles or a character array, then copy it to the output.  
	for(i=1;i<=numInArgs;i++){

		format = PsychGetArgType(i);

		switch(format){

			case PsychArgType_double:

				if(PsychGetArgM(i)==1 && PsychGetArgN(i)==1){

					tempValue=i;  //if 1x1 double then the default return value is the arg position.

					PsychCopyInDoubleArg(i, FALSE, &tempValue);

					PsychCopyOutDoubleArg(i, FALSE, tempValue);

				}else{

					PsychAllocInDoubleMatArg(i, FALSE, &m, &n, &p, &array);

					PsychCopyOutDoubleMatArg(i, FALSE, m, n, p, array);

				}

				break;

			case PsychArgType_char:

				str=NULL; //This tells PsychGetCharArg() to use its own (volatile) memory. 

				PsychAllocInCharArg(i, FALSE, &str); 

				PsychCopyOutCharArg(i, FALSE, str);

				break;

			case PsychArgType_default:

				PsychCopyOutCharArg(i, FALSE, defaultString);

				break;

		}

	}

	return(PsychError_none);
}
