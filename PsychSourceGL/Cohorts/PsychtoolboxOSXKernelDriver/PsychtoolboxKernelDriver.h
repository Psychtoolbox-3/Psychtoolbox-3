/*
	File:			PsychtoolboxKernelDriver.h
	
	Description:	This file implements the I/O Kit driver kernel extension for Psychtoolbox (KEXT).

	Copyright:		Copyright © 2008-2012 Mario Kleiner, derived from an Apple example code.

	Change History of original Apple sample code (most recent first):

            1.1			05/22/2007			User client performs endian swapping when called from a user process 
											running using Rosetta. Updated to produce a universal binary.
											Now requires Xcode 2.2.1 or later to build.
			
			1.0d3	 	01/14/2003			New sample.

*/


#include <IOKit/IOService.h>
#include <IOKit/pci/IOPCIDevice.h>
#include <IOKit/IOInterruptEventSource.h>
#include <IOKit/IOFilterInterruptEventSource.h>

// Definitions of GPU registers etc.:
#include "PsychGraphicsCardRegisterSpecs.h"

// Handy IOLog/printf format strings for dealing with types that have a different
// length on LP64.
#if __LP64__
#define UInt32_FORMAT       "%u"
#define UInt32_x_FORMAT     "0x%08x"
#define PhysAddr_FORMAT     "0x%016llx"
#define PhysLen_FORMAT      "%llu"
#define VirtAddr_FORMAT     "0x%016llx"
#define ByteCount_FORMAT    "%llu"
#else
#define UInt32_FORMAT       "%lu"
#define UInt32_x_FORMAT     "0x%08lx"
#define PhysAddr_FORMAT     "0x%08lx"
#define PhysLen_FORMAT      UInt32_FORMAT
#define VirtAddr_FORMAT     "0x%08x"
#define ByteCount_FORMAT    UInt32_FORMAT
#endif

// PTB driver revision:
#define PTBKDRevision 0

struct PsychKDCommandStruct;

class PsychtoolboxKernelDriver : public IOService
{
    OSDeclareDefaultStructors(PsychtoolboxKernelDriver)
	
private:
	UInt32							fDeviceType;
	UInt32							fCardType;
	UInt16							fPCIDeviceId;
	UInt16							fPCIVendorId;
	IOPCIDevice*					fPCIDevice;
	IOMemoryMap*					fRadeonMap;
	IOVirtualAddress				fRadeonRegs;
	UInt32							fRadeonSize;
	UInt32							fRadeonLowlimit;
    UInt32                          fNumDisplayHeads;
	IOFilterInterruptEventSource*	fInterruptSrc;
	
	UInt32							fInterruptCookie;
	UInt32							fInterruptCounter;
	UInt32							fVBLCounter[2];
    UInt32                          oldDither[(DCE4_MAXHEADID + 1)];

	// Initialize our own interrupt handler for snooping on gfx-card state changes:
	bool InitializeInterruptHandler(void);	

	// Fast-Path interrupt handler: Gets called in direct primary hardware interrupt context!
	bool fastPathInterruptHandler(OSObject* myself, IOFilterInterruptEventSource* mySource);

	// Slow-Path WorkLoop interrupt handler: Gets called if the fast-path handler returns "true".
	void workLoopInterruptHandler(OSObject* myself, IOInterruptEventSource* mySource, int pendingIRQs);
	
	// Handle VBLANK IRQ's for display head 'headId':
	void handleVBLIRQ(UInt32 headId);
	
	// Read 32 bit control register at 'offset':
	UInt32	ReadRegister(UInt32 offset);

	// Write 32 bit control register at 'offset' with 'value':
	void	WriteRegister(UInt32 offset, UInt32 value);

	// Return current vertical rasterbeam position of display head 'headId' (0=Primary CRTC1, 1=Secondary CRTC2):
	UInt32 GetBeamPosition(UInt32 headId);
	
	// Instantaneously resynchronize display heads of a Radeon dual-head gfx-card:
	SInt32	FastSynchronizeAllDisplayHeads(void);
	
	// Try to change hardware dither mode on GPU:
	void SetDitherMode(UInt32 headId, UInt32 ditherOn);
	
	// Perform instant state snapshot of interesting registers and return'em:
	void	GetStateSnapshot(PsychKDCommandStruct* outStruct);
	
	// Dump interesting register state to system log:
	void	DumpGfxState(void);
	
    // Returns multiple flags with info like PCI Vendor/device id, display engine type etc.
    void GetGPUInfo(UInt32 *inOutArgs);

    // Query if LUT for given headId is all-zero:
    UInt32 GetLUTState(UInt32 headId, UInt32 debug);

    // Load an identity LUT into display head 'headid':
    UInt32 LoadIdentityLUT(UInt32 headId);

	// Is a given ATI/AMD GPU a DCE4 type ASIC, i.e., with the new display engine?
	bool isDCE4(void);

	// Is a given ATI/AMD GPU a DCE4.1 type ASIC, i.e., with the new display engine?
	bool isDCE41(void);

	// Is a given ATI/AMD GPU a DCE5 type ASIC, i.e., with the new display engine?
	bool isDCE5(void);

	// Is a given ATI/AMD GPU a DCE6 type ASIC, i.e., with the new display engine?
	bool isDCE6(void);

	// Is a given ATI/AMD GPU a DCE6.1 type ASIC, i.e., with the new display engine?
	bool isDCE61(void);

	// Is a given ATI/AMD GPU a DCE6.4 type ASIC, i.e., with the new display engine?
	bool isDCE64(void);
	
public:
	// IOService methods
	virtual bool init(OSDictionary* dictionary = 0);
	virtual void free(void);
	
	virtual IOService* probe(IOService* provider, SInt32* score);
	
	virtual bool start(IOService* provider);
	virtual void stop(IOService* provider);
	
	virtual bool willTerminate(IOService* provider, IOOptionBits options);
	virtual bool didTerminate(IOService* provider, IOOptionBits options, bool* defer);
	
    virtual bool terminate(IOOptionBits options = 0);
    virtual bool finalize(IOOptionBits options);

	// SimpleDriver methods
//	virtual IOReturn ScalarIStructI(int inNumber, MySampleStruct* inStruct, IOByteCount inStructSize);
//	virtual IOReturn ScalarIStructO(int inNumber1, int inNumber2, MySampleStruct* outStruct, IOByteCount* outStructSize);
//	virtual IOReturn ScalarIScalarO(int inNumber1, int inNumber2, int* outNumber);
	virtual IOReturn PsychKDDispatchCommand(PsychKDCommandStruct* inStruct, PsychKDCommandStruct* outStruct, IOByteCount inStructSize, IOByteCount* outStructSize);
};