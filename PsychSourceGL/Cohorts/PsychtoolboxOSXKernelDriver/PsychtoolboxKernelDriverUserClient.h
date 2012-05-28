/*
	File:			PsychtoolboxKernelDriverUserClient.h
	
	Description:	This file implements a simple I/O Kit user client that is Rosetta-aware.

	Copyright:		Copyright © 2008 - 2012 Mario Kleiner, derived from an Apple example code.

	Change History of original Apple sample code (most recent first):
                    
            1.1			05/22/2007			User client performs endian swapping when called from a user process 
											running using Rosetta. Updated to produce a universal binary.
											Now requires Xcode 2.2.1 or later to build.
			
			1.0d3	 	01/14/2003			New sample.

*/


#include <IOKit/IOService.h>
#include <IOKit/IOUserClient.h>
#include "PsychUserKernelShared.h"

class PsychtoolboxKernelDriver;

class PsychtoolboxKernelDriverUserClient : public IOUserClient
{
    OSDeclareDefaultStructors(PsychtoolboxKernelDriverUserClient)
    
protected:
    PsychtoolboxKernelDriver*			fProvider;
    task_t								fTask;
	bool								fCrossEndian;
    static const IOExternalMethod		fMethods[kNumberOfMethods];
      
public:
    // IOUserClient methods
    virtual void stop(IOService* provider);
    virtual bool start(IOService* provider);
    
	virtual bool initWithTask(task_t owningTask, void* securityToken, UInt32 type, OSDictionary* properties);

    virtual IOReturn clientClose(void);
    virtual IOReturn clientDied(void);

	virtual bool willTerminate(IOService* provider, IOOptionBits options);
	virtual bool didTerminate(IOService* provider, IOOptionBits options, bool* defer);
	
    virtual bool terminate(IOOptionBits options = 0);
    virtual bool finalize(IOOptionBits options);
	
protected:	
    virtual IOExternalMethod* getTargetAndMethodForIndex(IOService** target, UInt32 index);

    // SimpleUserClient methods
	virtual IOReturn openUserClient(void);
    virtual IOReturn closeUserClient(void);
    
//	virtual IOReturn ScalarIStructI(int inNumber, MySampleStruct* inStruct, IOByteCount inStructSize);
//	virtual IOReturn ScalarIStructO(int inNumber1, int inNumber2, MySampleStruct* outStruct, IOByteCount* outStructSize);
	virtual IOReturn PsychKDDispatchCommand(PsychKDCommandStruct* inStruct, PsychKDCommandStruct* outStruct, IOByteCount inStructSize, IOByteCount* outStructSize);
};