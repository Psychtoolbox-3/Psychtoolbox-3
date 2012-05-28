/*
	File:			PsychtoolboxKernelDriverUserClient.cpp
	
	Description:	This file implements a simple I/O Kit user client that is Rosetta-aware.

	Copyright:		Copyright © 2008-2012 Mario Kleiner, derived from an Apple example code.

	Change History of original Apple sample code (most recent first):

            1.1			05/22/2007			User client performs endian swapping when called from a user process 
											running using Rosetta. Updated to produce a universal binary.
											Now requires Xcode 2.2.1 or later to build.
			
			1.0d3	 	01/14/2003			New sample.

*/


#include <IOKit/IOLib.h>
#include <IOKit/IOKitKeys.h>
#include <libkern/OSByteOrder.h>
#include "PsychtoolboxKernelDriverUserClient.h"
#include "PsychtoolboxKernelDriver.h"


#define super IOUserClient
OSDefineMetaClassAndStructors(PsychtoolboxKernelDriverUserClient, IOUserClient)


// The following I/O Registry property keys are currently undefined in the PowerPC architecture headers.
// We define them here so the same source will build on both architectures. 

#ifndef kIOUserClientCrossEndianKey
#define kIOUserClientCrossEndianKey				"IOUserClientCrossEndian"
#endif
#ifndef kIOUserClientCrossEndianCompatibleKey
#define kIOUserClientCrossEndianCompatibleKey	"IOUserClientCrossEndianCompatible"
#endif


// Sentinel values for the method dispatch table
enum {
    kMethodObjectThis = 0,
    kMethodObjectProvider
};


// User client method dispatch table.
//
// The user client mechanism is designed to allow calls from a user process to be dispatched to
// any IOService-based object in the kernel. Almost always this mechanism is used to dispatch calls to
// either member functions of the user client itself or of the user client's provider. The provider is
//  the driver which the user client is connecting to the user process.
//
// While this sample shows one case of dispatching calls directly to the driver (ScalarIScalarO),
// it is recommended that calls be dispatched to the user client. This allows the user client to perform
// error checking on the parameters before passing them to the driver. It also allows the user client to
// do any endian-swapping of parameters in the cross-endian case. (See ScalarIStructI below for further
// discussion of this subject.)
//
// The dispatch table makes use of the sentinel values kMethodObjectThis and kMethodObjectProvider to
// represent at compile time the values of the this pointer and fProvider respectively at run time.  
const IOExternalMethod PsychtoolboxKernelDriverUserClient::fMethods[kNumberOfMethods] = {
	{   // kMyUserClientOpen
		(IOService *) kMethodObjectThis,									// Target object is this user client.
		(IOMethod) &PsychtoolboxKernelDriverUserClient::openUserClient,	// Method pointer.
		kIOUCScalarIScalarO,												// Scalar Input, Scalar Output.
		0,																	// No scalar input values.
		0																	// No scalar output values.
	},
	{   // kMyUserClientClose
		(IOService *) kMethodObjectThis,									// Target object is this user client.
		(IOMethod) &PsychtoolboxKernelDriverUserClient::closeUserClient,	// Method pointer.
		kIOUCScalarIScalarO,												// Scalar Input, Scalar Output.
		0,																	// No scalar input values.
		0																	// No scalar output values.
	},
//	{   // kMyScalarIStructIMethod
//		(IOService *) kMethodObjectThis,									// Target object is this user client.
//		(IOMethod) &PsychtoolboxKernelDriverUserClient::ScalarIStructI,	// Method pointer.
//		kIOUCScalarIStructI,												// Scalar Input, Struct Input.
//		1,																	// One scalar input value.
//		sizeof(MySampleStruct)												// The size of the input struct.
//	},
//	{   // kMyScalarIStructOMethod
//		(IOService *) kMethodObjectThis,									// Target object is this user client.
//		(IOMethod) &PsychtoolboxKernelDriverUserClient::ScalarIStructO,	// Method pointer.
//		kIOUCScalarIStructO,												// Scalar Input, Struct Output.
//		2,																	// Two scalar input values.
//		sizeof(MySampleStruct)												// The size of the output struct.
//	},
//	{   // kMyScalarIScalarOMethod
//		(IOService *) kMethodObjectProvider,								// Target object is this user client's provider
//																			// (the driver).
//		(IOMethod) &PsychtoolboxKernelDriver::ScalarIScalarO,		// Method pointer.
//		kIOUCScalarIScalarO,												// Scalar Input, Scalar Output.
//		2,																	// Two scalar input values.
//		1																	// One scalar output value.
//	},
	{   // kPsychKDDispatchCommand
		(IOService *) kMethodObjectThis,									// Target object is this user client.
		(IOMethod) &PsychtoolboxKernelDriverUserClient::PsychKDDispatchCommand,	// Method pointer.
		kIOUCStructIStructO,												// Struct Input, Struct Output.
		sizeof(PsychKDCommandStruct),												// The size of the input struct.
		sizeof(PsychKDCommandStruct)												// The size of the output struct.
	}
};
    
    
IOExternalMethod* PsychtoolboxKernelDriverUserClient::getTargetAndMethodForIndex(IOService** target, UInt32 index)
{
    // IOLog("PsychtoolboxKernelDriverUserClient::getTargetAndMethodForIndex(index = %d)\n", (int)index);
    
    // Make sure that the index of the function we're calling actually exists in the function table.
    if (index < (UInt32) kNumberOfMethods) {
		if (fMethods[index].object == (IOService *) kMethodObjectThis) {
			*target = this;	   
        }
		else {
			*target = fProvider;	   
		}
		return (IOExternalMethod *) &fMethods[index];
    }
	else {
		*target = NULL;
		return NULL;
	}
}


// There are two forms of IOUserClient::initWithTask, the second of which accepts an additional OSDictionary* parameter.
// If your user client needs to modify its behavior when it's being used by a process running using Rosetta,
// you need to implement the form of initWithTask with this additional parameter.
//
// initWithTask is called as a result of the user process calling IOServiceOpen.
bool PsychtoolboxKernelDriverUserClient::initWithTask(task_t owningTask, void* securityToken, UInt32 type,
														 OSDictionary* properties)
{
    bool	success;

    // IOLog("PsychtoolboxKernelDriverUserClient::initWithTask()\n");
    
	success = super::initWithTask(owningTask, securityToken, type, properties);	    
	
	if (success) {
		// This code will do the right thing on both PowerPC- and Intel-based systems because the cross-endian
		// property will never be set on PowerPC-based Macs. 
		fCrossEndian = false;
	
		if (properties != NULL && properties->getObject(kIOUserClientCrossEndianKey)) {
			// A connection to this user client is being opened by a user process running using Rosetta.
			
			// Indicate that this user client can handle being called from cross-endian user processes by 
			// setting its IOUserClientCrossEndianCompatible property in the I/O Registry.
			if (setProperty(kIOUserClientCrossEndianCompatibleKey, kOSBooleanTrue)) {
				fCrossEndian = true;
				IOLog("PsychtoolboxKernelDriverUserClient::initWithTask(): fCrossEndian = true\n");
			}
		}
	}
	
    fTask = owningTask;
    fProvider = NULL;
        
    return success;
}


// start is called after initWithTask as a result of the user process calling IOServiceOpen.
bool PsychtoolboxKernelDriverUserClient::start(IOService* provider)
{
    bool	success;
	
	// IOLog("PsychtoolboxKernelDriverUserClient::start()\n");
    
    // Verify that this user client is being started with a provider that it knows
	// how to communicate with.
	fProvider = OSDynamicCast(PsychtoolboxKernelDriver, provider);
    success = (fProvider != NULL);
    
    if (success) {
		// It's important not to call super::start if some previous condition
		// (like an invalid provider) would cause this function to return false. 
		// I/O Kit won't call stop on an object if its start function returned false.
		success = super::start(provider);
	}
	
    return success;
}


// We override stop only to log that it has been called to make it easier to follow the user client's lifecycle.
void PsychtoolboxKernelDriverUserClient::stop(IOService* provider)
{
    // IOLog("PsychtoolboxKernelDriverUserClient::stop()\n");
    
    super::stop(provider);
}


// clientClose is called as a result of the user process calling IOServiceClose.
IOReturn PsychtoolboxKernelDriverUserClient::clientClose(void)
{
    // IOLog("PsychtoolboxKernelDriverUserClient::clientClose()\n");
    
    // Defensive coding in case the user process called IOServiceClose
	// without calling closeUserClient first.
    (void) closeUserClient();
    
	// Inform the user process that this user client is no longer available. This will also cause the
	// user client instance to be destroyed.
	//
	// terminate would return false if the user process still had this user client open.
	// This should never happen in our case because this code path is only reached if the user process
	// explicitly requests closing the connection to the user client.
	bool success = terminate();
	if (!success) {
		IOLog("PsychtoolboxKernelDriverUserClient::clientClose(): terminate() failed.\n");
	}

    // DON'T call super::clientClose, which just returns kIOReturnUnsupported.
    
    return kIOReturnSuccess;
}


// clientDied is called if the client user process terminates unexpectedly (crashes).
// We override clientDied only to log that it has been called to make it easier to follow the user client's lifecycle.
// Production user clients need to override clientDied only if they need to take some alternate action if the user process
// crashes instead of exiting normally.
IOReturn PsychtoolboxKernelDriverUserClient::clientDied(void)
{
    IOReturn result = kIOReturnSuccess;

    IOLog("PsychtoolboxKernelDriverUserClient::clientDied()\n");

    // The default implementation of clientDied just calls clientClose.
    result = super::clientDied();

    return result;
}


// willTerminate is called at the beginning of the termination process. It is a notification
// that a provider has been terminated, sent before recursing up the stack, in root-to-leaf order.
//
// This is where any pending I/O should be terminated. At this point the user client has been marked
// inactive and any further requests from the user process should be returned with an error.
bool PsychtoolboxKernelDriverUserClient::willTerminate(IOService* provider, IOOptionBits options)
{
	// IOLog("PsychtoolboxKernelDriverUserClient::willTerminate()\n");
	
	return super::willTerminate(provider, options);
}


// didTerminate is called at the end of the termination process. It is a notification
// that a provider has been terminated, sent after recursing up the stack, in leaf-to-root order.
bool PsychtoolboxKernelDriverUserClient::didTerminate(IOService* provider, IOOptionBits options, bool* defer)
{
	// IOLog("PsychtoolboxKernelDriverUserClient::didTerminate()\n");
	
	// If all pending I/O has been terminated, close our provider. If I/O is still outstanding, set defer to true
	// and the user client will not have stop called on it.
	closeUserClient();
	*defer = false;
	
	return super::didTerminate(provider, options, defer);
}


// We override terminate only to log that it has been called to make it easier to follow the user client's lifecycle.
// Production user clients will rarely need to override terminate. Termination processing should be done in
// willTerminate or didTerminate instead.
bool PsychtoolboxKernelDriverUserClient::terminate(IOOptionBits options)
{
    bool	success;
    
    // IOLog("PsychtoolboxKernelDriverUserClient::terminate()\n");

    success = super::terminate(options);
    
    return success;
}


// We override finalize only to log that it has been called to make it easier to follow the user client's lifecycle.
// Production user clients will rarely need to override finalize.
bool PsychtoolboxKernelDriverUserClient::finalize(IOOptionBits options)
{
    bool	success;
    
    // IOLog("PsychtoolboxKernelDriverUserClient::finalize()\n");
    
    success = super::finalize(options);
    
    return success;
}


IOReturn PsychtoolboxKernelDriverUserClient::openUserClient(void)
{
    IOReturn	result = kIOReturnSuccess;
	
	IOLog("PsychtoolboxKernelDriverUserClient::openUserClient(): Client connects to us :-)\n");
    
    if (fProvider == NULL || isInactive()) {
		// Return an error if we don't have a provider. This could happen if the user process
		// called openUserClient without calling IOServiceOpen first. Or, the user client could be
		// in the process of being terminated and is thus inactive.
        result = kIOReturnNotAttached;
	}
    else if (!fProvider->open(this)) {
		// The most common reason this open call will fail is because the provider is already open
		// and it doesn't support being opened by more than one client at a time.
		result = kIOReturnExclusiveAccess;
	}
        
    return result;
}


IOReturn PsychtoolboxKernelDriverUserClient::closeUserClient(void)
{
    IOReturn	result = kIOReturnSuccess;
	
    IOLog("PsychtoolboxKernelDriverUserClient::closeUserClient(): Client disconnects from us :-(\n");
            
    if (fProvider == NULL) {
		// Return an error if we don't have a provider. This could happen if the user process
		// called closeUserClient without calling IOServiceOpen first. 
		result = kIOReturnNotAttached;
		IOLog("PsychtoolboxKernelDriverUserClient::closeUserClient(): returning kIOReturnNotAttached.\n");
	}
	else if (fProvider->isOpen(this)) {
		// Make sure we're the one who opened our provider before we tell it to close.
		fProvider->close(this);
	}
	else {
		result = kIOReturnNotOpen;
		IOLog("PsychtoolboxKernelDriverUserClient::closeUserClient(): returning kIOReturnNotOpen.\n");
	}
	
    return result;
}

/*
IOReturn PsychtoolboxKernelDriverUserClient::ScalarIStructI(int inNumber, MySampleStruct* inStruct, IOByteCount inStructSize)
{
	IOReturn	result;

	IOLog("PsychtoolboxKernelDriverUserClient::ScalarIStructI(inNumber = %d, int16 = %d, int32 = %d, inStructSize = %d)\n",
		  inNumber, inStruct->int16, (int)inStruct->int32, (int)inStructSize);
    
	// Endian-swap structure parameters in the user client before passing them to the driver.
	//
	// This may require adding new functions to your user client and modifying the dispatch table in
	// getTargetAndMethodForIndex to point to these new functions.
	//
	// This approach is greatly preferable because it avoids the complexity of a driver which can be opened by multiple clients,
	// each of which may or may not be cross-endian. It also avoids having to change the driver to make it cross-endian-aware.
	//
	// Note that fCrossEndian will always be false if running on a PowerPC-based Mac.
	
	if (fProvider == NULL || isInactive()) {
		// Return an error if we don't have a provider. This could happen if the user process
		// called ScalarIStructI without calling IOServiceOpen first. Or, the user client could be
		// in the process of being terminated and is thus inactive.
		result = kIOReturnNotAttached;
	}
	else if (!fProvider->isOpen(this)) {
		// Return an error if we do not have the driver open. This could happen if the user process
		// did not call openUserClient before calling this function.
		result = kIOReturnNotOpen;
	}
	else {
		if (fCrossEndian) {
			// Structures aren't automatically swapped by the user client mechanism as it has no knowledge of how the fields
			// structure are laid out.
			
			// Swap the fields of the structure passed by the client user process before passing it to the driver.
			// Use the unconditional swap macros here as we know only at runtime if we're being called from a
			// cross-endian user process running using Rosetta.
			
			inStruct->int16 = OSSwapInt16(inStruct->int16);
			inStruct->int32 = OSSwapInt32(inStruct->int32);
			
			IOLog("PsychtoolboxKernelDriverUserClient::ScalarIStructI(after swap: inNumber = %d, int16 = %d, int32 = %d, inStructSize = %d)\n",
				  inNumber, inStruct->int16, (int) inStruct->int32, (int) inStructSize);
		}
		
		result = fProvider->ScalarIStructI(inNumber, inStruct, inStructSize);
	}
	
	return result;
}


IOReturn PsychtoolboxKernelDriverUserClient::ScalarIStructO(int inNumber1, int inNumber2,
															   MySampleStruct* outStruct, IOByteCount* outStructSize)
{
	IOReturn	result;

    IOLog("PsychtoolboxKernelDriverUserClient::ScalarIStructO(inNumber1 = %d, inNumber2 = %d)\n", inNumber1, inNumber2);

	if (fProvider == NULL || isInactive()) {
		// Return an error if we don't have a provider. This could happen if the user process
		// called ScalarIStructO without calling IOServiceOpen first. Or the user client could be
		// in the process of being terminated and is thus inactive.
		result = kIOReturnNotAttached;
	}
	else if (!fProvider->isOpen(this)) {
		// Return an error if we do not have the driver open. This could happen if the user process
		// did not call openUserClient before calling this function.
		result = kIOReturnNotOpen;
	}
	else {
		result = fProvider->ScalarIStructO(inNumber1, inNumber2, outStruct, outStructSize);

		// Note that fCrossEndian will always be false if running on a PowerPC-based Mac.
		if (fCrossEndian) {
			// Swap the fields of the structure returned by the driver before returning it to the client user process.
			// Use the unconditional swap macros here as we know only at runtime if we're being called from a
			// cross-endian user process running using Rosetta.

			outStruct->int16 = OSSwapInt16(outStruct->int16);
			outStruct->int32 = OSSwapInt32(outStruct->int32);

			IOLog("PsychtoolboxKernelDriverUserClient::ScalarIStructO(output after swap: int16 = %d, int32 = %d, outStructSize = %d)\n",
				  outStruct->int16, (int) outStruct->int32, (int) *outStructSize);
		}
	}

	return result;
}
*/

IOReturn PsychtoolboxKernelDriverUserClient::PsychKDDispatchCommand(PsychKDCommandStruct* inStruct, PsychKDCommandStruct* outStruct,
															   IOByteCount inStructSize, IOByteCount* outStructSize)
{
	IOReturn	result;

	if (fProvider == NULL || isInactive()) {
		// Return an error if we don't have a provider. This could happen if the user process
		// called StructIStructO without calling IOServiceOpen first. Or, the user client could be
		// in the process of being terminated and is thus inactive.
		result = kIOReturnNotAttached;
	}
	else if (!fProvider->isOpen(this)) {
		// Return an error if we do not have the driver open. This could happen if the user process
		// did not call openUserClient before calling this function.
		result = kIOReturnNotOpen;
	}
	else {
		// Note that fCrossEndian will always be false if running on a PowerPC-based Mac.
		if (fCrossEndian) {
			// Swap the fields of the structure passed to the driver:
			// Use the unconditional swap macros here as we know only at runtime if we're being called from a
			// cross-endian user process running using Rosetta.
			inStruct->command = OSSwapInt32(inStruct->command);
			for (UInt32 i=0; i < kPsychMaxDriverArgs; i++) inStruct->inOutArgs[i] = OSSwapInt32(inStruct->inOutArgs[i]);
		}
		
		result = fProvider->PsychKDDispatchCommand(inStruct, outStruct, inStructSize, outStructSize);

		if (fCrossEndian) {
			// Swap the results returned by the driver before returning them to the client user process.
			outStruct->command = OSSwapInt32(outStruct->command);
			for (UInt32 i=0; i < kPsychMaxDriverArgs; i++) outStruct->inOutArgs[i] = OSSwapInt32(outStruct->inOutArgs[i]);
		}
	}
    
    return result;
}
