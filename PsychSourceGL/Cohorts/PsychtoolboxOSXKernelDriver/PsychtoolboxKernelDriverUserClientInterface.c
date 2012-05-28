/*
	File:			PsychtoolboxKernelDriverUserClientInterface.c
	
	Description:	Implements an abstraction layer between client applications and the user client.

	Copyright:		Copyright © 2008 - 2012 Mario Kleiner, derived from an Apple example code.

	Change History of original Apple sample code (most recent first):

            1.1			05/22/2007			Split out user client interface functions from main test tool.

*/


#include <IOKit/IOKitLib.h>
#include "PsychtoolboxKernelDriverUserClientInterface.h"

kern_return_t MyOpenUserClient(io_connect_t connect)
{
    // This calls the openUserClient method in SimpleUserClient inside the kernel. Though not mandatory, it's good
    // practice to use open and close semantics in your driver to prevent multiple user space applications from
    // using your driver at the same time.
    
    kern_return_t kernResult = IOConnectMethodScalarIScalarO(connect, kMyUserClientOpen, 0, 0);
    
    return kernResult;
}


kern_return_t MyCloseUserClient(io_connect_t connect)
{
    // This calls the closeUserClient method in SimpleUserClient inside the kernel, which in turn closes
	// the driver.
    
	kern_return_t kernResult = IOConnectMethodScalarIScalarO(connect, kMyUserClientClose, 0, 0);

    return kernResult;
}

/*
kern_return_t MyScalarIStructureI(io_connect_t connect, const int scalarI, 
								  const IOByteCount structISize, const MySampleStruct* structI)
{
    // This calls the function ScalarIStructI in SimpleUserClient inside the kernel.
	    
	kern_return_t kernResult =
		IOConnectMethodScalarIStructureI(connect,					// an io_connect_t returned from IOServiceOpen().
										 kMyScalarIStructIMethod,	// Index to the function to be called via the user client.
										 1,							// the number of scalar input values.
										 structISize,				// the size of the input structure paramter.
										 scalarI,					// a scalar input parameter.
										 structI					// a pointer to a struct input parameter.
										 );
    
    return kernResult;
}


kern_return_t MyScalarIStructureO(io_connect_t connect, const int scalarI_1, const int scalarI_2,
								  IOByteCount* structOSize, MySampleStruct* structO)
{
    // This calls the function ScalarIStructureO in SimpleUserClient inside the kernel.

	kern_return_t kernResult =
		IOConnectMethodScalarIStructureO(connect,					// an io_connect_t returned from IOServiceOpen().
                                         kMyScalarIStructOMethod,	// an index to the function to be called via the user client.
										 2,							// the number of scalar input values.
										 structOSize,				// the size of the struct output paramter.
										 scalarI_1,					// a scalar input parameter.
										 scalarI_2,					// another scalar input parameter.
										 &structO					// a pointer to a struct output parameter.
										 );
        
    return kernResult;
}


kern_return_t MyScalarIScalarO(io_connect_t connect, const int scalarI_1, const int scalarI_2, int* scalarO)
{
    // This calls the function ScalarIScalarO in SimpleUserClient inside the kernel.
    
    kern_return_t kernResult =
		IOConnectMethodScalarIScalarO(connect,						// an io_connect_t returned from IOServiceOpen().
									  kMyScalarIScalarOMethod,		// an index to the function to be called via the user client.
									  2,							// the number of scalar input values.
									  1,							// the number of scalar output values.
									  scalarI_1,					// a scalar input parameter.
									  scalarI_2,					// another scalar input parameter.
									  scalarO						// a scalar output parameter.
									  );
    
    return kernResult;
}
*/

kern_return_t PsychKDDispatchCommand(io_connect_t connect, const IOByteCount structISize, const PsychKDCommandStruct* structI,
									 IOByteCount* structOSize, PsychKDCommandStruct* structO)
{
    // This calls the function StructureIStructureO in SimpleUserClient inside the kernel.
    
    kern_return_t kernResult =
		IOConnectMethodStructureIStructureO(connect,					// an io_connect_t returned from IOServiceOpen().
                                            kPsychKDDispatchCommand,	// an index to the function to be called via the user client.
                                            structISize,				// the size of the input struct paramter.
                                            structOSize,				// a pointer to the size of the output struct paramter.
                                            (PsychKDCommandStruct*) structI,	// a pointer to the input struct parameter.
                                            structO						// a pointer to the output struct parameter.
                                            );
    
    return kernResult;
}
