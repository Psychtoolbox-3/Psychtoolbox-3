/*
	File:			PsychtoolboxKernelDriverUserClientInterface.h
	
	Description:	Implements an abstraction layer between client applications and the user client.

	Copyright:		Copyright © 2008 - 2012 Mario Kleiner, derived from an Apple example code.
					
	Change History of original Apple sample code (most recent first):

            1.1			05/22/2007			Split out user client interface functions from main test tool.

*/

#include "PsychUserKernelShared.h"

// Open connection to driver:
kern_return_t MyOpenUserClient(io_connect_t connect);

// Close connection to driver:
kern_return_t MyCloseUserClient(io_connect_t connect);

// Dispatch a command to the driver and return its results:
kern_return_t PsychKDDispatchCommand(io_connect_t connect, const IOByteCount structISize, const PsychKDCommandStruct* structI,
									 IOByteCount* structOSize, PsychKDCommandStruct* structO);

/*
kern_return_t MyScalarIStructureI(io_connect_t connect, const int scalarI, const IOByteCount structISize, 
								  const MySampleStruct* structI);
kern_return_t MyScalarIStructureO(io_connect_t connect, const int scalarI_1, const int scalarI_2,
								  IOByteCount* structOSize, MySampleStruct* structO);
kern_return_t MyScalarIScalarO(io_connect_t connect, const int scalarI_1, const int scalarI_2, int* scalarO);
*/
