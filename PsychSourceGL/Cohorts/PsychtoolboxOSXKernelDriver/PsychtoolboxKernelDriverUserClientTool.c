/*
	File:			PsychtoolboxKernelDriverUserClientTool.c
	
	Description:	This file shows how to communicate with the I/O Kit user client of the PsychtoolboxKernelDriver.

	Copyright:		Copyright © 2008-2012 Mario Kleiner, derived from an Apple example code.
	
    This kernel driver and associated tools as a whole are licensed to you as follows: (MIT style license):

    * Permission is hereby granted, free of charge, to any person obtaining a
    * copy of this software and associated documentation files (the "Software"),
    * to deal in the Software without restriction, including without limitation
    * the rights to use, copy, modify, merge, publish, distribute, sublicense,
    * and/or sell copies of the Software, and to permit persons to whom the
    * Software is furnished to do so, subject to the following conditions:
    *
    * The above copyright notice and this permission notice shall be included in
    * all copies or substantial portions of the Software.
    *
    * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
    * THE COPYRIGHT HOLDER(S) OR AUTHOR(S) BE LIABLE FOR ANY CLAIM, DAMAGES OR
    * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
    * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
    * OTHER DEALINGS IN THE SOFTWARE.

				
	Change History of original Apple sample code (most recent first):

            1.1			05/22/2007			User client performs endian swapping when called from a user process 
											running using Rosetta. Updated to produce a universal binary.
											Now requires Xcode 2.2.1 or later to build.
			
			1.0d3	 	01/14/2003			New sample.

*/


#include <IOKit/IOKitLib.h>
#include <ApplicationServices/ApplicationServices.h>
#include "PsychtoolboxKernelDriverUserClientInterface.h"
#include <limits.h>
#include <unistd.h>

// Definitions of GPU registers etc.:
#include "PsychGraphicsCardRegisterSpecs.h"

#define kMyPathToSystemLog			"/var/log/system.log"

// IO_OBJECT_NULL was added in Mac OS X 10.4. If it's not defined in the headers, define it ourselves.
// This allows the same code to build using the 10.4.0 SDK for Intel-based systems and earlier SDKs for backwards compatibility
// with PowerPC-based systems.
#ifndef IO_OBJECT_NULL
#define	IO_OBJECT_NULL	((io_object_t) 0)
#endif

void G80SetCursorPosition(int headId, int ditherOn);

kern_return_t MyUserClientOpenExample(io_service_t service, io_connect_t *connect)
{
    // This call will cause the user client to be instantiated. It returns an io_connect_t handle
	// that is used for all subsequent calls to the user client.
    kern_return_t kernResult = IOServiceOpen(service, mach_task_self(), 0, connect);
	
    if (kernResult != KERN_SUCCESS) {
        fprintf(stderr, "IOServiceOpen returned 0x%08x\n", kernResult);
    }
	else {
		// This is an example of calling our user client's openUserClient method.
		kernResult = MyOpenUserClient(*connect);
			
		if (kernResult == KERN_SUCCESS) {
			//printf("MyOpenUserClient was successful.\n\n");
		}
		else {
			fprintf(stderr, "MyOpenUserClient returned 0x%08x.\n\n", kernResult);
		}
    }
		
	return kernResult;
}


void MyUserClientCloseExample(io_connect_t connect)
{
	kern_return_t kernResult = MyCloseUserClient(connect);
        
    if (kernResult == KERN_SUCCESS) {
        //printf("MyCloseUserClient was successful.\n\n");
    }
	else {
		fprintf(stderr, "MyCloseUserClient returned 0x%08x.\n\n", kernResult);
	}
    
    kernResult = IOServiceClose(connect);
    
    if (kernResult == KERN_SUCCESS) {
        //printf("IOServiceClose was successful.\n\n");
    }
    else {
	    fprintf(stderr, "IOServiceClose returned 0x%08x\n\n", kernResult);
    }
}

/*
void MyScalarIStructureIExample(io_connect_t connect)
{
    MySampleStruct	sampleStruct = { 586, 8756 };		// These are random numbers I picked.
    int				sampleNumber = 15;					// Another random number.
    
    kern_return_t kernResult = 
		MyScalarIStructureI(connect, sampleNumber, sizeof(MySampleStruct), &sampleStruct);
	    
    if (kernResult == KERN_SUCCESS) {
        printf("MyScalarIStructureI was successful.\n\n");
    }
	else {
		fprintf(stderr, "MyScalarIStructureI returned 0x%08x.\n\n", kernResult);
	}
}


void MyScalarIStructureOExample(io_connect_t connect)
{
    MySampleStruct	sampleStruct;
    int				sampleNumber1 = 154;	// This number is random.
    int				sampleNumber2 = 863;	// This number is random.
    IOByteCount		structSize = sizeof(MySampleStruct); 
    
	kern_return_t kernResult =
		MyScalarIStructureO(connect, sampleNumber1, sampleNumber2, &structSize, &sampleStruct);
        
    if (kernResult == KERN_SUCCESS) {
        printf("MyScalarIStructureO was successful.\n");
        printf("int16 = %d, int32 = %ld\n\n", sampleStruct.int16, sampleStruct.int32);
    }
	else {
		fprintf(stderr, "MyScalarIStructureO returned 0x%08x.\n\n", kernResult);
	}
}


void MyScalarIScalarOExample(io_connect_t connect)
{
    int				sampleNumber1 = 10;		// Random number with no meaning.
    int				sampleNumber2 = 32768;	// Another random number.
    int				resultNumber;
    
	kern_return_t kernResult =
		MyScalarIScalarO(connect, sampleNumber1, sampleNumber2, &resultNumber);
	
    if (kernResult == KERN_SUCCESS) {
        printf("MyScalarIScalarO was successful.\n");
        printf("resultNumber = %d\n\n", resultNumber);
    }
	else {
		fprintf(stderr, "MyScalarIScalarO returned 0x%08x.\n\n", kernResult);
	}
}
*/

void TestPsychKDDispatchCommand(io_connect_t connect)
{
    PsychKDCommandStruct	sampleStruct1 = { 0xdeadbeef, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };	// These are random numbers I picked.
    PsychKDCommandStruct	sampleStruct2;
    IOByteCount		structSize1 = sizeof(PsychKDCommandStruct);
    IOByteCount		structSize2 = sizeof(PsychKDCommandStruct);
    int i;
	
	kern_return_t kernResult = PsychKDDispatchCommand(connect, structSize1, &sampleStruct1, &structSize2, &sampleStruct2);    
    
	if (kernResult == KERN_SUCCESS) {
        printf("MyStructureIStructureO was successful: Command code %lx.\n", sampleStruct2.command);
		for (i=0; i < kPsychMaxDriverArgs; i++) printf("[%i]: %i, ", i, sampleStruct2.inOutArgs[i]);
		printf("\n");
    }
	else {
		fprintf(stderr, "MyStructureIStructureO returned 0x%08x.\n\n", kernResult);
	}
}


void MyLaunchConsoleApp()
{
    CFURLRef pathRef;

    pathRef = CFURLCreateWithFileSystemPath(kCFAllocatorDefault, CFSTR(kMyPathToSystemLog), kCFURLPOSIXPathStyle, false);
    
    if (pathRef) {
        LSOpenCFURLRef(pathRef, NULL);
        CFRelease(pathRef);
    }
}

int PsychOSSynchronizeAllDisplayHeads(io_connect_t connect) {
	PsychKDCommandStruct syncCommand;
    IOByteCount			 structSize1 = sizeof(PsychKDCommandStruct);
	
	// Set command code for display sync:
	syncCommand.command = kPsychKDFastSyncAllHeads;
	
	// Issue request:
	kern_return_t kernResult = PsychKDDispatchCommand(connect, structSize1, &syncCommand, &structSize1, &syncCommand);    
	if (kernResult == KERN_SUCCESS) {
		printf("PTB-INFO: Graphics display heads resynchronized. Residual vertical beamposition error is %ld scanlines.\n", syncCommand.inOutArgs[0]);
	}
	else {
		printf("PTB-ERROR: Graphics display head synchronization failed (Kernel error code: %lx).\n", kernResult);
		return(INT_MAX);
	}

	return((int) syncCommand.inOutArgs[0]);
}

int PsychOSKDGetBeamposition(io_connect_t connect, int headid)
{
	PsychKDCommandStruct syncCommand;
    IOByteCount			 structSize1 = sizeof(PsychKDCommandStruct);
	
	// Set command code for display sync:
	syncCommand.command = kPsychKDGetBeamposition;
	syncCommand.inOutArgs[0] = headid;
	
	// Issue request:
	kern_return_t kernResult = PsychKDDispatchCommand(connect, structSize1, &syncCommand, &structSize1, &syncCommand);    
	if (kernResult != KERN_SUCCESS) {
		printf("PTB-ERROR: Kernel driver beamposition failed (Kernel error code: %lx).\n", kernResult);
		// A value of -1 signals beamposition queries unsupported:
		return(-1);
	}
	
	// Return queried position:
	return((int) syncCommand.inOutArgs[0]);
}

unsigned int PsychOSKDReadRegister(io_connect_t connect, unsigned int offset)
{
	PsychKDCommandStruct syncCommand;
    IOByteCount			 structSize1 = sizeof(PsychKDCommandStruct);
	
	// Set command code for display sync:
	syncCommand.command = kPsychKDReadRegister;
	syncCommand.inOutArgs[0] = offset;
	
	// Issue request:
	kern_return_t kernResult = PsychKDDispatchCommand(connect, structSize1, &syncCommand, &structSize1, &syncCommand);    
	if (kernResult != KERN_SUCCESS) {
		printf("PTB-ERROR: Kernel driver register read failed (Kernel error code: %lx).\n", kernResult);
		// A value of 0xffffffff signals failure:
		return(0xffffffff);
	}
	
	// Return queried position:
	return((int) syncCommand.inOutArgs[0]);
}

unsigned int PsychOSKDWriteRegister(io_connect_t connect, unsigned int offset, unsigned int value)
{
	PsychKDCommandStruct syncCommand;
    IOByteCount			 structSize1 = sizeof(PsychKDCommandStruct);
	
	// Set command code for display sync:
	syncCommand.command = KPsychKDWriteRegister;
	syncCommand.inOutArgs[0] = offset;
	syncCommand.inOutArgs[1] = value;
	
	// Issue request:
	kern_return_t kernResult = PsychKDDispatchCommand(connect, structSize1, &syncCommand, &structSize1, &syncCommand);    
	if (kernResult != KERN_SUCCESS) {
		printf("PTB-ERROR: Kernel driver register write failed (Kernel error code: %lx).\n", kernResult);
		// A value of 1 signals failure:
		return(1);
	}
	
	// Return success:
	return(0);
}

unsigned int PsychOSKDSetDitherMode(io_connect_t connect, unsigned int head, unsigned int dither)
{
	PsychKDCommandStruct syncCommand;
    IOByteCount			 structSize1 = sizeof(PsychKDCommandStruct);
	
	// Set command code for display sync:
	syncCommand.command = kPsychKDSetDitherMode;
	syncCommand.inOutArgs[0] = head;
	syncCommand.inOutArgs[1] = dither;
	
	// Issue request:
	kern_return_t kernResult = PsychKDDispatchCommand(connect, structSize1, &syncCommand, &structSize1, &syncCommand);    
	if (kernResult != KERN_SUCCESS) {
		printf("PTB-ERROR: Kernel driver kPsychKDSetDitherMode failed (Kernel error code: %lx).\n", kernResult);
		// A value of 1 signals failure:
		return(1);
	}
	
	// Return success:
	return(0);
}

unsigned int PsychOSKDGetRevision(io_connect_t connect)
{
	PsychKDCommandStruct syncCommand;
    IOByteCount			 structSize1 = sizeof(PsychKDCommandStruct);
	
	// Set command code for display sync:
	syncCommand.command = kPsychKDGetRevision;
	
	// Issue request:
	kern_return_t kernResult = PsychKDDispatchCommand(connect, structSize1, &syncCommand, &structSize1, &syncCommand);    
	if (kernResult != KERN_SUCCESS) {
		printf("PTB-ERROR: Kernel driver revision read failed (Kernel error code: %lx).\n", kernResult);
		// A value of 0xffffffff signals failure:
		return(0xffffffff);
	}
	
	// Return revision:
	return((unsigned int) syncCommand.inOutArgs[0]);
}

void PsychOSKDGetGPUInfo(io_connect_t connect)
{
	PsychKDCommandStruct syncCommand;
    IOByteCount			 structSize1 = sizeof(PsychKDCommandStruct);
	
	// Set command code for display sync:
	syncCommand.command = kPsychKDGetGPUInfo;
	
	// Issue request:
	kern_return_t kernResult = PsychKDDispatchCommand(connect, structSize1, &syncCommand, &structSize1, &syncCommand);    
	if (kernResult != KERN_SUCCESS) {
		printf("PTB-ERROR: Kernel driver gpu info read failed (Kernel error code: %lx).\n", kernResult);
		// A value of 0xffffffff signals failure:
		return;
	}
	
	printf("PTB-INFO: GPU-Vendor %i, PCIId %x, GPU-Type %x.\n", syncCommand.inOutArgs[0], syncCommand.inOutArgs[1], syncCommand.inOutArgs[2]);
	return;
}

unsigned int PsychOSKDGetLUTState(io_connect_t connect, unsigned int head)
{
	PsychKDCommandStruct syncCommand;
    IOByteCount			 structSize1 = sizeof(PsychKDCommandStruct);
	
	// Set command code for display sync:
	syncCommand.command = kPsychKDGetLUTState;
	syncCommand.inOutArgs[0] = head;
	syncCommand.inOutArgs[1] = 1;
	
	// Issue request:
	kern_return_t kernResult = PsychKDDispatchCommand(connect, structSize1, &syncCommand, &structSize1, &syncCommand);    
	if (kernResult != KERN_SUCCESS) {
		printf("PTB-ERROR: Kernel driver lut state read failed (Kernel error code: %lx).\n", kernResult);
		// A value of 0xffffffff signals failure:
		return(0xffffffff);
	}
	
	// Return lut state:
	return((unsigned int) syncCommand.inOutArgs[0]);
}

unsigned int PsychOSKDLoadIdentityLUT(io_connect_t connect, unsigned int head)
{
	PsychKDCommandStruct syncCommand;
    IOByteCount			 structSize1 = sizeof(PsychKDCommandStruct);
	
	// Set command code for display sync:
	syncCommand.command = kPsychKDSetIdentityLUT;
	syncCommand.inOutArgs[0] = head;
	
	// Issue request:
	kern_return_t kernResult = PsychKDDispatchCommand(connect, structSize1, &syncCommand, &structSize1, &syncCommand);    
	if (kernResult != KERN_SUCCESS) {
		printf("PTB-ERROR: Kernel driver identity lut setup failed (Kernel error code: %lx).\n", kernResult);
		// A value of 0xffffffff signals failure:
		return(0xffffffff);
	}
	
	// Return lut setup rc:
	return((unsigned int) syncCommand.inOutArgs[0]);
}

// Define globally for subroutines:
io_connect_t	connect;

int main(int argc, char* argv[])
{
    kern_return_t	kernResult; 
    io_service_t	service;
    io_iterator_t 	iterator;
    CFDictionaryRef	classToMatch;
    int				i,j;
	unsigned int	ov;
	
    if (argc < 2) {
        printf("\n");
        printf("%s: Command line tool for controlling the MacOS/X PsychtoolboxKernelDriver:\n", (argc > 0) ? argv[0] : "unknown");
        printf("%s: Command line arguments missing. Usage:\n", (argc > 0) ? argv[0] : "unknown");
        printf("%s -setdither headId ditherEnable\n", (argc > 0) ? argv[0] : "unknown");
        printf("%s -resyncdisplays\n", (argc > 0) ? argv[0] : "unknown");
        printf("%s -revision\n", (argc > 0) ? argv[0] : "unknown");
        printf("%s -gpuinfo\n", (argc > 0) ? argv[0] : "unknown");
        printf("%s -lutstatus headId\n", (argc > 0) ? argv[0] : "unknown");
        printf("%s -identitylut headId\n", (argc > 0) ? argv[0] : "unknown");
        
        printf("%s -nvidiasetcursor x y\n", (argc > 0) ? argv[0] : "unknown");
        printf("%s -radeonswaptests\n", (argc > 0) ? argv[0] : "unknown");
        printf("\n\n");
        return(-2);
    }
    
    // This will launch the Console.app so you can see the IOLogs from the KEXT.
    // MyLaunchConsoleApp();

	classToMatch = IOServiceMatching(kMyDriversIOKitClassName);
    
    if (classToMatch == NULL) {
        fprintf(stderr, "IOServiceMatching returned a NULL dictionary.\n");
        return -1;
    }
    
    // This creates an io_iterator_t of all instances of our driver that exist in the I/O Registry.
    kernResult = IOServiceGetMatchingServices(kIOMasterPortDefault, classToMatch, &iterator);
    
    if (kernResult != KERN_SUCCESS) {
        fprintf(stderr, "IOServiceGetMatchingServices returned 0x%08x\n\n", kernResult);
        return -1;
    }
        
    // In a real product you would want to handle the case where more than one of your devices is attached.
	// The iterator would return multiple instances of your driver and you'd need to decide which one to connect to.
	// One way to do this would be to provide some user interface to allow the user to choose the device.
	// In the interest of simplicity this sample does not handle this case but instead just gets the first item from the iterator.
    service = IOIteratorNext(iterator);
    
    // Release the io_iterator_t now that we're done with it.
    IOObjectRelease(iterator);
    
    if (service == IO_OBJECT_NULL) {
        fprintf(stderr, "Couldn't find any matching drivers.\n");
    }
    else {
		// Instantiate a connection to the user client.
		kernResult = MyUserClientOpenExample(service, &connect);
		
		// Release the io_service_t now that we're done with it.
		IOObjectRelease(service);

		if (kernResult != KERN_SUCCESS) {
			fprintf(stderr, "Failed to open connection to driver!\n");
			if (kernResult == kIOReturnExclusiveAccess) fprintf(stderr, "Check if other app is using driver already [kIOReturnExclusiveAccess]!\n");
			return(-1);
		}
		
		if (connect != IO_OBJECT_NULL) {	

            if (argc > 1 && (strcmp(argv[1], "-setdither") == 0)) {
                if (argc < 4) {
                    printf("Error: You must provide the headId and ditherEnable flag!\n");
                    return(-1);
                }
                
                printf("Trying to %s dithering on display head %i.\n", (atoi(argv[3]) > 0) ? "enable" : "disable", atoi(argv[2]));
                PsychOSKDSetDitherMode(connect, atoi(argv[2]), atoi(argv[3]));
            }
            
            if (argc > 1 && (strcmp(argv[1], "-resyncdisplays") == 0)) {
                printf("Trying to resync all display heads.\n");
				printf("Residual offset between displays is %i scanlines.\n\n", PsychOSSynchronizeAllDisplayHeads(connect));
            }
            
            if (argc > 1 && (strcmp(argv[1], "-revision") == 0)) {
				printf("Revision number of kernel driver is %i.\n\n", PsychOSKDGetRevision(connect));
            }

            if (argc > 1 && (strcmp(argv[1], "-gpuinfo") == 0)) {
				printf("GPU info:\n\n");
                PsychOSKDGetGPUInfo(connect);
            }
            
            if (argc > 1 && (strcmp(argv[1], "-lutstatus") == 0)) {
                if (argc < 3) {
                    printf("Error: You must provide the headId!\n");
                    return(-1);
                }
                
                printf("LUT status of display head %i is %x\n", atoi(argv[2]), PsychOSKDGetLUTState(connect, atoi(argv[2])));
            }

            if (argc > 1 && (strcmp(argv[1], "-identitylut") == 0)) {
                if (argc < 3) {
                    printf("Error: You must provide the headId!\n");
                    return(-1);
                }
                
                printf("Setting identity LUT on display head %i. RC = %i\n", atoi(argv[2]), PsychOSKDLoadIdentityLUT(connect, atoi(argv[2])));
            }

			if (argc > 3 && (strcmp(argv[1], "-nvidiasetcursor") == 0)) {
				// NVidia test:
                printf("Trying to move NVidia G80 cursor to position %i x %i.\n", atoi(argv[2]), atoi(argv[3]));
				G80SetCursorPosition(atoi(argv[2]), atoi(argv[3]));
			}

            if (argc > 1 && (strcmp(argv[1], "-radeonswaptests") == 0)) {				
				printf("D1: Primary surface is %lx\n\n", PsychOSKDReadRegister(connect, RADEON_D1GRPH_PRIMARY_SURFACE_ADDRESS));
				printf("D1: Secondary surface is %lx\n\n", PsychOSKDReadRegister(connect, RADEON_D1GRPH_SECONDARY_SURFACE_ADDRESS));
				printf("D1: Pitch is %lx\n\n", PsychOSKDReadRegister(connect, RADEON_D1GRPH_PITCH));
				printf("D2: Primary surface is %lx\n\n", PsychOSKDReadRegister(connect, RADEON_D2GRPH_PRIMARY_SURFACE_ADDRESS));
				printf("D2: Secondary surface is %lx\n\n", PsychOSKDReadRegister(connect, RADEON_D2GRPH_SECONDARY_SURFACE_ADDRESS));
				printf("D2: Pitch is %lx\n\n", PsychOSKDReadRegister(connect, RADEON_D2GRPH_PITCH));
				
				while(1) {
					if (PsychOSKDReadRegister(connect, RADEON_D2GRPH_PRIMARY_SURFACE_ADDRESS)!=i) {
						i=PsychOSKDReadRegister(connect, RADEON_D2GRPH_PRIMARY_SURFACE_ADDRESS);
						printf("D2: New Primary surface is %lx\n\n", i);
						printf("D2: New Secondary surface is %lx\n\n", PsychOSKDReadRegister(connect, RADEON_D2GRPH_SECONDARY_SURFACE_ADDRESS));
					}
				}
				
				if (0) {
					printf("D1: RADEON_D1GRPH_CONTROL %lx\n\n", PsychOSKDReadRegister(connect, RADEON_D1GRPH_CONTROL));
					PsychOSKDWriteRegister(connect, RADEON_D1GRPH_LUT_SEL, 0x1 << 8);
					printf("D1: RADEON_D1GRPH_LUT_SEL %lx\n\n", PsychOSKDReadRegister(connect, RADEON_D1GRPH_LUT_SEL));
					
					printf("Switching to ARGB2101010...\n");
					ov = PsychOSKDReadRegister(connect, RADEON_D1GRPH_CONTROL);
					PsychOSKDWriteRegister(connect, RADEON_D1GRPH_CONTROL, ((0x1 << 8) | ov));
					for(j=1; j<10; j++) for(i=1; i< 1000*1000*1000; i++);
					printf("Switching back to ARGB8888...\n");
					PsychOSKDWriteRegister(connect, RADEON_D1GRPH_CONTROL, ov);
					PsychOSKDWriteRegister(connect, RADEON_D1GRPH_LUT_SEL, 0);
					
					printf("D1: RADEON_LVTMA_BIT_DEPTH_CONTROL %lx\n\n", PsychOSKDReadRegister(connect, RADEON_LVTMA_BIT_DEPTH_CONTROL));
				}
			}
			
			// Close the user client and tear down the connection.
			MyUserClientCloseExample(connect);
		}
	}
	
	return 0;
}

// This blob of code allows to change the G90 GPU's hardware cursor position.
// Proof that our MMIO mapping works correctly, but other than that pretty pointles...
void G80SetCursorPosition(int headId, int ditherOn)
{
    const int headOff = 0x1000 * headId;
	int x, y, xr, yr;
	
	for (yr=0; yr< 500; yr+=20) {
		for (xr=0; xr< 500; xr++) {
			x = xr & 0xffff;
			y = yr & 0xffff;
			PsychOSKDWriteRegister(connect, (0x00647084 + headOff), y << 16 | x);
			PsychOSKDWriteRegister(connect, (0x00647080 + headOff), 0);
			usleep(1000);
		}
	}
	return;
}
