/*
 *  PsychSourceGL/Source/OSX/PsychHID/PsychHIDGenericUSBSupport.c
 *
 *  PROJECTS: PsychHID
 *
 *  Platform: OS X
 *
 *  Authors:
 *
 *  Christopher Broussard <chrg@sas.upenn.edu>		cgb
 *
 *  HISTORY:
 *
 *	4.4.2009	Created.
 *
 *  DESCRIPTION:
 *
 *  OS/X specific support routines that implement generic USB device handling.
 *
 *  Each OS platform needs its own specific version of this file.
 *
 */

#include "PsychHID.h"

// Function Declarations
IOReturn ConfigureDevice(IOUSBDeviceInterface **dev, int configIdx);

// Perform device control transfer on USB device:
int PsychHIDOSControlTransfer(PsychUSBDeviceRecord* devRecord, psych_uint8 bmRequestType, psych_uint8 bRequest, psych_uint16 wValue, psych_uint16 wIndex, psych_uint16 wLength, void *pData)
{
	IOUSBDevRequest request;
	IOUSBDeviceInterface **dev;
	
	// Setup the USB request data structure.
	memset(&request, 0, sizeof(IOUSBDevRequest));
	request.bmRequestType = bmRequestType;
	request.bRequest = bRequest;
	request.wValue = wValue;
	request.wLength = wLength;
	request.wIndex = wIndex;
	request.pData = pData;

	dev = (IOUSBDeviceInterface **) devRecord->device;
	if (dev == NULL) {
		PsychErrorExitMsg(PsychError_internal, "IOUSBDeviceInterface** device points to NULL device!");
	}
	
	// Send the data across the USB bus by executing the device control request. Return status code.
	// On success, zero aka kIOReturnSuccess is returned, otherwise some kernel error return code.
	return((int) ((*dev)->DeviceRequest(dev, &request)) );
}

// Close USB device, mark device record as "free/invalid":
void PsychHIDOSCloseUSBDevice(PsychUSBDeviceRecord* devRecord)
{	
	(void)(*(devRecord->device))->USBDeviceClose(devRecord->device);
	(void)(*(devRecord->device))->Release(devRecord->device);
	devRecord->device = NULL;
	devRecord->valid = 0;
}


// Open first USB device that satisfies given matching critera, mark device record as "active/valid":
// errorcode would contain a diagnostic error code on failure, but is not yet used.
// spec contains the specification of the device to open and how to configure it at open time.
// Returns true on success, false on error or if no matching device could be found.
psych_bool PsychHIDOSOpenUSBDevice(PsychUSBDeviceRecord* devRecord, int* errorcode, PsychUSBSetupSpec* spec)
{
	kern_return_t           kr;
	CFMutableDictionaryRef  matchingDict;
	SInt32                  usbVendor = (SInt32) spec->vendorID;
	SInt32                  usbProduct = (SInt32) spec->deviceID;
	IOUSBDeviceInterface182 **dev = NULL;
	io_iterator_t           iterator;
	IOCFPlugInInterface     **plugInInterface = NULL;
	HRESULT                 result;
	io_service_t            usbDevice;
	SInt32                  score;
	UInt16                  vendor;
	UInt16                  product;
	UInt16                  release;
    int                     i;
	psych_bool					deviceFound = false;
	char *name="",*description="";

	// Init errorcode to failed:
	*errorcode = -2;
	
	// Set up matching dictionary for class IOUSBDevice and its subclasses
	matchingDict = IOServiceMatching(kIOUSBDeviceClassName);
	if (!matchingDict) {
		PsychErrorExitMsg(PsychError_system, "Couldn't create a USB matching dictionary.");
	}
	
	//Add the vendor and product IDs to the matching dictionary.
	//This is the second key in the table of device-matching keys of the
	//USB Common Class Specification
	CFDictionarySetValue(matchingDict, CFSTR(kUSBVendorName),
						 CFNumberCreate(kCFAllocatorDefault,
										kCFNumberSInt32Type, &usbVendor));
	
	CFDictionarySetValue(matchingDict, CFSTR(kUSBProductName),
						 CFNumberCreate(kCFAllocatorDefault,
										kCFNumberSInt32Type, &usbProduct));
	
	kr = IOServiceGetMatchingServices(kIOMasterPortDefault, matchingDict, &iterator);
	if (kr) {
		*errorcode = (int) kr;
		PsychHIDErrors(NULL, *errorcode, &name, &description); 
		printf("PsychHID: PsychHIDOSOpenUSBDevice: IO-KIT error: %08x = '%s' [%s]\n", *errorcode, name, description);
		PsychErrorExitMsg(PsychError_system, "Couldn't get matching services\n");
	}
	
	// Attempt to find the correct device.
	while (usbDevice = IOIteratorNext(iterator)) {
		// Create an intermediate plug-in
		kr = IOCreatePlugInInterfaceForService(usbDevice,
											   kIOUSBDeviceUserClientTypeID, kIOCFPlugInInterfaceID,
											   &plugInInterface, &score);
		
		// Don't need the device object after intermediate plug-in is created
		IOObjectRelease(usbDevice);

		if ((kIOReturnSuccess != kr) || !plugInInterface) {
			*errorcode = (int) kr;
			PsychHIDErrors(NULL, *errorcode, &name, &description); 
			printf("PsychHID: PsychHIDOSOpenUSBDevice: IO-KIT error: %08x = '%s' [%s]\n", *errorcode, name, description);			
			printf("PsychHID: PsychHIDOSOpenUSBDevice: WARNING! Unable to create a plug-in (%08x)\n", kr);
			continue;
		}
		
		// Now create the device interface
		result = (*plugInInterface)->QueryInterface(plugInInterface,
													CFUUIDGetUUIDBytes(kIOUSBDeviceInterfaceID182),
													(LPVOID *)&dev);
		
		// Don't need the intermediate plug-in after device interface is created.
		(*plugInInterface)->Release(plugInInterface);
		
		if (result || !dev) {
			*errorcode = (int) result;
			PsychHIDErrors(NULL, *errorcode, &name, &description); 
			printf("PsychHID: PsychHIDOSOpenUSBDevice: IO-KIT error: %08x = '%s' [%s]\n", *errorcode, name, description);
			printf("PsychHID: PsychHIDOSOpenUSBDevice: WARNING! Couldn't create a device interface (%08x)\n", (int) result);
			continue;
		}
		
		// Check these values for confirmation.
		kr = (*dev)->GetDeviceVendor(dev, &vendor);
		kr = (*dev)->GetDeviceProduct(dev, &product);
		kr = (*dev)->GetDeviceReleaseNumber(dev, &release);
		if (((int) vendor != spec->vendorID) || ((int) product != spec->deviceID)) {
			printf("PsychHID: PsychHIDOSOpenUSBDevice: WARNING! Found unwanted device (vendor = %d, device = %d)\n", vendor, product);
			(void) (*dev)->Release(dev);
			continue;
		}
		else {
			deviceFound = true;
			//printf("Vendor: 0x%x\nProduct: 0x%x\nRelease: 0x%x\n", vendor, product, release);
			break;
		}
	}

    // Release the iterator now that we're done with it.
    IOObjectRelease(iterator);

	// At this point, all allocated ressources and objects have been released and
	// either deviceFound == false on failure, or deviceFound == true and 'dev'
	// is the device interface to our device.
	if (deviceFound) {
		// Open the device to change its state. Retry opening up to 5 times if
        // it doesn't immediately succeed:
        for (i = 0; i < 5; i++) {
            // Do a seize-open: This asks other drivers or apps claiming the device to
            // release it for us if possible -- increases our chance of successful open:
            kr = (*dev)->USBDeviceOpenSeize(dev);
            
            // If we don't get a kIOReturnExclusiveAccess, we're done, either successfully
            // or with fatal error:
            if (kr != kIOReturnExclusiveAccess) break;
            
            // Got kIOReturnExclusiveAccess -- Device claimed by other app or driver.
            // Wait a second, then possibly retry, in the hope that the seize-request
            // has been honored after a small amount of time:
            PsychWaitIntervalSeconds(1);
        }
        
        // Final verdict on our open attempt?
		if (kr != kIOReturnSuccess) {
			(void) (*dev)->Release(dev);
			
			*errorcode = (int) kr;
			PsychHIDErrors(NULL, *errorcode, &name, &description); 
			printf("PsychHID: PsychHIDOSOpenUSBDevice: IO-KIT error: %08x = '%s' [%s]\n", *errorcode, name, description);
			PsychErrorExitMsg(PsychError_system, "Unable to open USB device.");
		}
		
		// Configure device
		kr = ConfigureDevice((IOUSBDeviceInterface **) dev, spec->configurationID);
		if (kr != kIOReturnSuccess) {
			(void) (*dev)->USBDeviceClose(dev);
			(void) (*dev)->Release(dev);

			*errorcode = (int) kr;
			PsychHIDErrors(NULL, *errorcode, &name, &description); 
			printf("PsychHID: PsychHIDOSOpenUSBDevice: IO-KIT error: %08x = '%s' [%s]\n", *errorcode, name, description);
			PsychErrorExitMsg(PsychError_system, "Unable to configure USB device.");
		}
		
		// Success! Assign device interface and mark device record as active/open/valid:
		devRecord->device = dev;
		devRecord->valid = 1;
		
		// Set errorcode to success:
		*errorcode = 0;
	}
	else {
		// No matching device found. NULL-out the record, we're done.
		// This is not strictly needed, as this NULL state is the initial
		// state of the record upon entering this function.
		devRecord->device = NULL;
		devRecord->valid = 0;
	}

	// Return the success status.
	return (deviceFound);
}


IOReturn ConfigureDevice(IOUSBDeviceInterface **dev, int configIdx)
{
	UInt8                           numConfig, current_bConfigurationValue;
	IOReturn                        kr;
	IOUSBConfigurationDescriptorPtr configDesc;
	
	// A configIdx == -1 means: Skip configuration.
	if (configIdx == -1) return(kIOReturnSuccess);
	
	// Get the number of configurations. The sample code always chooses
	// the first configuration (at index 0) but your code may need a
	// different one
	kr = (*dev)->GetNumberOfConfigurations(dev, &numConfig);	
	if (kr || (numConfig == 0)) {
		printf("PsychHID: USB ConfigureDevice: ERROR! Error getting number of configurations or no configurations available at all (err = %08x)\n", kr);
		return(kr);
	}
	
	if (configIdx < 0 || configIdx >= (int) numConfig) {
		printf("PsychHID: USB ConfigureDevice: ERROR! Provided configuration index %i outside support range 0 - %i for this device!\n", configIdx, (int) numConfig);
		return(-1);
	}
	
	// Get the configuration descriptor for index 'configIdx':
	kr = (*dev)->GetConfigurationDescriptorPtr(dev, (UInt8) configIdx, &configDesc);
	if (kr) {
		printf("PsychHID: USB ConfigureDevice: ERROR! Couldn't get configuration descriptor for index %d (err = %08x)\n", configIdx, kr);
		return(kr);
	}
	
	kr = (*dev)->GetConfiguration(dev, &current_bConfigurationValue);
	if (kr) {
		printf("PsychHID: USB ConfigureDevice: ERROR! Couldn't get configuration of device (err = %08x)\n", kr);
		return(kr);
	}

	// If current value is already identical to requested value, we're done:
	if (current_bConfigurationValue == configDesc->bConfigurationValue) return kIOReturnSuccess;

	// Set the device's configuration. The configuration value is found in
	// the bConfigurationValue field of the configuration descriptor
	kr = (*dev)->SetConfiguration(dev, configDesc->bConfigurationValue);
	if (kr) {
		printf("PsychHID: USB ConfigureDevice: ERROR! Couldn't set configuration to value %d (err = %08x)\n", (int) configDesc->bConfigurationValue, kr);
		return(kr);
	}
	
	return kIOReturnSuccess;
}
