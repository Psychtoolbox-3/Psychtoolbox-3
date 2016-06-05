/*
    File:             PsychtoolboxKernelDriver.cpp

    Description:    This file implements a I/O Kit driver kernel extension (KEXT) for Psychtoolbox-3.
                    The driver allows to augment the standard graphics driver with some useful features
                    for visual psychophysics applications. It is orthogonal to the system graphics driver in that
                    it doesn't change that drivers functionality, nor does it interact with that driver.

                    The driver currently works on AMD/ATI Radeon graphics cards of the X1000, HD2000 and later.
                    It may work on older ATI cards, but that is not tested or supported.

                    The driver also supports rasterposition queries on the majority of NVidia GPU's, but no
                    other functions.

                    The driver is experimental in nature, with some small but non-zero chance of causing
                    malfunctions or crashes of your system, so USE AT YOUR OWN RISK AND WITH CAUTION!

                    The following features are currently supported:

                    * Rasterbeamposition queries on Intel-based Macintosh computers.

                    * A function to quickly synchronize the video refresh cycles of the two output heads of
                      a dual-head graphics card, e.g., for binocular or stereo stimulation on dual-display setups.
                      or all six display heads of an AMD evergreen GPU (HD-5000 and later).

                    * Control of digital display bit depths truncation and dithering on AMD hardware.

                    * Control and testing of gamma tables and other transformations.

                    * Guided queries and logging of a few interesting gfx-state registers for vision science applications.

                    * A low level interface for reading- and writing from/to arbitrary gfx-hardware control registers.

                    * A dump function that dumps interesting settings to the system log.

                    The rasterbeamposition queries and fast dual-head syncing are used by Psychtoolbox if
                    the driver is loaded at Screen() startup time.

                    All other features can be accessed by the stand-alone UNIX command line tool
                    "PsychtoolboxKernelDriverUserClientTool", to be found in the PsychAlpha subfolder.

    Copyrights and license:

                    This drivers copyright:
                    Copyright © 2008-2015 Mario Kleiner.

                    This driver contains code for radeon DCE-4 and AVIVO which are derived from the free software radeon kms
                    driver for Linux (radeon_display.c, radeon_regs.h). The Radeon kms driver has the following copyright:

                    * Copyright 2007-8 Advanced Micro Devices, Inc.
                    * Copyright 2008 Red Hat Inc.
                    *
                    * Authors: Dave Airlie
                    *          Alex Deucher

                    This driver also contains small fragments of code for NVidia GPU model detection which are derived from the
                    free software Nouveau kms driver on Linux, specifically nouveau_state.c

                    The nouveau kms driver has the same license as radeon kms, except that copyright is:

                    Copyright 2005 Stephane Marchesin
                    Copyright 2008 Stuart Bennett

                    This driver as a whole is licensed to you as follows: (MIT style license):

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

*/


#include <IOKit/IOLib.h>
#include <IOKit/assert.h>
#include "PsychUserKernelShared.h"
#include "PsychtoolboxKernelDriver.h"

#define super IOService
OSDefineMetaClassAndStructors(PsychtoolboxKernelDriver, IOService)

/* Mappings up to date for May 2016 (last update e-mail patch / commit 2016-05-18). Would need updates for anything after start of June 2016 */

/* Is a given ATI/AMD GPU a DCE11 type ASIC, i.e., with the new display engine? */
bool PsychtoolboxKernelDriver::isDCE11(void)
{
    bool isDCE11 = false;

    // POLARIS10/11 are DCE11.2, but for our purpose we can so far
    // treat them as DCE11.0:

    // POLARIS10: 0x67C0 - 0x67DF
    if ((fPCIDeviceId & 0xFFF0) == 0x67C0) isDCE11 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x67D0) isDCE11 = true;

    // POLARIS11: 0x67E0 - 0x67FF
    if ((fPCIDeviceId & 0xFFF0) == 0x67E0) isDCE11 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x67F0) isDCE11 = true;

    // CARRIZO and STONEY are DCE11 -- This is part of the "Volcanic Islands" GPU family.

    // CARRIZO: 0x987x so far.
    if ((fPCIDeviceId & 0xFFF0) == 0x9870) isDCE11 = true;

    // STONEY: 0x98E4 so far.
    if ((fPCIDeviceId & 0xFFFF) == 0x98E4) isDCE11 = true;

    return(isDCE11);
}

/* Is a given ATI/AMD GPU a DCE10 type ASIC, i.e., with the new display engine? */
bool PsychtoolboxKernelDriver::isDCE10(void)
{
    bool isDCE10 = false;

    // TONGA and FIJI are DCE10 -- This is part of the "Volcanic Islands" GPU family.

    // TONGA: 0x692x - 0x693x so far.
    if ((fPCIDeviceId & 0xFFF0) == 0x6920) isDCE10 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x6930) isDCE10 = true;

    // FIJI in 0x7300 range:
    if ((fPCIDeviceId & 0xFF00) == 0x7300) isDCE10 = true;

    // All DCE11 are also DCE10, so far...
    if (isDCE11()) isDCE10 = true;

    return(isDCE10);
}

/* Is a given ATI/AMD GPU a DCE8 type ASIC, i.e., with the new display engine? */
bool PsychtoolboxKernelDriver::isDCE8(void)
{
    bool isDCE8 = false;

    // Everything >= BONAIRE is DCE8 -- This is part of the "Sea Islands" GPU family.

    // BONAIRE in 0x664x - 0x665x range:
    if ((fPCIDeviceId & 0xFFF0) == 0x6640) isDCE8 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x6650) isDCE8 = true;

    // KABINI in 0x983x range:
    if ((fPCIDeviceId & 0xFFF0) == 0x9830) isDCE8 = true;

    // KAVERI in 0x13xx range:
    if ((fPCIDeviceId & 0xFF00) == 0x1300) isDCE8 = true;

    // HAWAII in 0x67Ax - 0x67Bx range:
    if ((fPCIDeviceId & 0xFFF0) == 0x67A0) isDCE8 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x67B0) isDCE8 = true;

    // MULLINS in 0x985x range:
    if ((fPCIDeviceId & 0xFFF0) == 0x9850) isDCE8 = true;

    return(isDCE8);
}

/* Is a given ATI/AMD GPU a DCE6.4 type ASIC, i.e., with the new display engine? */
bool PsychtoolboxKernelDriver::isDCE64(void)
{
    bool isDCE64 = false;

    // Everything == OLAND is DCE6.4 -- This is part of the "Southern Islands" GPU family.

    // OLAND in 0x66xx range:
    if ((fPCIDeviceId & 0xFF00) == 0x6600) isDCE64 = true;

    return(isDCE64);
}

/* Is a given ATI/AMD GPU a DCE6.1 type ASIC, i.e., with the new display engine? */
bool PsychtoolboxKernelDriver::isDCE61(void)
{
    bool isDCE61 = false;

    // Everything >= ARUBA which is an IGP is DCE6.1 -- This is the "Trinity" GPU family.

    // ARUBA in 0x99xx range: This is the "Trinity" chip family.
    if ((fPCIDeviceId & 0xFF00) == 0x9900) isDCE61 = true;

    // KAVERI in 0x13xx range:
    if ((fPCIDeviceId & 0xFF00) == 0x1300) isDCE61 = true;

    return(isDCE61);
}

/* Is a given ATI/AMD GPU a DCE6 type ASIC, i.e., with the new display engine? */
bool PsychtoolboxKernelDriver::isDCE6(void)
{
    bool isDCE6 = false;

    // Everything >= ARUBA is DCE6 -- This is the "Southern Islands" GPU family.
    // TAHITI is first real DCE-6 core in 0x678x - 0x679x range:
    if ((fPCIDeviceId & 0xFFF0) == 0x6780) isDCE6 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x6790) isDCE6 = true;

    // PITCAIRN, VERDE in 0x6800 - 0x683x range:
    if ((fPCIDeviceId & 0xFFF0) == 0x6800) isDCE6 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x6810) isDCE6 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x6820) isDCE6 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x6830) isDCE6 = true;

    // And one outlier PITCAIRN:
    if ((fPCIDeviceId & 0xFFFF) == 0x684c) isDCE6 = true;

    // Then HAINAN in the 0x666x range:
    if ((fPCIDeviceId & 0xFFF0) == 0x6660) isDCE6 = true;

    // All DCE-6.1 engines are also DCE-6:
    if (isDCE61()) isDCE6 = true;

    // All DCE-6.4 engines are also DCE-6:
    if (isDCE64()) isDCE6 = true;

    // All DCE-8 engines are also DCE-6:
    if (isDCE8()) isDCE6 = true;

    return(isDCE6);
}

/* Is a given ATI/AMD GPU a DCE5 type ASIC, i.e., with the new display engine? */
bool PsychtoolboxKernelDriver::isDCE5(void)
{
    bool isDCE5 = false;

    // Everything after BARTS is DCE5 -- This is the "Northern Islands" GPU family.
    // Barts, Turks, Caicos, Cayman, Antilles in 0x67xx range:
    if ((fPCIDeviceId & 0xFF00) == 0x6700) isDCE5 = true;

    // More Turks ids:
    if ((fPCIDeviceId & 0xFFF0) == 0x6840) isDCE5 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x6850) isDCE5 = true;

    // All DCE-6 engines are also DCE-5:
    if (isDCE6()) isDCE5 = true;

    return(isDCE5);
}

/* Is a given ATI/AMD GPU a DCE-4.1 type ASIC, i.e., with the new display engine? */
bool PsychtoolboxKernelDriver::isDCE41(void)
{
    bool isDCE41 = false;

    // Everything after PALM which is an IGP is DCE-4.1
    // Currently these are Palm, Sumo and Sumo2.
    // DCE-4.1 is a real subset of DCE-4, with all its
    // functionality, except it only has 2 crtcs instead of 6.

    // Palm in 0x980x range:
    if ((fPCIDeviceId & 0xFFF0) == 0x9800) isDCE41 = true;

    // Sumo/Sumo2 in 0x964x range:
    if ((fPCIDeviceId & 0xFFF0) == 0x9640) isDCE41 = true;

    return(isDCE41);
}

/* Is a given ATI/AMD GPU a DCE4 type ASIC, i.e., with the new display engine? */
bool PsychtoolboxKernelDriver::isDCE4(void)
{
    bool isDCE4 = false;

    // Everything after CEDAR is DCE4. The Linux radeon kms driver defines
    // in radeon_family.h which chips are CEDAR or later, and the mapping to
    // these chip codes is done by matching against pci device id's in a
    // mapping table inside linux/include/drm/drm_pciids.h
    // Mapping of chip codes to DCE-generations is in drm/radeon/radeon.h
    // Maintaining a copy of that table is impractical for PTB, so we simply
    // check which range of PCI device id's is covered by the DCE-4 chips and
    // code up matching rules here. This should do for now...

    // Caiman, Cedar, Redwood, Juniper, Cypress, Hemlock in 0x6xxx range:
    if ((fPCIDeviceId & 0xF000) == 0x6000) isDCE4 = true;

    // All DCE-4.1 engines are also DCE-4, except for lower crtc count:
    if (isDCE41()) isDCE4 = true;

    // All DCE-5 engines are also DCE-4:
    if (isDCE5()) isDCE4 = true;

    return(isDCE4);
}

bool PsychtoolboxKernelDriver::isDCE3(void)
{
    bool isDCE3 = false;

    // RV620, RV635, RS780, RS880, RV770, RV710, RV730, RV740,
    // aka roughly HD4330 - HD5165, HD5xxV, and some HD4000 parts.

    if ((fPCIDeviceId & 0xFFF0) == 0x9440) isDCE3 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x9450) isDCE3 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x9460) isDCE3 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x9470) isDCE3 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x9480) isDCE3 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x9490) isDCE3 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x94A0) isDCE3 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x94B0) isDCE3 = true;

    if ((fPCIDeviceId & 0xFFF0) == 0x9540) isDCE3 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x9550) isDCE3 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x9590) isDCE3 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x95C0) isDCE3 = true;

    if ((fPCIDeviceId & 0xFFF0) == 0x9610) isDCE3 = true;

    if ((fPCIDeviceId & 0xFFF0) == 0x9710) isDCE3 = true;

    return(isDCE3);
}

/* The start() method is called at driver load time: It tries to find the MMIO register
 * range of the ATI Radeon X1000/HD2000/HD3000/HDxxx series graphics cards and then
 * memory-maps it into our own virtual address space, so we can read/write registers
 * by simple memory access. It will also query a few registers and print out their
 * values to the system-log.
 *
 */
bool PsychtoolboxKernelDriver::start(IOService* provider)
{
    bool                        success;
    UInt8                       pciBARReg;
    IOMemoryDescriptor *        mem;
    IOMemoryMap *               map;
    IOPhysicalAddress           candidate_phys_addr;
    IOPhysicalLength            candidate_size;
    UInt32                      candidate_count = 0;
    const char*                 providerName = NULL;
    UInt32                      chipset, reg0;

    // Say Hello ;-)
    IOLog("%s::start() called from IOKit. Starting up and trying to attach to GPU:\n", (getName()) ? getName() : "PTB-3");

    // Start our parent provider:
    success = super::start(provider);
    if (!success) return (false);

    // Our provider is start'ed...

    /*
     * Our provider class is specified in the driver property table
     * as IOPCIDevice, so the provider must be of that class.
     * The assert is just to make absolutely sure for debugging.
     */
    assert( OSDynamicCast( IOPCIDevice, provider ));

    // It is a PCI device, check if it's a display driver:
    // Our info.plist for this kext is set up, so it only matches for
    // PCI devices which are graphics cards, and whose vendor is AMD/ATI,
    // so we should be reasonably safe up to that point. Of course it could
    // be something else than the supported X1000, HD2000, HD3000, HDxxx series
    // cards...

    // Double check if we're connected to proper provider:
    // Actually, don't: Keeping track of all the different incompatible names that
    // Apple gives to the provider is too much of a hazzle for small benefit :(
    providerName = provider->getName();
    IOLog("%s: Our provider service is: %s\n", getName(), providerName);

    // ioreg -b -x -p IODeviceTree -n display

    // Ok, store internal reference to our provider:
    fPCIDevice = (IOPCIDevice *) provider;

    // Read PCI device id register with 16 bit device id:
    fPCIDeviceId = fPCIDevice->configRead16(2);

    // Read PCI configuration register 0, a 16 bit register with the
    // Vendor ID:
    fPCIVendorId = fPCIDevice->configRead16(0);

    IOLog("%s: PCI device id is %04x, vendor id is %04x\n", getName(), fPCIDeviceId, fPCIVendorId);

    // Assume two display heads by default:
    fNumDisplayHeads = 2;

    // NVidia card?
    if (PCI_VENDOR_ID_NVIDIA == fPCIVendorId) {
        // Assume it is a NVIDIA GeForce GPU: BAR 0 is the MMIO registers:
        pciBARReg = kIOPCIConfigBaseAddress0;
        fDeviceType = kPsychGeForce;
        IOLog("%s: This is a NVidia GPU, hopefully a compatible GeForce...\n", getName());
    }

    // AMD/ATI card?
    if ((PCI_VENDOR_ID_AMD == fPCIVendorId) || (PCI_VENDOR_ID_ATI == fPCIVendorId)) {
        // Assume it is a ATI Radeon GPU: BAR 2 is the MMIO registers for old GPU's, BAR 5 for DCE-8 "Sea Islands" and later:
        pciBARReg = (isDCE8() || isDCE10()) ? kIOPCIConfigBaseAddress5 : kIOPCIConfigBaseAddress2;
        fDeviceType = kPsychRadeon;
        IOLog("%s: This is a ATI/AMD GPU, hopefully a compatible Radeon...\n", getName());
        if (PCI_VENDOR_ID_ATI == fPCIDevice->configRead16(0)) IOLog("%s: Confirmed to have ATI's vendor id.\n", getName());
        if (PCI_VENDOR_ID_AMD == fPCIDevice->configRead16(0)) IOLog("%s: Confirmed to have AMD's vendor id.\n", getName());

        IOLog("%s: This is a GPU with %s display engine.\n", getName(), isDCE11() ? "DCE-11" : isDCE10() ? "DCE-10" : isDCE8() ? "DCE-8" : isDCE6() ? "DCE-6" : (isDCE5() ? "DCE-5" : (isDCE4() ? "DCE-4" : (isDCE3() ? "DCE-3" : "AVIVO"))));

        // Setup for DCE-4/5/6/8:
        if (isDCE4() || isDCE5() || isDCE6() || isDCE8()) {
            fRadeonLowlimit = 0;

            // Offset of crtc blocks of evergreen gpu's for each of the six possible crtc's:
            crtcoff[0] = EVERGREEN_CRTC0_REGISTER_OFFSET;
            crtcoff[1] = EVERGREEN_CRTC1_REGISTER_OFFSET;
            crtcoff[2] = EVERGREEN_CRTC2_REGISTER_OFFSET;
            crtcoff[3] = EVERGREEN_CRTC3_REGISTER_OFFSET;
            crtcoff[4] = EVERGREEN_CRTC4_REGISTER_OFFSET;
            crtcoff[5] = EVERGREEN_CRTC5_REGISTER_OFFSET;

            // Also, DCE-4 and DCE-5 and DCE-6, but not DCE-4.1 or DCE-6.4 (which have only 2) or DCE-6.1 (4 heads), supports up to six display heads:
            if (!isDCE41() && !isDCE61() && !isDCE64()) fNumDisplayHeads = 6;

            // DCE-6.1 "Trinity" chip family supports 4 display heads:
            if (!isDCE41() && isDCE61()) fNumDisplayHeads = 4;
        }

        // Setup for DCE-10/11:
        if (isDCE10() || isDCE11()) {
            // DCE-10/11 of the "Volcanic Islands" gpu family uses (mostly) the same register specs,
            // but the offsets for the different CRTC blocks are different wrt. to pre DCE-10. Therefore
            // need to initialize the offsets differently. Also, some of these parts seem to support up
            // to 7 display engines instead of the old limit of 6 engines:
            fRadeonLowlimit = 0;

            // Offset of crtc blocks of Volcanic Islands DCE-10/11 gpu's for each of the possible crtc's:
            crtcoff[0] = DCE10_CRTC0_REGISTER_OFFSET;
            crtcoff[1] = DCE10_CRTC1_REGISTER_OFFSET;
            crtcoff[2] = DCE10_CRTC2_REGISTER_OFFSET;
            crtcoff[3] = DCE10_CRTC3_REGISTER_OFFSET;
            crtcoff[4] = DCE10_CRTC4_REGISTER_OFFSET;
            crtcoff[5] = DCE10_CRTC5_REGISTER_OFFSET;
            crtcoff[6] = DCE10_CRTC6_REGISTER_OFFSET;

            // DCE-10 has 6 display controllers:
            if (isDCE10()) fNumDisplayHeads = 6;

            // DCE-11 has 3 display controllers:
            if (isDCE11()) fNumDisplayHeads = 3;
        }
    }

    // Intel card?
    if (PCI_VENDOR_ID_INTEL == fPCIVendorId) {
        // Assume it is an Intel IGP: BAR 0 is the MMIO registers:
        pciBARReg = kIOPCIConfigBaseAddress0;
        fDeviceType = kPsychIntelIGP;

        // GEN-7+ (IvyBridge and later) and maybe GEN-6 (SandyBridge) has 3 display
        // heads, older IGP's have 2. Let's be optimistic and assume 3, to safe us
        // from lots of new detection code:
        fNumDisplayHeads = 3;

        IOLog("%s: This is a Intel GPU, hopefully a compatible one...\n", getName());
    }

    // GPU type and vendor detected?
    if (kPsychUnknown == fDeviceType) {
        // Failed! Cleanup and exit:
        IOLog("%s: This is not a supported GPU from AMD, NVidia or Intel! Will abort here - this is not a safe ground for proceeding.\n", getName());

        // Detach our device handle, so our cleanup routine won't do anything later on:
        fPCIDevice = NULL;

        // We will exit with "success", but the fPCIDevice == NULL signals failure whenever any
        // of our methods gets called, so we are basically idle from now on...
        return(success);
    }

    /*
     * Enable memory response from the card
     */
    fPCIDevice->setMemoryEnable( true );

    /*
     * Log some info about the device
     */

    /*    // DEBUG CODE: Detect both parents and clients of our provider... Left here for documentation, not useful in current impl.
    OSIterator * clientiter = provider->getClientIterator();
    OSObject* client;

    for (client = clientiter->getNextObject(); client; client = clientiter->getNextObject()) {
        altProvider = (IOService*) client;
        IOLog("Client is %s\n", altProvider->getName());
        if (altProvider->compareName(OSString::withCStringNoCopy("ATY,Wormy"))) break;
    }

    clientiter->release();
    if (altProvider->compareName(OSString::withCStringNoCopy("ATY,Wormy"))) {
        IOLog("Radeon driver found!\n");
    }
    else {
        altProvider = NULL;
    }

    //    altProvider = provider->getClient();
    //    if (altProvider) altProvider = altProvider->getProvider();
    if (altProvider) {
        IOLog("Our providers provider is %s\n", altProvider->getName());


        int irq = 0;
        int irqtype;
        IOReturn rc;
        while(irq < 10 && kIOReturnNoInterrupt != (rc=altProvider->getInterruptType(irq, &irqtype))) {
            IOLog("IRQ %ld available: Type %ld Status %lx (%s)\n", irq, irqtype, rc, (irqtype == kIOInterruptTypeLevel) ? "Leveltriggered" : "Other"); irq++;
        }

        rc = altProvider->registerInterrupt(0, this, OSMemberFunctionCast(IOInterruptAction, this, &PsychtoolboxKernelDriver::fastPathInterruptHandler));
        IOLog("Resultcode for reg. interrupt 0: %lx (%s)\n", rc, (rc == kIOReturnNoResources) ? "No ressources" : "Other");
        if (rc == kIOReturnSuccess) {
            altProvider->enableInterrupt(0);
            IOSleep(1000);
            altProvider->disableInterrupt(0);
            altProvider->unregisterInterrupt(0);
        }

    }

*/

    // Perform complex probing and double-checking on Radeon GPU's:
    if (fDeviceType == kPsychRadeon) {
        /* print all the device's memory ranges */
        for( UInt32 index = 0;
            index < fPCIDevice->getDeviceMemoryCount();
            index++ ) {

            mem = fPCIDevice->getDeviceMemoryWithIndex( index );
            assert( mem );
            // Disabled to handle Apple braindamage: IOLog("%s: PCI Range[%ld] %08lx:%08lx\n", getName(), index, mem->getPhysicalAddress(), mem->getLength());
            IOLog("%s: PCI Range[%d] Size " ByteCount_FORMAT "\n", getName(), index, mem->getLength());

            // Find the MMIO range of expected size around 0x10000 - 0x20000: We find the one with size > 0x1000 and not bigger
            // than 0x40000, ie. in the range 4 Kb < size < 256 Kb. This likely represents the register space.

            // Much bigger than 0x40000 would be either 3D engine command FIFO's or VRAM framebuffer. Much smaller
            // would be the PCI config space registerset (or maybe VGA set?!?)
            if (mem->getLength() >= 0x1000 && mem->getLength() <= 0x40000) {
                // A possible candidate:
                // Disabled to handle Apple braindamage: candidate_phys_addr = mem->getPhysicalAddress();
                candidate_phys_addr = 0xdeadbeef;
                candidate_size = mem->getLength();
                candidate_count++;
                // Disabled to handle Apple braindamage: IOLog("%s: ==> AMD/ATI Radeon PCI Range[%ld] %08lx:%08lx is %ld. candidate for register block.\n", getName(), index, mem->getPhysicalAddress(), mem->getLength(), candidate_count);
                IOLog("%s: ==> AMD/ATI Radeon PCI Range[%d] Size " ByteCount_FORMAT " is %d. candidate for register block.\n", getName(), index, mem->getLength(), candidate_count);
            }
        }

        /* More than one candidate found? */
        if (candidate_count != 1) {
            IOLog("%s: Found %d candidates for Radeon register block.\n", getName(), candidate_count);
        }

        /* look up a range based on its config space base address register */
        mem = fPCIDevice->getDeviceMemoryWithRegister(pciBARReg);
        if( mem ) {
            // Disabled to handle Apple braindamage: IOLog("%s: Range@0x%x %08lx:%08lx\n", getName(), pciBARReg, mem->getPhysicalAddress(), mem->getLength());
            IOLog("%s: Range@0x%x Length " ByteCount_FORMAT "\n", getName(), pciBARReg, mem->getLength());
        }
        else {
            IOLog("%s: Could not find MMIO mapping for config base address register 0x%x!\n", getName(), pciBARReg);
        }

        /* Valid mem'ory address returned for mapping? */
        if (mem == NULL) {
            IOLog("%s: Could not resolve relevant MMIO register block!! Will abort here - this is not a safe ground for proceeding.\n", getName());

            // Detach our device handle, so our cleanup routine won't do anything later on:
            fPCIDevice = NULL;

            // We will exit with "success", but the fPCIDevice == NULL signals failure whenever any
            // of our methods gets called, so we are basically idle from now on...
            return(success);
        }

        // End of Radeon specific probing.
    }

    // Ok, we think we have a good candidate for the MMIO block in "mem". Let's try to map it
    // into our address space...

    /* map a range based on its config space base address register,
     * this is how the driver gets access to its memory mapped registers
     * the getVirtualAddress() method returns a kernel virtual address
     * for the register mapping */
    map = fPCIDevice->mapDeviceMemoryWithRegister( pciBARReg );
    if( map ) {
        // Disabled to handle Apple braindamage: IOLog("%s: GPU MMIO register block Range@0x%x (%08lx) mapped to kernel virtual address %08x\n", getName(), pciBARReg, map->getPhysicalAddress(), map->getVirtualAddress());
        IOLog("%s: GPU MMIO register block BAR@0x%x mapped to kernel virtual address " VirtAddr_FORMAT "\n", getName(), pciBARReg, map->getVirtualAddress());

        /* Assign the map object, and the mapping itself to our private members: */
        fRadeonMap = map;
        fRadeonRegs = map->getVirtualAddress();

        if (map->getLength() > 0xffffffff) {
            IOLog("%s: WARNING: MMIO register block size greater than 4 GB! Will clamp to 4GB aka 0xffffffff. Fix and recompile me for larger MMIO space!\n", getName());
            fRadeonSize = (UInt32) 0xffffffff;
        }
        else {
            fRadeonSize = (UInt32) map->getLength();
        }
    }
    else {
        // Failed! Cleanup and exit:
        IOLog("%s: Could not memory-map relevant MMIO register block into my kernel memory address space! Will abort here - this is not a safe ground for proceeding.\n", getName());

        // Detach our device handle, so our cleanup routine won't do anything later on:
        fPCIDevice = NULL;

        // We will exit with "success", but the fPCIDevice == NULL signals failure whenever any
        // of our methods gets called, so we are basically idle from now on...
        return(success);
    }

    // Execute an I/O memory barrier, so order of execution and CPU <-> GPU sync is guaranteed:
    OSSynchronizeIO();

    // NVidia GPU? If so, we can detect specific chipset family:
    if (fDeviceType == kPsychGeForce) {
        // Debug check: Report if it is a big-endian configured GPU:
        if(ReadRegister(NV03_PMC_BOOT_1)) IOLog("%s: Big Endian NVidia GPU detected.\n", getName());

        // Get hardware id code from gpu register:
        reg0 = ReadRegister(NV03_PMC_BOOT_0);

        /* We're dealing with >=NV10 */
        if ((reg0 & 0x0f000000) > 0) {
            /* Bit 27-20 contain the architecture in hex */
            chipset = (reg0 & 0xff00000) >> 20;
            /* NV04 or NV05 */
        } else if ((reg0 & 0xff00fff0) == 0x20004000) {
            if (reg0 & 0x00f00000)
                chipset = 0x05;
            else
                chipset = 0x04;
        } else {
            chipset = 0xff;
        }

        switch (chipset & 0xf0) {
            case 0x00:
                // NV_04/05: RivaTNT , RivaTNT2
                fCardType = 0x04;
                break;
            case 0x10:
            case 0x20:
            case 0x30:
                // NV30 or earlier: GeForce-5 / GeForceFX and earlier:
                fCardType = chipset & 0xf0;
                break;
            case 0x40:
            case 0x60:
                // NV40: GeForce6/7 series: "Curie"
                fCardType = 0x40;
                break;
            case 0x50:
            case 0x80:
            case 0x90:
            case 0xa0:
                // NV50: GeForce8/9/Gxxx: "Tesla"
                fCardType = 0x50;
                break;
            case 0xc0:
                // NVC0: GeForce-400/500: "Fermi"
                fCardType = 0xc0;
                break;
            case 0xd0:
            case 0xe0:
            case 0xf0:
                // NVE0: GeForce-600: "Kepler"
                fCardType = 0xe0;
                break;
            default:
                IOLog("%s: Unknown NVidia chipset 0x%08x. Optimistically assuming latest generation GPU too new to be known.\n", getName(), reg0);
                fCardType = 0x00;
        }

        // "Kepler" chip family and later supports 4 display heads:
        if ((fCardType == 0x0) || (fCardType >= 0xe0)) fNumDisplayHeads = 4;

        if (fCardType > 0x00) IOLog("%s: NV-%02x GPU with %d display heads detected.\n", getName(), fCardType, fNumDisplayHeads);
    }

    // The following code chunk if enabled, will detaching the Radeon driver IRQ handler and
    // instead attach our own fastPathInterruptHandler() to IRQ 1 of Radeon GPU's.
    //
    // Detaching the Radeon driver seems to do "limited damage" to OS/X operation:
    //
    // 1. Our VBL userspace timestamping doesn't work anymore, because it relies on VBL timestamps from the driver
    //    This gets detected, handled and reported by PTB's consistency checks -- So this is a good method of
    //    fault-injection for testing PTB's self-diagnostic.
    //
    // 2. The window-server behaves more "nervous", e.g., mouse movements become slightly more jerky,
    //    probably due to the server not getting any timestamps from the driver anymore, so it can't
    //    smooth out display updates in time for a smooth experience. However, OpenGL/PTB seems to operate
    //    normally and the general system/GUI is still very useable -- A noticeable, but not really problematic
    //    nuisance.
    //
    //      Of course, this behaviour may change on each different gfx-card and OS/X release, so it is
    //    not desirable for production use -- Tradeoffs all over the place...
    if (false && (fDeviceType == kPsychRadeon) && !isDCE4() && !isDCE10()) {
        // Setup our own interrupt handler for gfx-card interrupts:
        if (!InitializeInterruptHandler()) {
            // Failed!
            IOLog("%s: Could not install our snoop - interrupt handler! This is non-fatal, will continue...\n", getName());
        }
        else {
            // :-)
            IOLog("%s: Graphics card IRQ Snoop - interrupt handler installed. Cool :-)\n", getName());
        }
    }

    // We should be ready...
    IOLog("\n");
    IOLog("%s: Psychtoolbox-3 kernel-level support driver V1.10 (Revision %d) for ATI/AMD/NVidia/Intel GPU's ready for use!\n", getName(), PTBKDRevision);
    IOLog("%s: This driver is copyright 2008 - 2015 Mario Kleiner and the Psychtoolbox-3 project developers.\n", getName());
    IOLog("%s: The driver is licensed to you under the MIT free and open-source software license.\n", getName());
    IOLog("%s: See the file License.txt in the Psychtoolbox root installation folder for details.\n", getName());
    IOLog("%s: The driver contains bits of code derived from the free software nouveau and radeon kms drivers on Linux.\n", getName());
    IOLog("%s: See driver source code for specific copyright notices of the compatible licenses of those bits.\n", getName());

    // If everything worked, we publish ourselves and our services:
    if (success) {
        // Publish ourselves so clients can find us
        registerService();
    }

    // Return final startup state:
    return success;
}


/* The stop() method is called at driver shutdown time: It checks if we're attached to
 * a gfx-card. If so, we unmap the register control block, detach, and clean up after
 * ourselves.
 */
void PsychtoolboxKernelDriver::stop(IOService* provider)
{
    IOWorkLoop* myWorkLoop = getWorkLoop();

    IOLog("%s::stop() request from IOKit: Shutting down...\n", getName());

    // Perform cleanup:
    if (fPCIDevice) {
        // Disable, detach and delete our interrupt handler, if any:
        if (fInterruptSrc) {
            fInterruptSrc->disable();
            // TODO: Some kind of release() function or detach from Workloop needed here?
            myWorkLoop->removeEventSource(fInterruptSrc);
            fInterruptSrc->release();

            IOLog("%s::stop(): Final total interrupt count is %d.\n", getName(), fInterruptCounter);
            IOLog("%s::stop(): Final VBL interrupt counts are Head[0] = %d, Head[1] = %d \n", getName(), fVBLCounter[0], fVBLCounter[1]);
            fastPathInterruptHandler(NULL, NULL);

            fInterruptSrc = NULL;
            fInterruptCookie = 0xdeadbeef;
            fInterruptCounter = 0;
        }

        // Detach our device handle:
        fPCIDevice = NULL;

        // Release memory mapping, if any:
        if (fRadeonMap) {
            IOLog("%s::stop(): Releasing MMIO memory mapped GPU register block...\n", getName());
            fRadeonMap->release();
        }

        // NULL-Out our now stale references:
        fRadeonMap = NULL;
        fRadeonRegs = NULL;
        fRadeonSize = 0;
    }

    IOLog("%s::stop(): Driver is stopped and ready for jettisoning or restart...\n", getName());

    // Call our parents stop routine:
    super::stop(provider);
}


/* The init() method is called at the very beginning of the driver life-cycle: Set our member
 * variables to safe and clean defaults:
 */
bool PsychtoolboxKernelDriver::init(OSDictionary* dictionary)
{
    UInt32 i;

    IOLog("PTB::init() called at driver load time...\n");

    if (!super::init(dictionary)) {
        return false;
    }

    fDeviceType = kPsychUnknown;
    fCardType = 0x00;
    fPCIDeviceId = 0x0000;
    fPCIDevice = NULL;
    fRadeonMap = NULL;
    fRadeonRegs = NULL;
    fRadeonSize = 0;
    fRadeonLowlimit = 0;
    fInterruptSrc = NULL;
    fInterruptCookie = 0xdeadbeef;
    fInterruptCounter = 0;
    fVBLCounter[0] = 0;
    fVBLCounter[1] = 0;
    for (i = 0; i <= MAXHEADID; i++) {
        oldDither[i] = 0;
        crtcoff[i] = 0;
    }

    return true;
}


// We override free only to log that it has been called to make it easier to follow the driver's lifecycle.
void PsychtoolboxKernelDriver::free(void)
{
    IOLog("%s::free() called at module jettisoning time -- Shutdown complete. Bye!\n", getName());
    super::free();
}


// We override probe only to log that it has been called to make it easier to follow the driver's lifecycle.
// Production drivers can override probe if they need to make an active decision whether the driver is appropriate for
// the provider.
IOService* PsychtoolboxKernelDriver::probe(IOService* provider, SInt32* score)
{
    IOLog("PTB::probe() request from IOKit received...\n");
    IOService *res = super::probe(provider, score);
    return res;
}


// We override willTerminate only to log that it has been called to make it easier to follow the driver's lifecycle.
//
// willTerminate is called at the beginning of the termination process. It is a notification
// that a provider has been terminated, sent before recursing up the stack, in root-to-leaf order.
//
// This is where any pending I/O should be terminated. At this point the user client has been marked
// inactive and any further requests from the user process should be returned with an error.
bool PsychtoolboxKernelDriver::willTerminate(IOService* provider, IOOptionBits options)
{
    if (false) IOLog("%s::willTerminate()\n", getName());
    return super::willTerminate(provider, options);
}


// We override didTerminate only to log that it has been called to make it easier to follow the driver's lifecycle.
//
// didTerminate is called at the end of the termination process. It is a notification
// that a provider has been terminated, sent after recursing up the stack, in leaf-to-root order.
bool PsychtoolboxKernelDriver::didTerminate(IOService* provider, IOOptionBits options, bool* defer)
{
    if (false) IOLog("%s::didTerminate()\n", getName());
    return super::didTerminate(provider, options, defer);
}


// We override terminate only to log that it has been called to make it easier to follow the driver's lifecycle.
// Production drivers will rarely need to override terminate. Termination processing should be done in
// willTerminate or didTerminate instead.
bool PsychtoolboxKernelDriver::terminate(IOOptionBits options)
{
    bool    success;
    if (false) IOLog("%s::terminate()\n", getName());
    success = super::terminate(options);
    return success;
}


// We override finalize only to log that it has been called to make it easier to follow the driver's lifecycle.
// Production drivers will rarely need to override finalize.
bool PsychtoolboxKernelDriver::finalize(IOOptionBits options)
{
    bool    success;
    if (false) IOLog("%s::finalize()\n", getName());
    success = super::finalize(options);
    return success;
}

// Initialize our own interrupt handler for snooping on gfx-card state changes
// that are interesting to Psychtoolbox, e.g., buffer swaps. This handler can
// also trigger event-based actions, e.g., emit a hardware trigger signal in
// response to a completed double-buffer swap:
bool PsychtoolboxKernelDriver::InitializeInterruptHandler(void)
{
    // Reset our special debug-members:
    fInterruptCookie = 0xdeadbeef;
    fInterruptCounter = 0;

    // Get handle to our assigned IOKit workloop:
    IOWorkLoop * myWorkLoop = (IOWorkLoop *) getWorkLoop();
    if (!myWorkLoop) {
        // Failed!
        IOLog("%s: In IRQ handler setup: Failed to get my workloop!\n", getName());
        return(false);
    }

    // Dump current IRQ status of GPU:
    IOLog("%s: In IRQ handler setup: Initial: Current hw irqControl is %x.\n", getName(), ReadRegister(RADEON_R500_GEN_INT_CNTL));
    IOLog("%s: In IRQ handler setup: Initial: Current hw irqStatus is  %x.\n",    getName(), ReadRegister(RADEON_R500_GEN_INT_STATUS));

    // Disable currently attached ATI interrupt handler on IRQ1 and detach it:
    if (kIOReturnSuccess != fPCIDevice->disableInterrupt(1)) {
        IOLog("%s: In IRQ handler setup: Failed to disable IRQ1 on Display provider!\n", getName());
        return(false);
    }

    // Dump current IRQ status of GPU:
    IOLog("%s: In IRQ handler setup: After disable IRQ: Current hw irqControl is %x.\n", getName(), ReadRegister(RADEON_R500_GEN_INT_CNTL));
    IOLog("%s: In IRQ handler setup: After disable IRQ: Current hw irqStatus is  %x.\n",    getName(), ReadRegister(RADEON_R500_GEN_INT_STATUS));

    if (kIOReturnSuccess != fPCIDevice->unregisterInterrupt(1)) {
        fPCIDevice->enableInterrupt(1);
        IOLog("%s: In IRQ handler setup: Failed to detach IRQ1 handler on Display provider! Old handler reenabled...\n", getName());
        return(false);
    }

    IOLog("%s: In IRQ handler setup: ATI IRQ1 handler on Display provider permanently detached...\n", getName());

    // Dump current IRQ status of GPU:
    IOLog("%s: In IRQ handler setup: After detach IRQ: Current hw irqControl is %x.\n", getName(), ReadRegister(RADEON_R500_GEN_INT_CNTL));
    IOLog("%s: In IRQ handler setup: After detach IRQ: Current hw irqStatus is  %x.\n",    getName(), ReadRegister(RADEON_R500_GEN_INT_STATUS));

    // Option A: Attach a filter handler in the primary interrupt path, and a deferred workloop handler
    // in the workloop --> More flexible than option B, but both don't yield useable results at the end
    // of the day :-(
    //
    // Create and register an interrupt event source. Our PCI provider will
    // take care of the low-level interrupt registration stuff.
    fInterruptSrc = IOFilterInterruptEventSource::filterInterruptEventSource(this,
                    OSMemberFunctionCast(IOInterruptEventAction, this, &PsychtoolboxKernelDriver::workLoopInterruptHandler),
                    OSMemberFunctionCast(IOFilterInterruptAction, this, &PsychtoolboxKernelDriver::fastPathInterruptHandler),
                    fPCIDevice,
                    1);
    if (!fInterruptSrc) {
        // Failed!
        IOLog("%s: In IRQ handler setup: Failed to get and attach an interrupt source for IRQ 1!\n", getName());
        return(false);
    }

    // Attach our interrupt source/handlers to the workloop:
    if (myWorkLoop->addEventSource(fInterruptSrc) != kIOReturnSuccess) {
        // TODO: Some kind of release() function or detach from Workloop needed here?
        fInterruptSrc = NULL;
        IOLog("%s: In IRQ handler setup: Failed to attach our interrupt source to our workloop!\n", getName());
        return(false);
    }


/*  // Option B: Direct attachment to primary interrupt: Works as good as option A, but maybe A is a bit cleaner...
    if (kIOReturnSuccess != fPCIDevice->registerInterrupt(1, this, OSMemberFunctionCast(IOInterruptAction, this, &PsychtoolboxKernelDriver::fastPathInterruptHandler))) {
        IOLog("%s: In IRQ handler setup: Failed to attach our interrupt handler for Radeon IRQ1 on Display provider!\n", getName());
        return(false);
    }
*/

    if (kIOReturnSuccess != fPCIDevice->enableInterrupt(1)) {
        IOLog("%s: In IRQ handler setup: Failed to reenable IRQ1 on Display provider!\n", getName());
        return(false);
    }

    // Ok, ready to rock: Enable interrupt handling:
    fInterruptSrc->enable();
    myWorkLoop->enableAllInterrupts();

    // Dump current IRQ status of GPU:
    IOLog("%s: In IRQ handler setup: After setup IRQ: Current hw irqControl is %x.\n", getName(), ReadRegister(RADEON_R500_GEN_INT_CNTL));
    IOLog("%s: In IRQ handler setup: After setup IRQ: Current hw irqStatus is  %x.\n",    getName(), ReadRegister(RADEON_R500_GEN_INT_STATUS));

    // Acknowledge all pending IRQ's via a special-call to our fastPathInterruptHandler:
    fastPathInterruptHandler(NULL, (IOFilterInterruptEventSource*) 0x1);

    // Return success status:
    return(true);
}

// Fast-Path interrupt handler: Gets called in direct primary hardware interrupt context!
// --> CAUTION - Think twice before doing anything here!!!
bool PsychtoolboxKernelDriver::fastPathInterruptHandler(OSObject* myself, IOFilterInterruptEventSource* mySource)
{
    static UInt32 myCounter = 0;
    static UInt32 irqStatus = 0xdeadbeef;
    static UInt32 displayirqStatus = 0;

    // Special readout request from outside interrupt path?
    if (myself == NULL && mySource == NULL) {
        // Yes: Log current internal counters value to log:
        IOLog("Internal fastPath invocation count [total number of hardware interrupts on this line] is %d.\n", myCounter);
        IOLog("Internal fastPath irqStatus is %x.\n", irqStatus);
        IOLog("Current hw irqControl is %x.\n",    ReadRegister(RADEON_R500_GEN_INT_CNTL));
        IOLog("Current hw irqStatus is  %x.\n",    ReadRegister(RADEON_R500_GEN_INT_STATUS));

        return(false);
    }

    // Regular call from interrupt subsystem. Do your thing...

    // Increment internal counter for debug purpose: This one can be read out at
    // shutdown time by a fastPathInterruptHandler(NULL, NULL) call, *after* the
    // handler has been detached from actual interrupt dispatch system!
    myCounter++;

    // Check if we have an address space reference to our member variables. If the
    // cookie 0xdeadbeef is found, we can (hopefully) safely access our member
    // variables. If that isn't the case, we're unable to access our members and
    // therefore to do a single useful thing.
    if(fInterruptCookie == 0xdeadbeef) {
        // Ok, member access seems to be safe. Do our actual job:

        // Retrieve the current interrupt status of the Radeon card:
        irqStatus = 0;
        irqStatus = ReadRegister(RADEON_R500_GEN_INT_STATUS);

        // We are specifically interested in display interrupts: Any such?
        if (irqStatus & RADEON_R500_DISPLAY_INT_STATUS) {
            // Display interrupt! Which one?
            displayirqStatus = ReadRegister(RADEON_R500_DISP_INTERRUPT_STATUS);

            // Handle IRQ acknowledge and triggered actions per display pipe:
            if (displayirqStatus & RADEON_R500_D1_VBLANK_INTERRUPT) {
                // Display head 0 [Pipe 1]:
                WriteRegister(RADEON_R500_D1MODE_VBLANK_STATUS, RADEON_R500_VBLANK_ACK);

                // Perform whatever actions need to be performed:
                handleVBLIRQ(0);
            }

            if (displayirqStatus & RADEON_R500_D2_VBLANK_INTERRUPT) {
                // Display head 1 [Pipe 2]:
                WriteRegister(RADEON_R500_D2MODE_VBLANK_STATUS, RADEON_R500_VBLANK_ACK);

                // Perform whatever actions need to be performed:
                handleVBLIRQ(1);
            }
        }

        // Acknowledge IRQs, if any:
        if (irqStatus) WriteRegister(RADEON_R500_GEN_INT_STATUS, irqStatus);

        // Increment the gloabl fInterruptCounter which counts the graphics cards interrupts we're interested in:
        fInterruptCounter++;
    }

    // We return false, so no secondary workloop handler is triggered and hardware
    // interrupts don't get disabled:
    return(false);
}

// Slow-Path WorkLoop interrupt handler: Gets called if the fast-path handler returns "true".
// Will be scheduled/executed by the IOWorkloop with some (possibly significant, non-deterministic)
// latency after the real hardware event has occured. Can use more of IOKit API without disaster.
void PsychtoolboxKernelDriver::workLoopInterruptHandler(OSObject* myself, IOInterruptEventSource* mySource, int pendingIRQs)
{
    // Nothing to do yet: Just return.
    return;
}

// Handle VBLANK IRQ's for display head 'headId':
void PsychtoolboxKernelDriver::handleVBLIRQ(UInt32 headId)
{
    // Increment per-head counter:
    fVBLCounter[headId]++;

    return;
}

// Register access convenience functions:

// Read 32 bit control register at 'offset':
UInt32    PsychtoolboxKernelDriver::ReadRegister(UInt32 offset)
{
    // Safety check: Don't allow reads past devices MMIO range:
    // We don't return error codes and don't log the problem,
    // because we could be called from primary Interrupt path, so IOLog() is not
    // an option!
    if (fRadeonRegs == NULL || offset < fRadeonLowlimit || offset >= fRadeonSize-4) return(0);

    // Read and return value:

    // Radeon: Don't know endianity behaviour: Play save, stick to LE assumption for now:
    if (fDeviceType == kPsychRadeon) return(OSReadLittleInt32((void*) fRadeonRegs, offset));

    // Read the register in native byte order: At least NVidia GPU's adapt their
    // endianity to match the host systems endianity, so no need for conversion:
    if (fDeviceType == kPsychGeForce) return(_OSReadInt32((void*) fRadeonRegs, offset));
    if (fDeviceType == kPsychIntelIGP) return(_OSReadInt32((void*) fRadeonRegs, offset));

    return(0);
}

// Write 32 bit control register at 'offset' with 'value':
void PsychtoolboxKernelDriver::WriteRegister(UInt32 offset, UInt32 value)
{
    // Safety check: Don't allow writes past devices MMIO range:
    // We don't return error codes and don't log the problem,
    // because we could be called from primary Interrupt path, so IOLog() is not
    // an option!
    if (fRadeonRegs == NULL || offset < fRadeonLowlimit || offset >= fRadeonSize-4) return;

    // Write the register in native byte order: At least NVidia GPU's adapt their
    // endianity to match the host systems endianity, so no need for conversion:
    if (fDeviceType == kPsychGeForce) _OSWriteInt32((void*) fRadeonRegs, offset, value);
    if (fDeviceType == kPsychIntelIGP) _OSWriteInt32((void*) fRadeonRegs, offset, value);

    // Radeon: Don't know endianity behaviour: Play save, stick to LE assumption for now:
    if (fDeviceType == kPsychRadeon) OSWriteLittleInt32((void*) fRadeonRegs, offset, value);

    // Execute memory I/O barrier to make sure that writes happen in-order and properly synced between CPU and GPU/PCI card:
    OSSynchronizeIO();

    return;
}

// Try to change hardware dither mode on GPU:
void PsychtoolboxKernelDriver::SetDitherMode(UInt32 headId, UInt32 ditherOn)
{
    UInt32 reg;

    // AMD/ATI Radeon, FireGL or FirePro GPU?
    if (fDeviceType == kPsychRadeon) {
        IOLog("%s: SetDitherMode: Trying to %s digital display dithering on display head %d.\n", getName(), (ditherOn) ? "enable" : "disable", headId);

        // Map headId to proper hardware control register offset:
        if (isDCE4() || isDCE5() || isDCE6() || isDCE8() || isDCE10()) {
            // DCE-4 display engine (CEDAR and later afaik): Up to six crtc's. Map to proper
            // register offset for this headId:
            if (headId > fNumDisplayHeads - 1) {
                // Invalid head - bail:
                IOLog("%s: SetDitherMode: ERROR! Invalid headId %d provided. Must be between 0 and %d. Aborted.\n", getName(), headId, fNumDisplayHeads - 1);
                return;
            }

            // Map to dither format control register for head 'headId':
            reg = EVERGREEN_FMT_BIT_DEPTH_CONTROL + crtcoff[headId];
        } else if (isDCE3()) {
            // DCE-3 display engine for R700: HD4330 - HD5165, HD5xxV, and some R600's:
            reg = (headId == 0) ? DCE3_FMT_BIT_DEPTH_CONTROL : DCE3_FMT_BIT_DEPTH_CONTROL + 0x800;
        } else {
            // AVIVO / DCE-1 / DCE-2 display engine (R300 - R600 afaik): At most two display heads for dual-head gpu's.

            // These have a weird wiring of encoders/transmitters to output connectors with no simple 1:1 correspondence
            // between crtc's and encoders. We need to probe each encoder block if it is enabled and sourcing from our headId,
            // respective its corresponding crtc to find which encoder block needs to be configured wrt. dithering on this
            // display headId:
            reg = 0x0;

            // TMDSA block enabled, and driven by headId? Then we control its encoder:
            if ((ReadRegister(0x7880) & 0x1) && ((ReadRegister(0x7884) & 0x1) == headId)) reg = RADEON_TMDSA_BIT_DEPTH_CONTROL;

            // LVTMA block enabled, and driven by headId? Then we control its encoder:
            if ((ReadRegister(0x7A80) & 0x1) && ((ReadRegister(0x7A84) & 0x1) == headId)) reg = RADEON_LVTMA_BIT_DEPTH_CONTROL;

            // DVOA block enabled, and driven by headId? Then we control its encoder:
            if ((ReadRegister(0x7980) & 0x1) && ((ReadRegister(0x7984) & 0x1) == headId)) reg = RADEON_DVOA_BIT_DEPTH_CONTROL;

            // If no digital encoder block was assigned, then this likely means we're connected to a
            // analog VGA monitor driven by the DAC. The DAC doesn't have dithering ever, so we are
            // done with a simple no-op:
            if (reg == 0x0) {
                IOLog("%s: SetDitherMode: Head %i connected to analog VGA DAC. Dithering control skipped.\n", getName(), headId);
                return;
            }
            else {
                switch (reg) {
                    case RADEON_TMDSA_BIT_DEPTH_CONTROL:
                        IOLog("%s: SetDitherMode: Head %i connected to TMDSA block.\n", getName(), headId);
                        break;

                    case RADEON_LVTMA_BIT_DEPTH_CONTROL:
                        IOLog("%s: SetDitherMode: Head %i connected to LVTMA block.\n", getName(), headId);
                        break;

                    case RADEON_DVOA_BIT_DEPTH_CONTROL:
                        IOLog("%s: SetDitherMode: Head %i connected to DVOA block.\n", getName(), headId);
                        break;
                }
            }
        }

        // Perform actual enable/disable/reconfigure sequence for target encoder/head:

        // Enable dithering?
        if (ditherOn) {
            // Reenable dithering with old, previously stored settings, if it is disabled:

            // Dithering currently off (all zeros)?
            if (ReadRegister(reg) == 0) {
                // Dithering is currently off. Do we know the old setting from a previous
                // disable?
                if (oldDither[headId] > 0) {
                    // Yes: Restore old "factory settings":
                    IOLog("%s: SetDitherMode: Dithering previously disabled by us. Reenabling with old control setting 0x%x.\n", getName(), oldDither[headId]);
                    WriteRegister(reg, oldDither[headId]);
                }
                else {
                    // No: Dithering was disabled all the time, so we don't know the
                    // OS defaults. Use the numeric value of 'ditherOn' itself:
                    IOLog("%s: SetDitherMode: Dithering off. Enabling with userspace provided setting 0x%x. Cross your fingers!\n", getName(), ditherOn);
                    WriteRegister(reg, ditherOn);
                }
            }
            else {
                // Dithering currently on.

                // Specific value for control reg specified?
                if (ditherOn > 1) {
                    // Yes. Use it "as is":
                    IOLog("%s: SetDitherMode: Setting dither mode to userspace provided setting 0x%x. Cross your fingers!\n", getName(), ditherOn);
                    WriteRegister(reg, ditherOn);
                }
                else IOLog("%s: SetDitherMode: Dithering already enabled with current control value 0x%x. Skipped.\n", getName(), ReadRegister(reg));
            }
        }
        else {
            // Disable all dithering if it is enabled: Clearing the register to all zero bits does this.
            IOLog("%s: SetDitherMode: Current dither setting before our dither disable on head %d is 0x%x.\n", getName(), headId, ReadRegister(reg));
            if (ReadRegister(reg) > 0) {
                oldDither[headId] = ReadRegister(reg);
                IOLog("%s: SetDitherMode: Current dither setting before our dither disable on head %d is 0x%x. Disabling.\n", getName(), headId, oldDither[headId]);
                WriteRegister(reg, 0x0);
            }
            else {
                IOLog("%s: SetDitherMode: Dithering already disabled. Skipped.\n", getName());
            }
        }

        // End of Radeon et al. support code.
    }
    else {
        // Other unsupported GPU:
        IOLog("%s: SetDitherMode: Tried to call me on a non ATI/AMD GPU. Unsupported.\n", getName());
    }

    return;
}

// Return current vertical rasterbeam position of display head 'headId' (0=Primary CRTC1, 1=Secondary CRTC2):
UInt32 PsychtoolboxKernelDriver::GetBeamPosition(UInt32 headId)
{
    SInt32 beampos = 0;

    if (headId > fNumDisplayHeads - 1) {
        // Invalid head - bail:
        IOLog("%s: GetBeamPosition: ERROR! Invalid headId %d provided. Must be between 0 and %d. Aborted.\n", getName(), headId, fNumDisplayHeads - 1);
        return(0);
    }

    // Query code for ATI/AMD Radeon/FireGL/FirePro:
    if (fDeviceType == kPsychRadeon) {
        if (isDCE4() || isDCE5() || isDCE10()) {
            // DCE-4 display engine (CEDAR and later afaik): Up to six crtc's.

            // Read raw beampostion from GPU:
            beampos = (SInt32) (ReadRegister(EVERGREEN_CRTC_STATUS_POSITION + crtcoff[headId]) & RADEON_VBEAMPOSITION_BITMASK);

            // Query end-offset of VBLANK interval of this GPU and correct for it:
            beampos = beampos - (SInt32) ((ReadRegister(EVERGREEN_CRTC_V_BLANK_START_END + crtcoff[headId]) >> 16) & RADEON_VBEAMPOSITION_BITMASK);

            // Correction for in-VBLANK range:
            if (beampos < 0) beampos = ((SInt32) ReadRegister(EVERGREEN_CRTC_V_TOTAL + crtcoff[headId])) + beampos;

        } else {
            // AVIVO / DCE-2 / DCE-3 display engine (R300 - R600 afaik): At most two display heads for dual-head gpu's.

            // Read raw beampostion from GPU:
            beampos = (SInt32) (ReadRegister((headId == 0) ? RADEON_D1CRTC_STATUS_POSITION : RADEON_D2CRTC_STATUS_POSITION) & RADEON_VBEAMPOSITION_BITMASK);

            // Query end-offset of VBLANK interval of this GPU and correct for it:
            beampos = beampos - (SInt32) ((ReadRegister((headId == 0) ? AVIVO_D1CRTC_V_BLANK_START_END : AVIVO_D2CRTC_V_BLANK_START_END) >> 16) & RADEON_VBEAMPOSITION_BITMASK);

            // Correction for in-VBLANK range:
            if (beampos < 0) beampos = ((SInt32) ReadRegister((headId == 0) ? AVIVO_D1CRTC_V_TOTAL : AVIVO_D2CRTC_V_TOTAL)) + beampos;
        }
    }

    // Query code for NVidia GeForce/Quadro:
    if (fDeviceType == kPsychGeForce) {
        // Pre NV-50 GPU? [Anything before GeForce-8 series]
        if ((fCardType > 0x00) && (fCardType < 0x50)) {
            // Pre NV-50, e.g., RivaTNT-1/2 and all GeForce 256/2/3/4/5/FX/6/7:

            // Lower 12 bits are vertical scanout position (scanline), bit 16 is "in vblank" indicator.
            // Offset between crtc's is 0x2000, we're only interested in scanline, not "in vblank" status:
            // beampos = (SInt32) (ReadRegister((headId == 0) ? 0x600808 : 0x600808 + 0x2000) & 0xFFF);

            // NV-47: Lower 16 bits are vertical scanout position (scanline), upper 16 bits are horizontal
            // scanout position. Offset between crtc's is 0x2000. We only use the lower 16 bits and
            // ignore horizontal scanout position for now:
            beampos = (SInt32) (ReadRegister((headId == 0) ? 0x600868 : 0x600868 + 0x2000) & 0xFFFF);
        } else {
            // NV-50 (GeForce-8) and later:

            // Lower 16 bits are vertical scanout position (scanline), upper 16 bits are vblank counter.
            // Offset between crtc's is 0x800, we're only interested in scanline, not vblank counter:
            beampos = (SInt32) (ReadRegister(0x616340 + (0x800 * headId)) & 0xFFFF);
        }
    }

    // Query code for Intel IGP's:
    if (fDeviceType == kPsychIntelIGP) {
        beampos = (SInt32) (ReadRegister(0x70000 + (headId * 0x1000)) & 0x1FFF);
    }

    // Safety measure: Cap to zero if something went wrong -> This will trigger proper high level error handling in PTB:
    if (beampos < 0) beampos = 0;

    return((UInt32) beampos);
}


// Instantaneously resynchronize display heads of a Radeon dual-head gfx-card:
SInt32 PsychtoolboxKernelDriver::FastSynchronizeAllDisplayHeads(void)
{
    SInt32      deltabeampos;
    UInt32      beampos0, beampos1;
    UInt32      old_crtc_master_enable = 0;
    UInt32      new_crtc_master_enable = 0;
    UInt32      i;

    // AMD GPU's only:
    if (fDeviceType != kPsychRadeon) {
        IOLog("%s: FastSynchronizeAllDisplayHeads(): This function is not supported on non-ATI/AMD GPU's! Aborted.\n", getName());
        return(0xffff);
    }

    IOLog("%s: FastSynchronizeAllDisplayHeads(): About to resynchronize all display heads by use of a 1 second CRTC stop->start cycle:\n", getName());

    // DCE-4 needs a different strategy for now:
    if (isDCE4() || isDCE5() || isDCE10()) {
        // Shut down heads one by one:

        // Detect enabled heads:
        old_crtc_master_enable = 0;
        for (i = 0; i < 6; i++) {
            // Bit 16 "CRTC_CURRENT_MASTER_EN_STATE" allows read-only polling
            // of current activation state of crtc:
            if (ReadRegister(EVERGREEN_CRTC_CONTROL + crtcoff[i]) & (0x1 << 16)) old_crtc_master_enable |= (0x1 << i);
        }

        // Shut down heads, one after each other, wait for each one to settle at its defined resting position:
        for (i = 0; i < 6; i++) {
            IOLog("%s: Head %d ...  ", getName(), i);
            if (old_crtc_master_enable & (0x1 << i)) {
                IOLog("active -> Shutdown. ");

                // Shut down this CRTC by clearing its master enable bit (bit 0):
                WriteRegister(EVERGREEN_CRTC_CONTROL + crtcoff[i], ReadRegister(EVERGREEN_CRTC_CONTROL + crtcoff[i]) & ~(0x1 << 0));

                // Wait 50 msecs, so CRTC has enough time to settle and disable at its
                // programmed resting position:
                IOSleep(50);

                // Double check - Poll until crtc is offline:
                while(ReadRegister(EVERGREEN_CRTC_CONTROL + crtcoff[i]) & (0x1 << 16));
                IOLog("-> Offline.\n");
            }
            else {
                IOLog("already offline.\n");
            }
        }

        // Sleep for 1 second: This is a blocking call, ie. the driver goes to sleep and may wakeup a bit later:
        IOSleep(1000);

        // Reenable all now disabled, but previously enabled display heads.
        // This must be a tight loop, as every microsecond counts for a good sync...
        for (i = 0; i < 6; i++) {
            if (old_crtc_master_enable & (0x1 << i)) {
                // Restart this CRTC by setting its master enable bit (bit 0):
                WriteRegister(EVERGREEN_CRTC_CONTROL + crtcoff[i], ReadRegister(EVERGREEN_CRTC_CONTROL + crtcoff[i]) | (0x1 << 0));
            }
        }

        // We don't have meaningful residual info in the multihead case. Just assume we succeeded:
        deltabeampos = 0;

        IOLog("\n%s: Display head resync operation finished.\n\n", getName());

        // Done.
        return(deltabeampos);
    }

    // AVIVO case:

    // A little pretest...
    IOLog("%s: Pretest...\n", getName());
    for (UInt32 i = 0; i<10; i++) {
        beampos0 = GetBeamPosition(0);
        beampos1 = GetBeamPosition(1);
        IOLog("%s: Sample %d: Beampositions are %d vs. %d . Offset %d\n", getName(), i, beampos0, beampos1, (SInt32) beampos1 - (SInt32) beampos0);
    }

    // Query the CRTC scan-converter master enable state: Bit 0 (value 0x1) controls Pipeline 1,
    // whereas Bit 1(value 0x2) controls Pipeline 2:
    old_crtc_master_enable = ReadRegister(RADEON_DC_CRTC_MASTER_ENABLE);
    IOLog("%s: Current CRTC Master enable state is %d . Trying to stop and reset all display heads.\n", getName(), old_crtc_master_enable);
    IOLog("%s: Will wait individually for each head to get close to scanline 0, then disable it.\n", getName());

    // Shut down heads, one after each other, each one at the start of a new refresh cycle:
    for (UInt32 i = 0; i <= 1; i++) {
        // Wait for head i to start a new display cycle (scanline 0), then shut it down - well if it is active at all:
        IOLog("%s: Head %d ...  ", getName(), i);
        if (old_crtc_master_enable & (0x1 << i)) {
            IOLog("active -> Shutdown. ");
            // Wait for beam going above scanline 240: We choose 240, because even at the lowest conceivable
            // useful display resolution of 640 x 480, 240 will be in the middle of the frame, aka far away
            // from both, the start and the end of a frame:
            while (GetBeamPosition(i) <= 240);

            // Beam is heading for the end of the frame + VBL area. Wait for wrap-around, ie. until
            // reaching a scanline value smaller than 100 --> Until it wraps back to zero or at least
            // a value close to zero:
            while (GetBeamPosition(i) > 240);

            // Start of new refresh interval! Shut down this heads CRTC!
            // We do so by clearing enable bit for this head:
            WriteRegister(RADEON_DC_CRTC_MASTER_ENABLE, ReadRegister(RADEON_DC_CRTC_MASTER_ENABLE) & ~(0x1 << i));
            IOLog("New state is %d.\n", ReadRegister(RADEON_DC_CRTC_MASTER_ENABLE));

            // Head should be down, close to scanline 0.
            IOSleep(50);
        }
        else {
            IOLog("already offline.\n");
        }
    }

    // Disable all display heads:
    // WriteRegister(RADEON_DC_CRTC_MASTER_ENABLE, 0);
    IOSleep(20);

    // Query current beamposition and check state:
    beampos0 = GetBeamPosition(0);
    beampos1 = GetBeamPosition(1);

    new_crtc_master_enable = ReadRegister(RADEON_DC_CRTC_MASTER_ENABLE);

    if (new_crtc_master_enable == 0) {
        IOLog("%s: CRTC's down (state %d): Beampositions are [0]=%d and [1]=%d. Synchronized restart in 1 second...\n", getName(), new_crtc_master_enable, beampos0, beampos1);
    }
    else {
        IOLog("%s: CRTC's shutdown failed!! (state %d): Beamposition are [0]=%d and [1]=%d. Will try to restart in 1 second...\n", getName(), new_crtc_master_enable, beampos0, beampos1);
    }

    // Sleep for 1 secs == 1000 milliseconds: This is a blocking call, ie. the thread goes to sleep and may wakeup a bit later:
    IOSleep(1000);

    // Reset all display heads enable state to original setting:
    WriteRegister(RADEON_DC_CRTC_MASTER_ENABLE, old_crtc_master_enable);

    // Query position and state after restart:
    beampos0 = GetBeamPosition(0);
    beampos1 = GetBeamPosition(1);
    new_crtc_master_enable = ReadRegister(RADEON_DC_CRTC_MASTER_ENABLE);
    if (new_crtc_master_enable == old_crtc_master_enable) {
        IOLog("%s: CRTC's restarted in sync: Master enable state is %d. Beampositions after restart: [0]=%d and [1]=%d.\n", getName(), new_crtc_master_enable, beampos0, beampos1);
    }
    else {
        IOLog("%s: CRTC's restart FAILED!!: Master enable state is %d. Beampositions: [0]=%d and [1]=%d.\n", getName(), new_crtc_master_enable, beampos0, beampos1);
    }

    deltabeampos = (SInt32) beampos1 - (SInt32) beampos0;
    IOLog("%s: Residual beam offset after display sync: %d.\n\n", getName(), deltabeampos);

    // A little posttest...
    IOLog("%s: Posttest...\n", getName());
    for (UInt32 i = 0; i<10; i++) {
        beampos0 = GetBeamPosition(0);
        beampos1 = GetBeamPosition(1);
        IOLog("%s: Sample %d: Beampositions are %d vs. %d . Offset %d\n", getName(), i, beampos0, beampos1, (SInt32) beampos1 - (SInt32) beampos0);
    }

    IOLog("\n%s: Display head resync operation finished.\n\n", getName());

    // Return offset after (re-)sync:
    return(deltabeampos);
}

// Returns multiple flags with info like PCI Vendor/device id, display engine type etc.
void PsychtoolboxKernelDriver::GetGPUInfo(UInt32 *inOutArgs)
{
    IOLog("%s: GetGPUInfo(): Returning GPU info.\n", getName());

    // First return arg is detected GPU vendor type:
    inOutArgs[0] = fDeviceType;

    // 2nd = PCI device id:
    inOutArgs[1] = (UInt32) fPCIDeviceId;

    // 3rd = Type of display engine:

    // Default to "don't know".
    inOutArgs[2] = 0;

    // On Radeons we distinguish between Avivo / DCE-2 (10), DCE-3 (30), or DCE-4 style (40) or DCE-5 (50) or DCE-6 (60) or DCE-8 (80), DCE-10 (100) or DCE-11 (110) for now.
    if (fDeviceType == kPsychRadeon) inOutArgs[2] = isDCE11() ? 110 : isDCE10() ? 100 : isDCE8() ? 80 : (isDCE6() ? 60 : (isDCE5() ? 50 : (isDCE4() ? 40 : (isDCE3() ? 30 : 10))));

    // On NVidia's we distinguish between chip family, e.g., 0x40 for the NV-40 family.
    if (fDeviceType == kPsychGeForce) inOutArgs[2] = fCardType;

    // 4th = Maximum number of crtc's:
    inOutArgs[3] = fNumDisplayHeads;

    return;
}

// Query if LUT for given headId is all-zero: 0 = Something else, 1 = Zero-LUT, 2 = It's an identity LUT,
// 3 = Not-quite-identity mapping, 0xffffffff = don't know.
UInt32 PsychtoolboxKernelDriver::GetLUTState(UInt32 headId, UInt32 debug)
{
    UInt32 i, v, r, m, bo, wo, offset, reg;
    UInt32 isZero = 1;
    UInt32 isIdentity = 1;

    // AMD GPU's:
    if (fDeviceType == kPsychRadeon) {
        IOLog("%s: GetLUTState(): Checking LUT and bias values on GPU for headId %d.\n", getName(), headId);

        if (isDCE4() || isDCE5() || isDCE10()) {
            // DCE-4.0 and later: Up to (so far) six display heads:
            if (headId > fNumDisplayHeads - 1) {
                // Invalid head - bail:
                IOLog("%s: GetLUTState: ERROR! Invalid headId %d provided. Must be between 0 and %d. Aborted.\n", getName(), headId, fNumDisplayHeads - 1);
                return(0xffffffff);
            }

            offset = crtcoff[headId];
            WriteRegister(EVERGREEN_DC_LUT_RW_MODE + offset, 0);
            WriteRegister(EVERGREEN_DC_LUT_RW_INDEX + offset, 0);
            reg = EVERGREEN_DC_LUT_30_COLOR + offset;

            // Find out if there are non-zero black offsets:
            bo = 0x0;
            bo|= ReadRegister(EVERGREEN_DC_LUT_BLACK_OFFSET_BLUE + offset);
            bo|= ReadRegister(EVERGREEN_DC_LUT_BLACK_OFFSET_GREEN + offset);
            bo|= ReadRegister(EVERGREEN_DC_LUT_BLACK_OFFSET_RED + offset);

            // Find out if there are non-0xffff white offsets:
            wo = 0x0;
            wo|= 0xffff - ReadRegister(EVERGREEN_DC_LUT_WHITE_OFFSET_BLUE + offset);
            wo|= 0xffff - ReadRegister(EVERGREEN_DC_LUT_WHITE_OFFSET_GREEN + offset);
            wo|= 0xffff - ReadRegister(EVERGREEN_DC_LUT_WHITE_OFFSET_RED + offset);
        }
        else {
            // AVIVO: Dualhead.
            offset = (headId > 0) ? 0x800 : 0x0;
            WriteRegister(AVIVO_DC_LUT_RW_SELECT, headId & 0x1);
            WriteRegister(AVIVO_DC_LUT_RW_MODE, 0);
            WriteRegister(AVIVO_DC_LUT_RW_INDEX, 0);
            reg = AVIVO_DC_LUT_30_COLOR;

            // Find out if there are non-zero black offsets:
            bo = 0x0;
            bo|= ReadRegister(AVIVO_DC_LUTA_BLACK_OFFSET_BLUE + offset);
            bo|= ReadRegister(AVIVO_DC_LUTA_BLACK_OFFSET_GREEN + offset);
            bo|= ReadRegister(AVIVO_DC_LUTA_BLACK_OFFSET_RED + offset);

            // Find out if there are non-0xffff white offsets:
            wo = 0x0;
            wo|= 0xffff - ReadRegister(AVIVO_DC_LUTA_WHITE_OFFSET_BLUE + offset);
            wo|= 0xffff - ReadRegister(AVIVO_DC_LUTA_WHITE_OFFSET_GREEN + offset);
            wo|= 0xffff - ReadRegister(AVIVO_DC_LUTA_WHITE_OFFSET_RED + offset);
        }

        if (debug) IOLog("%s: Offsets: Black %d : White %d.\n", getName(), bo, wo);

        for (i = 0; i < 256; i++) {
            // Read 32 bit value of this slot, mask out upper 2 bits,
            // so the least significant 30 bits are left, as these
            // contain the 3 * 10 bits for the 10 bit R,G,B channels:
            v = ReadRegister(reg) & (0xffffffff >> 2);

            // All zero as they should be for a all-zero LUT?
            if (v > 0) isZero = 0;

            // Compare with expected value in slot i for a perfect 10 bit identity LUT
            // intended for a 8 bit output encoder, i.e., 2 least significant bits
            // zero to avoid dithering and similar stuff:
            r = i << 2;
            m = (r << 20) | (r << 10) | (r << 0);

            // Mismatch? Not a perfect identity LUT:
            if (v != m) isIdentity = 0;

            if (debug) {
                IOLog("%d:%d,%d,%d\n", i, (v >> 20) & 0x3ff, (v >> 10) & 0x3ff, (v >> 0) & 0x3ff);
                IOSleep(1);
            }
        }

        if (isZero) return(1);  // All zero LUT.

        if (isIdentity) {
            // If wo or bo is non-zero then it is not quite an identity
            // mapping, as the black and white offset are not neutral.
            // Return 3 in this case:
            if ((wo | bo) > 0) return(3);

            // Perfect identity LUT:
            return(2);
        }

        // Regular LUT:
        return(0);
    }

    // Unhandled:
    IOLog("%s: GetLUTState(): This function is not supported on this GPU. Returning 0xffffffff.\n", getName());
    return(0xffffffff);
}

// Load an identity LUT into display head 'headid': Return 1 on success, 0 on failure or if unsupported for this GPU:
UInt32 PsychtoolboxKernelDriver::LoadIdentityLUT(UInt32 headId)
{
    UInt32 i, r, m, offset, reg;

    // AMD GPU's:
    if (fDeviceType == kPsychRadeon) {
        IOLog("%s: LoadIdentityLUT(): Uploading identity LUT and bias values into GPU for headId %d.\n", getName(), headId);

        if (isDCE4() || isDCE5() || isDCE10()) {
            // DCE-4.0 and later: Up to (so far) six display heads:
            if (headId > fNumDisplayHeads - 1) {
                // Invalid head - bail:
                IOLog("%s: LoadIdentityLUT: ERROR! Invalid headId %d provided. Must be between 0 and %d. Aborted.\n", getName(), headId, fNumDisplayHeads - 1);
                return(0);
            }

            offset = crtcoff[headId];
            reg = EVERGREEN_DC_LUT_30_COLOR + offset;

            WriteRegister(EVERGREEN_DC_LUT_CONTROL + offset, 0);

            if (isDCE5() || isDCE10()) {
                WriteRegister(NI_INPUT_CSC_CONTROL + offset,
                              (NI_INPUT_CSC_GRPH_MODE(NI_INPUT_CSC_BYPASS) |
                               NI_INPUT_CSC_OVL_MODE(NI_INPUT_CSC_BYPASS)));
                WriteRegister(NI_PRESCALE_GRPH_CONTROL + offset,
                              NI_GRPH_PRESCALE_BYPASS);
                WriteRegister(NI_PRESCALE_OVL_CONTROL + offset,
                              NI_OVL_PRESCALE_BYPASS);
                WriteRegister(NI_INPUT_GAMMA_CONTROL + offset,
                              (NI_GRPH_INPUT_GAMMA_MODE(NI_INPUT_GAMMA_USE_LUT) |
                               NI_OVL_INPUT_GAMMA_MODE(NI_INPUT_GAMMA_USE_LUT)));
            }

            // Set zero black offsets:
            WriteRegister(EVERGREEN_DC_LUT_BLACK_OFFSET_BLUE  + offset, 0x0);
            WriteRegister(EVERGREEN_DC_LUT_BLACK_OFFSET_GREEN + offset, 0x0);
            WriteRegister(EVERGREEN_DC_LUT_BLACK_OFFSET_RED   + offset, 0x0);

            // Set 0xffff white offsets:
            WriteRegister(EVERGREEN_DC_LUT_WHITE_OFFSET_BLUE  + offset, 0xffff);
            WriteRegister(EVERGREEN_DC_LUT_WHITE_OFFSET_GREEN + offset, 0xffff);
            WriteRegister(EVERGREEN_DC_LUT_WHITE_OFFSET_RED   + offset, 0xffff);

            WriteRegister(EVERGREEN_DC_LUT_RW_MODE + offset, 0);
            WriteRegister(EVERGREEN_DC_LUT_WRITE_EN_MASK + offset, 0x00000007);

            WriteRegister(EVERGREEN_DC_LUT_RW_INDEX + offset, 0);
        }
        else {
            // AVIVO: Dualhead.
            offset = (headId > 0) ? 0x800 : 0x0;
            reg = AVIVO_DC_LUT_30_COLOR;

            WriteRegister(AVIVO_DC_LUTA_CONTROL + offset, 0);

            // Set zero black offsets:
            WriteRegister(AVIVO_DC_LUTA_BLACK_OFFSET_BLUE  + offset, 0x0);
            WriteRegister(AVIVO_DC_LUTA_BLACK_OFFSET_GREEN + offset, 0x0);
            WriteRegister(AVIVO_DC_LUTA_BLACK_OFFSET_RED   + offset, 0x0);

            // Set 0xffff white offsets:
            WriteRegister(AVIVO_DC_LUTA_WHITE_OFFSET_BLUE  + offset, 0xffff);
            WriteRegister(AVIVO_DC_LUTA_WHITE_OFFSET_GREEN + offset, 0xffff);
            WriteRegister(AVIVO_DC_LUTA_WHITE_OFFSET_RED   + offset, 0xffff);

            WriteRegister(AVIVO_DC_LUT_RW_SELECT, headId & 0x1);
            WriteRegister(AVIVO_DC_LUT_RW_MODE, 0);
            WriteRegister(AVIVO_DC_LUT_WRITE_EN_MASK, 0x0000003f);

            WriteRegister(AVIVO_DC_LUT_RW_INDEX, 0);
        }

        for (i = 0; i < 256; i++) {
            // Compute perfect value for slot i for a perfect 10 bit identity LUT
            // intended for a 8 bit output encoder, i.e., 2 least significant bits
            // zero to avoid dithering and similar stuff, the 8 most significant
            // bits for each 10 bit color channel linearly increasing one unit
            // per slot:
            r = i << 2;
            m = (r << 20) | (r << 10) | (r << 0);

            // Write 32 bit value of this slot:
            WriteRegister(reg, m);
        }

        if (isDCE5() || isDCE10()) {
            WriteRegister(NI_DEGAMMA_CONTROL + offset,
                          (NI_GRPH_DEGAMMA_MODE(NI_DEGAMMA_BYPASS) |
                           NI_OVL_DEGAMMA_MODE(NI_DEGAMMA_BYPASS) |
                           NI_ICON_DEGAMMA_MODE(NI_DEGAMMA_BYPASS) |
                           NI_CURSOR_DEGAMMA_MODE(NI_DEGAMMA_BYPASS)));
            WriteRegister(NI_GAMUT_REMAP_CONTROL + offset,
                          (NI_GRPH_GAMUT_REMAP_MODE(NI_GAMUT_REMAP_BYPASS) |
                           NI_OVL_GAMUT_REMAP_MODE(NI_GAMUT_REMAP_BYPASS)));
            WriteRegister(NI_REGAMMA_CONTROL + offset,
                          (NI_GRPH_REGAMMA_MODE(NI_REGAMMA_BYPASS) |
                           NI_OVL_REGAMMA_MODE(NI_REGAMMA_BYPASS)));
            WriteRegister(NI_OUTPUT_CSC_CONTROL + offset,
                          (NI_OUTPUT_CSC_GRPH_MODE(NI_OUTPUT_CSC_BYPASS) |
                           NI_OUTPUT_CSC_OVL_MODE(NI_OUTPUT_CSC_BYPASS)));
            /* XXX match this to the depth of the crtc fmt block, move to modeset? */
            WriteRegister(0x6940 + offset, 0);
        }

        // Done.
        return(1);
    }

    // Unhandled:
    IOLog("%s: LoadIdentityLUT(): This function is not supported on this GPU. Returning 0.\n", getName());
    return(0);
}

// Perform instant state snapshot of interesting registers and return'em:
void PsychtoolboxKernelDriver::GetStateSnapshot(PsychKDCommandStruct* outStruct)
{
    // TODO...
    return;
}

// Dump interesting register state to system log:
void PsychtoolboxKernelDriver::DumpGfxState(void)
{
    UInt32 regidx;
    UInt32 col = 0;

    // Output some range of registers to IOLog:
    for (regidx = 0; regidx <= 0x200; regidx+=4) {
        IOLog("Reg[%x] = %x :: ", regidx, ReadRegister(regidx));
        if ((col++) % 5 == 0) IOLog("\n");
    }

    return;
}

/*
IOReturn PsychtoolboxKernelDriver::ScalarIStructI(int inNumber, MySampleStruct* inStruct, IOByteCount inStructSize)
{
    IOLog("SimpleDriver::ScalarIStructI(inNumber = %d, int16 = %d, int32 = %d, inStructSize = %d)\n",
          inNumber, inStruct->int16, (int)inStruct->int32, (int)inStructSize);

    return kIOReturnSuccess;
}


IOReturn PsychtoolboxKernelDriver::ScalarIStructO(int inNumber1, int inNumber2, MySampleStruct* outStruct, IOByteCount* outStructSize)
{
    IOLog("SimpleDriver::ScalarIStructO(inNumber1 = %d, inNumber2 = %d)\n", inNumber1, inNumber2);

    // The scalars inNumber1 and inNumber2 are automatically endian-swapped by the user client mechanism.
    // The user client will swap the fields in outStruct if necessary.

    // The output is simply the two input scalars copied to the output struct.
    outStruct->int16 = inNumber1;
    outStruct->int32 = inNumber2;
    *outStructSize = sizeof(MySampleStruct);

    return kIOReturnSuccess;
}


IOReturn PsychtoolboxKernelDriver::ScalarIScalarO(int inNumber1, int inNumber2, int* outNumber)
{
    IOLog("SimpleDriver::ScalarIScalarO(inNumber1 = %d, inNumber2 = %d)\n", inNumber1, inNumber2);

    // The output is the product of the two input scalars.
    *outNumber = inNumber1 * inNumber2;

    return kIOReturnSuccess;
}
*/

IOReturn PsychtoolboxKernelDriver::PsychKDDispatchCommand(PsychKDCommandStruct* inStruct, PsychKDCommandStruct* outStruct, IOByteCount inStructSize, IOByteCount* outStructSize)
{
    // Default return value is success:
    IOReturn ret = kIOReturnSuccess;

    // Fetch command code:
    UInt32 cmd = inStruct->command;

    // Validate it:
    if (cmd >= kPsychKDMaxCommands) {
        IOLog("%s: Invalid/Unknown command code received - Out of range.\n", getName());
        return(kIOReturnBadArgument);
    }

    // Dispatch depending on command code:
    switch(cmd) {
        case kPsychKDFastSyncAllHeads:    // Fast synchronize all display heads:
            // One return argument, the residual offset after sync:
            outStruct->inOutArgs[0] = FastSynchronizeAllDisplayHeads();
        break;

        case kPsychKDDumpMiscInfo:        // Dump miscellaneous info about gfx-state to System log:
            // No arguments, just do it:
            DumpGfxState();
        break;

        case kPsychKDGetBeamposition:    // Query current beamposition for a specific display head:
            // One input, the headid. One output, the beamposition:
            outStruct->inOutArgs[0] = GetBeamPosition(inStruct->inOutArgs[0]);
        break;

        case kPsychKDGetStateSnapshot:    // Perform snapshot of many relevant registers and return it timestamped:
            GetStateSnapshot(outStruct);
        break;

        case kPsychKDReadRegister:        // Read specific register and return its value:
            outStruct->inOutArgs[0] = ReadRegister(inStruct->inOutArgs[0]);
        break;

        case KPsychKDWriteRegister:        // Write value to specific register:
            WriteRegister(inStruct->inOutArgs[0], inStruct->inOutArgs[1]);
        break;

        case kPsychKDSetDitherMode:        // Write value to specific register:
            SetDitherMode(inStruct->inOutArgs[0], inStruct->inOutArgs[1]);
        break;

        case kPsychKDGetRevision:    // Query current driver revision:
            // One output, the revision:
            outStruct->inOutArgs[0] = PTBKDRevision;
        break;

        case kPsychKDGetGPUInfo:    // Query info about installed GPU:
            // Returns multiple flags with info like PCI Vendor/device id, display engine type etc.
            GetGPUInfo(outStruct->inOutArgs);
        break;

        case kPsychKDGetLUTState:    // Query if LUT for given headId is all-zero:
            // One output, the binary flag (1 = All-zero LUT, 0 = Non-Zero LUT)
            outStruct->inOutArgs[0] = GetLUTState(inStruct->inOutArgs[0], inStruct->inOutArgs[1]);
        break;

        case kPsychKDSetIdentityLUT:    // Load an identity LUT into display head 'headid':
            // One return argument, the binary flag (1 = Success, 0 = Not supported)
            outStruct->inOutArgs[0] = LoadIdentityLUT(inStruct->inOutArgs[0]);
        break;

        default:
            ret = kIOReturnBadArgument;
            IOLog("%s: Invalid/Unknown command code received - No match in dispatch table.\n", getName());
    }

    *outStructSize = sizeof(PsychKDCommandStruct);

    return kIOReturnSuccess;
}
