/*
    PsychToolbox3/Source/Linux/Screen/PsychScreenGlue.c

    PLATFORMS:

        This is the Linux version.
        It contains all code shared across Linux backends and
        the code for the non-Wayland backend - specifically for
        the X11 backend and the GBM+DRM/KMS backend, although
        the GBM backend is just a shim which doesn't do anything
        but prevent the beast from crashing.

    AUTHORS:

        Mario Kleiner       mk      mario.kleiner.de@gmail.com

    HISTORY:

        2/20/06             mk      Wrote it. Derived from Windows version.

    DESCRIPTION:

        Functions in this file comprise an abstraction layer for probing and controlling screen state.

        Each C function which implements a particular Screen subcommand should be platform neutral.  For example, the source to SCREENPixelSizes()
        should be platform-neutral, despite that the calls in OS X and Linux to detect available pixel sizes are
        different.  The platform specificity is abstracted out in C files which end it "Glue", for example PsychScreenGlue, PsychWindowGlue,
        PsychWindowTextClue.

        In addition to glue functions for windows and screen there are functions which implement shared functionality between between Screen commands,
        such as ScreenTypes.c and WindowBank.c.

*/

#include "Screen.h"

/* These are needed for our GPU specific beamposition query implementation: */
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

// To use libpciaccess for GPU device detection and mmaping():
#include "pciaccess.h"
#define PCI_CLASS_DISPLAY    0x03

// Maximum number of slots in a gamma table to set/query: This should be plenty.
#define MAX_GAMMALUT_SIZE 16384

// Defined in PsychGraphicsHardwareHALSupport.c, but accessed and initialized here:
extern unsigned int crtcoff[kPsychMaxPossibleCrtcs];

// If non-zero entries, then we are or have been running on Mesa
// and mesaversion encodes major.minor.patchlevel:
extern int mesaversion[3];

// Event and error base for XRandR extension:
int xr_event, xr_error;
psych_bool has_xrandr_1_2 = FALSE;
psych_bool has_xrandr_1_3 = FALSE;

/* Following structures are needed by our ATI/AMD/NVIDIA beamposition query implementation: */
/* Location and format of the relevant hardware registers of the ATI R500/R600 chips
 * was taken from the official register spec for that chips which was released to
 * the public by AMD/ATI end of 2007 and is available for download at XOrg.
 *
 * http://www.x.org/docs/AMD/
 *
 * Register spec's for DCE-4 hardware are from Linux kms driver and Alex Deucher.
 * This should work on any AVIVO or DCE4/5 display hardware chip, i.e., R300 and
 * later. It won't work on ancient pre-AVIVO hardware.
 */

#include "PsychGraphicsCardRegisterSpecs.h"
#include <endian.h>

// gfx_cntl_mem is mapped to the actual device's memory mapped control area.
// Not the address but what it points to is volatile.
struct pci_device *gpu = NULL;
unsigned char * volatile gfx_cntl_mem = NULL;
unsigned long gfx_length = 0;
unsigned long gfx_lowlimit = 0;
unsigned int  fDeviceType = 0;
unsigned int  fCardType = 0;
unsigned int  fPCIDeviceId = 0;
unsigned int  fNumDisplayHeads = 0;

// Minimum allowed physical crtc id for assignment to X-Screens. Used
// for the (x-screen, output) -> physical crtc id mapping heuristic
// for multi-x-screen ZaphodHead display setups:
static int    minimum_crtcid = 0;

// Count of kernel drivers:
static int    numKernelDrivers = 0;

// Internal helper function prototype:
void PsychInitNonX11(void);

/* Mappings up to date for May 2016 (last update e-mail patch / commit 2016-05-18). Would need updates for anything after start of June 2016 */

/* Is a given ATI/AMD GPU a DCE11 type ASIC, i.e., with the new display engine? */
static psych_bool isDCE11(int screenId)
{
    psych_bool isDCE11 = false;

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
static psych_bool isDCE10(int screenId)
{
    psych_bool isDCE10 = false;

    // TONGA and FIJI are DCE10 -- This is part of the "Volcanic Islands" GPU family.

    // TONGA: 0x692x - 0x693x so far.
    if ((fPCIDeviceId & 0xFFF0) == 0x6920) isDCE10 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x6930) isDCE10 = true;

    // FIJI in 0x7300 range:
    if ((fPCIDeviceId & 0xFF00) == 0x7300) isDCE10 = true;

    // All DCE11 are also DCE10, so far...
    if (isDCE11(screenId)) isDCE10 = true;

    return(isDCE10);
}

/* Is a given ATI/AMD GPU a DCE8 type ASIC, i.e., with the new display engine? */
static psych_bool isDCE8(int screenId)
{
    psych_bool isDCE8 = false;

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

    // CAUTION: DCE 10 and higher are *not* DCE8 as well!
    // These new parts have a different register layout, so
    // need separate code!

    return(isDCE8);
}

/* Is a given ATI/AMD GPU a DCE6.4 type ASIC, i.e., with the new display engine? */
static psych_bool isDCE64(int screenId)
{
    psych_bool isDCE64 = false;

    // Everything == OLAND is DCE6.4 -- This is part of the "Southern Islands" GPU family.

    // OLAND in 0x66xx range:
    if ((fPCIDeviceId & 0xFF00) == 0x6600) isDCE64 = true;

    return(isDCE64);
}

/* Is a given ATI/AMD GPU a DCE6.1 type ASIC, i.e., with the new display engine? */
static psych_bool isDCE61(int screenId)
{
    psych_bool isDCE61 = false;

    // Everything >= ARUBA which is an IGP is DCE6.1 -- This is the "Trinity" GPU family.

    // ARUBA in 0x99xx range: This is the "Trinity" chip family.
    if ((fPCIDeviceId & 0xFF00) == 0x9900) isDCE61 = true;

    // KAVERI in 0x13xx range:
    if ((fPCIDeviceId & 0xFF00) == 0x1300) isDCE61 = true;

    return(isDCE61);
}

/* Is a given ATI/AMD GPU a DCE6 type ASIC, i.e., with the new display engine? */
static psych_bool isDCE6(int screenId)
{
    psych_bool isDCE6 = false;

    // Everything >= ARUBA is DCE6 -- This is the "Southern Islands" GPU family.
    // First real DCE-6.0 is TAHITI in 0x678x - 0x679x range:
    if ((fPCIDeviceId & 0xFFF0) == 0x6780) isDCE6 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x6790) isDCE6 = true;

    // Then PITCAIRN, VERDE in 0x6800 - 0x683x range:
    if ((fPCIDeviceId & 0xFFF0) == 0x6800) isDCE6 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x6810) isDCE6 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x6820) isDCE6 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x6830) isDCE6 = true;

    // And one outlier PITCAIRN:
    if ((fPCIDeviceId & 0xFFFF) == 0x684c) isDCE6 = true;

    // Then HAINAN in the 0x666x range:
    if ((fPCIDeviceId & 0xFFF0) == 0x6660) isDCE6 = true;

    // All DCE-6.1 engines are also DCE-6:
    if (isDCE61(screenId)) isDCE6 = true;

    // All DCE-6.4 engines are also DCE-6:
    if (isDCE64(screenId)) isDCE6 = true;

    // All DCE-8 engines are also DCE-6:
    if (isDCE8(screenId)) isDCE6 = true;

    return(isDCE6);
}

/* Is a given ATI/AMD GPU a DCE5 type ASIC, i.e., with the new display engine? */
static psych_bool isDCE5(int screenId)
{
    psych_bool isDCE5 = false;

    // Everything after BARTS is DCE5 -- This is the "Northern Islands" GPU family.
    // Barts, Turks, Caicos, Cayman, Antilles in 0x67xx range:
    if ((fPCIDeviceId & 0xFF00) == 0x6700) isDCE5 = true;

    // More Turks ids:
    if ((fPCIDeviceId & 0xFFF0) == 0x6840) isDCE5 = true;
    if ((fPCIDeviceId & 0xFFF0) == 0x6850) isDCE5 = true;

    // All DCE-6 engines are also DCE-5:
    if (isDCE6(screenId)) isDCE5 = true;

    return(isDCE5);
}

/* Is a given ATI/AMD GPU a DCE-4.1 type ASIC, i.e., with the new display engine? */
static psych_bool isDCE41(int screenId)
{
    psych_bool isDCE41 = false;

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
static psych_bool isDCE4(int screenId)
{
    psych_bool isDCE4 = false;

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
    if (isDCE41(screenId)) isDCE4 = true;

    // All DCE-5 engines are also DCE-4:
    if (isDCE5(screenId)) isDCE4 = true;

    return(isDCE4);
}

static psych_bool isDCE3(int screenId)
{
    psych_bool isDCE3 = false;

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

// Helper routine: Read a single 32 bit unsigned int hardware register at
// offset 'offset' and return its value:
static unsigned int ReadRegister(unsigned long offset)
{
    unsigned int value;

    // Safety check: Don't allow reads past devices MMIO range:
    // We don't return error codes and don't log the problem,
    // because we could be called from primary Interrupt path, so IOLog() is not
    // an option!
    if (gfx_cntl_mem == NULL || offset > gfx_length-4 || offset < gfx_lowlimit) {
        printf("PTB-ERROR: In GPU ReadRegister(): MMIO not mapped or reg offset %p out of range [%p - %p]! Nop zero return!\n", offset, gfx_lowlimit, gfx_length-4);
        return(0);
    }

    // Read and return value:
    value = *(unsigned int * volatile)(gfx_cntl_mem + offset);

    // Enforce a full memory barrier: This is a gcc intrinsic:
    __sync_synchronize();

    // Radeon: Don't know endianity behaviour: Play save, stick to LE assumption for now:
    if (fDeviceType == kPsychRadeon) return(le32toh(value));

    // Read the register in native byte order: At least NVidia GPU's adapt their
    // endianity to match the host systems endianity, so no need for conversion:
    if (fDeviceType == kPsychGeForce)  return(value);
    if (fDeviceType == kPsychIntelIGP) return(value);

    // No-Op return:
    printf("PTB-ERROR: In GPU ReadRegister(): UNKNOWN fDeviceType of GPU! NO OPERATION!\n");
    return(0);
}

// Helper routine: Write a single 32 bit unsigned int hardware register at
// offset 'offset':
static void WriteRegister(unsigned long offset, unsigned int value)
{
    // Safety check: Don't allow reads past devices MMIO range:
    // We don't return error codes and don't log the problem,
    // because we could be called from primary Interrupt path, so IOLog() is not
    // an option!
    if (gfx_cntl_mem == NULL || offset > gfx_length-4 || offset < gfx_lowlimit) {
        printf("PTB-ERROR: In GPU WriteRegister(): MMIO not mapped or reg offset %p out of range [%p - %p]! Nop zero return!\n", offset, gfx_lowlimit, gfx_length-4);
        return;
    }

    // Write the register in native byte order: At least NVidia GPU's adapt their
    // endianity to match the host systems endianity, so no need for conversion:
    if (fDeviceType == kPsychGeForce)  value = value;
    if (fDeviceType == kPsychIntelIGP) value = value;

    // Radeon: Don't know endianity behaviour: Play save, stick to LE assumption for now:
    if (fDeviceType == kPsychRadeon) value = htole32(value);

    *(unsigned int* volatile)(gfx_cntl_mem + offset) = value;

    // Enforce a full memory barrier: This is a gcc intrinsic:
    __sync_synchronize();
}

void PsychScreenUnmapDeviceMemory(void)
{
    // Any mapped?
    if (gfx_cntl_mem) {
        // Unmap:
        pci_device_unmap_range(gpu, (void*) gfx_cntl_mem, gfx_length);
        gfx_cntl_mem = NULL;
        gfx_length = 0;
        gpu = NULL;
        numKernelDrivers = 0;
    }

    // Shutdown PCI access library, release all resources:
    pci_system_cleanup();

    return;
}

// Helper routine: Check if a supported GPU is installed, and mmap() its MMIO register
// control block into our address space for direct register access:
psych_bool PsychScreenMapRadeonCntlMemory(void)
{
    struct pci_device_iterator *iter;
    struct pci_device *dev;
    struct pci_mem_region *region = NULL;
    int ret = 0;
    int screenId = 0;
    int currentgpuidx = 0, targetgpuidx = -1;

    // A bit of a hack for now: Allow user to select which gpu in a multi-gpu
    // system should be used for low-level mmio based features. If the environment
    // variable PSYCH_USE_GPUIDX is set to a number, it will try to use that GPU:
    // TODO: Replace this by true multi-gpu support and - far in the future? -
    // automatic mapping of screens to gpu's:
    if (getenv("PSYCH_USE_GPUIDX")) {
        targetgpuidx = atoi(getenv("PSYCH_USE_GPUIDX"));
        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Will try to use GPU number %i for low-level access during this session.\n", targetgpuidx);
    }

    // Safe-guard:
    if (gfx_cntl_mem || gpu) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Redundant call to PsychScreenMapRadeonCntlMemory()! Ignored for now. This should not happen!\n");
        return(TRUE);
    }

    // Start with default setting: No low-level access possible.
    gfx_cntl_mem = NULL;
    gfx_length = 0;
    gpu = NULL;

    // On ARM architecture system? If so, we assume it is a SoC without PCI bus, ergo no low-level PCI MMIO mapping/access:
    #if defined(__arm__) || defined(__thumb__) || defined(__aarch64__)
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Not using any low-level access to GPU, as this is an ARM SoC.\n");
        return(FALSE);
    #endif

    // Initialize libpciaccess:
    ret = pci_system_init();
    if (ret) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Could not establish low-level access to GPU for screenId %i - Could not initialize PCI system.\n", screenId);
        return(FALSE);
    }

    // Enumerate them:
    iter = pci_id_match_iterator_create(NULL);
    while ((dev = pci_device_next(iter)) != NULL) {
        if (PsychPrefStateGet_Verbosity() > 4) {
            printf("PTB-DEBUG: Checking PCI device [%s %s] with class x%08x ...\n", pci_device_get_vendor_name(dev), pci_device_get_device_name(dev), dev->device_class);
        }

        // "Upgrade" pre PCI 2.0 class device to PCI 2.0 class equivalent
        // to simplify matching:
        if (dev->device_class == 0x00000101) dev->device_class = 0x00030000;

        // GPU aka display device class?
        if ((dev->device_class & 0x00ff0000) == (PCI_CLASS_DISPLAY << 16)) {
            // dev is our current candidate gpu. Matching vendor?
            if (dev->vendor_id == PCI_VENDOR_ID_NVIDIA || dev->vendor_id == PCI_VENDOR_ID_ATI || dev->vendor_id == PCI_VENDOR_ID_AMD || dev->vendor_id == PCI_VENDOR_ID_INTEL) {
                // Yes. This is our baby from NVidia or ATI/AMD or Intel:

                // From the land of bad hacks: Keep track early if we enumerated an Intel GPU, because
                // some KDE-specific hacks for override_redirect handling in window setup need to know
                // about running on an Intel IGP, even if we don't actually mmap() the gpu, and the hacks
                // need to know about this before we have an easy to probe OpenGL context online:
                if (dev->vendor_id == PCI_VENDOR_ID_INTEL) fDeviceType = kPsychIntelIGP;

                // Skip intel gpu's, unless the PSYCH_ALLOW_DANGEROUS env variable is set:
                // Intel IGP's have a design defect which can cause machine hard lockup if multiple
                // regs are accessed simultaneously! As we can't serialize our MMIO reads with the
                // kms-driver, using our MMIO code on Intel is unsafe. Horrible crashes are reported
                // against Haswell on the freedesktop bug tracker for this issue.
                if ((dev->vendor_id == PCI_VENDOR_ID_INTEL) && !getenv("PSYCH_ALLOW_DANGEROUS")) {
                    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Skipping detected Intel GPU for safety reasons. setenv('PSYCH_ALLOW_DANGEROUS', '1') to override.\n");
                    continue;
                }

                // Select the targetgpuidx'th detected gpu:
                // TODO: Replace this hack by true multi-gpu support and - far in the future? -
                // automatic mapping of screens to gpu's:
                if (currentgpuidx >= targetgpuidx) {
                    if ((PsychPrefStateGet_Verbosity() > 2) && (targetgpuidx >= 0)) printf("PTB-INFO: Choosing GPU number %i for low-level access during this session.\n", currentgpuidx);

                    // Assign as gpu:
                    gpu = dev;
                    break;
                }

                currentgpuidx++;
            }
        }
    }

    // Enumeration finished - Release iterator:
    pci_iterator_destroy(iter);

    // Found matching GPU?
    if (gpu) {
        // Yes!
        if (PsychPrefStateGet_Verbosity() > 2) {
            printf("PTB-INFO: %s - %s GPU found. Trying to establish low-level access...\n", pci_device_get_vendor_name(gpu), pci_device_get_device_name(gpu));
            fflush(NULL);
        }

        // Need to zero-out errno to work around a bug in shipping libpciaccess.so versions prior June 2011 which would cause
        // pci_device_probe(gpu) to report failure even on success at every invocation after the 1st invocation.
        // See <http://cgit.freedesktop.org/xorg/lib/libpciaccess/commit/src/linux_sysfs.c?id=f9159b97834ba4b4e42a07953a33866e7ac90dbd>
        errno = 0;

        // Pull in remaining info about gpu:
        ret = pci_device_probe(gpu);
        if (ret) {
            if (PsychPrefStateGet_Verbosity() > 1) {
                printf("PTB-INFO: Could not probe properties of GPU for screenId %i [%s]\n", screenId, strerror(ret));
                printf("PTB-INFO: Beamposition timestamping and other special features disabled.\n");
                fflush(NULL);
            }

            gpu = NULL;

            // Cleanup:
            pci_system_cleanup();

            return(FALSE);
        }

        // Store PCI device id:
        fPCIDeviceId = gpu->device_id;

        // Find out which BAR to use for mapping MMIO registers. Depends on GPU vendor:
        if (gpu->vendor_id == PCI_VENDOR_ID_NVIDIA) {
            // BAR 0 is MMIO:
            region = &gpu->regions[0];
            fDeviceType = kPsychGeForce;
            fNumDisplayHeads = 2;
        }

        if (gpu->vendor_id == PCI_VENDOR_ID_ATI || gpu->vendor_id == PCI_VENDOR_ID_AMD) {
            // BAR 2 is MMIO on old AMD gpus, BAR 5 is MMIO on DCE-8/10/11/... "Sea Islands" gpus and later models:
            region = &gpu->regions[(isDCE8(screenId) || isDCE10(screenId)) ? 5 : 2];
            fDeviceType = kPsychRadeon;
            fNumDisplayHeads = 2;
        }

        if (gpu->vendor_id == PCI_VENDOR_ID_INTEL) {
            // On non GEN-2 hardware, BAR 0 is MMIO:
            region = &gpu->regions[0];
            fCardType = 0;

            // On GEN-2 hardware, BAR 1 is MMIO: Detect known IGP's of GEN-2.
            if ((fPCIDeviceId == 0x3577) || (fPCIDeviceId == 0x2562) || (fPCIDeviceId == 0x3582) || (fPCIDeviceId == 0x358e) || (fPCIDeviceId == 0x2572)) {
                region = &gpu->regions[1];
                fCardType = 2;
            }

            fDeviceType = kPsychIntelIGP;

            // GEN-7+ (IvyBridge and later) and maybe GEN-6 (SandyBridge) has 3 display
            // heads, older IGP's have 2. Let's be optimistic and assume 3, to safe us
            // from lots of new detection code:
            fNumDisplayHeads = 3;
        }

        if (region) {
            // Try to MMAP MMIO registers with write access, assign their base address to gfx_cntl_mem on success:
            if (PsychPrefStateGet_Verbosity() > 4) {
                printf("PTB-DEBUG: Mapping GPU BAR address %p ...\n", region->base_addr);
                printf("PTB-DEBUG: Mapping %p bytes...\n", region->size);
                fflush(NULL);
            }

            ret = pci_device_map_range(gpu, region->base_addr, region->size, PCI_DEV_MAP_FLAG_WRITABLE, (void**) &gfx_cntl_mem);
            // Mapping MMIO for write access is a nono on Intel with latest kernels, so retry a readonly mapping:
            if ((ret == EAGAIN) && (fDeviceType == kPsychIntelIGP)) {
                ret = pci_device_map_range(gpu, region->base_addr, region->size, 0, (void**) &gfx_cntl_mem);
            }
        }
        else {
            // Unsupported GPU type:
            gfx_cntl_mem = NULL;
            ret = 0;
        }

        if (ret || (NULL == gfx_cntl_mem)) {
            if (PsychPrefStateGet_Verbosity() > 1) {
                printf("PTB-INFO: Failed to map GPU low-level control registers for screenId %i [%s].\n", screenId, strerror(ret));
                printf("PTB-INFO: Beamposition timestamping and other special functions disabled.\n");
                printf("PTB-INFO: You need to run Matlab/Octave with root-privileges, or run the script PsychLinuxConfiguration once for this to work.\n");
                printf("PTB-INFO: However, if you are using the free graphics drivers, there isn't any need for this.\n");
                fflush(NULL);
            }

            // Failed:
            gpu = NULL;

            // Cleanup:
            pci_system_cleanup();

            return(FALSE);
        }

        // Success! Identify GPU:
        gfx_length = region->size;

        // Lowest allowable MMIO register offset for given GPU:
        gfx_lowlimit = 0;

        if (fDeviceType == kPsychGeForce) {
            fCardType = PsychGetNVidiaGPUType(NULL);

            // GPU powered down and therefore offline? Signalled by special return code 0xffffffff.
            if (fCardType == 0xffffffff) {
                // Yes. This can happen in Optimus setups or other hybrid graphics setups.
                // Bail out of this GPU and cleanup:
                PsychScreenUnmapDeviceMemory();

                if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: This NVidia GPU is powered down, probably in a Optimus setup. Skipping its use for beamposition timestamping.\n");

                // We are done.
                return(FALSE);
            }

            // NV-E0 "Kepler" and later have 4 display heads:
            if ((fCardType == 0x0) || (fCardType >= 0xe0)) fNumDisplayHeads = 4;

            if (PsychPrefStateGet_Verbosity() > 2) {
                printf("PTB-INFO: Connected to NVidia %s GPU of NV-%03x family with %i display heads. Beamposition timestamping enabled.\n", pci_device_get_device_name(gpu), fCardType, fNumDisplayHeads);
                fflush(NULL);
            }
        }

        if (fDeviceType == kPsychRadeon) {
            // On Radeons we distinguish between Avivo / DCE-2 (10), DCE-3 (30), or DCE-4 style (40) or DCE-5 (50) or DCE-6 (60), DCE-8 (80), DCE-10 (100), DCE-11 (110) for now.
            fCardType = isDCE11(screenId) ? 110 : isDCE10(screenId) ? 100 : isDCE8(screenId) ? 80 : (isDCE6(screenId) ? 60 : (isDCE5(screenId) ? 50 : (isDCE4(screenId) ? 40 : (isDCE3(screenId) ? 30 : 10))));

            // Setup for DCE-4/5/6/8:
            if (isDCE4(screenId) || isDCE5(screenId) || isDCE6(screenId) || isDCE8(screenId)) {
                gfx_lowlimit = 0;

                // Offset of crtc blocks of evergreen gpu's for each of the six possible crtc's:
                crtcoff[0] = EVERGREEN_CRTC0_REGISTER_OFFSET;
                crtcoff[1] = EVERGREEN_CRTC1_REGISTER_OFFSET;
                crtcoff[2] = EVERGREEN_CRTC2_REGISTER_OFFSET;
                crtcoff[3] = EVERGREEN_CRTC3_REGISTER_OFFSET;
                crtcoff[4] = EVERGREEN_CRTC4_REGISTER_OFFSET;
                crtcoff[5] = EVERGREEN_CRTC5_REGISTER_OFFSET;

                // Also, DCE-4 and DCE-5 and DCE-6, but not DCE-4.1 or DCE-6.4 (which have only 2) or DCE-6.1 (4 heads), supports up to six display heads:
                if (!isDCE41(screenId) && !isDCE61(screenId) && !isDCE64(screenId)) fNumDisplayHeads = 6;

                // DCE-6.1 "Trinity" chip family supports 4 display heads:
                if (!isDCE41(screenId) && isDCE61(screenId)) fNumDisplayHeads = 4;
            }

            // Setup for DCE-10/11:
            if (isDCE10(screenId) || isDCE11(screenId)) {
                // DCE-10/11 of the "Volcanic Islands" gpu family uses (mostly) the same register specs,
                // but the offsets for the different CRTC blocks are different wrt. to pre DCE-10. Therefore
                // need to initialize the offsets differently. Also, some of these parts seem to support up
                // to 7 display engines instead of the old limit of 6 engines:
                gfx_lowlimit = 0;

                // Offset of crtc blocks of Volcanic Islands DCE-10/11 gpu's for each of the possible crtc's:
                crtcoff[0] = DCE10_CRTC0_REGISTER_OFFSET;
                crtcoff[1] = DCE10_CRTC1_REGISTER_OFFSET;
                crtcoff[2] = DCE10_CRTC2_REGISTER_OFFSET;
                crtcoff[3] = DCE10_CRTC3_REGISTER_OFFSET;
                crtcoff[4] = DCE10_CRTC4_REGISTER_OFFSET;
                crtcoff[5] = DCE10_CRTC5_REGISTER_OFFSET;
                crtcoff[6] = DCE10_CRTC6_REGISTER_OFFSET;

                // DCE-10 has 6 display controllers:
                if (isDCE10(screenId)) fNumDisplayHeads = 6;

                // DCE-11 has 3 display controllers:
                if (isDCE11(screenId)) fNumDisplayHeads = 3;
            }

            if (PsychPrefStateGet_Verbosity() > 2) {
                printf("PTB-INFO: Connected to %s %s GPU with DCE-%.1f display engine [%i heads]. Beamposition timestamping enabled.\n", pci_device_get_vendor_name(gpu), pci_device_get_device_name(gpu), (float) fCardType / 10, fNumDisplayHeads);
                fflush(NULL);
            }
        }

        if (fDeviceType == kPsychIntelIGP) {
            if (PsychPrefStateGet_Verbosity() > 2) {
                printf("PTB-INFO: Connected to Intel %s GPU%s. Beamposition timestamping enabled.\n", pci_device_get_device_name(gpu), (fCardType == 2) ? " of GEN-2 type" : "");
                fflush(NULL);
            }
        }

        // Perform auto-detection of screen to head mappings: This will no-op if users script
        // has already manually specified mappings via Screen('Preference','ScreenToHead', ...):
        PsychAutoDetectScreenToHeadMappings(fNumDisplayHeads);

        // Ready to rock!
    } else {
        // No candidate.
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: No low-level controllable GPU on screenId %i. Beamposition timestamping and other special functions disabled.\n", screenId);
        fflush(NULL);

        // Cleanup:
        pci_system_cleanup();
    }

    // Keep track if something is mapped:
    if (gfx_cntl_mem) numKernelDrivers++;

    // Return final success or failure status:
    return((gfx_cntl_mem) ? TRUE : FALSE);
}

/*
 * Return identifying information about GPU for a given screen screenNumber:
 *
 * Returns TRUE on success, and the actual info in int variables, FALSE if info
 * not available:
 * Input: screenNumber of the screen for which to query GPU.
 *
 * Output: All optional - NULL == Don't return info.
 *
 * gpuMaintype = Basically what vendor.
 * gpuMinortype = Vendor specific id meaningful to us to define a certain class or generation of hardware.
 * pciDeviceId = The PCI device id.
 * numDisplayHeads = Maximum number of crtc's.
 *
 */
psych_bool PsychGetGPUSpecs(int screenNumber, int* gpuMaintype, int* gpuMinortype, int* pciDeviceId, int* numDisplayHeads)
{
  // Provide the basic device type, ie., unknown, intel, amd, ...
  if (gpuMaintype) *gpuMaintype = fDeviceType;

  // Remaining info is only available for mapped gpu's:
  if (!PsychOSIsKernelDriverAvailable(screenNumber)) return(FALSE);

  if (gpuMinortype) *gpuMinortype = fCardType;
  if (pciDeviceId) *pciDeviceId = fPCIDeviceId;
  if (numDisplayHeads) *numDisplayHeads = fNumDisplayHeads;

  return(TRUE);
}

// These are shared with display backend specific code, e.g., PsychScreenGlueWayland:
PsychScreenSettingsType     displayOriginalCGSettings[kPsychMaxPossibleDisplays];            //these track the original video state before the Psychtoolbox changed it.
psych_bool                  displayOriginalCGSettingsValid[kPsychMaxPossibleDisplays];
psych_bool                  displayCursorHidden[kPsychMaxPossibleDisplays];
CGDisplayCount              numDisplays;

// displayCGIDs stores the X11 Display* handles to the display connections of each PTB logical screen:
CGDirectDisplayID           displayCGIDs[kPsychMaxPossibleDisplays];

// Only used in this file:
static psych_bool           displayBeampositionHealthy[kPsychMaxPossibleDisplays];
static psych_bool           displayLockSettingsFlags[kPsychMaxPossibleDisplays];
// displayX11Screens stores the mapping of PTB screenNumber's to corresponding X11 screen numbers:
static int                  displayX11Screens[kPsychMaxPossibleDisplays];
static XRRScreenResources*  displayX11ScreenResources[kPsychMaxPossibleDisplays];
static Atom                 displayX11ScreenCompositionAtom[kPsychMaxPossibleDisplays];

// XInput-2 extension data per display:
static int                  xi_opcode = 0, xi_event = 0, xi_error = 0;
static int                  xinput_ndevices[kPsychMaxPossibleDisplays];
static XIDeviceInfo*        xinput_info[kPsychMaxPossibleDisplays];

// File local functions:
void PsychLockScreenSettings(int screenNumber);
void PsychUnlockScreenSettings(int screenNumber);
psych_bool PsychCheckScreenSettingsLock(int screenNumber);
void InitPsychtoolboxKernelDriverInterface(void);

// X11 has a different - and much more powerful and flexible - concept of displays than OS-X or Windows:
// One can have multiple X11 connections to different logical displays. A logical display corresponds
// to a specific X-Server. This X-Server could run on the same machine as Matlab+PTB or on a different
// machine connected via network somewhere in the building or the world. A single machine can even run
// multiple X-Servers. Each display itself can consist of multiple screens. Each screen represents
// a single physical display device. E.g., a dual-head gfx-adaptor could be driven by a single X-Server and have
// two screens for each physical output. A single X-Server could also drive multiple different gfx-cards
// and therefore have many screens. A Linux render-cluster could consist of multiple independent machines,
// each with multiple screens aka gfx heads connected to each machine (aka X11 display).
//
// By default, PTB just connects to the same display as the one that Matlab is running on and tries to
// detect and enumerate all physical screens connected to that display. The default display is set either
// via Matlab command option "-display" or via the Shell environment variable $DISPLAY. Typically, it
// is simply $DISPLAY=:0.0, which means the local gfx-adaptor attached to the machine the user is logged into.
//
// If a user wants to make use of other displays than the one Matlab is running on, (s)he can set the
// environment variable $PSYCHTOOLBOX_DISPLAYS to a list of all requested displays. PTB will then try
// to connect to each of the listed displays, enumerate all attached screens and build its list of
// available screens as a merge of all screens of all displays.
// E.g., export PSYCHTOOLBOX_DISPLAYS=":0.0,kiwi.kyb.local:0.0,coriander.kyb.local:0.0" would enumerate
// all screens of all gfx-adaptors on the local machine ":0.0", and the network connected machines
// "kiwi.kyb.local" and "coriander.kyb.local".
//
// Possible applications: Multi-display setups on Linux, possibly across machines, e.g., render-clusters
// Weird experiments with special setups. Show stimulus on display 1, query mouse or keyboard from
// different machine... 

static int x11_errorval = 0;
static int x11_errorbase = 0;
static int (*x11_olderrorhandler)(Display*, XErrorEvent*);

//file local functions
void InitCGDisplayIDList(void);

/* Lock graphics access
 *
 * This is mostly to protect access to XLib functions (and also GLX functions, as they
 * are using XLib/X11 as transport), but could also be used to protect access to other
 * backend libraries for display configuration and stimulus onset control.
 *
 * We use this as we don't want to rely on XLib's builtin locking and thread-safety.
 * Why? Because XLib's multi-threading protection must be initialized at application
 * startup time *before any other XLib calls* via a call to XInitThreads(), and we
 * do not have any control if or when that XInitThreads() call happens, so it is better
 * to not rely on Octave, Matlab or other future hosting environments doing the right
 * thing and implement our own locking around potentially racy XLib calls.
 *
 * XLib itself is reentrant, but concurrent access by multiple threads to the same xdisplay
 * connection handle, and thereby to the same x-windows connection queue is not safe.
 * We use one shared xdisplay handle for all onscreen windows and screens, across all of
 * Screen's threads (main thread and async flipper/frame-sequential stereo threads). We
 * need to use one shared handle because otherwise OpenGL context resource sharing and the
 * OML_sync_control timing and timestamping functions won't work properly. Therefore we
 * must lock-protect all calls by all threads to XLib. N.b. PsychHID uses its own private
 * xdisplay handle for its KbQueue thread, but that handle is really only used exclusively
 * by one thread, so no need to worry about locking that one.
 *
 * For now we only maintain a global lock and extend this into finer-grained locking if
 * neccessary. Using XCB would be another option, but that would be a huge rewrite with
 * rather sparse documentation available...
 */
psych_mutex displayLock;
double tLockDisplay;

void PsychLockDisplay(void)
{
    PsychLockMutex(&displayLock);

    if (PsychPrefStateGet_Verbosity() > 15) {
        printf("PTB-DEBUG: PsychLockDisplay(): Display locked!\n");
        fflush(NULL);
        PsychGetAdjustedPrecisionTimerSeconds(&tLockDisplay);
    }
}

void PsychUnlockDisplay(void)
{
    if (PsychPrefStateGet_Verbosity() > 15) {
        double tUnlockDisplay;
        PsychGetAdjustedPrecisionTimerSeconds(&tUnlockDisplay);
        printf("PTB-DEBUG: PsychUnlockDisplay(): Display unlocked! Lock hold time was %f msecs.\n", 1000 * (tUnlockDisplay - tLockDisplay));
        fflush(NULL);
    }

    PsychUnlockMutex(&displayLock);
}

int PsychGetXScreenIdForScreen(int screenNumber)
{
    if ((screenNumber >= numDisplays) || (screenNumber < 0)) PsychErrorExit(PsychError_invalidScumber);
    return(displayX11Screens[screenNumber]);
}

void PsychGetCGDisplayIDFromScreenNumber(CGDirectDisplayID *displayID, int screenNumber)
{
    if ((screenNumber >= numDisplays) || (screenNumber < 0)) PsychErrorExit(PsychError_invalidScumber);
    *displayID=displayCGIDs[screenNumber];
}

/*  About locking display settings:
 *
 *    SCREENOpenWindow and SCREENOpenOffscreenWindow  set the lock when opening  windows and
 *    SCREENCloseWindow unsets upon the close of the last of a screen's windows. PsychSetVideoSettings checks for a lock
 *    before changing the settings.  Anything (SCREENOpenWindow or SCREENResolutions) which attemps to change
 *    the display settings should report that attempts to change a dipslay's settings are not allowed when its windows are open.
 *
 *    PsychSetVideoSettings() may be called by either SCREENOpenWindow or by Resolutions().  If called by Resolutions it both sets the video settings
 *    and caches the video settings so that subsequent calls to OpenWindow can use the cached mode regardless of whether interceding calls to OpenWindow
 *    have changed the display settings or reverted to the virgin display settings by closing.  SCREENOpenWindow should thus invoke the cached mode
 *    settings if they have been specified and not current actual display settings.
 *
 */
void PsychLockScreenSettings(int screenNumber)
{
    displayLockSettingsFlags[screenNumber]=TRUE;
}

void PsychUnlockScreenSettings(int screenNumber)
{
    displayLockSettingsFlags[screenNumber]=FALSE;
}

psych_bool PsychCheckScreenSettingsLock(int screenNumber)
{
    return(displayLockSettingsFlags[screenNumber]);
}

/* Because capture and lock will always be used in conjuction, capture calls lock, and SCREENOpenWindow must only call Capture and Release */
void PsychCaptureScreen(int screenNumber)
{
    if(screenNumber>=numDisplays) PsychErrorExit(PsychError_invalidScumber);
    PsychLockScreenSettings(screenNumber);
}

/*
 *    PsychReleaseScreen()
 */
void PsychReleaseScreen(int screenNumber)
{
    if(screenNumber>=numDisplays) PsychErrorExit(PsychError_invalidScumber);
    PsychUnlockScreenSettings(screenNumber);
}

psych_bool PsychIsScreenCaptured(int screenNumber)
{
    return(PsychCheckScreenSettingsLock(screenNumber));
}

//Read display parameters.
/*
 *    PsychGetNumDisplays()
 *    Get the number of video displays connected to the system.
 */
int PsychGetNumDisplays(void)
{
    return((int) numDisplays);
}

//Initialization functions
void InitializePsychDisplayGlue(void)
{
    static psych_bool firstTime = TRUE;
    int i;

    //init the display settings flags.
    for(i=0;i<kPsychMaxPossibleDisplays;i++){
        displayLockSettingsFlags[i]=FALSE;
        displayOriginalCGSettingsValid[i]=FALSE;
        displayCursorHidden[i]=FALSE;
        displayBeampositionHealthy[i]=TRUE;
        displayX11ScreenResources[i] = NULL;
        displayX11ScreenCompositionAtom[i] = None;
        xinput_ndevices[i]=0;
        xinput_info[i]=NULL;
    }

    has_xrandr_1_2 = FALSE;
    has_xrandr_1_3 = FALSE;

    // Set Mesa version to undefined:
    mesaversion[0] = mesaversion[1] = mesaversion[2] = 0;

    // Initialize our mutex for locking of display function access, e.g., XLib/GLX calls:
    PsychInitMutex(&displayLock);

    if (firstTime) {
        firstTime = FALSE;

        // We must initialize XLib for multithreading-safe operations / access on first
        // call if usercode explicitely requests this via environment variable PSYCH_XINITTHREADS.
        //
        // We can only do this as opt-in, as XInitThreads() must be called as very first
        // XLib function after process startup or bad things will happen! We don't have control
        // over this wrt. Matlab or Octave (especially future Octave 3.7+ with its QT based GUI),
        // so we implemented our own locking in Screen() and don't need it in PsychHID, as PsychHID's
        // x-connection is exclusively used by PsychHID's Xinput processing thread. However, there
        // may be some cases when our own locking is insufficient, due to deficiencies in the
        // DRI2 XOrg FOSS Mesa graphics driver stack, so some users may want to opt-into use
        // XLib's threading protection as a work-around if they can guarantee Octave or Matlab
        // hasn't called any XLib calls already during its running session:
        if (getenv("PSYCH_XINITTHREADS")) XInitThreads();
    }

    //init the list of Core Graphics display IDs.
    InitCGDisplayIDList();

    // Attach to kernel-level Psychtoolbox graphics card interface driver if possible
    // *and* allowed by settings, setup all relevant mappings:
    InitPsychtoolboxKernelDriverInterface();
}

int PsychGetScreenDepthValue(int screenNumber)
{
    PsychDepthType    depthStruct;

    PsychInitDepthStruct(&depthStruct);
    PsychGetScreenDepth(screenNumber, &depthStruct);
    return(PsychGetValueFromDepthStruct(0,&depthStruct));
}

/*
 *    PsychGetScreenSettings()
 *
 *    Fills a structure describing the screen settings such as x, y, depth, frequency, etc.
 *
 *    Consider inverting the calling sequence so that this function is at the bottom of call hierarchy.
 */
void PsychGetScreenSettings(int screenNumber, PsychScreenSettingsType *settings)
{
    settings->screenNumber=screenNumber;
    PsychGetScreenRect(screenNumber, settings->rect);
    PsychInitDepthStruct(&(settings->depth));
    PsychGetScreenDepth(screenNumber, &(settings->depth));
    settings->mode=PsychGetColorModeFromDepthStruct(&(settings->depth));
    settings->nominalFrameRate= (int) (PsychGetNominalFramerate(screenNumber) + 0.5);
    //settings->dacbits=PsychGetDacBits(screenNumber);
}

PsychColorModeType PsychGetScreenMode(int screenNumber)
{
    PsychDepthType depth;

    PsychInitDepthStruct(&depth);
    PsychGetScreenDepth(screenNumber, &depth);
    return(PsychGetColorModeFromDepthStruct(&depth));
}

// TODO: Hack to keep running under Wayland:
XIDeviceInfo* PsychGetInputDevicesForScreen(int screenNumber, int* nDevices)
{
    if(screenNumber >= numDisplays) PsychErrorExit(PsychError_invalidScumber);
    if (nDevices) *nDevices = xinput_ndevices[screenNumber];
    return(xinput_info[screenNumber]);
}

// Wayland native builds do have their own screen glue for the X11 specific bits:
#ifndef PTB_USE_WAYLAND

psych_bool PsychGetCGModeFromVideoSetting(CFDictionaryRef *cgMode, PsychScreenSettingsType *setting);

/*
 *    PsychCheckVideoSettings()
 *
 *    Check all available video display modes for the specified screen number and return true if the
 *    settings are valid and false otherwise.
 */
psych_bool PsychCheckVideoSettings(PsychScreenSettingsType *setting)
{
    CFDictionaryRef cgMode;
    return(PsychGetCGModeFromVideoSetting(&cgMode, setting));
}

// Init XInput extension: Called under display lock protection:
static void InitXInputExtensionForDisplay(CGDirectDisplayID dpy, int idx)
{
    int major, minor;
    int rc, i;

    // XInputExtension supported? If so do basic init:
    if (!XQueryExtension(dpy, "XInputExtension", &xi_opcode, &xi_event, &xi_error)) {
        printf("PTB-WARNING: XINPUT1/2 extension unsupported. Will only be able to handle one mouse and mouse cursor.\n");
        goto out;
    }

    // XInput V 2.0 or later supported?
    major = 2;
    minor = 0;
    rc = XIQueryVersion(dpy, &major, &minor);
    if (rc == BadRequest) {
        printf("PTB-WARNING: No XInput-2 support. Server supports version %d.%d only.\n", major, minor);
        printf("PTB-WARNING: XINPUT1/2 extension unsupported. Will only be able to handle one mouse and mouse cursor.\n");
        goto out;
    } else if (rc != Success) {
        printf("PTB-ERROR: Internal error during XInput-2 extension init sequence! This is a bug in Xlib!\n");
        printf("PTB-WARNING: XINPUT1/2 extension unsupported. Will only be able to handle one mouse and mouse cursor.\n");
        goto out;
    }

    // printf("PsychHID: INFO: XI2 supported. Server provides version %d.%d.\n", major, minor);

    // Enumerate all XI2 input devices for this x-display:
    xinput_info[idx] = XIQueryDevice(dpy, XIAllDevices, &xinput_ndevices[idx]);

    out:
    return;
}

// ProcessRandREvents: Must be called called under display lock protection!
static void ProcessRandREvents(int screenNumber)
{
    XEvent evt;

    if (!has_xrandr_1_2) return;

    // Check for screen config change events and dispatch them:
    while (XCheckTypedWindowEvent(displayCGIDs[screenNumber], RootWindow(displayCGIDs[screenNumber], displayX11Screens[screenNumber]), xr_event + RRScreenChangeNotify, &evt)) {
        // Screen changed: Dispatch new configuration to X-Lib:
        XRRUpdateConfiguration(&evt);
    }
}

// GetRandRScreenConfig: Must be called called under display lock protection!
static void GetRandRScreenConfig(CGDirectDisplayID dpy, int idx)
{
    int major, minor;
    int o, m, num_crtcs, isPrimary, crtcid, crtccount;
    int primaryOutput = -1, primaryCRTC = -1, primaryCRTCIdx = -1;
    int crtcs[100];

    // Preinit to "undefined":
    displayX11ScreenResources[idx] = NULL;

    // XRandR extension supported? If so do basic init:
    if (!XRRQueryExtension(dpy, &xr_event, &xr_error) ||
        !XRRQueryVersion(dpy, &major, &minor)) {
        printf("PTB-WARNING: XRandR extension unsupported. Display infos and configuration functions will be very limited!\n");
        return;
    }

    // Detect version of XRandR:
    if (major > 1 || (major == 1 && minor >= 2)) has_xrandr_1_2 = TRUE;
    if (major > 1 || (major == 1 && minor >= 3)) has_xrandr_1_3 = TRUE;

    // Select screen configuration notify events to get delivered to us:
    Window root = RootWindow(dpy, displayX11Screens[idx]);
    XRRSelectInput(dpy, root, RRScreenChangeNotifyMask);

    // Fetch current screen configuration info for this screen and display:
    XRRScreenResources* res = XRRGetScreenResourcesCurrent(dpy, root);
    displayX11ScreenResources[idx] = res;
    if (NULL == res) {
        printf("PTB-WARNING: Could not query configuration of x-screen %i on display %s. Display infos and configuration will be very limited.\n",
                displayX11Screens[idx], DisplayString(dpy));
        return;
    }

    if (!has_xrandr_1_2) {
        printf("PTB-WARNING: XRandR version 1.2 unsupported! Could not query useful info for x-screen %i on display %s. Infos and configuration will be very limited.\n",
                displayX11Screens[idx], DisplayString(dpy));
        return;
    }

    // Total number of assigned crtc's for this screen:
    crtccount = 0;

    // Iterate over all outputs for this screen:
    for (o = 0; o < res->noutput; o++) {
        XRROutputInfo *output_info = XRRGetOutputInfo(dpy, res, res->outputs[o]);
        if (!output_info) {
            printf("PTB-WARNING: Could not get output info for %i'th output of screen %i [display %s]!\n", o, displayX11Screens[idx], DisplayString(dpy));
            continue;
        }

        // Get info about this output:
        if (has_xrandr_1_3 && (XRRGetOutputPrimary(dpy, root) > 0)) {
            isPrimary = (XRRGetOutputPrimary(dpy, root) == res->outputs[o]) ? 1 : 0;
        }
        else {
            isPrimary = -1;
        }

        for (crtcid = 0; crtcid < res->ncrtc; crtcid++) {
            if (res->crtcs[crtcid] == output_info->crtc) break;
        }
        if (crtcid == res->ncrtc) crtcid = -1;

        // Store crtc for this output:
        crtcs[o] = crtcid;

        if (PsychPrefStateGet_Verbosity() > 3) {
            printf("PTB-INFO: Display '%s' : X-Screen %i : Output %i [%s]: %s : ", DisplayString(dpy), displayX11Screens[idx], o, (const char*) output_info->name, (isPrimary > -1) ? ((isPrimary == 1) ? "Primary output" : "Secondary output") : "Unknown output priority");
            printf("%s : CRTC %i [XID %i]\n", (output_info->connection == RR_Connected) ? "Connected" : "Offline", crtcid, (int) output_info->crtc);
        }

        if ((isPrimary > 0) && (crtcid >= 0)) {
            primaryOutput = o;
            primaryCRTC = crtcid;
        }

        // Is this output - and its crtc - really enabled for this screen?
        if (crtcid >=0) {
            // Yes: Add its crtcid to the list of crtc's for this screen:
            PsychSetScreenToHead(idx, crtcid, crtccount);
            PsychSetScreenToCrtcId(idx, minimum_crtcid, crtccount);

            // Increment id of next allocated crtc scanout engine on GPU:
            // We assume they are allocated in the order of activated outputs,
            // e.g., first output of first X-Screen -> crtc 0, 2nd output of first
            // X-Screen -> crtc 1, first output of 2nd X-Screen -> crtc 2 etc.
            //
            // This is as heuristic as previous approach, but it should continue
            // work as well or as bad as previous one, except it should fix the
            // problem reported in forum message #14200 for AMD Catalyst driver.
            // It should work for bog-standard ZaphodHead setups. It will work in
            // any case on single-display setups or multi-display setups where a single
            // X-Screen spans multiple display outputs aka multiple crtcs.
            //
            // The working assumption is that the user of a ZaphodHead config assigned the different
            // GPU outputs, and thereby their associated physical crtc's, in an ascending order to
            // X-Screens. This is a reasonable assumption, but in no way guaranteed by the system.
            // Therefore this heuristic can go wrong on non-standard ZaphodHead Multi-X-Screen setups.
            // In such cases the user can always use the Screen('Preference', 'ScreenToHead', ...);
            // command or the PSYCHTOOLBOX_PIPEMAPPINGS environment variable to override the wrong
            // mapping, although it would be a pita for complex setups.
            minimum_crtcid++;

            // Increment running count of active outputs (displays) attached to
            // the currently processed X-Screend idx:
            crtccount++;
        }

        // Release info for this output:
        XRRFreeOutputInfo(output_info);
    }

    // Found a defined primary output?
    if (primaryOutput == -1) {
        // Could not find primary output -- none defined. Use first connected
        // output as primary output:
        for (o = 0; o < res->noutput; o++) {
            XRROutputInfo *output_info = XRRGetOutputInfo(dpy, res, res->outputs[o]);
            if (output_info && (output_info->connection == RR_Connected) && (crtcs[o] >= 0)) {
                primaryOutput = o;
                primaryCRTC = crtcs[o];
                XRRFreeOutputInfo(output_info);
                break;
            }

            if (output_info) XRRFreeOutputInfo(output_info);
        }

        // Still undefined? Use first output as primary output:
        if (primaryOutput == -1) {
            primaryOutput = 0;
            primaryCRTC = (crtcs[0] >= 0) ? crtcs[0] : 0;
        }
    }

    // Assign primary crtc of primary output - index 0 - as default display head for this screen:
    // We swap the contents of slot 0 - the primary crtc slot - and whatever slot currently
    // contains the crtcid of our primaryCRTC. This way we shuffle the primary crtc into the
    // 1st slot (zero):
    for (o = 0; o < crtccount; o++) {
        if (PsychScreenToHead(idx, o) == primaryCRTC) {
            PsychSetScreenToHead(idx, PsychScreenToHead(idx, 0), o);
            primaryCRTCIdx = PsychScreenToCrtcId(idx, o);
            PsychSetScreenToCrtcId(idx, PsychScreenToCrtcId(idx, 0), o);
        }
    }

    PsychSetScreenToHead(idx, primaryCRTC, 0);
    PsychSetScreenToCrtcId(idx, primaryCRTCIdx, 0);

    if (PsychPrefStateGet_Verbosity() > 2) {
        printf("PTB-INFO: Display '%s' : X-Screen %i : Assigning primary output as %i with RandR-CRTC %i and GPU-CRTC %i.\n", DisplayString(dpy), displayX11Screens[idx], primaryOutput, primaryCRTC, primaryCRTCIdx);
    }

    return;
}

// Linux only: Retrieve modeline and crtc_info for a specific output on a specific screen:
// Caution: If crtc is non-NULL and receives a valid XRRCrtcInfo*, then this pointer must
//          be released by the caller via XRRFreeCrtcInfo(crtc), or resources will leak!
// Must be called under display lock protection!
XRRModeInfo* PsychOSGetModeLine(int screenId, int outputIdx, XRRCrtcInfo **crtc)
{
    int m;
    XRRModeInfo *mode = NULL;
    XRRCrtcInfo *crtc_info = NULL;

    // Query info about video modeline and crtc of output 'outputIdx':
    XRRScreenResources *res = displayX11ScreenResources[screenId];
    if (has_xrandr_1_2 && (PsychScreenToHead(screenId, outputIdx) >= 0)) {
        crtc_info = XRRGetCrtcInfo(displayCGIDs[screenId], res, res->crtcs[PsychScreenToHead(screenId, outputIdx)]);

        for (m = 0; (m < res->nmode) && crtc_info; m++) {
            if (res->modes[m].id == crtc_info->mode) {
                mode = &res->modes[m];
                break;
            }
        }
    }

    // Optionally return crtc_info in *crtc:
    if (crtc) {
        // Return crtc_info, if any - NULL otherwise:
        *crtc = crtc_info;
    }
    else {
        // crtc_info not required by caller. We release it:
        if (crtc_info) XRRFreeCrtcInfo(crtc_info);
    }

    return(mode);
}

const char* PsychOSGetOutputProps(int screenId, int outputIdx, unsigned long *mm_width, unsigned long *mm_height)
{
    static char outputName[100];
    int o;
    XRROutputInfo *output_info = NULL;
    XRRScreenResources *res = displayX11ScreenResources[screenId];
    RRCrtc crtc = res->crtcs[PsychScreenToHead(screenId, outputIdx)];

    // Find output associated with the crtc for this outputIdx on this screen:
    PsychLockDisplay();
    for (o = 0; o < res->noutput; o++) {
        output_info = XRRGetOutputInfo(displayCGIDs[screenId], res, res->outputs[o]);
        if (output_info->crtc == crtc) break;
        XRRFreeOutputInfo(output_info);
    }
    PsychUnlockDisplay();

    // Got it?
    if (o == res->noutput) PsychErrorExitMsg(PsychError_user, "Invalid output index provided! No such output for this screen!");

    // Store output name to return:
    sprintf(outputName, "%s", output_info->name);

    // And width / height in mm:
    if (mm_width) *mm_width = output_info->mm_width;
    if (mm_height) *mm_height = output_info->mm_height;

    XRRFreeOutputInfo(output_info);

    // Return output name:
    return(&outputName[0]);
}

void InitCGDisplayIDList(void)
{
    int major, minor;
    int rc, i, j, k, count, scrnid;
    char* ptbdisplays = NULL;
    char displayname[1000];
    CGDirectDisplayID x11_dpy = NULL;

    // NULL-out array of displays:
    for(i=0;i<kPsychMaxPossibleDisplays;i++) displayCGIDs[i]=NULL;

    // Preinit screen to head mappings to identity default:
    PsychInitScreenToHeadMappings(0);

    // Initial count of screens is zero:
    numDisplays = 0;

    // Initial minimum allowed crtc id is zero:
    minimum_crtcid = 0;

    // Multiple X11 display specifier strings provided in the environment variable
    // $PSYCHTOOLBOX_DISPLAYS? If so, we connect to all of them and enumerate all
    // available screens on them.
    ptbdisplays = getenv("PSYCHTOOLBOX_DISPLAYS");
    if (ptbdisplays) {
        // Displays explicitely specified. Parse the string and connect to all of them:
        j=0;
        for (i = 0; i <= (int) strlen(ptbdisplays) && j < 1000; i++) {
            // Accepted separators are ',', '"', white-space and end of string...
            if (ptbdisplays[i]==',' || ptbdisplays[i]=='"' || ptbdisplays[i]==' ' || i == (int) strlen(ptbdisplays)) {
                // Separator or end of string detected. Try to connect to display:
                displayname[j]=0;
                printf("PTB-INFO: Trying to connect to X-Display %s ...", displayname);

                x11_dpy = XOpenDisplay(displayname);
                if (x11_dpy == NULL) {
                    // Failed.
                    printf(" ...Failed! Skipping this display...\n");
                }
                else {
                    // Query number of available screens on this X11 display:
                    count=ScreenCount(x11_dpy);
                    scrnid=0;

                    // Set the screenNumber --> X11 display mappings up:
                    for (k=numDisplays; (k<numDisplays + count) && (k<kPsychMaxPossibleDisplays); k++) {
                        if (k == numDisplays) {
                            // 1st entry for this x-display: Init XInput2 extension for it:
                            InitXInputExtensionForDisplay(x11_dpy, numDisplays);
                        } else {
                            // Successive entry. Copy info from 1st entry:
                            xinput_info[k] = xinput_info[numDisplays];
                            xinput_ndevices[k] = xinput_ndevices[numDisplays];
                        }

                        // Mapping of logical screenNumber to X11 Display:
                        displayCGIDs[k]= x11_dpy;
                        // Mapping of logical screenNumber to X11 screenNumber for X11 Display:
                        displayX11Screens[k]=scrnid++;

                        // Get all relevant screen config info and cache it internally:
                        GetRandRScreenConfig(x11_dpy, k);
                    }

                    printf(" ...success! Added %i new physical display screens of %s as PTB screens %i to %i.\n",
                        scrnid, displayname, numDisplays, k-1);

                    // Update total count:
                    numDisplays = k;
                }

                // Reset idx:
                j=0;
            }
            else {
                // Add character to display name:
                displayname[j++]=ptbdisplays[i];
            }
        }

        // At least one screen enumerated?
        if (numDisplays < 1) {
            // We're screwed :(
            PsychErrorExitMsg(PsychError_system, "FATAL ERROR: Couldn't open any X11 display connection to any X-Server!!!");
        }
    }
    else {
        // User didn't setup env-variable with any special displays. We just use
        // the default $DISPLAY or -display of Matlab:
        x11_dpy = XOpenDisplay(NULL);
        if (x11_dpy == NULL) {
            #ifndef PTB_USE_WAFFLE
                // We're screwed :(
                PsychErrorExitMsg(PsychError_system, "FATAL ERROR: Couldn't open default X11 display connection to X-Server!!!");
            #endif

            // No X-Display available, but we are configured with waffle support, so
            // probably user wants to use a non-X11 based display backend.
            printf("PTB-INFO: Could not open any X11/X-Windows system based display connection. Trying other display backends.\n");
            PsychInitNonX11();
            return;
        }

        // Query number of available screens on this X11 display:
        count=ScreenCount(x11_dpy);

        InitXInputExtensionForDisplay(x11_dpy, 0);

        // Set the screenNumber --> X11 display mappings up:
        for (i=0; i<count && i<kPsychMaxPossibleDisplays; i++) {
            displayCGIDs[i]= x11_dpy;
            displayX11Screens[i]=i;
            xinput_info[i] = xinput_info[0];
            xinput_ndevices[i] = xinput_ndevices[0];

            // Get all relevant screen config info and cache it internally:
            GetRandRScreenConfig(x11_dpy, i);
        }
        numDisplays=i;
    }

    if (numDisplays>1) printf("PTB-INFO: A total of %i X-Windows display screens is available for use.\n", numDisplays);

    // Initialize screenId -> GPU headId mapping to identity mappings,
    // unless already setup by XRandR setup code:
    if (!has_xrandr_1_2) PsychInitScreenToHeadMappings(numDisplays);

    // Mark screen+output mappings as not yet overriden by user. We used PsychScreenToCrtcId() above which
    // set that flag, but we only want this flag set if triggered by true override from usercode by use of
    // Screen('Preference', 'ScreenToHead', screenId, ...):
    PsychResetCrtcIdUserOverride();

    // Prepare atoms for "Desktop composition active?" queries:
    // Each atom corresponds to one X-Screen. It is selection-owned by the
    // desktop compositor for that screen, if a compositor is actually active.
    // It is not owned by anybody if desktop composition is off or suspended for
    // that screen, so checking if such an atom has an owner allows to check if
    // the corresponding x-screen is subject to desktop composition atm.
    for (i = 0; i < numDisplays; i++) {
        CGDirectDisplayID dpy;
        char cmatomname[100];

        // Retrieve and cache an atom for this screen on this display:
        PsychGetCGDisplayIDFromScreenNumber(&dpy, i);
        sprintf(cmatomname, "_NET_WM_CM_S%d", PsychGetXScreenIdForScreen(i));
        displayX11ScreenCompositionAtom[i] = XInternAtom(dpy, cmatomname, False);
    }

    return;
}

void PsychInitNonX11(void)
{
    int i;

    // Set the screenNumber --> X11 display mappings up:
    for (i=0; i < kPsychMaxPossibleDisplays; i++) {
        displayCGIDs[i] = NULL;
        displayX11Screens[i] = i;
        xinput_info[i] = NULL;
        xinput_ndevices[i] = 0;
    }

    // Just make something up;
    numDisplays = 1;

    return;
}

void PsychCleanupDisplayGlue(void)
{
    CGDirectDisplayID dpy, last_dpy;
    int i;

    // Make sure the mmio mapping is shut down:
    PsychOSShutdownPsychtoolboxKernelDriverInterface();

    PsychLockDisplay();

    last_dpy = NULL;
    // Go trough full screen list:
    for (i=0; i < PsychGetNumDisplays(); i++) {
        // Get display-ptr for this screen:
        PsychGetCGDisplayIDFromScreenNumber(&dpy, i);

        // No X11 display associated with this screen? Skip it.
        if (!dpy) continue;

        // Did we close this connection already (dpy==last_dpy)?
        if (dpy != last_dpy) {
            // Nope. Keep track of it...
            last_dpy=dpy;
            // ...and close display connection to X-Server:
            if ((PsychGetNumDisplays() > 1) && (mesaversion[0] > 0)) {
                // This is a multi-x-screen (probably ZaphodHead style) setup, we are using the
                // OSS Mesa OpenGL library, and Mesa has been already initialized by us for the
                // session which is just about to end. Current Mesa versions have a bug in both
                // the nouveau and radeon winsys drivers which would cause a crash if we'd close
                // the X connection here and later Screen() gets reloaded again within the same
                // Octave/Matlab session. By avoiding to close the connection here we can avoid
                // this crash bug, at the minor expense of leaking a bit of memory and X resources.
                //
                // The bug has been already fixed by myself upstream for nouveau and radeon,
                // so we just use this workaround on Mesa versions known to contain the bug.
                // These would be current stable Mesa 10.6.1 and earlier stable versions, but
                // backports of these fixes are already available in Mesa 10.5.9 and expected
                // for Mesa 10.6.2+.
                if ((mesaversion[0] > 10) ||
                    ((mesaversion[0] == 10) && (mesaversion[1] >= 7)) ||
                    ((mesaversion[0] == 10) && (mesaversion[1] == 6) && (mesaversion[2] >= 2)) ||
                    ((mesaversion[0] == 10) && (mesaversion[1] == 5) && (mesaversion[2] >= 9))) {
                    // This Mesa version is safe to close the connection:
                    XCloseDisplay(dpy);
                }
                else {
                    if (PsychPrefStateGet_Verbosity() > 3)
                        printf("PTB-INFO: Skipping XCloseDisplay(screen %i) on multi-x-screen setup, as workaround for buggy Mesa version %i.%i.%i.\n",
                               i, mesaversion[0], mesaversion[1], mesaversion[2]);
                }
            }
            else {
                // Always close display connection on single x-screen setups, or if we aren't
                // running under Mesa OpenGL, or if we haven't initialized our Mesa instance
                // yet, as closing the connection is then safe in any case:
                XCloseDisplay(dpy);
            }

            // Release actual xinput info list for this x11 display connection:
            if (xinput_info[i]) {
                XIFreeDeviceInfo(xinput_info[i]);
            }
        }

        // NULL-Out xinput extension data:
        xinput_info[i] = NULL;
        xinput_ndevices[i] = 0;

        // Release per-screen RandR info structures:
        if (displayX11ScreenResources[i]) XRRFreeScreenResources(displayX11ScreenResources[i]);
        displayX11ScreenResources[i] = NULL;
    }

    PsychUnlockDisplay();

    // Destroy the display lock mutex, now that we're done with it for this Screen() session instance:
    PsychDestroyMutex(&displayLock);

    // All connections should be closed now. We can't NULL-out the display list, but
    // Matlab will flush the Screen - Mexfile anyway...
    return;
}

void PsychGetScreenDepths(int screenNumber, PsychDepthType *depths)
{
    int* x11_depths = NULL;
    int  i, count = 0;

    if(screenNumber>=numDisplays || screenNumber < 0) PsychErrorExitMsg(PsychError_internal, "screenNumber is out of range"); //also checked within SCREENPixelSizes

    // Update XLib's view of this screens configuration:
    PsychLockDisplay();
    ProcessRandREvents(screenNumber);

    if (displayCGIDs[screenNumber]) {
        x11_depths = XListDepths(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber), &count);
    }

    PsychUnlockDisplay();

    if (x11_depths && depths && count > 0) {
        // Query successful: Add all values to depth struct:
        for(i=0; i<count; i++) PsychAddValueToDepthStruct(x11_depths[i], depths);
        XFree(x11_depths);
    }
    else {
        // Query failed: Assume at least 32 bits is available.
        printf("PTB-WARNING: Couldn't query available display depths values! Returning a made up list...\n");
        fflush(NULL);
        PsychAddValueToDepthStruct(32, depths);
        PsychAddValueToDepthStruct(24, depths);
        PsychAddValueToDepthStruct(16, depths);
    }
}

double PsychOSVRefreshFromMode(XRRModeInfo *mode)
{
    double dot_clock = (double) mode->dotClock / 1000.0;
    double vrefresh;

    // Sanity check:
    if ((dot_clock <= 0) || (mode->hTotal < 1) || (mode->vTotal < 1)) return(0);

    vrefresh = (((dot_clock * 1000.0) / mode->hTotal) * 1000.0) / mode->vTotal;

    // Divide vrefresh by 1000 to get real Hz - value:
    vrefresh = vrefresh / 1000.0;

    // Doublescan mode? If so, divide vrefresh by 2:
    if (mode->modeFlags & RR_DoubleScan) vrefresh /= 2.0;

    // Interlaced mode? If so, multiply vrefresh by 2:
    if (mode->modeFlags & RR_Interlace) vrefresh *= 2.0;

    return(vrefresh);
}

/*   PsychGetAllSupportedScreenSettings()
 *
 *     Queries the display system for a list of all supported display modes, ie. all valid combinations
 *     of resolution, pixeldepth and refresh rate. Allocates temporary arrays for storage of this list
 *     and returns it to the calling routine. This function is basically only used by Screen('Resolutions').
 */
int PsychGetAllSupportedScreenSettings(int screenNumber, int outputId, long** widths, long** heights, long** hz, long** bpp)
{
    int i, j, o, nsizes, nrates, numPossibleModes;
    XRRModeInfo *mode = NULL;
    XRROutputInfo *output_info = NULL;

    if(screenNumber >= numDisplays || screenNumber < 0) PsychErrorExit(PsychError_invalidScumber);

    // Only supported with RandR 1.2 or later:
    if (!has_xrandr_1_2) return(0);

    if (outputId < 0) {
        PsychLockDisplay();

        // Iterate over all screen sizes and count number of size x refresh rate combos:
        numPossibleModes = 0;
        XRRScreenSize *scs = XRRSizes(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber), &nsizes);
        for (i = 0; i < nsizes; i++) {
            XRRRates(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber), i, &nrates);
            numPossibleModes += nrates;
        }

        PsychUnlockDisplay();

        // Allocate output arrays: These will get auto-released at exit from Screen():
        *widths  = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));
        *heights = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));
        *hz      = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));
        *bpp     = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));

        // Reiterate and fill all slots:
        numPossibleModes = 0;
        for (i = 0; i < nsizes; i++) {
            PsychLockDisplay();
            short* rates = XRRRates(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber), i, &nrates);
            PsychUnlockDisplay();
            for (j = 0; j < nrates; j++) {
                (*widths)[numPossibleModes]  = (long) scs[i].width;
                (*heights)[numPossibleModes] = (long) scs[i].height;
                (*hz)[numPossibleModes]      = (long) rates[j];
                (*bpp)[numPossibleModes]     = (long) PsychGetScreenDepthValue(screenNumber);
                numPossibleModes++;
            }
        }

        // Done:
        return(numPossibleModes);
    }

    // Find crtc for given outputid and screen:
    XRRScreenResources *res = displayX11ScreenResources[screenNumber];
    if (outputId >= kPsychMaxPossibleCrtcs) PsychErrorExitMsg(PsychError_user, "Invalid output index provided! No such output for this screen!");

    outputId = PsychScreenToHead(screenNumber, outputId);
    if (outputId >= res->ncrtc || outputId < 0) PsychErrorExitMsg(PsychError_user, "Invalid output index provided! No such output for this screen!");

    RRCrtc crtc = res->crtcs[outputId];

    // Find output associated with the crtc for this outputId:
    PsychLockDisplay();
    for (o = 0; o < res->noutput; o++) {
        output_info = XRRGetOutputInfo(displayCGIDs[screenNumber], res, res->outputs[o]);
        if (output_info->crtc == crtc) break;
        XRRFreeOutputInfo(output_info);
    }
    PsychUnlockDisplay();

    // Got it?
    if (o == res->noutput) PsychErrorExitMsg(PsychError_user, "Invalid output index provided! No such output for this screen!");

    // Got it: output_info contains a list of all modes (XID's) supported by this
    // display output / crtc combo: Iterate over all of them and return them.
    numPossibleModes = output_info->nmode;

    // Allocate output arrays: These will get auto-released at exit from Screen():
    *widths  = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));
    *heights = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));
    *hz      = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));
    *bpp     = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));

    for (i = 0; i < numPossibleModes; i++) {
        // Fetch modeline for i'th mode:
        for (j = 0; j < res->nmode; j++) {
            if (res->modes[j].id == output_info->modes[i]) break;
        }

        (*widths)[i] = (long) res->modes[j].width;
        (*heights)[i] = (long) res->modes[j].height;
        (*hz)[i] = (long) (PsychOSVRefreshFromMode(&res->modes[j]) + 0.5);
        (*bpp)[i] = (long) 32;
    }

    // Free output info:
    XRRFreeOutputInfo(output_info);

    // Done:
    return(numPossibleModes);
}

/*
 * PsychGetCGModeFromVideoSettings()
 */
psych_bool PsychGetCGModeFromVideoSetting(CFDictionaryRef *cgMode, PsychScreenSettingsType *setting)
{
    int i, j, nsizes = 0, nrates = 0;

    // No op on system without RandR:
    if (!has_xrandr_1_2) {
        // Dummy assignment:
        *cgMode = -1;

        // Also cannot restore display settings at Window / Screen / Runtime close time, so disable it:
        displayOriginalCGSettingsValid[setting->screenNumber] = FALSE;

        // Some info for the reader:
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Getting or setting display video modes unsupported on this graphics driver due to lack of RandR v1.2 support.\n");

        // Return success in good faith that its ok.
        return(TRUE);
    }

    // Extract parameters from setting struct:
    CGDirectDisplayID dpy = displayCGIDs[setting->screenNumber];
    int width  = (int) PsychGetWidthFromRect(setting->rect);
    int height = (int) PsychGetHeightFromRect(setting->rect);
    int fps    = (int) (setting->nominalFrameRate + 0.5);

    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Trying to validate/find closest video mode for requested spec: width = %i x height = %i, rate %i Hz.\n", width, height, fps);

    // Find matching mode:
    int size_index = -1;

    PsychLockDisplay();

    XRRScreenSize *scs = XRRSizes(dpy, PsychGetXScreenIdForScreen(setting->screenNumber), &nsizes);
    for (i = 0; i < nsizes; i++) {
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Testing against mode of resolution w x h = %i x %i with refresh rates: ", scs[i].width, scs[i].height);
        if ((width == scs[i].width) && (height == scs[i].height)) {
            short *rates = XRRRates(dpy, PsychGetXScreenIdForScreen(setting->screenNumber), i, &nrates);
            for (j = 0; j < nrates; j++) {
                if (PsychPrefStateGet_Verbosity() > 3) printf("%i ", (int) rates[j]);
                if (rates[j] == (short) fps) {
                    // Our requested size x fps combo is supported:
                    size_index = i;
                    if (PsychPrefStateGet_Verbosity() > 3) printf("--> Got it! Mode id %i. ", size_index);
                }
            }
        }
        if (PsychPrefStateGet_Verbosity() > 3) printf("\n");
    }

    PsychUnlockDisplay();

    if ((nsizes == 0 || nrates == 0) && (PsychPrefStateGet_Verbosity() > 1)) printf("PTB-WARNING: Getting or setting display video modes unsupported on this graphics driver despite advertised RandR v1.2 support.\n");

    // Found valid settings?
    if (size_index == -1) return(FALSE);

    *cgMode = size_index;
    return(TRUE);
}

/*
    PsychGetScreenDepth()

    The caller must allocate and initialize the depth struct. 
*/
void PsychGetScreenDepth(int screenNumber, PsychDepthType *depth)
{
    if (screenNumber>=numDisplays || screenNumber < 0) PsychErrorExitMsg(PsychError_internal, "screenNumber is out of range"); //also checked within SCREENPixelSizes

    // Update XLib's view of this screens configuration:
    PsychLockDisplay();
    ProcessRandREvents(screenNumber);

    if (displayCGIDs[screenNumber]) {
        PsychAddValueToDepthStruct(DefaultDepth(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber)), depth);
    }
    else {
        PsychAddValueToDepthStruct(32, depth);
    }

    PsychUnlockDisplay();
}

// We build with VidModeExtension support unless forcefully disabled at compile time via a -DNO_VIDMODEEXTS
#ifndef NO_VIDMODEEXTS
#define USE_VIDMODEEXTS 1
#endif

#ifdef USE_VIDMODEEXTS
// Functions for setup and query of hw gamma CLUTS and for monitor refresh rate query:
#include <X11/extensions/xf86vmode.h>
#else
#define XF86VidModeNumberErrors 0
#endif

float PsychGetNominalFramerate(int screenNumber)
{
    if (PsychPrefStateGet_ConserveVRAM() & kPsychIgnoreNominalFramerate) return(0);

    if(screenNumber >= numDisplays || screenNumber < 0)
        PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetNominalFramerate() is out of range");

    // No-Op on non-X11:
    if (!displayCGIDs[screenNumber]) return(0);

#ifdef USE_VIDMODEEXTS

    // Information returned by the XF86VidModeExtension:
    XF86VidModeModeLine mode_line;  // The mode line of the current video mode.
    int dot_clock;                  // The RAMDAC / TDMS pixel clock frequency.

    // We start with a default vrefresh of zero, which means "couldn't query refresh from OS":
    float vrefresh = 0;

    // First we try to get modeline of primary crtc from RandR:
    PsychLockDisplay();
    XRRModeInfo *mode = PsychOSGetModeLine(screenNumber, 0, NULL);
    PsychUnlockDisplay();

    // Modeline with plausible values returned by RandR?
    if (mode && (mode->hTotal > mode->width) && (mode->vTotal > mode->height)) {
        if (PsychPrefStateGet_Verbosity() > 4) {
            printf ("RandR: %s (0x%x) %6.1fMHz\n",
                    mode->name, (int)mode->id,
                    (double)mode->dotClock / 1000000.0);
            printf ("        h: width  %4d start %4d end %4d total %4d skew %4d\n",
                    mode->width, mode->hSyncStart, mode->hSyncEnd,
                    mode->hTotal, mode->hSkew);
            printf ("        v: height %4d start %4d end %4d total %4d\n",
                    mode->height, mode->vSyncStart, mode->vSyncEnd, mode->vTotal);
        }

        dot_clock = (int) ((double) mode->dotClock / 1000.0);
        mode_line.htotal = mode->hTotal;
        mode_line.vtotal = mode->vTotal;
        mode_line.flags = 0;
        mode_line.flags |= (mode->modeFlags & RR_DoubleScan) ? 0x0020 : 0x0;
        mode_line.flags |= (mode->modeFlags & RR_Interlace) ? 0x0010 : 0x0;
    }
    else {
        // No modeline from RandR or invalid modeline. Retry with vidmode extensions:
        if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: PsychGetNominalFramerate: No mode or invalid mode from RandR. Using XF86VidModeExt fallback path...\n");

        PsychLockDisplay();
        if (!XF86VidModeSetClientVersion(displayCGIDs[screenNumber])) {
            // Failed to use VidMode-Extension. We just return a vrefresh of zero.
            if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: PsychGetNominalFramerate: XF86VidModeExt fallback path failed in init.\n");
            PsychUnlockDisplay();
            return(0);
        }

        if (!XF86VidModeGetModeLine(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber), &dot_clock, &mode_line)) {
            // Failed to use VidMode-Extension. We just return a vrefresh of zero.
            if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: PsychGetNominalFramerate: XF86VidModeExt fallback path failed in modeline query.\n");
            PsychUnlockDisplay();
            return(0);
        }

        PsychUnlockDisplay();
    }

    // More child-protection: (utterly needed!)
    if ((dot_clock <= 0) || (mode_line.htotal < 1) || (mode_line.vtotal < 1)) {
        if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: PsychGetNominalFramerate: Invalid modeline retrieved from RandR/VidModeExt. Giving up!\n");
        return(0);
    }

    // Query vertical refresh rate. If it fails we default to the last known good value...
    // Vertical refresh rate is: RAMDAC pixel clock / width of a scanline in clockcylces /
    // number of scanlines per videoframe.
    vrefresh = (((dot_clock * 1000) / mode_line.htotal) * 1000) / mode_line.vtotal;

    // Divide vrefresh by 1000 to get real Hz - value:
    vrefresh = vrefresh / 1000.0f;

    // Definitions from xserver's hw/xfree86/common/xf86str.h
    // V_INTERLACE    = 0x0010,
    // V_DBLSCAN    = 0x0020,

    // Doublescan mode? If so, divide vrefresh by 2:
    if (mode_line.flags & 0x0020) vrefresh /= 2;

    // Interlaced mode? If so, multiply vrefresh by 2:
    if (mode_line.flags & 0x0010) vrefresh *= 2;

    // Done.
    return(vrefresh);
#else
    return(0);
#endif
}

// Error callback handler for X11 errors:
static int x11VidModeErrorHandler(Display* dis, XErrorEvent* err)
{
    // If x11_errorbase not yet setup, simply return and ignore this error:
    if (x11_errorbase == 0) return(0);

    // Setup: Check if its an XVidMode-Error - the only one we do handle.
    if (((err->error_code >=x11_errorbase) && (err->error_code < x11_errorbase + XF86VidModeNumberErrors)) ||
        (err->error_code == BadValue)) {
        // We caused some error. Set error flag:
        x11_errorval = 1;
    }

    // Done.
    return(0);
}

float PsychSetNominalFramerate(int screenNumber, float requestedHz)
{
#ifdef USE_VIDMODEEXTS
    // Information returned by/sent to the XF86VidModeExtension:
    XF86VidModeModeLine mode_line;  // The mode line of the current video mode.
    int dot_clock;                  // The RAMDAC / TDMS pixel clock frequency.
    int rc;
    int event_base;

    // We start with a default vrefresh of zero, which means "couldn't query refresh from OS":
    float vrefresh = 0;

    if(screenNumber>=numDisplays || screenNumber < 0)
        PsychErrorExitMsg(PsychError_internal, "screenNumber is out of range");

    // Not available on non-X11:
    if (!displayCGIDs[screenNumber]) return(0);

    PsychLockDisplay();

    if (!XF86VidModeSetClientVersion(displayCGIDs[screenNumber])) {
        // Failed to use VidMode-Extension. We just return a vrefresh of zero.
        PsychUnlockDisplay();
        return(0);
    }

    if (!XF86VidModeQueryExtension(displayCGIDs[screenNumber], &event_base, &x11_errorbase)) {
        // Failed to use VidMode-Extension. We just return a vrefresh of zero.
        PsychUnlockDisplay();
        return(0);
    }

    // Attach our error callback handler and reset error-state:
    x11_errorval = 0;
    x11_olderrorhandler = XSetErrorHandler(x11VidModeErrorHandler);

    // Step 1: Query current dotclock and modeline:
    if (!XF86VidModeGetModeLine(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber), &dot_clock, &mode_line)) {
        // Restore default error handler:
        XSetErrorHandler(x11_olderrorhandler);
        PsychUnlockDisplay();

        PsychErrorExitMsg(PsychError_internal, "Failed to query video dotclock and modeline!");
    }

    // Step 2: Calculate updated modeline:
    if (requestedHz > 10) {
        // Step 2-a: Given current dot-clock and modeline and requested vrefresh, compute
        // modeline for closest possible match:
        requestedHz*=1000.0f;
        vrefresh = (((dot_clock * 1000) / mode_line.htotal) * 1000) / requestedHz;

        // Assign it to closest modeline setting:
        mode_line.vtotal = (int)(vrefresh + 0.5f);
    }
    else {
        // Step 2-b: Delta mode. requestedHz represents a direct integral offset
        // to add or subtract from current modeline setting:
        mode_line.vtotal+=(int) requestedHz;
    }

    // Step 3: Try to set new modeline:
    if (!XF86VidModeModModeLine(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber), &mode_line)) {
        // Restore default error handler:
        XSetErrorHandler(x11_olderrorhandler);
        PsychUnlockDisplay();

        // Invalid modeline? Signal this:
        return(-1);
    }

    // We synchronize and wait for X-Request completion. If the modeline was invalid,
    // this will trigger an invocation of our errorhandler, which in turn will
    // set the x11_errorval to a non-zero value:
    XSync(displayCGIDs[screenNumber], FALSE);

    // Restore default error handler:
    XSetErrorHandler(x11_olderrorhandler);

    PsychUnlockDisplay();

    // Check for error:
    if (x11_errorval) {
        // Failed to set new mode! Must be invalid. We return -1 to signal this:
        return(-1);
    }

    // No error...

    // Step 4: Query new settings and return them:
    vrefresh = PsychGetNominalFramerate(screenNumber);

    // Done.
    return(vrefresh);
#else
    return(0);
#endif
}

/* Returns the physical display size as reported by X11: */
void PsychGetDisplaySize(int screenNumber, int *width, int *height)
{
    if(screenNumber>=numDisplays || screenNumber < 0)
        PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetDisplaySize() is out of range");

    // Not available on non-X11:
    if (!displayCGIDs[screenNumber]) { *width = 0; *height = 0; return; }

    // Update XLib's view of this screens configuration:
    PsychLockDisplay();
    ProcessRandREvents(screenNumber);

    *width = (int) XDisplayWidthMM(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber));
    *height = (int) XDisplayHeightMM(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber));

    PsychUnlockDisplay();
}

void PsychGetScreenPixelSize(int screenNumber, long *width, long *height)
{
    // For now points == pixels, so just a dumb wrapper, as long as we
    // don't have special "Retina Display" / HiDPI handling in place on X11:
    PsychGetScreenSize(screenNumber, width, height);
}

void PsychGetScreenSize(int screenNumber, long *width, long *height)
{
    if(screenNumber>=numDisplays || screenNumber < 0) PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetScreenDepths() is out of range");

    // Not available on non-X11: MK TODO FIXME - How to get real values?
    if (!displayCGIDs[screenNumber]) { *width = 1680; *height = 1050; return; }

    // Update XLib's view of this screens configuration:
    PsychLockDisplay();
    ProcessRandREvents(screenNumber);

    *width=XDisplayWidth(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber));
    *height=XDisplayHeight(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber));

    PsychUnlockDisplay();
}

void PsychGetGlobalScreenRect(int screenNumber, double *rect)
{
    // Create an empty rect:
    PsychMakeRect(rect, 0, 0, 1, 1);
    // Fill it with meaning by PsychGetScreenRect():
    PsychGetScreenRect(screenNumber, rect);
}

void PsychGetScreenRect(int screenNumber, double *rect)
{
    long width, height; 

    PsychGetScreenSize(screenNumber, &width, &height);
    rect[kPsychLeft]=0;
    rect[kPsychTop]=0;
    rect[kPsychRight]=(int)width;
    rect[kPsychBottom]=(int)height; 
} 

/*
    This is a place holder for a function which uncovers the number of dacbits.  To be filled in at a later date.
    If you know that your card supports >8 then you can fill that in the PsychtPreferences and the psychtoolbox
    will act accordingly.

    There seems to be no way to uncover the dacbits programatically.  According to apple CoreGraphics
    sends a 16-bit word and the driver throws out whatever it chooses not to use.

    For now we just use 8 to avoid false precision.

    If we can uncover the video card model then  we can implement a table lookup of video card model to number of dacbits.
*/
int PsychGetDacBitsFromDisplay(int screenNumber)
{
    return(8);
}

//Set display parameters

/* Linux only: PsychOSSetOutputConfig()
 * Set a video mode and other settings for a specific crtc of a specific output 'outputId'
 * for a specific screen 'screenNumber'.
 *
 * Returns true on success, false on failure.
 */
int PsychOSSetOutputConfig(int screenNumber, int outputId, int newWidth, int newHeight, int newHz, int newX, int newY)
{
    int modeid, maxw, maxh, output, widthMM, heightMM;
    XRRCrtcInfo *crtc_info = NULL, *crtc_info2;
    CGDirectDisplayID dpy = displayCGIDs[screenNumber];
    XRRScreenResources *res = displayX11ScreenResources[screenNumber];

    // Need this later:
    PsychGetDisplaySize(screenNumber, &widthMM, &heightMM);

    if (has_xrandr_1_2 && (PsychScreenToHead(screenNumber, outputId) >= 0)) {
        PsychLockDisplay();
        crtc_info = XRRGetCrtcInfo(dpy, res, res->crtcs[PsychScreenToHead(screenNumber, outputId)]);
        PsychUnlockDisplay();
    }
    else {
        // Failed!
        return(FALSE);
    }

    // Disable auto-restore of screen settings - It would end badly...
    displayOriginalCGSettingsValid[screenNumber] = FALSE;

    // Find matching mode:
    for (modeid = 0; modeid < res->nmode; modeid++) {
        if (((int) res->modes[modeid].width == newWidth) && ((int) res->modes[modeid].height == newHeight) &&
            (newHz == (int)(PsychOSVRefreshFromMode(&res->modes[modeid]) + 0.5))) {
            break;
        }
    }

    // Matching mode found for modesetting?
    if (modeid < res->nmode) {
        PsychLockDisplay();

        // Assign default panning:
        if (newX < 0) newX = crtc_info->x;
        if (newY < 0) newY = crtc_info->y;

        // Iterate over all outputs and compute the new screen bounding box:
        maxw = maxh = 0;
        for (output = 0; (PsychScreenToHead(screenNumber, output) >= 0) && (output < kPsychMaxPossibleCrtcs); output++) {
            if (output == outputId) continue;

            crtc_info2 = XRRGetCrtcInfo(dpy, res, res->crtcs[PsychScreenToHead(screenNumber, output)]);
            if (crtc_info2->x + (int) crtc_info2->width > maxw) maxw = crtc_info2->x + (int) crtc_info2->width;
            if (crtc_info2->y + (int) crtc_info2->height > maxh) maxh = crtc_info2->y + (int) crtc_info2->height;
            XRRFreeCrtcInfo(crtc_info2);
        }

        // Incorporate our soon reconfigured crtc:
        if (newX + newWidth  > maxw) maxw = newX + newWidth;
        if (newY + newHeight > maxh) maxh = newY + newHeight;

        // [0, 0, maxw, maxh] is the new bounding rectangle of the scanned out framebuffer. Set screen size accordingly:

        // Prevent clients from getting confused by our config sequence:
        // XGrabServer(dpy);

        // Disable target crtc:
        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Disabling crtc %i.\n", outputId);
        Status rc = XRRSetCrtcConfig(dpy, res, res->crtcs[PsychScreenToHead(screenNumber, outputId)], crtc_info->timestamp,
                                    0, 0, None, RR_Rotate_0, NULL, 0);

        // Resize screen: MK Don't! Skip this for now, use PsychSetScreenSettings() aka Screen('Resolution') to resize
        // the screen without changing the crtc / output settings. More flexible...
        // if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Resizing screen %i to %i x %i pixels.\n", screenNumber, maxw, maxh);
        // XRRSetScreenSize(dpy, RootWindow(dpy, PsychGetXScreenIdForScreen(screenNumber)), maxw, maxh, widthMM, heightMM);

        // Switch mode of target crtc and reenable it:
        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Enabling crtc %i.\n", outputId);

        crtc_info2 = XRRGetCrtcInfo(dpy, res, res->crtcs[PsychScreenToHead(screenNumber, outputId)]);
        rc = XRRSetCrtcConfig(dpy, res, res->crtcs[PsychScreenToHead(screenNumber, outputId)], crtc_info2->timestamp,
                                newX, newY, res->modes[modeid].id, crtc_info->rotation,
                                crtc_info->outputs, crtc_info->noutput);
        XRRFreeCrtcInfo(crtc_info);
        XRRFreeCrtcInfo(crtc_info2);

        // XUngrabServer(dpy);

        // Make sure the screen change gets noticed by XLib:
        ProcessRandREvents(screenNumber);

        PsychUnlockDisplay();

        return(TRUE);
    } else {
        XRRFreeCrtcInfo(crtc_info);
        return(FALSE);
    }
}

/*
    PsychSetScreenSettings()

    Accept a PsychScreenSettingsType structure holding a video mode and set the display mode accordingly.

    If we can not change the display settings because of a lock (set by open window or close window) then return false.

    SCREENOpenWindow should capture the display before it sets the video mode.  If it doesn't, then PsychSetVideoSettings will
    detect that and exit with an error.  SCREENClose should uncapture the display. 

    The duties of SCREENOpenWindow are:
    -Lock the screen which serves the purpose of preventing changes in video setting with open Windows.
    -Capture the display which gives the application synchronous control of display parameters.

    TO DO: for 8-bit palletized mode there is probably more work to do.

*/

psych_bool PsychSetScreenSettings(psych_bool cacheSettings, PsychScreenSettingsType *settings)
{
    CFDictionaryRef cgMode;
    psych_bool      isValid, isCaptured;
    Rotation        rotation;
    short           rate;
    Time            cfg_timestamp;
    CGDirectDisplayID dpy;

    if (settings->screenNumber >= numDisplays || settings->screenNumber < 0) PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychSetScreenSettings() is out of range");
    dpy = displayCGIDs[settings->screenNumber];

    // Not available on non-X11:
    if (!dpy) return(false);

    //Check for a lock which means onscreen or offscreen windows tied to this screen are currently open.
    // MK Disabled: if(PsychCheckScreenSettingsLock(settings->screenNumber)) return(false);  //calling function should issue an error for attempt to change display settings while windows were open.

    //Check to make sure that this display is captured, which OpenWindow should have done.  If it has not been done, then exit with an error.  
    isCaptured=PsychIsScreenCaptured(settings->screenNumber);
    if(!isCaptured) PsychErrorExitMsg(PsychError_internal, "Attempt to change video settings without capturing the display");

    // Store the original display mode if this is the first time we have called this function.  The psychtoolbox will disregard changes in 
    // the screen state made through the control panel after the Psychtoolbox was launched. That is, OpenWindow will by default continue to 
    // open windows with finder settings which were in place at the first call of OpenWindow.  That's not intuitive, but not much of a problem
    // either. 
    if(!displayOriginalCGSettingsValid[settings->screenNumber]) {
        PsychGetScreenSettings(settings->screenNumber, &displayOriginalCGSettings[settings->screenNumber]);
        displayOriginalCGSettingsValid[settings->screenNumber] = TRUE;
    }

    // Multi-Display configuration?
    if (PsychScreenToHead(settings->screenNumber, 1) != -1) {
        // Yes: At least two display heads attached. We can't use the XRRSetScreenConfigAndRate() method,
        // it is only suitable for single-display setups. In this case, we only set the screen size, aka
        // framebuffer size. User scripts can use the 'ConfigureDisplay' function to setup the crtc's:

        // Also cannot restore display settings at Window / Screen / Runtime close time, so disable it:
        displayOriginalCGSettingsValid[settings->screenNumber] = FALSE;

        // Resize screen:
        int widthMM, heightMM;
        PsychGetDisplaySize(settings->screenNumber, &widthMM, &heightMM);
        int width  = (int) PsychGetWidthFromRect(settings->rect);
        int height = (int) PsychGetHeightFromRect(settings->rect);

        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Resizing screen %i to %i x %i pixels.\n", settings->screenNumber, width, height);

        PsychLockDisplay();

        XRRSetScreenSize(dpy, RootWindow(dpy, PsychGetXScreenIdForScreen(settings->screenNumber)), width, height, widthMM, heightMM);

        // Make sure the screen change gets noticed by XLib:
        ProcessRandREvents(settings->screenNumber);

        PsychUnlockDisplay();

        // Done.
        return(true);
    }

    // Single display configuration, go ahead:

    //Find core graphics video settings which correspond to settings as specified withing by an abstracted psychsettings structure.  
    isValid = PsychGetCGModeFromVideoSetting(&cgMode, settings);
    if (!isValid || (int) cgMode < 0) {
        // This is an internal error because the caller is expected to check first. 
        PsychErrorExitMsg(PsychError_user, "Attempt to set invalid video settings or function unsupported with this graphics-driver.");
    }

    // Change the display mode.
    PsychLockDisplay();
    XRRScreenConfiguration *sc = XRRGetScreenInfo(dpy, RootWindow(dpy, PsychGetXScreenIdForScreen(settings->screenNumber)));

    // Extract parameters from settings struct:
    rate = (short) (settings->nominalFrameRate + 0.5);

    // Fetch current rotation, so we can (re)apply it -- We don't support changing rotation yet:
    XRRConfigCurrentConfiguration(sc, &rotation);

    // Fetch config timestamp so we can prove to the server we're trustworthy:
    Time timestamp = XRRConfigTimes(sc, &cfg_timestamp);

    // Apply new configuration - combo of old rotation with new size (encoded in cgMode) and refresh rate:
    Status rc = XRRSetScreenConfigAndRate(dpy, sc, RootWindow(dpy, PsychGetXScreenIdForScreen(settings->screenNumber)), cgMode, rotation, rate, timestamp);

    // Cleanup:
    XRRFreeScreenConfigInfo(sc);

    // Make sure the screen change gets noticed by XLib:
    ProcessRandREvents(settings->screenNumber);

    PsychUnlockDisplay();

    // Done:
    return((rc != BadValue) ? true : false);
}

/*
    PsychRestoreVideoSettings()

    Restores video settings to the state set by the finder.  Returns TRUE if the settings can be restored or false if they 
    can not be restored because a lock is in effect, which would mean that there are still open windows.    

*/
psych_bool PsychRestoreScreenSettings(int screenNumber)
{
    CFDictionaryRef             cgMode;
    psych_bool                  isValid, isCaptured;
    Rotation                    rotation;
    short                       rate;
    Time                        cfg_timestamp;
    CGDirectDisplayID           dpy;
    PsychScreenSettingsType     *settings;

    if(screenNumber >= numDisplays || screenNumber < 0)
        PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetScreenDepths() is out of range"); //also checked within SCREENPixelSizes

    //Check to make sure that the original graphics settings were cached.  If not, it means that the settings were never changed, so we can just
    //return true. 
    if(!displayOriginalCGSettingsValid[screenNumber]) return(true);

    //Check to make sure that this display is captured, which OpenWindow should have done.  If it has not been done, then exit with an error.  
    isCaptured=PsychIsScreenCaptured(screenNumber);
    if(!isCaptured) PsychErrorExitMsg(PsychError_internal, "Attempt to change video settings without capturing the display");

    // Retrieve original screen settings which we should restore for this screen:
    settings = &displayOriginalCGSettings[screenNumber];

    // Invalidate settings - we want a fresh game after restoring the resolution:
    displayOriginalCGSettingsValid[screenNumber] = FALSE;

    //Find core graphics video settings which correspond to settings as specified withing by an abstracted psychsettings structure.  
    isValid = PsychGetCGModeFromVideoSetting(&cgMode, settings);
    if (!isValid || (int) cgMode < 0){
        // This is an internal error because the caller is expected to check first. 
        PsychErrorExitMsg(PsychError_user, "Attempt to restore invalid video settings or function unsupported with this graphics-driver.");
    }

    //Change the display mode.
    dpy = displayCGIDs[settings->screenNumber];

    PsychLockDisplay();

    XRRScreenConfiguration *sc = XRRGetScreenInfo(dpy, RootWindow(dpy, PsychGetXScreenIdForScreen(settings->screenNumber)));

    // Extract parameters from settings struct:
    rate = (short) (settings->nominalFrameRate + 0.5);

    // Fetch current rotation, so we can (re)apply it -- We don't support changing rotation yet:
    XRRConfigCurrentConfiguration (sc, &rotation);

    // Fetch config timestamp so we can prove to the server we're trustworthy:
    Time timestamp = XRRConfigTimes(sc, &cfg_timestamp);

    // Apply new configuration - combo of old rotation with new size (encoded in cgMode) and refresh rate:
    Status rc = XRRSetScreenConfigAndRate(dpy, sc, RootWindow(dpy, PsychGetXScreenIdForScreen(settings->screenNumber)), cgMode, rotation, rate, timestamp);

    // Cleanup:
    XRRFreeScreenConfigInfo(sc);

    // Make sure the screen change gets noticed by XLib:
    ProcessRandREvents(settings->screenNumber);

    PsychUnlockDisplay();

    // Done:
    return ((rc != BadValue) ? true : false);
}

void PsychOSDefineX11Cursor(int screenNumber, int deviceId, Cursor cursor)
{
    PsychWindowRecordType **windowRecordArray;
    int i, numWindows;

    // Not available on non-X11:
    if (!displayCGIDs[screenNumber]) return;

    // Iterate over all open onscreen windows associated with this screenNumber and
    // apply new X11 cursor definition to each of them:
    PsychCreateVolatileWindowRecordPointerList(&numWindows, &windowRecordArray);

    PsychLockDisplay();
    for(i = 0; i < numWindows; i++) {
        if (PsychIsOnscreenWindow(windowRecordArray[i]) && (windowRecordArray[i]->screenNumber == screenNumber) &&
            windowRecordArray[i]->targetSpecific.xwindowHandle) {
            // Candidate.
            if (deviceId >= 0) {
                // XInput extension for per-device settings:
                XIDefineCursor(displayCGIDs[screenNumber], deviceId, windowRecordArray[i]->targetSpecific.xwindowHandle, cursor);
            }
            else {
                // Old-School global settings:
                XDefineCursor(displayCGIDs[screenNumber], windowRecordArray[i]->targetSpecific.xwindowHandle, cursor);
            }
        }
    }
    PsychUnlockDisplay();
    PsychDestroyVolatileWindowRecordPointerList(windowRecordArray);

    return;
}

void PsychHideCursor(int screenNumber, int deviceIdx)
{
    // Static "Cursor" object which defines a completely transparent - and therefore invisible
    // X11 cursor for the mouse-pointer.
    static Cursor nullCursor = -1;

    // Check for valid screenNumber:
    if((screenNumber >= numDisplays) || (screenNumber < 0)) PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychHideCursor() is out of range"); //also checked within SCREENPixelSizes

    // Not available on non-X11:
    if (!displayCGIDs[screenNumber]) return;

    // Cursor already hidden on screen? If so, nothing to do:
    if ((deviceIdx < 0) && displayCursorHidden[screenNumber]) return;

    // nullCursor already ready?
    if( nullCursor == (Cursor) -1 ) {
        // Create one:
        Pixmap cursormask;
        XGCValues xgc;
        GC gc;
        XColor dummycolour;

        PsychLockDisplay();

        cursormask = XCreatePixmap(displayCGIDs[screenNumber], RootWindow(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber)), 1, 1, 1/*depth*/);
        xgc.function = GXclear;
        gc = XCreateGC(displayCGIDs[screenNumber], cursormask, GCFunction, &xgc );
        XFillRectangle(displayCGIDs[screenNumber], cursormask, gc, 0, 0, 1, 1 );
        dummycolour.pixel = 0;
        dummycolour.red   = 0;
        dummycolour.flags = 04;
        nullCursor = XCreatePixmapCursor(displayCGIDs[screenNumber], cursormask, cursormask, &dummycolour, &dummycolour, 0, 0 );
        XFreePixmap(displayCGIDs[screenNumber], cursormask );
        XFreeGC(displayCGIDs[screenNumber], gc );

        PsychUnlockDisplay();
    }

    if (deviceIdx < 0) {
        // Attach nullCursor to our onscreen window:
        PsychOSDefineX11Cursor(screenNumber, deviceIdx, nullCursor);
        PsychLockDisplay();
        XFlush(displayCGIDs[screenNumber]);
        PsychUnlockDisplay();
        displayCursorHidden[screenNumber]=TRUE;
    } else {
        // XInput cursor: Master pointers only.
        int nDevices;
        XIDeviceInfo* indevs = PsychGetInputDevicesForScreen(screenNumber, &nDevices);

        // Sanity check:
        if (NULL == indevs) PsychErrorExitMsg(PsychError_user, "Sorry, your system does not support individual mouse pointers.");
        if (deviceIdx >= nDevices) PsychErrorExitMsg(PsychError_user, "Invalid 'mouseIndex' provided. No such cursor pointer.");
        if (indevs[deviceIdx].use != XIMasterPointer) PsychErrorExitMsg(PsychError_user, "Invalid 'mouseIndex' provided. No such master cursor pointer.");

        // Attach nullCursor to our onscreen window:
        PsychOSDefineX11Cursor(screenNumber, indevs[deviceIdx].deviceid, nullCursor);
        PsychLockDisplay();
        XFlush(displayCGIDs[screenNumber]);
        PsychUnlockDisplay();
    }

    return;
}

void PsychShowCursor(int screenNumber, int deviceIdx)
{
    Cursor arrowCursor;

    // Check for valid screenNumber:
    if(screenNumber >= numDisplays || screenNumber < 0) PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychHideCursor() is out of range"); //also checked within SCREENPixelSizes

    // Not available on non-X11:
    if (!displayCGIDs[screenNumber]) return;

    if (deviceIdx < 0) {
        // Cursor not hidden on screen? If so, nothing to do:
        if(!displayCursorHidden[screenNumber]) return;

        // Reset to standard Arrow-Type cursor, which is a visible one.
        PsychLockDisplay();
        arrowCursor = XCreateFontCursor(displayCGIDs[screenNumber], 2);
        PsychUnlockDisplay();

        PsychOSDefineX11Cursor(screenNumber, deviceIdx, arrowCursor);

        PsychLockDisplay();
        XFlush(displayCGIDs[screenNumber]);
        PsychUnlockDisplay();

        displayCursorHidden[screenNumber]=FALSE;
    } else {
        // XInput cursor: Master pointers only.
        int nDevices;
        XIDeviceInfo* indevs = PsychGetInputDevicesForScreen(screenNumber, &nDevices);

        // Sanity check:
        if (NULL == indevs) PsychErrorExitMsg(PsychError_user, "Sorry, your system does not support individual mouse pointers.");
        if (deviceIdx >= nDevices) PsychErrorExitMsg(PsychError_user, "Invalid 'mouseIndex' provided. No such cursor pointer.");
        if (indevs[deviceIdx].use != XIMasterPointer) PsychErrorExitMsg(PsychError_user, "Invalid 'mouseIndex' provided. No such master cursor pointer.");

        // Reset to standard Arrow-Type cursor, which is a visible one.
        PsychLockDisplay();
        arrowCursor = XCreateFontCursor(displayCGIDs[screenNumber], 2);
        PsychUnlockDisplay();

        PsychOSDefineX11Cursor(screenNumber, indevs[deviceIdx].deviceid, arrowCursor);

        PsychLockDisplay();
        XFlush(displayCGIDs[screenNumber]);
        PsychUnlockDisplay();
    }
}

void PsychPositionCursor(int screenNumber, int x, int y, int deviceIdx)
{
    // Not available on non-X11:
    if (!displayCGIDs[screenNumber]) return;

    // Reposition the mouse cursor:
    if (deviceIdx < 0) {
        // Core protocol cursor:
        PsychLockDisplay();
        if (XWarpPointer(displayCGIDs[screenNumber], None, RootWindow(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber)), 0, 0, 0, 0, x, y)==BadWindow) {
            PsychUnlockDisplay();
            PsychErrorExitMsg(PsychError_internal, "Couldn't position the mouse cursor! (XWarpPointer() failed).");
        }
    } else {
        // XInput cursor: Master pointers only.
        int nDevices;
        XIDeviceInfo* indevs = PsychGetInputDevicesForScreen(screenNumber, &nDevices);

        // Sanity check:
        if (NULL == indevs) PsychErrorExitMsg(PsychError_user, "Sorry, your system does not support individual mouse pointers.");
        if (deviceIdx >= nDevices) PsychErrorExitMsg(PsychError_user, "Invalid 'mouseIndex' provided. No such cursor pointer.");
        if (indevs[deviceIdx].use != XIMasterPointer) PsychErrorExitMsg(PsychError_user, "Invalid 'mouseIndex' provided. No such master cursor pointer.");

        PsychLockDisplay();
        if (XIWarpPointer(displayCGIDs[screenNumber], indevs[deviceIdx].deviceid, None, RootWindow(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber)), 0, 0, 0, 0, x, y)) {
            PsychUnlockDisplay();
            PsychErrorExitMsg(PsychError_internal, "Couldn't position the mouse cursor! (XIWarpPointer() failed).");
        }
    }

    XFlush(displayCGIDs[screenNumber]);
    PsychUnlockDisplay();
}

/*
    PsychReadNormalizedGammaTable()
*/
void PsychReadNormalizedGammaTable(int screenNumber, int outputId, int *numEntries, float **redTable, float **greenTable, float **blueTable)
{
    CGDirectDisplayID cgDisplayID;
    static float localRed[MAX_GAMMALUT_SIZE], localGreen[MAX_GAMMALUT_SIZE], localBlue[MAX_GAMMALUT_SIZE];

    // The X-Windows hardware LUT has 3 tables for R,G,B.
    // Each entry is a 16-bit word with the n most significant bits used for an n-bit DAC or display encoder:
    psych_uint16    RTable[MAX_GAMMALUT_SIZE];
    psych_uint16    GTable[MAX_GAMMALUT_SIZE];
    psych_uint16    BTable[MAX_GAMMALUT_SIZE];
    int i, n;

    // Initial assumption: Failed.
    n = 0;

    // Not available on non-X11:
    if (!displayCGIDs[screenNumber]) { *numEntries = 0; return; }

    // Query OS for gamma table:
    PsychGetCGDisplayIDFromScreenNumber(&cgDisplayID, screenNumber);

    if (has_xrandr_1_2) {
        // Use RandR V 1.2 for per-crtc query:
        XRRScreenResources *res = displayX11ScreenResources[screenNumber];

        if (outputId >= kPsychMaxPossibleCrtcs) PsychErrorExitMsg(PsychError_user, "Invalid output index provided! No such output for this screen!");
        outputId = PsychScreenToHead(screenNumber, ((outputId < 0) ? 0 : outputId));
        if (outputId >= res->ncrtc || outputId < 0) PsychErrorExitMsg(PsychError_user, "Invalid output index provided! No such output for this screen!");

        RRCrtc crtc = res->crtcs[outputId];
        PsychLockDisplay();
        XRRCrtcGamma *lut = XRRGetCrtcGamma(cgDisplayID, crtc);
        PsychUnlockDisplay();

        n = (lut) ? lut->size : 0;

        if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: PsychReadNormalizedGammaTable: Provided RandR HW-LUT size is n=%i.\n", n);

        // Gamma lut query successfull?
        if (n > 0) {
            if ((n > MAX_GAMMALUT_SIZE) && (PsychPrefStateGet_Verbosity() > 1)) {
                printf("PTB-WARNING: ReadNormalizedGammatable: Hardware gamma table size of %i slots exceeds our maximum of %i slots. Clamping returned size to %i slots. Returned LUT's may be wrong!\n", n, MAX_GAMMALUT_SIZE, MAX_GAMMALUT_SIZE);
            }

            // Clamp for safety:
            n = (n <= MAX_GAMMALUT_SIZE) ? n : MAX_GAMMALUT_SIZE;

            // Convert tables: Map 16-bit values into 0-1 normalized floats:
            for (i = 0; i < n; i++) localRed[i]   = ((float) lut->red[i]) / 65535.0f;
            for (i = 0; i < n; i++) localGreen[i] = ((float) lut->green[i]) / 65535.0f;
            for (i = 0; i < n; i++) localBlue[i]  = ((float) lut->blue[i]) / 65535.0f;
        }

        if (lut) XRRFreeGamma(lut);
    }

    // Do we need the fallback path with XVidmodeExtension on systems which don't support RandR-1.2?
    // This applies to, e.g., the NVidia binary blob in the year 2011 (Sadly, not a joke):
    if (n <= 0) {
        // Use old-fashioned VidmodeExt path: No control over which output is queried on multi-display setups,
        // except in a ZaphodHead configuration where each display corresponds to a separate x-screen:
        #ifdef USE_VIDMODEEXTS
            // Query size of to-be-returned gamma table:
            PsychLockDisplay();
            XF86VidModeGetGammaRampSize(cgDisplayID, PsychGetXScreenIdForScreen(screenNumber), &n);
            PsychUnlockDisplay();

            if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: PsychReadNormalizedGammaTable: Provided XF86VidMode HW-LUT size is n=%i.\n", n);

            // Make sure we stay within our limits:
            if (n > MAX_GAMMALUT_SIZE) {
                printf("PTB-ERROR: ReadNormalizedGammatable: Hardware gamma table size of %i slots exceeds our maximum of %i slots! Gamma table query failed!\n", n, MAX_GAMMALUT_SIZE);
                PsychErrorExitMsg(PsychError_user, "Gamma table query failed due to size mismatch. Please report this on the Psychtoolbox forum.");
            }

            // Make sure this works at all:
            if (n <= 0) PsychErrorExitMsg(PsychError_user, "Gamma table query failed while trying XF86VidModeExtension fallback path.");

            // Retrieve gamma table with n slots:
            PsychLockDisplay();
            XF86VidModeGetGammaRamp(cgDisplayID, PsychGetXScreenIdForScreen(screenNumber), n, (unsigned short*) RTable, (unsigned short*) GTable, (unsigned short*) BTable);
            PsychUnlockDisplay();
        #else
            PsychErrorExitMsg(PsychError_user, "Sorry, this graphics card and driver does not support gamma table queries!");
        #endif

        // Convert tables: Map 16-bit values into 0-1 normalized floats:
        for (i = 0; i < n; i++) localRed[i]   = ((float) RTable[i]) / 65535.0f;
        for (i = 0; i < n; i++) localGreen[i] = ((float) GTable[i]) / 65535.0f;
        for (i = 0; i < n; i++) localBlue[i]  = ((float) BTable[i]) / 65535.0f;
    }

    // Assign output lut's:
    *redTable=localRed; *greenTable=localGreen; *blueTable=localBlue;

    // Assign size of LUT's::
    *numEntries = n;

    return;
}

/* Copy provided input LUT into target output LUT. If input is smaller than target LUT, but
 * fits nicely because output size is an integral multiple of input, then oversample input
 * to create proper output. If size doesn't match or output is smaller than input, abort
 * with error.
 *
 * Rationale: LUT's of standard GPUs are 256 slots, LUT's of high-end GPU's can be bigger
 *            powers-of-two sizes, e.g., 512, 1024, 2048, 4096 for some NVidia QuadroFX
 *            parts. For reasons of backwards compatibility we always want to be able to
 *            accept a good'ol 256 slots input LUT, even if the GPU expects something bigger.
 *            -> This simple oversampling via replication allows us to do that without obvious
 *               bad consequences for precision.
 *
 */
static void ConvertLUTToHwLUT(int nOut, psych_uint16* Rout, psych_uint16* Gout, psych_uint16* Bout, int nIn, float *redTable, float *greenTable, float *blueTable)
{
    int i, s;

    // Can't handle too big input tables for GPU:
    if (nOut < nIn) {
        printf("PTB-ERROR: Provided gamma table has %i slots, but graphics card accepts at most %i slots!\n", nIn, nOut);
        PsychErrorExitMsg(PsychError_user, "Provided gamma table has too many slots!");
    }

    // Can't handle tables which don't fit:
    if ((nOut % nIn) != 0) {
        printf("PTB-ERROR: Provided gamma table has %i slots, but graphics card expects %i slots.\n", nIn, nOut);
        printf("PTB-ERROR: Unfortunately, graphics card LUT size is not a integral multiple of provided gamma table size.\n");
        printf("PTB-ERROR: I can not handle this case, as it could cause ugly distortions in the displayed color range.\n");
        PsychErrorExitMsg(PsychError_user, "Provided gamma table has wrong number of slots!");
    }

    // Compute oversampling factor:
    s = nOut / nIn;
    if (PsychPrefStateGet_Verbosity() > 5) {
        printf("PTB-DEBUG: PsychLoadNormalizedGammaTable: LUT size nIn %i <= nOut %i --> Oversample %i times.\n", nIn, nOut, s);
    }

    // Copy and oversample: Each slot in red/green/blueTable is replicated
    // into s consecutive slots of R/G/Bout:
    for (i = 0; i < nOut; i++) {
        Rout[i] = (psych_uint16) (redTable[i/s]   * 65535.0f + 0.5f);
        Gout[i] = (psych_uint16) (greenTable[i/s] * 65535.0f + 0.5f);
        Bout[i] = (psych_uint16) (blueTable[i/s]  * 65535.0f + 0.5f);
    }
}

unsigned int PsychLoadNormalizedGammaTable(int screenNumber, int outputId, int numEntries, float *redTable, float *greenTable, float *blueTable)
{
    CGDirectDisplayID cgDisplayID;
    int i, j, n;
    RRCrtc crtc;

    static psych_uint16    RTable[MAX_GAMMALUT_SIZE];
    static psych_uint16    GTable[MAX_GAMMALUT_SIZE];
    static psych_uint16    BTable[MAX_GAMMALUT_SIZE];

    // The X-Windows hardware LUT has 3 tables for R,G,B.
    // Each entry is a 16-bit word with the n most significant bits used for an n-bit DAC or display encoder.

    // Set new gammaTable:
    PsychGetCGDisplayIDFromScreenNumber(&cgDisplayID, screenNumber);

    // Initial assumption: Failure.
    n = 0;

    // Not available on non-X11:
    if (!displayCGIDs[screenNumber]) return(0);

    if (has_xrandr_1_2) {
        // Use RandR V 1.2 for per-crtc setup:

        // Setup of all crtc's with this gamma table requested?
        if (outputId < 0) {
            // Yes: Iterate over all outputs, set via recursive call:
            j = 1;
            for (i = 0; (j > 0) && (i < kPsychMaxPossibleCrtcs) && (PsychScreenToHead(screenNumber, i) > -1); i++) {
                j = PsychLoadNormalizedGammaTable(screenNumber, i, numEntries, redTable, greenTable, blueTable);
            }

            // Done trying to set all crtc's. Return status:
            return((unsigned int) j);
        }

        // No, or recursive self-call: Load a specific crtc for output 'outputId':
        XRRScreenResources *res = displayX11ScreenResources[screenNumber];

        if (outputId >= kPsychMaxPossibleCrtcs) PsychErrorExitMsg(PsychError_user, "Invalid output index provided! No such output for this screen!");
        outputId = PsychScreenToHead(screenNumber, outputId);
        if (outputId >= res->ncrtc || outputId < 0) PsychErrorExitMsg(PsychError_user, "Invalid output index provided! No such output for this screen!");

        crtc = res->crtcs[outputId];

        // Get required size of gamma table:
        PsychLockDisplay();
        n = XRRGetCrtcGammaSize(cgDisplayID, crtc);
        PsychUnlockDisplay();
        if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: PsychLoadNormalizedGammaTable: Required RandR HW-LUT size is n=%i.\n", n);
    }

    // RandR 1.2 supported and has ability to set LUT's?
    if (has_xrandr_1_2 && (n > 0)) {
        // Allocate table of appropriate size:
        XRRCrtcGamma *lut = XRRAllocGamma(n);
        n = lut->size;

        // Convert tables: Map 16-bit values into 0-1 normalized floats:
        ConvertLUTToHwLUT(n, lut->red, lut->green, lut->blue, numEntries, redTable, greenTable, blueTable);

        // Assign to crtc:
        PsychLockDisplay();
        XRRSetCrtcGamma(cgDisplayID, crtc, lut);
        PsychUnlockDisplay();

        // Release lut:
        XRRFreeGamma(lut);
    }

    // RandR unsupported or failed?
    if (n <= 0) {
        // Use old-fashioned VidmodeExt fallback-path: No control over which output is setup on multi-display setups,
        // except in a ZaphodHead configuration where each display corresponds to a separate x-screen:
        if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: PsychLoadNormalizedGammaTable: Using XF86VidModeExt fallback path...\n");

        #ifdef USE_VIDMODEEXTS
            // Query size of to-be-set hw-gamma table:
            PsychLockDisplay();
            XF86VidModeGetGammaRampSize(cgDisplayID, PsychGetXScreenIdForScreen(screenNumber), &n);
            PsychUnlockDisplay();
            if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: PsychLoadNormalizedGammaTable: Required HW-LUT size is n=%i.\n", n);

            // Make sure we stay within our limits:
            if (n > MAX_GAMMALUT_SIZE) {
                printf("PTB-ERROR: LoadNormalizedGammatable: Hardware gamma table size of %i slots exceeds our maximum of %i slots! Gamma table setup failed!\n", n, MAX_GAMMALUT_SIZE);
                PsychErrorExitMsg(PsychError_user, "Gamma table setup failed due to size mismatch. Please report this on the Psychtoolbox forum.");
            }

            // Make sure this works at all:
            if (n <= 0) PsychErrorExitMsg(PsychError_user, "Gamma table setup failed while trying XF86VidModeExtension fallback path.");

            // Convert input table to X11 specific gammaTable:
            ConvertLUTToHwLUT(n, RTable, GTable, BTable, numEntries, redTable, greenTable, blueTable);

            // Assign to X-Screen:
            PsychLockDisplay();
            XF86VidModeSetGammaRamp(cgDisplayID, PsychGetXScreenIdForScreen(screenNumber), n, (unsigned short*) RTable, (unsigned short*) GTable, (unsigned short*) BTable);
            PsychUnlockDisplay();
        #else
            PsychErrorExitMsg(PsychError_user, "Sorry, this graphics card and driver does not support gamma table setup!");
        #endif
    }

    PsychLockDisplay();
    XFlush(cgDisplayID);
    PsychUnlockDisplay();

    // Return "success":
    return(1);
}

// Return true (non-zero) if a desktop compositor is likely active on screen 'screenNumber':
int PsychOSIsDWMEnabled(int screenNumber)
{
    int rc;
    CGDirectDisplayID dpy;
    PsychGetCGDisplayIDFromScreenNumber(&dpy, screenNumber);

    // Not available on non-X11: Assume no desktop compositor active.
    if (!dpy) return(0);

    // According to ICCCM spec, a compositing window manager who does composition on a
    // specific X-Screen must aquire "selection ownership" of the atom specified in our
    // displayX11ScreenCompositionAtom[screenNumber] for the given screenNumber. Therefore,
    // if the atom corresponding to 'screenNumber' does have a XGetSelectionOwner (!= None),
    // then that owner is the compositor, ergo desktop composition for that screenNumber is
    // active. Ref: http://standards.freedesktop.org/wm-spec/wm-spec-latest.html#id2579173
    //
    PsychLockDisplay();
    rc = (XGetSelectionOwner(dpy, displayX11ScreenCompositionAtom[screenNumber]) != None) ? 1 : 0;
    PsychUnlockDisplay();

    return(rc);
}

// !PTB_USE_WAYLAND
#endif

// PsychGetDisplayBeamPosition() contains the implementation of display beamposition queries.
// It requires both, a cgDisplayID handle, and a logical screenNumber and uses one of both for
// deciding which display pipe to query, whatever of both is more efficient or suitable for the
// host platform -- This is ugly, but neccessary, because the mapping with only one of these
// specifiers would be either ambigous (wrong results!) or usage would be inefficient and slow
// (bad for such a time critical low level call!). On some systems it may even ignore the arguments,
// because it's not capable of querying different pipes - ie., it will always query a hard-coded pipe.
//
int PsychGetDisplayBeamPosition(CGDirectDisplayID cgDisplayId, int screenNumber)
{
    // Beamposition queries aren't supported by the X11 graphics system.
    // However, for gfx-hardware where we have reliable register specs, we
    // can do it ourselves, bypassing the X server.

    // On systems that we can't handle, we return -1 as an indicator
    // to high-level routines that we don't know the rasterbeam position.
    double tdeadline, tnow;
    int vblbias, vbltotal;
    int beampos = -1;

    // Get beamposition from low-level driver code:
    if (PsychOSIsKernelDriverAvailable(screenNumber) && displayBeampositionHealthy[screenNumber]) {
        // Is application of the beamposition workaround requested by high-level code?
        // Or is this a NVidia GPU? In the latter case we always use the workaround,
        // because many NVidia GPU's (especially pre NV-50 hardware) need this in many
        // setups. It helps if needed, and doesn't hurt if not needed - burns at most
        // 25 insignificant microseconds of time.
        if ((PsychPrefStateGet_ConserveVRAM() & kPsychUseBeampositionQueryWorkaround) ||
            (fDeviceType == kPsychGeForce)) {
            // Yes: Avoid queries that return zero -- If query result is zero, retry
            // until it becomes non-zero: Some hardware may needs this to resolve...
            // We use a timeout of 100 msecs though to prevent hangs if we try to
            // query a disabled crtc's scanout position or similar bad things happen...
            PsychGetPrecisionTimerSeconds(&tdeadline);
            tdeadline += 0.1;
            while (0 == (beampos = PsychOSKDGetBeamposition(screenNumber))) {
                PsychGetPrecisionTimerSeconds(&tnow);
                if (tnow > tdeadline) {
                    // Trouble: Hanging here for more than 100 msecs?
                    // This display head is dead. Output a info to the user
                    // and disable it for further beamposition queries.
                    displayBeampositionHealthy[screenNumber] = FALSE;
                    beampos = -1;

                    if (PsychPrefStateGet_Verbosity() > 1) {
                        printf("PTB-WARNING: Hang in beamposition query detected! Seems my mapping of screen numbers to GPU's and display outputs is wrong?\n");
                        printf("PTB-WARNING: In a single GPU system you can resolve this by plugging in your monitors in a different order, changing the\n");
                        printf("PTB-WARNING: display arrangement in the control panel, or using the Screen('Preference', 'ScreenToHead', ...);\n");
                        printf("PTB-WARNING: command at the top of your scripts to set the mapping manually. See 'help DisplayOutputMappings' for more info.\n");
                        printf("PTB-WARNING: \n");
                        printf("PTB-WARNING: I am not yet able to handle multi-GPU systems reliably at all. If you have such a system it may work if\n");
                        printf("PTB-WARNING: you plug your monitor(s) into one of the other GPU's output connectors, trying different combinations.\n");
                        printf("PTB-WARNING: Or you simply live without high precision stimulus onset timestamping for now. Or you use the free and open-source\n");
                        printf("PTB-WARNING: graphics drivers (intel, radeon, or nouveau) instead of the proprietary Catalyst or NVidia binary drivers.\n");
                        printf("PTB-WARNING: I've disabled high precision timestamping for this screen for the remainder of the session.\n\n");
                        fflush(NULL);
                    }

                    break;
                }
            }
        } else {
            // Read final beampos:
            beampos = PsychOSKDGetBeamposition(screenNumber);
        }
    }

    // Return failure, if indicated:
    if (beampos == -1) return(-1);

    // Apply corrective offsets if any (i.e., if non-zero):
    // Note: In case of Radeon's, these are zero, because the code above already has applied proper corrections.
    PsychGetBeamposCorrection(screenNumber, &vblbias, &vbltotal);
    beampos = beampos - vblbias;
    if (beampos < 0) beampos = vbltotal + beampos;

    // Return our result or non-result:
    return(beampos);
}

// Try to attach to kernel level ptb support driver and setup everything, if it works:
void InitPsychtoolboxKernelDriverInterface(void)
{
    // This is currently a no-op on Linux, as most low-level stuff is done via mmapped() MMIO access...
    return;
}

// Try to detach to kernel level ptb support driver and tear down everything:
void PsychOSShutdownPsychtoolboxKernelDriverInterface(void)
{
    if (numKernelDrivers > 0) {
        // Make sure we're really detached from any MMIO memory before we shutdown:
        PsychScreenUnmapDeviceMemory();
    }

    // Ok, whatever happened, we're detached (for good or bad):
    numKernelDrivers = 0;

    return;
}

psych_bool PsychOSIsKernelDriverAvailable(int screenId)
{
    // Currently our "kernel driver" is available if MMIO mem could be mapped:
    // A real driver would indicate its presence via numKernelDrivers > 0 (see init/teardown code just above this routine):
    return((gfx_cntl_mem) ? TRUE : FALSE);
}

int PsychOSCheckKDAvailable(int screenId, unsigned int * status)
{
    // This doesn't make much sense on Linux yet. 'connect' should be something like a handle
    // to a kernel driver connection, e.g., the filedescriptor fd of the devicefile for ioctl()s
    // but we don't have such a thing yet.  Could be also a pointer to a little struct with all
    // relevant info...
    // Currently we do a dummy assignment...
    int connect = PsychScreenToCrtcId(screenId, 0);

    if ((numKernelDrivers<=0) && (gfx_cntl_mem == NULL)) {
        if (status) *status = ENODEV;
        return(0);
    }

    if (connect < 0) {
        if (status) *status = ENODEV;
        if (PsychPrefStateGet_Verbosity() > 6) printf("PTB-DEBUGINFO: Could not access kernel driver connection for screenId %i - No such connection.\n", screenId);
        return(0);
    }

    if (status) *status = 0;

    // Force this to '1', so the truth value is non-zero aka TRUE.
    connect = 1;
    return(connect);
}

unsigned int PsychOSKDReadRegister(int screenId, unsigned int offset, unsigned int* status)
{
    // Check availability of connection:
    int connect;
    if (!(connect = PsychOSCheckKDAvailable(screenId, status))) {
        if (status) *status = ENODEV;
        return(0xffffffff);
    }

    if (status) *status = 0;

    // Return readback register value:
    return(ReadRegister(offset));
}

unsigned int PsychOSKDWriteRegister(int screenId, unsigned int offset, unsigned int value, unsigned int* status)
{
    // Check availability of connection:
    int connect;
    if (!(connect = PsychOSCheckKDAvailable(screenId, status))) {
        if (status) *status = ENODEV;
        return(0xffffffff);
    }

    if (status) *status = 0;

    // Write the register:
    WriteRegister(offset, value);

    // Return success:
    return(0);
}

// Synchronize display screens video refresh cycle of DCE-4 (and later) GPU's, aka Evergreen. See PsychSynchronizeDisplayScreens() for help and details...
static PsychError PsychOSSynchronizeDisplayScreensDCE4(int *numScreens, int* screenIds, int* residuals, unsigned int syncMethod, double syncTimeOut, int allowedResidual)
{
    int             screenId = 0;
    double          abortTimeOut, now;
    int             residual;
    int             i, iter;
    unsigned int    old_crtc_master_enable = 0;

    // Check availability of connection:
    int             connect;
    unsigned int    status;

    // No support for other methods than fast hard sync:
    if (syncMethod > 1) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Could not execute display resync operation with requested non hard sync method. Not supported for this setup and settings.\n");
        return(PsychError_unimplemented);
    }

    // The current implementation only supports syncing the heads of a single card
    if (*numScreens <= 0) {
        // Resync all displays requested: Choose screenID zero for connect handle:
        screenId = 0;
    }
    else {
        // Resync of specific display requested: We only support resync of all heads of a single multi-head card,
        // therefore just choose the screenId of the passed master-screen for resync handle:
        screenId = screenIds[0];
    }

    if (!(connect = PsychOSCheckKDAvailable(screenId, &status))) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Could not execute display resync operation for master screenId %i. Not supported for this setup and settings.\n", screenId);
        return(PsychError_unimplemented);
    }

    if (fDeviceType != kPsychRadeon) {
        printf("PTB-INFO: PsychOSSynchronizeDisplayScreens(): This function is not supported on non-ATI/AMD GPU's! Aborted.\n");
        return(PsychError_unimplemented);
    }

    // Setup deadline for abortion or repeated retries:
    PsychGetAdjustedPrecisionTimerSeconds(&abortTimeOut);
    abortTimeOut+=syncTimeOut;
    residual = INT_MAX;

    // Repeat until timeout or good enough result:
    do {
        // If this isn't the first try, wait 0.5 secs before retry:
        if (residual != INT_MAX) PsychWaitIntervalSeconds(0.5);

        residual = INT_MAX;

        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: PsychOSSynchronizeDisplayScreens(): About to resynchronize all DCE-4+ display heads by use of a 1 second CRTC stop->start cycle:\n");

        if (PsychPrefStateGet_Verbosity() > 3) {
            printf("Trying to stop and reset all display heads by disabling them one by one.\n");
            printf("Will wait individually for each head to reach its defined resting position.\n");
        }

        // Detect enabled heads:
        old_crtc_master_enable = 0;
        for (iter = 0; iter < kPsychMaxPossibleCrtcs; iter++) {
            // Map 'iter'th head for this screenId to crtc index 'i'. Iterate over all crtc's for screen:
            if ((i = PsychScreenToCrtcId(screenId, iter)) < 0) break;

            // Sanity check crtc id i:
            if (i > (int) fNumDisplayHeads - 1) {
                printf("PTB-ERROR: PsychOSSynchronizeDisplayScreens(): Invalid headId %i provided! Must be between 0 and %i. Aborted.\n", i, (fNumDisplayHeads - 1));
                return(PsychError_user);
            }

            // Bit 16 "CRTC_CURRENT_MASTER_EN_STATE" allows read-only polling
            // of current activation state of crtc:
            if (ReadRegister(EVERGREEN_CRTC_CONTROL + crtcoff[i]) & (0x1 << 16)) old_crtc_master_enable |= (0x1 << i);
        }

        // Shut down heads, one after each other, wait for each one to settle at its defined resting position:
        for (iter = 0; iter < kPsychMaxPossibleCrtcs; iter++) {
            // Map 'iter'th head for this screenId to crtc index 'i'. Iterate over all crtc's for screen:
            if ((i = PsychScreenToCrtcId(screenId, iter)) < 0) break;

            if (PsychPrefStateGet_Verbosity() > 3) printf("Head %ld ...  ", i);
            if (old_crtc_master_enable & (0x1 << i)) {
                if (PsychPrefStateGet_Verbosity() > 3) printf("active -> Shutdown. ");

                // Shut down this CRTC by clearing its master enable bit (bit 0):
                WriteRegister(EVERGREEN_CRTC_CONTROL + crtcoff[i], ReadRegister(EVERGREEN_CRTC_CONTROL + crtcoff[i]) & ~(0x1 << 0));

                // Wait 50 msecs, so CRTC has enough time to settle and disable at its
                // programmed resting position:
                PsychWaitIntervalSeconds(0.050);

                // Double check - Poll until crtc is offline:
                while(ReadRegister(EVERGREEN_CRTC_CONTROL + crtcoff[i]) & (0x1 << 16));
                if (PsychPrefStateGet_Verbosity() > 3) printf("-> Offline.\n");
            }
            else {
                if (PsychPrefStateGet_Verbosity() > 3) printf("already offline.\n");
            }
        }

        // Need realtime priority for following synchronized start to minimize delays:
        PsychRealtimePriority(true);

        // Sleep for 1 second: This is a blocking call, ie. the thread goes to sleep and may wakeup a bit later:
        PsychWaitIntervalSeconds(1);

        // Reenable all now disabled, but previously enabled display heads.
        // This must be a tight loop, as every microsecond counts for a good sync...
        for (iter = 0; iter < kPsychMaxPossibleCrtcs; iter++) {
            // Map 'iter'th head for this screenId to crtc index 'i'. Iterate over all crtc's for screen:
            if ((i = PsychScreenToCrtcId(screenId, iter)) < 0) break;

            if (old_crtc_master_enable & (0x1 << i)) {
                // Restart this CRTC by setting its master enable bit (bit 0):
                WriteRegister(EVERGREEN_CRTC_CONTROL + crtcoff[i], ReadRegister(EVERGREEN_CRTC_CONTROL + crtcoff[i]) | (0x1 << 0));
            }
        }

        // Done with realtime bits:
        PsychRealtimePriority(false);

        // We don't have meaningful residual info. Just assume we succeeded:
        residual = 0;
        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Graphics display heads hopefully resynchronized.\n");

        // Timestamp:
        PsychGetAdjustedPrecisionTimerSeconds(&now);
    } while ((now < abortTimeOut) && (abs(residual) > allowedResidual));

    // Return residual value if wanted:
    if (residuals) {
        residuals[0] = residual;
    }

    if (abs(residual) > allowedResidual) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Failed to synchronize heads down to the allowable residual of +/- %i scanlines. Final residual %i lines.\n", allowedResidual, residual);
    }

    // TODO: Error handling not really worked out...
    if (residual == INT_MAX) return(PsychError_system);

    // Done.
    return(PsychError_none);
}

// Helper function for PsychOSSynchronizeDisplayScreens().
static unsigned int GetBeamPosition(int headId)
{
    return((unsigned int) ReadRegister((headId == 0) ? RADEON_D1CRTC_STATUS_POSITION : RADEON_D2CRTC_STATUS_POSITION) & RADEON_VBEAMPOSITION_BITMASK);
}

// Synchronize display screens video refresh cycle. See PsychSynchronizeDisplayScreens() for help and details...
PsychError PsychOSSynchronizeDisplayScreens(int *numScreens, int* screenIds, int* residuals, unsigned int syncMethod, double syncTimeOut, int allowedResidual)
{
    int             screenId = 0;
    double          abortTimeOut, now;
    int             residual;
    int             deltabeampos;
    unsigned int    beampos0, beampos1, i;
    unsigned int    old_crtc_master_enable = 0;
    unsigned int    new_crtc_master_enable = 0;

    // Check availability of connection:
    int             connect;
    unsigned int    status;

    // No support for other methods than fast hard sync:
    if (syncMethod > 1) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Could not execute display resync operation with requested non hard sync method. Not supported for this setup and settings.\n");
        return(PsychError_unimplemented);
    }

    // The current implementation only supports syncing all heads of a single card
    if (*numScreens <= 0) {
        // Resync all displays requested: Choose screenID zero for connect handle:
        screenId = 0;
    }
    else {
        // Resync of specific display requested: We only support resync of all heads of a single multi-head card,
        // therefore just choose the screenId of the passed master-screen for resync handle:
        screenId = screenIds[0];
    }

    if (!(connect = PsychOSCheckKDAvailable(screenId, &status))) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Could not execute display resync operation for master screenId %i. Not supported for this setup and settings.\n", screenId);
        return(PsychError_unimplemented);
    }

    if (fDeviceType != kPsychRadeon) {
        printf("PTB-INFO: PsychOSSynchronizeDisplayScreens(): This function is not supported on non-ATI/AMD GPU's! Aborted.\n");
        return(PsychError_unimplemented);
    }

    // DCE-4 display engine of Evergreen or later?
    if (isDCE4(screenId) || isDCE5(screenId) || isDCE10(screenId) || isDCE11(screenId)) {
        // Yes. Use DCE-4+ specific sync routine:
        return(PsychOSSynchronizeDisplayScreensDCE4(numScreens, screenIds, residuals, syncMethod, syncTimeOut, allowedResidual));
    }

    // Setup deadline for abortion or repeated retries:
    PsychGetAdjustedPrecisionTimerSeconds(&abortTimeOut);
    abortTimeOut+=syncTimeOut;
    residual = INT_MAX;

    // Repeat until timeout or good enough result:
    do {
        // If this isn't the first try, wait 0.5 secs before retry:
        if (residual != INT_MAX) PsychWaitIntervalSeconds(0.5);

        residual = INT_MAX;

        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: PsychOSSynchronizeDisplayScreens(): About to resynchronize all display heads by use of a 1 second CRTC stop->start cycle:\n");

        // A little pretest...
        if (PsychPrefStateGet_Verbosity() > 3) printf("Pretest...\n");
        for (i = 0; i < 10; i++) {
            beampos0 = GetBeamPosition(0);
            beampos1 = GetBeamPosition(1);
            if (PsychPrefStateGet_Verbosity() > 3) printf("Sample %ld: Beampositions are %ld vs. %ld . Offset %ld\n", i, beampos0, beampos1, (int) beampos1 - (int) beampos0);
        }

        // Query the CRTC scan-converter master enable state: Bit 0 (value 0x1) controls Pipeline 1,
        // whereas Bit 1(value 0x2) controls Pipeline 2:
        old_crtc_master_enable = ReadRegister(RADEON_DC_CRTC_MASTER_ENABLE);
        if (PsychPrefStateGet_Verbosity() > 3) {
            printf("Current CRTC Master enable state is %ld . Trying to stop and reset all display heads.\n", old_crtc_master_enable);
            printf("Will wait individually for each head to get close to scanline 0, then disable it.\n");
        }

        // Shut down heads, one after each other, each one at the start of a new refresh cycle:
        for (i = 0; i <= 1; i++) {
            // Wait for head i to start a new display cycle (scanline 0), then shut it down - well if it is active at all:
            if (PsychPrefStateGet_Verbosity() > 3) printf("Head %ld ...  ", i);
            if (old_crtc_master_enable & (0x1 << i)) {
                if (PsychPrefStateGet_Verbosity() > 3) printf("active -> Shutdown. ");
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
                if (PsychPrefStateGet_Verbosity() > 3) printf("New state is %ld.\n", ReadRegister(RADEON_DC_CRTC_MASTER_ENABLE));

                // Head should be down, close to scanline 0.
                PsychWaitIntervalSeconds(0.050);
            }
            else {
                if (PsychPrefStateGet_Verbosity() > 3) printf("already offline.\n");
            }
        }

        // All display heads should be disabled now.
        PsychWaitIntervalSeconds(0.100);

        // Query current beamposition and check state:
        beampos0 = GetBeamPosition(0);
        beampos1 = GetBeamPosition(1);

        new_crtc_master_enable = ReadRegister(RADEON_DC_CRTC_MASTER_ENABLE);

        if (new_crtc_master_enable == 0) {
            if (PsychPrefStateGet_Verbosity() > 3) printf("CRTC's down (state %ld): Beampositions are [0]=%ld and [1]=%ld. Synchronized restart in 1 second...\n", new_crtc_master_enable, beampos0, beampos1);
        }
        else {
            if (PsychPrefStateGet_Verbosity() > 3) printf("CRTC's shutdown failed!! (state %ld): Beamposition are [0]=%ld and [1]=%ld. Will try to restart in 1 second...\n", new_crtc_master_enable, beampos0, beampos1);
        }

        // Sleep for 1 second: This is a blocking call, ie. the thread goes to sleep and may wakeup a bit later:
        PsychWaitIntervalSeconds(1);

        // Reset all display heads enable state to original setting:
        WriteRegister(RADEON_DC_CRTC_MASTER_ENABLE, old_crtc_master_enable);

        // Query position and state after restart:
        beampos0 = GetBeamPosition(0);
        beampos1 = GetBeamPosition(1);
        new_crtc_master_enable = ReadRegister(RADEON_DC_CRTC_MASTER_ENABLE);
        if (new_crtc_master_enable == old_crtc_master_enable) {
            if (PsychPrefStateGet_Verbosity() > 3) printf("CRTC's restarted in sync: Master enable state is %ld. Beampositions after restart: [0]=%ld and [1]=%ld.\n", new_crtc_master_enable, beampos0, beampos1);
        }
        else {
            if (PsychPrefStateGet_Verbosity() > 3) printf("CRTC's restart FAILED!!: Master enable state is %ld. Beampositions: [0]=%ld and [1]=%ld.\n", new_crtc_master_enable, beampos0, beampos1);
        }

        deltabeampos = (int) beampos1 - (int) beampos0;
        if (PsychPrefStateGet_Verbosity() > 3) printf("Residual beam offset after display sync: %ld.\n\n", deltabeampos);

        // A little posttest...
        if (PsychPrefStateGet_Verbosity() > 3) printf("Posttest...\n");
        for (i = 0; i < 10; i++) {
            beampos0 = GetBeamPosition(0);
            beampos1 = GetBeamPosition(1);
            if (PsychPrefStateGet_Verbosity() > 3) printf("Sample %ld: Beampositions are %ld vs. %ld . Offset %ld\n", i, beampos0, beampos1, (int) beampos1 - (int) beampos0);
        }

        if (PsychPrefStateGet_Verbosity() > 3) printf("Display head resync operation finished.\n\n");

        // Assign residual for this iteration:
        residual = deltabeampos;

        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Graphics display heads resynchronized. Residual vertical beamposition error is %ld scanlines.\n", residual);

        // Timestamp:
        PsychGetAdjustedPrecisionTimerSeconds(&now);
    } while ((now < abortTimeOut) && (abs(residual) > allowedResidual));

    // Return residual value if wanted:
    if (residuals) {
        residuals[0] = residual;
    }

    if (abs(residual) > allowedResidual) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: Failed to synchronize heads down to the allowable residual of +/- %i scanlines. Final residual %i lines.\n", allowedResidual, residual);
    }

    // TODO: Error handling not really worked out...
    if (residual == INT_MAX) return(PsychError_system);

    // Done.
    return(PsychError_none);
}

int PsychOSKDGetBeamposition(int screenId)
{
    int beampos = -1;
    int headId  = PsychScreenToCrtcId(screenId, 0);
    static psych_bool firstTime = TRUE;

    if (headId < 0 || headId > ((int) fNumDisplayHeads - 1)) {
        printf("PTB-ERROR: PsychOSKDGetBeamposition: Invalid headId %i provided! Must be between 0 and %i. Aborted.\n", headId, (fNumDisplayHeads - 1));
        return(beampos);
    }

    // MMIO registers mapped?
    if (gfx_cntl_mem) {
        // Query code for ATI/AMD Radeon/FireGL/FirePro:
        if (fDeviceType == kPsychRadeon) {
            if (isDCE4(screenId) || isDCE5(screenId) || isDCE10(screenId) || isDCE11(screenId)) {
                // DCE-4+ display engine (CEDAR and later afaik): Up to six or seven crtc's.

                // Read raw beampostion from GPU:
                beampos = (int) (ReadRegister(EVERGREEN_CRTC_STATUS_POSITION + crtcoff[headId]) & RADEON_VBEAMPOSITION_BITMASK);

                // Query end-offset of VBLANK interval of this GPU and correct for it:
                beampos = beampos - (int) ((ReadRegister(EVERGREEN_CRTC_V_BLANK_START_END + crtcoff[headId]) >> 16) & RADEON_VBEAMPOSITION_BITMASK);

                // Correction for in-VBLANK range:
                if (beampos < 0) beampos = ((int) ReadRegister(EVERGREEN_CRTC_V_TOTAL + crtcoff[headId])) + beampos;

            } else {
                // AVIVO / DCE-2 / DCE-3 display engine (R300 - R700 afaik): At most two display heads for dual-head gpu's.

                // Read raw beampostion from GPU:
                beampos = (int) (ReadRegister((headId == 0) ? RADEON_D1CRTC_STATUS_POSITION : RADEON_D2CRTC_STATUS_POSITION) & RADEON_VBEAMPOSITION_BITMASK);

                // Query end-offset of VBLANK interval of this GPU and correct for it:
                beampos = beampos - (int) ((ReadRegister((headId == 0) ? AVIVO_D1CRTC_V_BLANK_START_END : AVIVO_D2CRTC_V_BLANK_START_END) >> 16) & RADEON_VBEAMPOSITION_BITMASK);

                // Correction for in-VBLANK range:
                if (beampos < 0) beampos = ((int) ReadRegister((headId == 0) ? AVIVO_D1CRTC_V_TOTAL : AVIVO_D2CRTC_V_TOTAL)) + beampos;
            }
        }

        // Query code for NVidia GeForce/Quadro:
        if (fDeviceType == kPsychGeForce) {
            // Pre NV-50 GPU? [Anything before GeForce-8 series]
            if ((fCardType > 0x0) && (fCardType < 0x50)) {
                // Pre NV-50, e.g., RivaTNT-1/2 and all GeForce 256/2/3/4/5/FX/6/7:

                // Lower 12 bits are vertical scanout position (scanline), bit 16 is "in vblank" indicator.
                // Offset between crtc's is 0x2000, we're only interested in scanline, not "in vblank" status:
                // beampos = (int) (ReadRegister((headId == 0) ? 0x600808 : 0x600808 + 0x2000) & 0xFFF);

                // NV-47: Lower 16 bits are vertical scanout position (scanline), upper 16 bits are horizontal
                // scanout position. Offset between crtc's is 0x2000. We only use the lower 16 bits and
                // ignore horizontal scanout position for now:
                beampos = (int) (ReadRegister((headId == 0) ? 0x600868 : 0x600868 + 0x2000) & 0xFFFF);
            } else {
                // NV-50 (GeForce-8) and later, also 4-CRTC NV-E0 and later:

                // Lower 16 bits are vertical scanout position (scanline), upper 16 bits are vblank counter.
                // Offset between crtc's is 0x800, we're only interested in scanline, not vblank counter:
                beampos = (int) (ReadRegister(0x616340 + (0x800 * headId)) & 0xFFFF);
                if (PsychPrefStateGet_Verbosity() > 11) printf("PTB-DEBUG: Head %i, HW-Vblankcount: %i\n", headId, (int) ((ReadRegister(0x616340 + (0x800 * headId)) >> 16) & 0xFFFF));

                if (FALSE && firstTime) {
                    firstTime = FALSE;
                    int newcount = -1;
                    int oldcount = (int) ((ReadRegister(0x616340 + (0x800 * headId)) >> 16) & 0xFFFF);
                    unsigned int foo;
                    while (newcount <= oldcount) {
                        foo = ReadRegister(0x616340 + (0x800 * headId));
                        newcount = (int) ((foo >> 16) & 0xFFFF);
                        beampos = (int) (foo & 0xFFFF);
                    }
                    printf("VBLIncrement %i -> %i at scanline %i\n", oldcount, newcount, beampos);
                }
            }
        }

        // Query code for Intel IGP's:
        if (fDeviceType == kPsychIntelIGP) {
            beampos = (int) (ReadRegister(0x70000 + (headId * 0x1000)) & 0x1FFF);
        }

        // Safety measure: Cap to zero if something went wrong -> This will trigger proper high level error handling in PTB:
        if (beampos < 0) beampos = -1;
    }

    return(beampos);
}

// Try to change hardware dither mode on GPU:
void PsychOSKDSetDitherMode(int screenId, unsigned int ditherOn)
{
    static unsigned int oldDither[kPsychMaxPossibleCrtcs] = { 0, 0, 0, 0, 0, 0 };
    unsigned int reg;
    int headId, iter;

    // MMIO registers mapped?
    if (!gfx_cntl_mem) return;

    // Check if the method is supported for this GPU type:
    // Currently ATI/AMD GPU's only...
    if (fDeviceType != kPsychRadeon) {
        // Other unsupported GPU:
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: SetDitherMode: Tried to call me on a non ATI/AMD GPU. Unsupported.\n");
        return;
    }

    // Start with headId undefined:
    headId = -1;

    for (iter = 0; iter < kPsychMaxPossibleCrtcs; iter++) {
        if (screenId >= 0) {
            // Positive screenId: Apply to all crtc's for this screenId:

            // Is there an iter'th crtc assigned to this screen?
            headId = PsychScreenToCrtcId(screenId, iter);

            // If end of list of associated crtc's for this screenId reached, then we're done:
            if (headId < 0) break;
        }
        else {
            // Negative screenId -> Only affect one head defined by screenId:
            if (headId < 0) {
                // Setup single target head in this iteration:
                headId = -screenId;
            }
            else {
                // Single target head already set up: We're done:
                break;
            }
        }

        // AMD/ATI Radeon, FireGL or FirePro GPU?
        if (fDeviceType == kPsychRadeon) {
            if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: SetDitherMode: Trying to %s digital display dithering on display head %d.\n", (ditherOn) ? "enable" : "disable", headId);

            // Map headId to proper hardware control register offset:
            if (isDCE4(screenId) || isDCE5(screenId) || isDCE6(screenId) || isDCE8(screenId) || isDCE10(screenId) || isDCE11(screenId)) {
                // DCE-4+ display engine (CEDAR and later afaik): Map to proper register offset for this headId:
                if (headId > ((int) fNumDisplayHeads - 1)) {
                    // Invalid head - bail:
                    if (PsychPrefStateGet_Verbosity() > 0) printf("SetDitherMode: ERROR! Invalid headId %d provided. Must be between 0 and %i. Aborted.\n", headId, (fNumDisplayHeads - 1));
                    continue;
                }

                // Map to dither format control register for head 'headId':
                reg = EVERGREEN_FMT_BIT_DEPTH_CONTROL + crtcoff[headId];
            } else if (isDCE3(screenId)) {
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
                if ((ReadRegister(0x7880) & 0x1) && ((ReadRegister(0x7884) & 0x1) == (unsigned int) headId)) reg = RADEON_TMDSA_BIT_DEPTH_CONTROL;

                // LVTMA block enabled, and driven by headId? Then we control its encoder:
                if ((ReadRegister(0x7A80) & 0x1) && ((ReadRegister(0x7A84) & 0x1) == (unsigned int) headId)) reg = RADEON_LVTMA_BIT_DEPTH_CONTROL;

                // DVOA block enabled, and driven by headId? Then we control its encoder:
                if ((ReadRegister(0x7980) & 0x1) && ((ReadRegister(0x7984) & 0x1) == (unsigned int) headId)) reg = RADEON_DVOA_BIT_DEPTH_CONTROL;

                // If no digital encoder block was assigned, then this likely means we're connected to a
                // analog VGA monitor driven by the DAC. The DAC doesn't have dithering ever, so we are
                // done with a simple no-op:
                if (reg == 0x0) {
                    if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: SetDitherMode: Screen %i, head %i connected to analog VGA DAC. Dithering control skipped.\n", screenId, headId);
                    continue;
                }
                else if (PsychPrefStateGet_Verbosity() > 3) {
                    switch (reg) {
                        case RADEON_TMDSA_BIT_DEPTH_CONTROL:
                            printf("PTB-INFO: SetDitherMode: Screen %i, head %i connected to TMDSA block.\n", screenId, headId);
                        break;

                        case RADEON_LVTMA_BIT_DEPTH_CONTROL:
                            printf("PTB-INFO: SetDitherMode: Screen %i, head %i connected to LVTMA block.\n", screenId, headId);
                        break;

                        case RADEON_DVOA_BIT_DEPTH_CONTROL:
                            printf("PTB-INFO: SetDitherMode: Screen %i, head %i connected to DVOA block.\n", screenId, headId);
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
                        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: SetDitherMode: Dithering previously disabled by us. Reenabling with old control setting %x.\n", oldDither[headId]);
                        WriteRegister(reg, oldDither[headId]);
                    }
                    else {
                        // No: Dithering was disabled all the time, so we don't know the
                        // OS defaults. Use the numeric value of 'ditherOn' itself:
                        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: SetDitherMode: Dithering off. Enabling with userspace provided setting %x. Cross your fingers!\n", ditherOn);
                        WriteRegister(reg, ditherOn);
                    }
                }
                else {
                    // Dithering currently on.

                    // Specific value for control reg specified?
                    if (ditherOn > 1) {
                        // Yes. Use it "as is":
                        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: SetDitherMode: Setting dithering mode to userspace provided setting %x. Cross your fingers!\n", ditherOn);
                        WriteRegister(reg, ditherOn);
                    }
                    else if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: SetDitherMode: Dithering already enabled with current control value %x. Skipped.\n", ReadRegister(reg));
                }
            }
            else {
                // Disable all dithering if it is enabled: Clearing the register to all zero bits does this.
                if (ReadRegister(reg) > 0) {
                    oldDither[headId] = ReadRegister(reg);
                    if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: SetDitherMode: Current dither setting before our dither disable on head %d is %x. Disabling.\n", headId, oldDither[headId]);
                    WriteRegister(reg, 0x0);
                }
                else {
                    if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: SetDitherMode: Dithering already disabled. Skipped.\n");
                }
            }
            // End of Radeon et al. support code.
        }
        // Next head for this screenId, if any...
    }

    return;
}

// Query if LUT for given headId is all-zero: 0 = Something else, 1 = Zero-LUT, 2 = It's an identity LUT,
// 3 = Not-quite-identity mapping, 0xffffffff = don't know.
unsigned int PsychOSKDGetLUTState(int screenId, unsigned int headId, unsigned int debug)
{
    unsigned int i, v, r, m, bo, wo, offset, reg;
    unsigned int isZero = 1;
    unsigned int isIdentity = 1;

    // AMD GPU's:
    if (fDeviceType == kPsychRadeon) {
        if (PsychPrefStateGet_Verbosity() > 3) printf("PsychOSKDGetLUTState(): Checking LUT and bias values on GPU for headId %d.\n", headId);

        if (isDCE4(screenId) || isDCE5(screenId) || isDCE10(screenId) || isDCE11(screenId)) {
            // DCE-4.0 and later: Up to (so far) six display heads:
            if (headId > (fNumDisplayHeads - 1)) {
                // Invalid head - bail:
                if (PsychPrefStateGet_Verbosity() > 2) printf("PsychOSKDGetLUTState: ERROR! Invalid headId %d provided. Must be between 0 and %i. Aborted.\n", headId, (fNumDisplayHeads - 1));
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

        if (debug) if (PsychPrefStateGet_Verbosity() > 3) printf("PsychOSKDOffsets: Black %d : White %d.\n", bo, wo);

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

            if (PsychPrefStateGet_Verbosity() > 4) {
                printf("%d:%d,%d,%d\n", i, (v >> 20) & 0x3ff, (v >> 10) & 0x3ff, (v >> 0) & 0x3ff);
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
    if (PsychPrefStateGet_Verbosity() > 3) printf("PsychOSKDGetLUTState(): This function is not supported on this GPU. Returning 0xffffffff.\n");
    return(0xffffffff);
}

// Load an identity LUT into display head 'headid': Return 1 on success, 0 on failure or if unsupported for this GPU:
unsigned int PsychOSKDLoadIdentityLUT(int screenId, unsigned int headId)
{
    unsigned int i, r, m, offset, reg;

    // AMD GPU's:
    if (fDeviceType == kPsychRadeon) {
        if (PsychPrefStateGet_Verbosity() > 3) printf("PsychOSKDLoadIdentityLUT(): Uploading identity LUT and bias values into GPU for headId %d.\n", headId);

        if (isDCE4(screenId) || isDCE5(screenId) || isDCE10(screenId) || isDCE11(screenId)) {
            // DCE-4.0+ and later: Up to (so far) six display heads:
            if (headId > (fNumDisplayHeads - 1)) {
                // Invalid head - bail:
                if (PsychPrefStateGet_Verbosity() > 3) printf("PsychOSKDLoadIdentityLUT: ERROR! Invalid headId %d provided. Must be between 0 and %i. Aborted.\n", headId, (fNumDisplayHeads - 1));
                return(0);
            }

            offset = crtcoff[headId];
            reg = EVERGREEN_DC_LUT_30_COLOR + offset;

            WriteRegister(EVERGREEN_DC_LUT_CONTROL + offset, 0);

            if (isDCE5(screenId) || isDCE10(screenId)) {
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

        if (isDCE5(screenId) || isDCE10(screenId)) {
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
    if (PsychPrefStateGet_Verbosity() > 3) printf("PsychOSKDLoadIdentityLUT(): This function is not supported on this GPU. Returning 0.\n");
    return(0);
}
