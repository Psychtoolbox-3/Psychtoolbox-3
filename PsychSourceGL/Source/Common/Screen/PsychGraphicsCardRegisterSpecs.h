/*
    PsychToolbox3/Source/Common/Screen/PsychGraphicsCardRegisterSpecs.h

    AUTHORS:

        mario.kleiner.de@gmail.com      mk
        Alex Deucher (indirectly - by providing useful pointers and hints)

    PLATFORMS:

        All.

    HISTORY:

    01/12/2008    mk        Created.

    DESCRIPTION:

    This file contains specifications of the low-level registers of different graphics
    cards.

    Currently it contains register offsets for recent ATI GPU's and some NVIDIA GPU's.

    Copyright:

    Parts of this file are direct copies of register defines from the Linux Radeon KMS
    driver. The file is copyright 2008-2011 Mario Kleiner and also
    Copyright 2010 Advanced Micro Devices, Inc.
    Authors: Alex Deucher

    License:

    *
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
    *

*/

//include once
#ifndef PSYCH_IS_INCLUDED_PsychGraphicsCardRegisterSpecs
#define PSYCH_IS_INCLUDED_PsychGraphicsCardRegisterSpecs

// PCI vendor ids:
#define PCI_VENDOR_ID_NVIDIA    0x10de
#define PCI_VENDOR_ID_AMD       0x1022
#define PCI_VENDOR_ID_ATI       0x1002
#define PCI_VENDOR_ID_INTEL     0x8086

// Settings for member fDeviceType: Decoded Vendor id's.
#define kPsychUnknown  0
#define kPsychGeForce  1
#define kPsychRadeon   2
#define kPsychIntelIGP 3

// The following register offsets and specs are taken from the official AMD/ATI
// specs, downloadable from http://www.x.org/docs/AMD/ the website of the X-ORG
// foundation, home of the open source Linux X11 server. These specs are official
// specs from AMD/ATI, released in autumn 2007 to the open-source community.

// R500 master interrupt control and status registers:
#define RADEON_R500_GEN_INT_CNTL   0x040
#define RADEON_R500_GEN_INT_STATUS 0x044

// Display controller interrupt occured? [Selection mask]:
#define RADEON_R500_DISPLAY_INT_STATUS    (1 << 0)

// Interrupt status register for reporting of display controller interrupts:
#define RADEON_R500_DISP_INTERRUPT_STATUS 0x7edc

// Interrupt mask for reporting display head 1 or 2 VBLANK interrupts:
#define RADEON_R500_D1_VBLANK_INTERRUPT (1 << 4)
#define RADEON_R500_D2_VBLANK_INTERRUPT (1 << 5)

// Acknowledge target registers to acknowledge display controller interrupts for...
// Display head 1
#define RADEON_R500_D1MODE_VBLANK_STATUS 0x6534
// Display head 2
#define RADEON_R500_D2MODE_VBLANK_STATUS 0x6d34
// Acknowledge bit for acknowledging display VBLANK interrupts:
#define RADEON_R500_VBLANK_ACK     (1<<4)

// The D1CRTC_STATUS_POSITION register (32 bits) encodes vertical beam position in
// bits 0:12 (the least significant 13 bits), and horizontal beam position in
// bits 16-28. D2 is the secondary display pipeline (e.g., the external video port
// on a laptop), whereas D1 is the primary pipeline (e.g., internal panel of a laptop).
// The addresses and encoding is the same for the X1000 series and the newer HD2000/3000
// series chips...
#define RADEON_D1CRTC_STATUS_POSITION 0x60a0
#define RADEON_D2CRTC_STATUS_POSITION 0x68a0
#define RADEON_VBEAMPOSITION_BITMASK  0x1fff
#define RADEON_HBEAMPOSITION_BITSHIFT 16

// Encodes true start- and endline of VBLANK interval, because vblank doesn't
// neccessarily start at vactive and doesn't stop at vtotal:
// Uppermost 16 bits contain end, lowermost 16 bits contain start, mask is
// RADEON_VBEAMPOSITION_BITMASK
#define AVIVO_D1CRTC_V_BLANK_START_END 0x6024
#define AVIVO_D2CRTC_V_BLANK_START_END 0x6824

#define AVIVO_D1CRTC_V_TOTAL           0x6020
#define AVIVO_D2CRTC_V_TOTAL           0x6820

// This (if we would use it) would give access to on-chip frame counters. These increment
// once per video refresh cycle - at the beginning of a new cycle (scanline zero) and
// can be software reset, but normally start at system bootup with zero. Not yet sure
// if we should/would wanna use 'em but we'll see...
#define RADEON_D1CRTC_STATUS_FRAME_COUNT 0x60a4
#define RADEON_D2CRTC_STATUS_FRAME_COUNT 0x68a4

// CRTC display scan converter master enable register: Bit 0 enables/disables CRTC1 head,
// whereas Bit 1 enables/disables CRTC2 head. Writing a zero into the reg. will disable and
// blank both displays, Writing a 0x1 | 0x2 == 0x3 will simulataneously enable both heads.
// --> As writing one register is atomic and both heads reset their beamposition on restart,
// this is a simple way to instantaneously and perfectly sync the two display heads.
#define RADEON_DC_CRTC_MASTER_ENABLE 0x60f8

// Registers which point to the start addresses of the onscreen framebuffers - the
// "OpenGL frontbuffers".
// In mono mode, primary surface is the frontbuffer memory pointer, in frame-sequential
// stereo mode, primary surface is the front-left buffer, secondary surface is the front
// right buffer. They don't change in OS/X windowed mode, ie., when the normal user interface
// is active, as the Quartz compositor seems to use fast VRAM->VRAM blits to recomposit
// the user interface, not double-buffering via page-flipping. In fullscreen mode
// (fullscreen onscreen windows open), they change at each Screen('Flip') execution, ie.
// at each completed double-buffer swap. They toggle between the VRAM addresses of the two
// (or four in stereo) buffers when they are exchanged:
#define RADEON_D1GRPH_PRIMARY_SURFACE_ADDRESS    0x6110
#define RADEON_D1GRPH_SECONDARY_SURFACE_ADDRESS  0x6118
#define RADEON_D2GRPH_PRIMARY_SURFACE_ADDRESS    0x6910
#define RADEON_D2GRPH_SECONDARY_SURFACE_ADDRESS  0x6918

#define RADEON_D1GRPH_UPDATE                    0x6144
#define RADEON_D2GRPH_UPDATE                    0x6944
#define RADEON_SURFACE_UPDATE_PENDING           4
#define RADEON_SURFACE_UPDATE_TAKEN             8

// Bit 0 of these enables (Setting of 0) or disables (Setting of 1 aka 0x1) internal routing
// the VSYNC signal (aka V_UPDATE), so this can prevent/hold back the double-buffered
// update of CRTC registers, esp. upate of the primary/secondary surface address.
// TODO TEST: Can this be used to delay flips, e.g., until some deadline is reached
// or some user defined trigger condition occurs, e.g., kernel-level sync swaps to
// trigger condition?
#define RADEON_D1MODE_MASTER_UPDATE_LOCK        0x60E0
#define RADEON_D2MODE_MASTER_UPDATE_LOCK        0x68E0

// Encodes pitch (length of a single scanline in pixels) as stored in the framebuffer:
#define RADEON_D1GRPH_PITCH 0x6120
#define RADEON_D2GRPH_PITCH 0x6920

// These registers control operation of the color/gamma hardware lookup tables:
// Bit 0 selects if hardware LUT A or B should be used for the given display pipeline.
// Bit 8 controls if the LUT's should be bypassed in ABGR2101010 operation mode, ie.,
// if a 10 bit per color channel framebuffer is used, should the 8->10 LUT's be ignored
// and the 10 bit values directly output to the RAMDAC's/TDMS/LVDS emitters? Usually
// you want to enable this bypass in 10 bpc mode for best results:
#define RADEON_D1GRPH_LUT_SEL    0x6108
#define RADEON_D2GRPH_LUT_SEL    0x6908

// These registers control the color depth and pixel formatting in the framebuffer -
// or better: The interpretation of the data by the CRTC's. They don't affect the
// rasterizer backends - the GPU will always draw ABGR8888 formatted data in 32 bit
// color mode. Bits 8-10 control data format. Especially interesting: If you choose
// 32 bit color depths (Bits 1:0 == 0x2) and you set Bits 10:8 to 0x1,ie. set bit 8,
// then the 32 bit pixel data will be interpreted as ABGR2101010 data, ie. 2 bits
// alpha, 10 bits for red, green and blue. Enabling LUT bypass (see above) will give
// you a real linear 10 bpc framebuffer from input to output! However, as the GPU
// always renders in 8888 format, you'll experience corrupted screen colors, unless
// you use the PTB imaging pipeline for automatic reformatting of floating point
// HDR images to RGBA8 formatted data, suitable for 2101010 mode. :-)
#define RADEON_D1GRPH_CONTROL    0x6104
#define RADEON_D2GRPH_CONTROL    0x6904

// These registers control output bit depths over TMDS, LVTM, DVO emitters,
// as well as the method to reduce the bit depths in a visually pleasing way
// if reduction is enabled.
// Writing all zero (0x0) will disable reduction, ie., output a 24bpp, 8bpc
// stream. Available methods of reduction are:
// a) Truncation: Cutting off least significant bits (Enable by setting bit 0),
// b) Spatial dithering: Set bit 8.
// c) Temporal dithering/modulation: Set bit 16.
// Bit 4 for a), 12 for b) and 20 for c) Control target bpp:
// either 18 bpp (bits cleared) or 24 bpp (bits sets).
// There are more higher order bits for controlling exact operation of temporal
// dithering...
#define RADEON_LVTMA_BIT_DEPTH_CONTROL  0x7A94
#define RADEON_TMDSA_BIT_DEPTH_CONTROL  0x7894
#define RADEON_DVOA_BIT_DEPTH_CONTROL   0x7988

#define AVIVO_DC_LUTA_CONTROL                   0x64C0
#define AVIVO_DC_LUTA_BLACK_OFFSET_BLUE         0x64c4
#define AVIVO_DC_LUTA_BLACK_OFFSET_GREEN        0x64c8
#define AVIVO_DC_LUTA_BLACK_OFFSET_RED          0x64cc
#define AVIVO_DC_LUTA_WHITE_OFFSET_BLUE         0x64d0
#define AVIVO_DC_LUTA_WHITE_OFFSET_GREEN        0x64d4
#define AVIVO_DC_LUTA_WHITE_OFFSET_RED          0x64d8

#define AVIVO_DC_LUT_RW_SELECT                  0x6480
#define AVIVO_DC_LUT_RW_MODE                    0x6484
#define AVIVO_DC_LUT_RW_INDEX                   0x6488
#define AVIVO_DC_LUT_SEQ_COLOR                  0x648c
#define AVIVO_DC_LUT_PWL_DATA                   0x6490
#define AVIVO_DC_LUT_30_COLOR                   0x6494
#define AVIVO_DC_LUT_READ_PIPE_SELECT           0x6498
#define AVIVO_DC_LUT_WRITE_EN_MASK              0x649c
#define AVIVO_DC_LUT_AUTOFILL                   0x64a0

// DCE-3 display engine (R700, roughly HD-4000 series):
#define DCE3_FMT_BIT_DEPTH_CONTROL              0x6710
#define DCE3_FMT_CLAMP_CONTROL                  0x672c

// Evergreen class hardware (DCE-4 display engine):
#define EVERGREEN_DC_LUT_RW_MODE                        0x69e0
#define EVERGREEN_DC_LUT_RW_INDEX                       0x69e4
#define EVERGREEN_DC_LUT_SEQ_COLOR                      0x69e8
#define EVERGREEN_DC_LUT_PWL_DATA                       0x69ec
#define EVERGREEN_DC_LUT_30_COLOR                       0x69f0
#define EVERGREEN_DC_LUT_VGA_ACCESS_ENABLE              0x69f4
#define EVERGREEN_DC_LUT_WRITE_EN_MASK                  0x69f8
#define EVERGREEN_DC_LUT_AUTOFILL                       0x69fc
#define EVERGREEN_DC_LUT_CONTROL                        0x6a00
#define EVERGREEN_DC_LUT_BLACK_OFFSET_BLUE              0x6a04
#define EVERGREEN_DC_LUT_BLACK_OFFSET_GREEN             0x6a08
#define EVERGREEN_DC_LUT_BLACK_OFFSET_RED               0x6a0c
#define EVERGREEN_DC_LUT_WHITE_OFFSET_BLUE              0x6a10
#define EVERGREEN_DC_LUT_WHITE_OFFSET_GREEN             0x6a14
#define EVERGREEN_DC_LUT_WHITE_OFFSET_RED               0x6a18
#define EVERGREEN_DC_LUT_10BIT_BYPASS                   0x6808
#define EVERGREEN_GRPH_CONTROL                          0x6804
#define EVERGREEN_GRPH_ENABLE                           0x6800

/* display controller offsets used for crtc/cur/lut/grph/viewport/etc. */
#define EVERGREEN_CRTC0_REGISTER_OFFSET                 (0x6df0 - 0x6df0)
#define EVERGREEN_CRTC1_REGISTER_OFFSET                 (0x79f0 - 0x6df0)
#define EVERGREEN_CRTC2_REGISTER_OFFSET                 (0x105f0 - 0x6df0)
#define EVERGREEN_CRTC3_REGISTER_OFFSET                 (0x111f0 - 0x6df0)
#define EVERGREEN_CRTC4_REGISTER_OFFSET                 (0x11df0 - 0x6df0)
#define EVERGREEN_CRTC5_REGISTER_OFFSET                 (0x129f0 - 0x6df0)

/* CRTC blocks at 0x6df0, 0x79f0, 0x105f0, 0x111f0, 0x11df0, 0x129f0 */
#define EVERGREEN_CRTC_V_TOTAL                          0x6e1c
#define EVERGREEN_CRTC_V_BLANK_START_END                0x6e34
#define EVERGREEN_CRTC_CONTROL                          0x6e70
#       define EVERGREEN_CRTC_MASTER_EN                 (1 << 0)
#       define EVERGREEN_CRTC_DISP_READ_REQUEST_DISABLE (1 << 24)
#define EVERGREEN_CRTC_BLANK_CONTROL                    0x6e74
#       define EVERGREEN_CRTC_BLANK_DATA_EN             (1 << 8)
#define EVERGREEN_CRTC_STATUS                           0x6e8c
#define EVERGREEN_CRTC_STATUS_POSITION                  0x6e90
#define EVERGREEN_CRTC_UPDATE_LOCK                      0x6ed4

// Primary/Secondary scanout buffer start addresses - low 32-Bits:
#define EVERGREEN_GRPH_PRIMARY_SURFACE_ADDRESS          0x6810
#define EVERGREEN_GRPH_SECONDARY_SURFACE_ADDRESS        0x6814

// upper 32-Bits:
#define EVERGREEN_GRPH_PRIMARY_SURFACE_ADDRESS_HIGH     0x681c
#define EVERGREEN_GRPH_SECONDARY_SURFACE_ADDRESS_HIGH   0x6820

// Flip status and control:
#define EVERGREEN_GRPH_UPDATE                           0x6844
#       define EVERGREEN_GRPH_SURFACE_UPDATE_PENDING    (1 << 2)

// Evergreen DCE-4 dithering control registers:
// Basics like on AVIVO: All zero disables dithering/bit depths truncation,
// Bits for enable are like cases a), b), c) above, but the meaning
// of other bits is different and there's more control bits to change
// dithering strategy.
#define EVERGREEN_FMT_BIT_DEPTH_CONTROL                 0x6fc8
#define EVERGREEN_FMT_CLAMP_CONTROL                     0x6fe4

/* northern islands - Following DCE5 specs direct excerpts from Linux Radeon KMS ni_reg.h */

#define NI_INPUT_GAMMA_CONTROL                         0x6840
#       define NI_GRPH_INPUT_GAMMA_MODE(x)             (((x) & 0x3) << 0)
#       define NI_INPUT_GAMMA_USE_LUT                  0
#       define NI_INPUT_GAMMA_BYPASS                   1
#       define NI_INPUT_GAMMA_SRGB_24                  2
#       define NI_INPUT_GAMMA_XVYCC_222                3
#       define NI_OVL_INPUT_GAMMA_MODE(x)              (((x) & 0x3) << 4)

#define NI_PRESCALE_GRPH_CONTROL                       0x68b4
#       define NI_GRPH_PRESCALE_BYPASS                 (1 << 4)

#define NI_PRESCALE_OVL_CONTROL                        0x68c4
#       define NI_OVL_PRESCALE_BYPASS                  (1 << 4)

#define NI_INPUT_CSC_CONTROL                           0x68d4
#       define NI_INPUT_CSC_GRPH_MODE(x)               (((x) & 0x3) << 0)
#       define NI_INPUT_CSC_BYPASS                     0
#       define NI_INPUT_CSC_PROG_COEFF                 1
#       define NI_INPUT_CSC_PROG_SHARED_MATRIXA        2
#       define NI_INPUT_CSC_OVL_MODE(x)                (((x) & 0x3) << 4)

#define NI_OUTPUT_CSC_CONTROL                          0x68f0
#       define NI_OUTPUT_CSC_GRPH_MODE(x)              (((x) & 0x7) << 0)
#       define NI_OUTPUT_CSC_BYPASS                    0
#       define NI_OUTPUT_CSC_TV_RGB                    1
#       define NI_OUTPUT_CSC_YCBCR_601                 2
#       define NI_OUTPUT_CSC_YCBCR_709                 3
#       define NI_OUTPUT_CSC_PROG_COEFF                4
#       define NI_OUTPUT_CSC_PROG_SHARED_MATRIXB       5
#       define NI_OUTPUT_CSC_OVL_MODE(x)               (((x) & 0x7) << 4)

#define NI_DEGAMMA_CONTROL                             0x6960
#       define NI_GRPH_DEGAMMA_MODE(x)                 (((x) & 0x3) << 0)
#       define NI_DEGAMMA_BYPASS                       0
#       define NI_DEGAMMA_SRGB_24                      1
#       define NI_DEGAMMA_XVYCC_222                    2
#       define NI_OVL_DEGAMMA_MODE(x)                  (((x) & 0x3) << 4)
#       define NI_ICON_DEGAMMA_MODE(x)                 (((x) & 0x3) << 8)
#       define NI_CURSOR_DEGAMMA_MODE(x)               (((x) & 0x3) << 12)

#define NI_GAMUT_REMAP_CONTROL                         0x6964
#       define NI_GRPH_GAMUT_REMAP_MODE(x)             (((x) & 0x3) << 0)
#       define NI_GAMUT_REMAP_BYPASS                   0
#       define NI_GAMUT_REMAP_PROG_COEFF               1
#       define NI_GAMUT_REMAP_PROG_SHARED_MATRIXA      2
#       define NI_GAMUT_REMAP_PROG_SHARED_MATRIXB      3
#       define NI_OVL_GAMUT_REMAP_MODE(x)              (((x) & 0x3) << 4)

#define NI_REGAMMA_CONTROL                             0x6a80
#       define NI_GRPH_REGAMMA_MODE(x)                 (((x) & 0x7) << 0)
#       define NI_REGAMMA_BYPASS                       0
#       define NI_REGAMMA_SRGB_24                      1
#       define NI_REGAMMA_XVYCC_222                    2
#       define NI_REGAMMA_PROG_A                       3
#       define NI_REGAMMA_PROG_B                       4
#       define NI_OVL_REGAMMA_MODE(x)                  (((x) & 0x7) << 4)

#define NI_REGAMMA_LUT_INDEX                               0x6a84
#define NI_REGAMMA_LUT_DATA                                0x6a88

/* "Volcanic Islands" DCE10/11 display controller offsets used for crtc/cur/lut/grph/viewport/etc. */
/* These are from amdgpu's vid.h defines. Multiply by 4 as we expect byte offsets, but vid.h defines
 * dword (= 4 bytes) aligned offsets. */
#define DCE10_CRTC0_REGISTER_OFFSET                 ((0x1b9c - 0x1b9c) * 4)
#define DCE10_CRTC1_REGISTER_OFFSET                 ((0x1d9c - 0x1b9c) * 4)
#define DCE10_CRTC2_REGISTER_OFFSET                 ((0x1f9c - 0x1b9c) * 4)
#define DCE10_CRTC3_REGISTER_OFFSET                 ((0x419c - 0x1b9c) * 4)
#define DCE10_CRTC4_REGISTER_OFFSET                 ((0x439c - 0x1b9c) * 4)
#define DCE10_CRTC5_REGISTER_OFFSET                 ((0x459c - 0x1b9c) * 4)
#define DCE10_CRTC6_REGISTER_OFFSET                 ((0x479c - 0x1b9c) * 4)

// NVIDIA REGISTERS:
// -----------------

// Card/Core id, e.g., NV20,30,40,... encoded in some bits:
#define NV03_PMC_BOOT_0        0x00000000
// Encodes or sets byte order (endianity) of card, starting with NV10:
// As OS/X only supports NV10 or later, we can use it unconditionally.
// On read: 0 == LE, non-0 == BE. On write, just write 0x00000001 to switch
// card to host system byte order:
#define NV03_PMC_BOOT_1        0x00000004

#endif // include once.
