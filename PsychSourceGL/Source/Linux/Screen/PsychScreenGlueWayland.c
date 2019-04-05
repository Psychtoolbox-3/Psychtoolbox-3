/*
    PsychToolbox3/Source/Linux/Screen/PsychScreenGlueWayland.c

    PLATFORMS:

        This is the Linux version with the Wayland specific backend
        code. PsychScreenGlue.c contains shared code which is shared/
        used from this file as well to avoid code duplication.

    AUTHORS:

        Mario Kleiner       mk      mario.kleiner.de@gmail.com

    HISTORY:

        4-Jan-2015          mk      Written.

    DESCRIPTION:

        Functions in this file comprise an abstraction layer for probing and controlling screen state on Wayland.

    TODO: Implement the following functions:

    PsychGetScreenDepth
    PsychSetNominalFramerate
    PsychOSSetOutputConfig

*/

/*
 * Some of the helper code for wl_output and wl_seat enumeration is copied
 * with modifications from the Weston example client "weston-info.c"
 *
 * That code carries the following copyright as of December 2014 (Weston 1.6+):
 *
 * Copyright © 2012 Philipp Brüschweiler
 *
 * Additionally we use bits and pieces from clients/window.c, with the following
 * copyright holders:
 *
 * Copyright © 2008 Kristian Høgsberg
 * Copyright © 2012-2013 Collabora, Ltd.
 *
 * All files carry the following license:
 *
 * Permission to use, copy, modify, distribute, and sell this software and
 * its documentation for any purpose is hereby granted without fee, provided
 * that the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation, and that the name of the copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific, written prior permission.  The copyright holders make
 * no representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 * THE COPYRIGHT HOLDERS DISCLAIM ALL WARRANTIES WITH REGARD TO THIS
 * SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
 * FITNESS, IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
 * RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
 * CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

/*
 * Some helper code for converting wl_output display identity info into
 * display info that is compatible with the colord color management daemon
 * is copied and adapted from Wayland/Westons cms-colord.c colord plugin.
 * This code has the following copyright compatible with our MIT/BSD-3 style
 * license:
 *
 * Copyright © 2013 Richard Hughes
 *
 * Permission to use, copy, modify, distribute, and sell this software and
 * its documentation for any purpose is hereby granted without fee, provided
 * that the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation, and that the name of the copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific, written prior permission.  The copyright holders make
 * no representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 * THE COPYRIGHT HOLDERS DISCLAIM ALL WARRANTIES WITH REGARD TO THIS
 * SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
 * FITNESS, IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
 * RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
 * CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#ifdef PTB_USE_WAYLAND

/* At the moment we still depend on Waffle even for "native" Wayland support.
 * Waffle is mostly used for initial Wayland connection setup, window creation,
 * and interacting with EGL to get the rendering api's setup and attached,
 * OpenGL contexts created and managed, GL extensions bound etc.
 *
 * Plan is to mix Wayland specific code with "minimal" Waffle here, and then incrementally
 * replace Waffle bits with Wayland bits as far as this is neccessary or useful.
 */
#ifndef PTB_USE_WAFFLE
#error PTB_USE_WAFFLE not defined for PTB_USE_WAYLAND build! We depend on Waffle at the moment, so this is a no-go!
#endif

#include "Screen.h"

#include <waffle-1/waffle_wayland.h>

#include <wayland-client.h>
#include <wayland-cursor.h>

// This header file and corresponding implementation file currently included in our
// source tree, as permitted by license. There's probably a better way to do this?
#include "presentation_timing-client-protocol.h"

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

// EvDev input style definition: Button and key codes.
#include <linux/input.h>
// libxkbcommon to map key codes to our KbName style stuff:
#include <xkbcommon/xkbcommon.h>

// colord client library for gamma table handling/color management:
#include <dlfcn.h>
#include <colord-1/colord.h>
#include <lcms2.h>

// Maximum number of slots in a gamma table to set/query: This should be plenty.
#define MAX_GAMMALUT_SIZE 16384

// These are defined in/shared with PsychScreenGlue.c:
extern PsychScreenSettingsType     displayOriginalCGSettings[kPsychMaxPossibleDisplays];            //these track the original video state before the Psychtoolbox changed it.
extern psych_bool                  displayOriginalCGSettingsValid[kPsychMaxPossibleDisplays];
extern psych_bool                  displayCursorHidden[kPsychMaxPossibleDisplays];
extern CGDisplayCount              numDisplays;

// displayCGIDs stores the wl_display handles to the display connections of each PTB logical screen.
extern CGDirectDisplayID           displayCGIDs[kPsychMaxPossibleDisplays];
extern psych_mutex displayLock;

// colord profiling inhibition state for each screen:
static psych_bool displayProfilingInhibit[kPsychMaxPossibleDisplays];

// Functions defined in PsychScreenGlue.c:
void PsychLockScreenSettings(int screenNumber);
void PsychUnlockScreenSettings(int screenNumber);
psych_bool PsychCheckScreenSettingsLock(int screenNumber);

// Hack:
extern psych_bool has_xrandr_1_2;

// Shared waffle display connection handle for whole session:
struct waffle_display *waffle_display = NULL;

// Shared corresponding wl_display handle for session:
struct wl_display* wl_display = NULL;

struct wl_compositor* wl_compositor = NULL;

// Also share native underlying EGL display:
EGLDisplay egl_display = NULL;

// Handle to the presentation extension:
struct wp_presentation *wayland_pres = NULL;

// And our presentation reference clock:
uint32_t wayland_presentation_clock_id;

static struct wl_registry *wl_registry = NULL;
static psych_bool wayland_roundtrip_needed = FALSE;

// Global shared memory object, used for cursor support:
static struct wl_shm *wl_shm = NULL;
static struct wl_cursor_theme *wayland_cursor_theme = NULL;

// Context and state for libxkbcommon - our keyboard keycode -> KbName handling:
static struct xkb_context *xkb_context = NULL;
static struct xkb_keymap *xkb_keymap = NULL;
static struct xkb_state *xkb_state = NULL;
static xkb_mod_mask_t xkb_control_mask;
static xkb_mod_mask_t xkb_alt_mask;
static xkb_mod_mask_t xkb_shift_mask;

// colord color daemon client connection object:
CdClient *colord_client = NULL;

/* Mapping of EDID style display manufacturer name to color style
 * names. This code is copied under a xxx license from Wayland/Westons
 * cms-colord.c implementation to make sure we use exactly the same
 * mapping a colord and can associate wl_output display ids with
 * colord display names. It is a modified version, adapted to our
 * execution environment.
 */
static GHashTable *pnp_ids = NULL; /* key = pnp-id, value = vendor */
static gchar *pnp_ids_data = NULL;

static bool edid_value_valid(const char *str)
{
    if (str == NULL)
        return false;
    if (str[0] == '\0')
        return false;
    if (strcmp(str, "unknown") == 0)
        return false;
    return true;
}

static gchar* get_output_id(const char* make, const char* model, const char* serial_number, int id)
{
    const gchar *tmp;
    GString *device_id;

    /* see https://github.com/hughsie/colord/blob/master/doc/device-and-profile-naming-spec.txt
     * for format and allowed values */
    device_id = g_string_new("xrandr");
    if (edid_value_valid(make)) {
        tmp = g_hash_table_lookup(pnp_ids, make);
        if (tmp == NULL)
            tmp = make;
        g_string_append_printf(device_id, "-%s", tmp);
    }
    if (edid_value_valid(model))
        g_string_append_printf(device_id, "-%s", model);
    if (edid_value_valid(serial_number))
        g_string_append_printf(device_id, "-%s", serial_number);

    /* no EDID data, so use fallback */
    if (strcmp(device_id->str, "xrandr") == 0)
        g_string_append_printf(device_id, "-drm-%i", id);

    return g_string_free(device_id, FALSE);
}

static void load_pnp_ids(void)
{
    gboolean ret = FALSE;
    gchar *tmp;
    GError *error = NULL;
    guint i;
    const gchar *pnp_ids_fn[] = { "/usr/share/hwdata/pnp.ids",
                                  "/usr/share/misc/pnp.ids",
                                  NULL };

    pnp_ids = g_hash_table_new_full(g_str_hash, g_str_equal, NULL, NULL);

    /* find and load file */
    for (i = 0; pnp_ids_fn[i] != NULL; i++) {
        if (!g_file_test(pnp_ids_fn[i], G_FILE_TEST_EXISTS))
            continue;
        ret = g_file_get_contents(pnp_ids_fn[i],
                                    &pnp_ids_data,
                                    NULL,
                                    &error);
        if (!ret) {
            printf("PTB-INFO: Failed to load PnP ids database %s: %s\n", pnp_ids_fn[i], error->message);
            g_error_free(error);
            return;
        }
        break;
    }
    if (!ret) {
        printf("PTB-INFO: No pnp.ids found!\n");
        return;
    }

    /* parse fixed offsets into lines */
    tmp = pnp_ids_data;
    for (i = 0; pnp_ids_data[i] != '\0'; i++) {
        if (pnp_ids_data[i] != '\n')
            continue;
        pnp_ids_data[i] = '\0';
        if (tmp[0] && tmp[1] && tmp[2] && tmp[3] == '\t' && tmp[4]) {
            tmp[3] = '\0';
            g_hash_table_insert(pnp_ids, tmp, tmp+4);
            tmp = &pnp_ids_data[i+1];
        }
    }
}

static void delete_pnp_ids(void)
{
    g_free(pnp_ids_data);
    g_hash_table_unref(pnp_ids);
}

// Helpers:
static void ProcessWaylandEvents(int screenNumber);

// Helpers so we can easily use/include/upgrade to Weston client sample code:
static void *
xmalloc(size_t s)
{
    return(malloc(s));
}

static void *
xzalloc(size_t s)
{
    return(calloc(1, s));
}

static char *
xstrdup(const char *s)
{
    return(strdup(s));
}

struct output_mode {
    struct wl_list link;

    uint32_t flags;
    int32_t width, height;
    int32_t refresh;
};

struct output_info {
    // struct global_info global;

    struct wl_output *output;

    struct {
        int32_t x, y;
        int32_t physical_width, physical_height;
        enum wl_output_subpixel subpixel;
        enum wl_output_transform output_transform;
        char *make;
        char *model;
    } geometry;

    struct wl_list modes;

    // Associated colord display device:
    CdDevice* colord_device;
};

// Array of information about all available Wayland outputs:
static struct output_info* displayOutputs[kPsychMaxPossibleDisplays];
// Same as above, but only stores a pointer to the wl_output:
// Currently shared with PsychWindoeGlueWayland.c:
struct wl_output* displayWaylandOutputs[kPsychMaxPossibleDisplays];

static void
print_output_info(void *data)
{
    struct output_info *output = data;
    struct output_mode *mode;
    const char *subpixel_orientation;
    const char *transform;

    //print_global_info(data);

    switch (output->geometry.subpixel) {
        case WL_OUTPUT_SUBPIXEL_UNKNOWN:
            subpixel_orientation = "unknown";
            break;
        case WL_OUTPUT_SUBPIXEL_NONE:
            subpixel_orientation = "none";
            break;
        case WL_OUTPUT_SUBPIXEL_HORIZONTAL_RGB:
            subpixel_orientation = "horizontal rgb";
            break;
        case WL_OUTPUT_SUBPIXEL_HORIZONTAL_BGR:
            subpixel_orientation = "horizontal bgr";
            break;
        case WL_OUTPUT_SUBPIXEL_VERTICAL_RGB:
            subpixel_orientation = "vertical rgb";
            break;
        case WL_OUTPUT_SUBPIXEL_VERTICAL_BGR:
            subpixel_orientation = "vertical bgr";
            break;
        default:
            fprintf(stderr, "unknown subpixel orientation %u\n",
                    output->geometry.subpixel);
            subpixel_orientation = "unexpected value";
            break;
    }

    switch (output->geometry.output_transform) {
        case WL_OUTPUT_TRANSFORM_NORMAL:
            transform = "normal";
            break;
        case WL_OUTPUT_TRANSFORM_90:
            transform = "90°";
            break;
        case WL_OUTPUT_TRANSFORM_180:
            transform = "180°";
            break;
        case WL_OUTPUT_TRANSFORM_270:
            transform = "270°";
            break;
        case WL_OUTPUT_TRANSFORM_FLIPPED:
            transform = "flipped";
            break;
        case WL_OUTPUT_TRANSFORM_FLIPPED_90:
            transform = "flipped 90°";
            break;
        case WL_OUTPUT_TRANSFORM_FLIPPED_180:
            transform = "flipped 180°";
            break;
        case WL_OUTPUT_TRANSFORM_FLIPPED_270:
            transform = "flipped 270°";
            break;
        default:
            fprintf(stderr, "unknown output transform %u\n",
                    output->geometry.output_transform);
            transform = "unexpected value";
            break;
    }

    printf("\tx: %d, y: %d,\n",
           output->geometry.x, output->geometry.y);
    printf("\tphysical_width: %d mm, physical_height: %d mm,\n",
           output->geometry.physical_width,
           output->geometry.physical_height);
    printf("\tmake: '%s', model: '%s',\n",
           output->geometry.make, output->geometry.model);
    printf("\tsubpixel_orientation: %s, output_transform: %s,\n",
           subpixel_orientation, transform);

    wl_list_for_each(mode, &output->modes, link) {
        printf("\tmode:\n");

        printf("\t\twidth: %d px, height: %d px, refresh: %.f Hz,\n",
               mode->width, mode->height,
               (float) mode->refresh / 1000);

        printf("\t\tflags:");
        if (mode->flags & WL_OUTPUT_MODE_CURRENT)
            printf(" current");
        if (mode->flags & WL_OUTPUT_MODE_PREFERRED)
            printf(" preferred");
        printf("\n");
    }
}

static void
output_handle_geometry(void *data, struct wl_output *wl_output,
                       int32_t x, int32_t y,
                       int32_t physical_width, int32_t physical_height,
                       int32_t subpixel,
                       const char *make, const char *model,
                       int32_t output_transform)
{
    struct output_info *output = data;
    GError *error = NULL;
    gchar *device_id = NULL;

    output->geometry.x = x;
    output->geometry.y = y;
    output->geometry.physical_width = physical_width;
    output->geometry.physical_height = physical_height;
    output->geometry.subpixel = subpixel;
    output->geometry.make = xstrdup(make);
    output->geometry.model = xstrdup(model);
    output->geometry.output_transform = output_transform;

    // Try to get a connection to the corresponding colord
    // device for color management of this output:
    if (colord_client) {
        // Map display manufacturer "make", display "model", (currently unknown) display serial number ""
        // and currently unknown output id "0" to a colord device_id:
        device_id = get_output_id(make, model, "", 0);
        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Mapping output via '%s' and '%s' to device_id '%s'\n", make, model, device_id);
        output->colord_device = cd_client_find_device_sync(colord_client, device_id, NULL, &error);
        if (!output->colord_device) {
            if (PsychPrefStateGet_Verbosity() > 1)
                printf("PTB-WARNING: Failed to find display device by device_id: %s. Retrying with modelname only.\n", error->message);
            g_error_free(error);
            error = NULL;
            // Failed. Retry with matching for model only - Better than nothing, but will end badly if multiple identical displays are connected:
            output->colord_device = cd_client_find_device_by_property_sync(colord_client, CD_DEVICE_PROPERTY_MODEL, model, NULL, &error);
        }
        if (!output->colord_device) {
            if (PsychPrefStateGet_Verbosity() > 0)
                printf("PTB-ERROR: On retry, failed to find display device for modelname '%s': %s. Gamma table control disabled for this output.\n", model, error->message);
            g_error_free(error);
        }
        else {
            if (!cd_device_connect_sync(output->colord_device, NULL, &error)) {
                if (PsychPrefStateGet_Verbosity() > 0)
                    printf("PTB-ERROR: Failed to connect to display device for device_id '%s': %s. Gamma table control disabled for this output.\n",
                           device_id, error->message);
                g_error_free(error);
                g_object_unref(output->colord_device);
                output->colord_device = NULL;
            }

            if (PsychPrefStateGet_Verbosity() > 2)
                printf("PTB-INFO: Found colord display '%s' - Vendor: '%s' Model: '%s' SerialNo: '%s' for output.\n",
                       device_id,
                       cd_device_get_vendor(output->colord_device),
                       cd_device_get_model(output->colord_device),
                       cd_device_get_serial(output->colord_device));
        }
    }
    else {
        output->colord_device = NULL;
    }

    if (device_id) g_free(device_id);
}

static void
output_handle_mode(void *data, struct wl_output *wl_output,
                   uint32_t flags, int32_t width, int32_t height,
                   int32_t refresh)
{
    struct output_info *output = data;
    struct output_mode *mode = xmalloc(sizeof *mode);

    mode->flags = flags;
    mode->width = width;
    mode->height = height;
    mode->refresh = refresh;

    wl_list_insert(output->modes.prev, &mode->link);
}

static const struct wl_output_listener output_listener = {
    output_handle_geometry,
    output_handle_mode,
    0,
    0
};

static void
destroy_output_info(void *data)
{
    struct output_info *output = data;
    struct output_mode *mode, *tmp;

    wl_output_destroy(output->output);

    if (output->geometry.make != NULL)
        free(output->geometry.make);
    if (output->geometry.model != NULL)
        free(output->geometry.model);

    wl_list_for_each_safe(mode, tmp, &output->modes, link) {
        wl_list_remove(&mode->link);
        free(mode);
    }

    if (output->colord_device) g_object_unref(output->colord_device);
    output->colord_device = NULL;
}

static void add_output_info(struct output_info** outputSlot, uint32_t id, uint32_t version)
{
    struct output_info *output = xzalloc(sizeof *output);

    // init_global_info(info, &output->global, id, "wl_output", version);
    // output->global.print = print_output_info;
    // output->global.destroy = destroy_output_info;
    *outputSlot = output;

    wl_list_init(&output->modes);

    output->output = wl_registry_bind(wl_registry, id,
                                      &wl_output_interface, 1);
    wl_output_add_listener(output->output, &output_listener,
                           output);
    wayland_roundtrip_needed = TRUE;
}

static void
wayland_set_presentation_clock_id(void *data, struct wp_presentation *presentation,
                                  uint32_t clk_id)
{
    struct wl_display *self = data;
    wayland_presentation_clock_id = clk_id;
    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-DEBUG: Wayland presentation clock set to id %i.\n", (int) clk_id);
}

static const struct wp_presentation_listener wayland_presentation_listener = {
    .clock_id = wayland_set_presentation_clock_id
};

struct seat_info {
//    struct global_info global;
    struct wl_seat *seat;
//    struct weston_info *info;

    uint32_t capabilities;
    char *name;

    int32_t repeat_rate;
    int32_t repeat_delay;

    // Keyboard state:
    PsychNativeBooleanType keyState[256];
    uint32_t keyTimestamp;

    // Pointer state:
    uint32_t last_serial;
    double buttonState[256];
    double posX;
    double posY;
    struct wl_surface *pointerFocusWindow;
    struct wl_surface *cursor_surface;
    struct wl_cursor *current_cursor;
    psych_bool cursor_hidden;
    int hotspot_x, hotspot_y;
};

#define kPsychMaxWaylandSeats 100
static struct seat_info* waylandSeats[kPsychMaxWaylandSeats];
static struct seat_info* waylandInputDevices[kPsychMaxWaylandSeats * 3];
static int numSeats = 0;
static int numInputDevices = 0;

static void
print_seat_info(void *data)
{
    struct seat_info *seat = data;

    //print_global_info(data);
    printf("\tNew Wayland input seat detected:\n");
    printf("\tname: %s\n", seat->name);
    printf("\tcapabilities:");

    if (seat->capabilities & WL_SEAT_CAPABILITY_POINTER)
        printf(" pointer");
    if (seat->capabilities & WL_SEAT_CAPABILITY_KEYBOARD)
        printf(" keyboard");
    if (seat->capabilities & WL_SEAT_CAPABILITY_TOUCH)
        printf(" touch");

    printf("\n");

    if (seat->repeat_rate > 0)
        printf("\tkeyboard repeat rate: %d\n", seat->repeat_rate);
    if (seat->repeat_delay > 0)
        printf("\tkeyboard repeat delay: %d\n", seat->repeat_delay);
}

static void
keyboard_handle_keymap(void *data, struct wl_keyboard *keyboard,
                       uint32_t format, int fd, uint32_t size)
{
    struct seat_info *seat = data;
    struct xkb_keymap *keymap;
    struct xkb_state *state;
    char *map_str;

    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Keyboard keymap received for keyboard of seat %p.\n", seat);

    if (!seat) {
        close(fd);
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Can't init keymap for seat %p: No seat!\n", seat);
        return;
    }

    if (format != WL_KEYBOARD_KEYMAP_FORMAT_XKB_V1) {
        close(fd);
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Can't init keymap for seat %p: Unsupported keymap format - not XKB-V1!\n", seat);
        return;
    }

    map_str = mmap(NULL, size, PROT_READ, MAP_SHARED, fd, 0);
    if (map_str == MAP_FAILED) {
        close(fd);
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Can't init keymap for seat %p: mmap() failed!\n", seat);
        return;
    }

    // Parse handed mmap()ed keymap description into a xkb keymap:
    keymap = xkb_keymap_new_from_string(xkb_context, map_str, XKB_KEYMAP_FORMAT_TEXT_V1, 0);
    munmap(map_str, size);
    close(fd);

    if (!keymap) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Can't init keymap for seat %p: Failed to compile keymap!\n", seat);
        return;
    }

    // Initialize our own state object for this keymap:
    state = xkb_state_new(keymap);
    if (!state) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Can't init keymap for seat %p: Failed to create xkb state!\n", seat);
        xkb_keymap_unref(keymap);
        return;
    }

    // Replace our old keymap and state with these new specimens:
    xkb_keymap_unref(xkb_keymap);
    xkb_state_unref(xkb_state);
    xkb_keymap = keymap;
    xkb_state = state;

    xkb_control_mask = 1 << xkb_keymap_mod_get_index(xkb_keymap, "Control");
    xkb_alt_mask = 1 << xkb_keymap_mod_get_index(xkb_keymap, "Mod1");
    xkb_shift_mask = 1 << xkb_keymap_mod_get_index(xkb_keymap, "Shift");

    return;
}

static void
keyboard_handle_enter(void *data, struct wl_keyboard *keyboard,
                      uint32_t serial, struct wl_surface *surface,
                      struct wl_array *keys)
{
    struct seat_info *seat = data;
    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Keyboard focus entered surface %p on keyboard of seat %p.\n", surface, seat);

    // Reset keyState to "nothing pressed":
    memset(&seat->keyState, 0, sizeof(seat->keyState));
    seat->keyTimestamp = 0;
}

static void
keyboard_handle_leave(void *data, struct wl_keyboard *keyboard,
                      uint32_t serial, struct wl_surface *surface)
{
    struct seat_info *seat = data;
    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Keyboard focus lost for surface %p on keyboard of seat %p.\n", surface, seat);

    // Reset keyState to "nothing pressed":
    memset(&seat->keyState, 0, sizeof(seat->keyState));
    seat->keyTimestamp = 0;
}

static void
keyboard_handle_key(void *data, struct wl_keyboard *keyboard,
                    uint32_t serial, uint32_t time, uint32_t key,
                    uint32_t state)
{
    struct seat_info *seat = data;
    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: %i: [%i msecs]: KEY %i -> %i\n", serial, time, key, state);

    // Update keyState for this key:
    if (key < 256) {
        seat->keyState[key] = (PsychNativeBooleanType) state;
        seat->keyTimestamp = time;
    }
}

static void
keyboard_handle_modifiers(void *data, struct wl_keyboard *keyboard,
                          uint32_t serial, uint32_t mods_depressed,
                          uint32_t mods_latched, uint32_t mods_locked,
                          uint32_t group)
{
    struct seat_info *seat = data;
    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Keyboard modifiers updated for keyboard of seat %p.\n", seat);
}

// Repeat info unused, and unsupported as long as we want protocol v3:
static void
keyboard_handle_repeat_info(void *data, struct wl_keyboard *keyboard,
                            int32_t rate, int32_t delay)
{
//     struct seat_info *seat = data;
//
//     seat->repeat_rate = rate;
//     seat->repeat_delay = delay;
}

static const struct wl_keyboard_listener keyboard_listener = {
    keyboard_handle_keymap,
    keyboard_handle_enter,
    keyboard_handle_leave,
    keyboard_handle_key,
    keyboard_handle_modifiers,
    keyboard_handle_repeat_info,
};

static void
pointer_handle_enter(void *data,
                      struct wl_pointer *wl_pointer,
                      uint32_t serial,
                      struct wl_surface *surface,
                      wl_fixed_t surface_x,
                      wl_fixed_t surface_y)
{
    struct wl_buffer *buffer;
    struct wl_cursor_image *image;
    struct seat_info *seat = data;

    if (PsychPrefStateGet_Verbosity() > 5) printf("PTB-DEBUG: Pointer focus entered surface %p for pointer of seat %p.\n", surface, seat);

    seat->last_serial = serial;
    seat->posX = wl_fixed_to_double(surface_x);
    seat->posY = wl_fixed_to_double(surface_y);
    seat->pointerFocusWindow = surface;

    // Get associated windowRecord for currently focused window / wl_surface for this
    // seats pointer. We can get the windows globalrect this way to convert the relative
    // surface x,y position to global x,y position as needed for proper multi-window
    // operation:
    if (seat->pointerFocusWindow) {
        PsychWindowRecordType *windowRecord = (PsychWindowRecordType*) wl_surface_get_user_data(seat->pointerFocusWindow);
        if (windowRecord) {
            seat->posX += windowRecord->globalrect[kPsychLeft];
            seat->posY += windowRecord->globalrect[kPsychTop];
        }
    }

    // Cursor setup for this pointer: Any cursor assigned and visible at the moment?
    if (seat->cursor_surface && seat->current_cursor) {
        // Yes, set it. Get first cursor image for current cursor:
        // There could be multiple, e.g., for animated cursors, but we only
        // care about one image - We don't do animated cursors:
        image = seat->current_cursor->images[0];

        // Assign cursor surface to this pointer:
        wl_pointer_set_cursor(wl_pointer, serial, seat->cursor_surface, image->hotspot_x, image->hotspot_y);
        seat->hotspot_x = image->hotspot_x;
        seat->hotspot_y = image->hotspot_y;
    }
    else {
        // No. Hide the cursor by assigning a NULL surface:
        wl_pointer_set_cursor(wl_pointer, serial, NULL, 0, 0);
    }

    return;
}

static void
pointer_handle_leave(void *data,
                      struct wl_pointer *wl_pointer,
                      uint32_t serial,
                      struct wl_surface *surface)
{
    struct seat_info *seat = data;
    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Pointer focus lost for surface %p for pointer of seat %p.\n", surface, seat);

    seat->pointerFocusWindow = NULL;
}

static void
pointer_handle_motion(void *data,
                      struct wl_pointer *wl_pointer,
                      uint32_t time,
                      wl_fixed_t surface_x,
                      wl_fixed_t surface_y)
{
    struct seat_info *seat = data;
    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Pointer motion on seat %p: time %i, x = %lf y = %lf\n", seat, time, wl_fixed_to_double(surface_x), wl_fixed_to_double(surface_y));
    seat->posX = wl_fixed_to_double(surface_x);
    seat->posY = wl_fixed_to_double(surface_y);

    // Get associated windowRecord for currently focused window / wl_surface for this
    // seats pointer. We can get the windows globalrect this way to convert the relative
    // surface x,y position to global x,y position as needed for proper multi-window
    // operation:
    if (seat->pointerFocusWindow) {
        PsychWindowRecordType *windowRecord = (PsychWindowRecordType*) wl_surface_get_user_data(seat->pointerFocusWindow);
        if (windowRecord) {
            seat->posX += windowRecord->globalrect[kPsychLeft];
            seat->posY += windowRecord->globalrect[kPsychTop];
        }
    }
}

static void
pointer_handle_button(void *data,
                      struct wl_pointer *wl_pointer,
                      uint32_t serial,
                      uint32_t time,
                      uint32_t button,
                      uint32_t state)
{
    struct seat_info *seat = data;
    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Pointer button changed on seat %p: time %i, button = %i, state = %i\n", seat, time, button, state);

    // Remap standard mouse buttons of 3-5 Button mice to our standard 0-4:
    switch (button) {
        case BTN_LEFT:
            button = 0;
        break;

        case BTN_MIDDLE:
            button = 1;
        break;

        case BTN_RIGHT:
            button = 2;
        break;

        case BTN_SIDE:
            button = 3;
        break;

        case BTN_EXTRA:
            button = 4;
        break;

        default:
            button = button - BTN_MISC + 5;
    }

    // No-Op for buttons higher than 255, as we don't keep stat for more buttons:
    if (button > 255) return;

    // Update button state:
    seat->buttonState[button] = (double) state;
}

static void
pointer_handle_axis(void *data,
                    struct wl_pointer *wl_pointer,
                    uint32_t time,
                    uint32_t axis,
                    wl_fixed_t value)
{
    struct seat_info *seat = data;
    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Pointer axis event on seat %p: time %i, axis = %i, value = %lf\n", seat, time, axis, wl_fixed_to_double(value));

    // TODO: No-Op for now, possibly map to our valuators...
    return;
}

static const struct wl_pointer_listener pointer_listener = {
    pointer_handle_enter,
    pointer_handle_leave,
    pointer_handle_motion,
    pointer_handle_button,
    pointer_handle_axis,
};

static void
seat_handle_capabilities(void *data, struct wl_seat *wl_seat,
                         enum wl_seat_capability caps)
{
    struct seat_info *seat = data;

    seat->capabilities = caps;

    /* we want listen for repeat_info from wl_keyboard, but only
     * do so if the seat info is >= 4 and if we actually have a
     * keyboard */
    //if (seat->global.version < 4)
    //    return;

    if (caps & WL_SEAT_CAPABILITY_KEYBOARD) {
        struct wl_keyboard *keyboard;

        keyboard = wl_seat_get_keyboard(seat->seat);
        wl_keyboard_add_listener(keyboard, &keyboard_listener,
                                 seat);

        //seat->info->roundtrip_needed = true;
        wayland_roundtrip_needed = TRUE;

        // Add backpointer to our input device array to
        // this seat, so we can access its keyboard state:
        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: New keyboard enumerated under index %i for seat %p.\n", numInputDevices, seat);
        waylandInputDevices[numInputDevices++] = seat;
    }

    if (caps & WL_SEAT_CAPABILITY_POINTER) {
        struct wl_pointer *pointer;

        pointer = wl_seat_get_pointer(seat->seat);
        wl_pointer_add_listener(pointer, &pointer_listener, seat);

        // Add backpointer to our input device array to
        // this seat, so we can access its pointer state:
        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: New pointer enumerated under index %i for seat %p.\n", numInputDevices, seat);
        waylandInputDevices[numInputDevices++] = seat;

        // Create cursor surface for this pointer:
        seat->cursor_surface = wl_compositor_create_surface(wl_compositor);
        if (!seat->cursor_surface) {
            if (PsychPrefStateGet_Verbosity() > 0) {
                printf("PTB-ERROR: Creating cursor surface for pointer device %i [seat %p] failed. Cursor support unavailable on this seat!\n", numInputDevices, seat);
            }
        }

        // First time load of cursor theme needed?
        if (!wayland_cursor_theme) {
            wayland_cursor_theme = wl_cursor_theme_load(NULL, 32, wl_shm);
            if (!wayland_cursor_theme) {
                if (PsychPrefStateGet_Verbosity() > 0) {
                    printf("PTB-ERROR: Loading Wayland cursor theme failed. Cursor support unavailable!\n");
                }
            }
        }

        if (wayland_cursor_theme) {
            // Assign initial default cursor: The classic left-tilted arrow pointer:
            seat->current_cursor = wl_cursor_theme_get_cursor(wayland_cursor_theme, "left_ptr");
            if (!seat->current_cursor) {
                if (PsychPrefStateGet_Verbosity() > 0) {
                    printf("PTB-ERROR: Loading Wayland cursor from theme for pointer device %i [seat %p] failed. Cursor support unavailable on this seat!\n", numInputDevices, seat);
                }
            }
        }

        // Cursor hidden by default:
        seat->cursor_hidden = TRUE;

        //seat->info->roundtrip_needed = true;
        wayland_roundtrip_needed = TRUE;
    }
}

static void
seat_handle_name(void *data, struct wl_seat *wl_seat,
                 const char *name)
{
    struct seat_info *seat = data;
    seat->name = xstrdup(name);
}

static const struct wl_seat_listener seat_listener = {
    seat_handle_capabilities,
    seat_handle_name,
};

static void
destroy_seat_info(void *data)
{
    struct seat_info *seat = data;

    wl_seat_destroy(seat->seat);

    if (seat->name != NULL)
        free(seat->name);

    // Delete cursor surface with the cursor image for this seat:
    seat->current_cursor = NULL;
    seat->cursor_hidden = TRUE;

    if (seat->cursor_surface) {
        wl_surface_destroy(seat->cursor_surface);
        seat->cursor_surface = NULL;
    }
}

static void
add_seat_info(struct seat_info **seatSlot, uint32_t id, uint32_t version)
{
    struct seat_info *seat = xzalloc(sizeof *seat);

    /* required to set roundtrip_needed to true in capabilities
     * handler */
    // seat->info = info;

    //init_global_info(info, &seat->global, id, "wl_seat", version);
    //seat->global.print = print_seat_info;
    //seat->global.destroy = destroy_seat_info;

    *seatSlot = seat;

    // Only bind to wl_seat_interface version 3 for now, so we can work
    // with our old wayland client library:
    // seat->seat = wl_registry_bind(wl_registry, id, &wl_seat_interface, MIN(version, 4));
    seat->seat = wl_registry_bind(wl_registry, id, &wl_seat_interface, 3);

    wl_seat_add_listener(seat->seat, &seat_listener, seat);

    seat->repeat_rate = seat->repeat_delay = -1;

    //info->roundtrip_needed = true;
    wayland_roundtrip_needed = TRUE;
}

static void
wayland_registry_listener_global(void *data,
                                 struct wl_registry *registry,
                                 uint32_t name,
                                 const char *interface,
                                 uint32_t version)
{
    struct wl_display *self = data;

    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Wayland registry extension candidate: %s\n", interface);

    // Look for stable presentation interface v1+:
    if (!strcmp(interface, "wp_presentation") && (version >= 1)) {
        wayland_pres = wl_registry_bind(registry, name, &wp_presentation_interface, 1);
        if (!wayland_pres) {
            if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: wl_registry_bind for wp_presentation_interface failed!\n");
            return;
        }

        wp_presentation_add_listener(wayland_pres, &wayland_presentation_listener, self);
        if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-DEBUG: Wayland wp_presentation_interface bound!\n");
    }

    // Look for Wayland outputs ~ video outputs ~ displays ~ our PTB screens:
    // Not yet sure if wl_output ~ PTB screen is the optimal abstraction/mapping,
    // but as a starter...
    if (!strcmp(interface, "wl_output") && (version >= 1)) {
        add_output_info(&displayOutputs[numDisplays], name, version);
        if (PsychPrefStateGet_Verbosity() > 4) {
            printf("PTB-DEBUG: New output display screen %i enumerated.\n", numDisplays);
        }
        numDisplays++;
    }

    // Look for Wayland wl_seat ~ collections of input devices:
    if (!strcmp(interface, "wl_seat") && (version >= 1)) {
        if (numSeats >= kPsychMaxWaylandSeats) {
            printf("PTB-WARNING: Could not enumerate all Wayland seats, as maximum number of %i seats reached.\n", (int) kPsychMaxWaylandSeats);
            return;
        }

        add_seat_info(&waylandSeats[numSeats], name, version);
        if (PsychPrefStateGet_Verbosity() > 4) {
            printf("PTB-DEBUG: New wayland seat %i enumerated.\n", numSeats);
        }
        numSeats++;
    }

    if (!strcmp(interface, "wl_shm")) {
        wl_shm = wl_registry_bind(registry, name, &wl_shm_interface, 1);
        // wl_shm_add_listener(shm, &shm_listener, NULL);
    }

    return;
}

static void
wayland_registry_listener_global_remove(void *data, struct wl_registry *registry, uint32_t name)
{
    return;
}

static const struct wl_registry_listener wayland_registry_listener = {
    .global = wayland_registry_listener_global,
    .global_remove = wayland_registry_listener_global_remove
};

struct wp_presentation *get_wayland_presentation_extension(PsychWindowRecordType* windowRecord)
{
    // Already have a cached presentation_interface? If so, just return it:
    return(wayland_pres);
}

psych_bool PsychWaylandGetKbNames(PsychGenericScriptType *kbNames)
{
    xkb_keycode_t keyCode;
    xkb_keysym_t keySym;
    char symName[64];

    // Ready?
    if (!xkb_state) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychWaylandGetKbNames(): No keymap available for KbName() initialisation!\n");
        return(FALSE);
    }

    // Iterate over all keyCodes and convert them to key names and return
    // them in the kbNames vector which is used by KbName() and friends:
    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Full keymap dump:\n");
    for (keyCode = 0; keyCode < 256; keyCode++) {
        keySym = xkb_state_key_get_one_sym(xkb_state, keyCode + 8);
        if (0 < xkb_keysym_get_name(keySym, &symName[0], 64)) {
            // Character found: Return its ASCII name string:
            PsychSetCellVectorStringElement((int) keyCode, symName, kbNames);

            if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: keyCode %i --> KeySym %i --> %s\n", (int) keyCode, (int) keySym, symName);
        }
        else {
            // No luck: Return empty string:
            PsychSetCellVectorStringElement((int) keyCode, "", kbNames);

            if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: keyCode %i --> KeySym %i --> No mapping\n", (int) keyCode, (int) keySym);
        }
    }
    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Done.\n\n");

    return(TRUE);
}

/* PsychWaylandEventTimeToRefTime - Convert from waylaned 32-Bit input
 * event timestamp in CLOCK_MONOTONIC time to PTB reference time.
 *
 * Unfortunately libinput as used by Wayland only provides millisecond
 * resolution timestamps.
 */
double PsychWaylandEventTimeToRefTime(uint32_t event_lo)
{
    struct timespec tp;
    psych_uint64 now, eventTimeMsecs;
    uint32_t now_lo, now_hi;
    double eventSecs;

    // Get current CLOCK_MONOTONIC time as baseline. Wayland uses libinput
    // to get input event timestamps in kernel reported CLOCK_MONOTONIC
    // time. Unfortunately libinput converts those timestamps from nanosecond
    // resolution to millisecond resolution and then only keeps the low 32-Bits.
    // We need to invert the process, recovering full 64-Bit milliseconds.

    // Get current time, convert to milliseconds, split in high and low 32-Bits:
    clock_gettime(CLOCK_MONOTONIC, &tp);
    now = tp.tv_sec * 1000 + tp.tv_nsec / 1000000;
    now_lo = now & 0xffffffff;
    now_hi = (now >> 32) & 0xffffffff;

    // As the event can't be far in the past, our upper 32-Bits of current time
    // are likely identical with the thrown away upper 32-Bits of the libinput
    // event, whereas the our lower 32-Bits of current time should be higher as
    // the input event was created in the past. Check this assumption. If this
    // isn't the case then system time must just recently have crossed a 32-Bit
    // time segment boundary, iow. the upper 32-Bits have incremented by 1, whereas
    // the lower 32-Bits have wrapped around to a zero or low value, therefore our
    // low timestamp appears to be smaller than that of a event from the past.
    // Correct for this by reducing now_hi by 1, so now_hi represents the 32-Bit
    // time segment corresponding to the input event. The only way this could fail
    // is if the input event was more than 2^32 milliseconds in the past, iow. over
    // 49 days have passed since the last user input -- very unlikely for a running
    // experiment.
    if (event_lo > now_lo) {
        now_hi--;
    }

    // Reconstruct 64-Bit msecs event time from our now_hi upper 32-Bits and
    // the events reported lower 32-Bits:
    eventTimeMsecs = ((psych_uint64) now_hi << 32) | event_lo;

    // Convert from integer milliseconds to double seconds in CLOCK_MONOTONIC time:
    eventSecs = ((double) eventTimeMsecs) / 1000.0;

    // Convert from monotonic time to Psychtoolbox "GetSecs" time, if needed:
    eventSecs = PsychOSMonotonicToRefTime(eventSecs);

    return(eventSecs);
}

psych_bool PsychWaylandGetKeyboardState(int deviceId, int numKeys, PsychNativeBooleanType *buttonStates, double *timeStamp)
{
    ProcessWaylandEvents(0);

    // deviceId -1 means "auto-detected default keyboard". Simply use the first found keyboard device:
    if (deviceId == -1) {
        for (deviceId = 0; deviceId < numInputDevices; deviceId++) {
            if (waylandInputDevices[deviceId] && (waylandInputDevices[deviceId]->capabilities & WL_SEAT_CAPABILITY_KEYBOARD))
                break;
        }
    }

    // Outside valid range for input devices?
    if (deviceId < 0 || deviceId >= numInputDevices) return(FALSE);

    // No device under that id, or device isn't a keyboard?
    if (!waylandInputDevices[deviceId] || !(waylandInputDevices[deviceId]->capabilities & WL_SEAT_CAPABILITY_KEYBOARD)) return(FALSE);

    // Copy current keyState to return vector:
    memcpy(buttonStates, &(waylandInputDevices[deviceId]->keyState[0]), MIN(numKeys, 256) * sizeof(*buttonStates));

    // Remap 32-Bit millisecond timestamp to GetSecs time:
    *timeStamp = PsychWaylandEventTimeToRefTime(waylandInputDevices[deviceId]->keyTimestamp);

    return(TRUE);
}

psych_bool PsychWaylandGetMouseState(int deviceId, int *mouseX, int *mouseY, int numButtons, double *buttonArray, void **focusWindow)
{
    ProcessWaylandEvents(0);

    // deviceId -1 means "auto-detected default pointer". Simply use the first found pointer device:
    if (deviceId == -1) {
        for (deviceId = 0; deviceId < numInputDevices; deviceId++) {
            if (waylandInputDevices[deviceId] && (waylandInputDevices[deviceId]->capabilities & WL_SEAT_CAPABILITY_POINTER))
                break;
        }
    }

    // Outside valid range for input devices?
    if (deviceId < 0 || deviceId >= numInputDevices) return(FALSE);

    // No device under that id, or device isn't a pointer?
    if (!waylandInputDevices[deviceId] || !(waylandInputDevices[deviceId]->capabilities & WL_SEAT_CAPABILITY_POINTER)) return(FALSE);

    // Copy current keyState to return vector:
    memcpy(buttonArray, &(waylandInputDevices[deviceId]->buttonState[0]), MIN(numButtons, 256) * sizeof(*buttonArray));

    *mouseX = (int) waylandInputDevices[deviceId]->posX;
    *mouseY = (int) waylandInputDevices[deviceId]->posY;
    *focusWindow = (void*) waylandInputDevices[deviceId]->pointerFocusWindow;

    return(TRUE);
}

static void ProcessWaylandEvents(int screenNumber)
{
    // Block until the Wayland server has processed all pending requests and
    // has sent out pending events on all event queues. This should ensure
    // that the registry listener has received announcement of the shell and
    // compositor.
    do {
        wayland_roundtrip_needed = FALSE;
        if (wl_display_roundtrip(wl_display) == -1) {
            if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: ProcessWaylandEvents(): wl_display_roundtrip failed!\n");
        }
        // Repeat until everything is enumerated.
    } while (wayland_roundtrip_needed);
}

struct output_mode* PsychWaylandGetCurrentMode(int screenId)
{
    struct output_info* output;
    struct output_mode* mode;

    if (screenId >= numDisplays || screenId < 0) PsychErrorExitMsg(PsychError_internal, "screenNumber is out of range");
    output = displayOutputs[screenId];

    wl_list_for_each(mode, &output->modes, link) {
        if (mode->flags & WL_OUTPUT_MODE_CURRENT)
            return(mode);
    }

    return(NULL);
}

// Linux only: Retrieve modeline and crtc_info for a specific output on a specific screen:
// Caution: If crtc is non-NULL and receives a valid XRRCrtcInfo*, then this pointer must
//          be released by the caller via XRRFreeCrtcInfo(crtc), or resources will leak!
// Must be called under display lock protection!
XRRModeInfo* PsychOSGetModeLine(int screenId, int outputIdx, XRRCrtcInfo **crtc)
{
    int m;
    static XRRModeInfo rrmode; // XRRAllocModeInfo (char *name, int nameLength);
    struct output_info* output;
    struct output_mode* mode;

    if (screenId >= numDisplays || screenId < 0) PsychErrorExitMsg(PsychError_internal, "screenNumber is out of range");
    output = displayOutputs[screenId];

    memset(&rrmode, 0, sizeof(rrmode));
    XRRCrtcInfo *crtc_info = calloc(1, sizeof(XRRCrtcInfo));

    // Query info about video modeline and crtc of output 'outputIdx':
    if (PsychScreenToHead(screenId, outputIdx) >= 0) {
        // Get current video mode for output:
        mode = PsychWaylandGetCurrentMode(screenId);

        // Assign minimal needed info to satisfy our callers:
        rrmode.width = mode->width;
        rrmode.height = mode->height;

        // Only fill (x,y) start position in compositor space:
        crtc_info->x = output->geometry.x;
        crtc_info->y = output->geometry.y;

        // Not quite true, as this should be viewport size, not
        // size of video mode:
        crtc_info->width = mode->width;
        crtc_info->height = mode->height;

        // Could do rotation if we wanted...
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

    // Return a pointer to our static rrmode. Welcome to hackistan!
    return(&rrmode);
}

const char* PsychOSGetOutputProps(int screenId, int outputIdx, unsigned long *mm_width, unsigned long *mm_height, unsigned long *rrOutputPrimary)
{
    static char outputName[100];
    struct output_info* output;

    if (screenId >= numDisplays || screenId < 0) PsychErrorExitMsg(PsychError_internal, "screenNumber is out of range");
    if (PsychScreenToHead(screenId, outputIdx) < 0) PsychErrorExitMsg(PsychError_internal, "outputIdx is out of range");

    output = displayOutputs[screenId];

    // Store output name to return:
    sprintf(outputName, "%s %s", output->geometry.make, output->geometry.model);

    // And width / height in mm:
    if (mm_width) *mm_width = (unsigned long) output->geometry.physical_width;
    if (mm_height) *mm_height = (unsigned long) output->geometry.physical_height;
    if (rrOutputPrimary) *rrOutputPrimary = 0;

    // Return output name:
    return(&outputName[0]);
}

void InitCGDisplayIDList(void)
{
    Dl_info libcolord_info;
    GError *error;
    int i;

    // Define waffle window system backend to use: Wayland only, obviously.
    static int32_t init_attrs[3] = {
        WAFFLE_PLATFORM,
        WAFFLE_PLATFORM_WAYLAND,
        0,
    };

    // NULL-out array of displays:
    for (i = 0; i < kPsychMaxPossibleDisplays; i++) {
        displayCGIDs[i] = NULL;
        displayWaylandOutputs[i] = NULL;
        displayProfilingInhibit[i] = FALSE;
    }

    // Clear all input devices and seats:
    for (i = 0; i< kPsychMaxWaylandSeats; i++) waylandSeats[i] = NULL;
    for (i = 0; i< 3* kPsychMaxWaylandSeats; i++) waylandInputDevices[i] = NULL;

    // Preinit screen to head mappings to identity default:
    PsychInitScreenToHeadMappings(0);

    // Initial count of screens is zero:
    numDisplays = 0;

    // Create libxkbcommon context:
    xkb_context = xkb_context_new(0);

    // Initialize Waffle for Wayland display system backend:
    if (PsychPrefStateGet_Verbosity() > 2) {
        printf("PTB-INFO: Using FOSS Waffle display backend library, written by Chad Versace, Copyright 2012 - 2015 Intel.\n");
    }

    if (!waffle_init(init_attrs) && (waffle_error_get_code() != WAFFLE_ERROR_ALREADY_INITIALIZED)) {
        // Game over:
        if (PsychPrefStateGet_Verbosity() > 0) {
            printf("PTB-ERROR: Could not initialize Waffle Wayland/EGL backend - Error: %s.\n", waffle_error_to_string(waffle_error_get_code()));
            printf("PTB-ERROR: Try to fix the reason for the error, then restart Octave/Matlab and retry.\n");
        }
        PsychErrorExitMsg(PsychError_system, "FATAL ERROR: Couldn't initialize Waffle's Wayland backend! Game over!");
    }

    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Trying to connect Wayland Waffle to display '%s'.\n", getenv("PSYCH_WAFFLE_DISPLAY"));
    waffle_display = waffle_display_connect(getenv("PSYCH_WAFFLE_DISPLAY"));

    if (!waffle_display) {
        if (PsychPrefStateGet_Verbosity() > 0) {
            printf("PTB-ERROR: Could not connect Waffle to display: %s.\n", waffle_error_to_string(waffle_error_get_code()));
            if (!getenv("WAYLAND_DISPLAY")) printf("PTB-ERROR: Seems Screen() is not running on a Wayland display server? That's a no-go for this Wayland-only Screen() mex file!\n");
        }
        PsychErrorExitMsg(PsychError_system, "FATAL ERROR: Couldn't open display connection to Wayland server! Game over!");
    }

    // Extract EGL_Display for backends which use EGL:
    union waffle_native_display* wafflenatdis = waffle_display_get_native(waffle_display);
    egl_display = wafflenatdis->wayland->egl_display;
    wl_display = wafflenatdis->wayland->wl_display;
    wl_compositor = wafflenatdis->wayland->wl_compositor;

    // Release the waffle_native_display:
    free(wafflenatdis);
    wafflenatdis = NULL;

    if (!getenv("PSYCH_DISABLE_COLORD")) {
        // Initialize connection to colord - the color management daemon:
        colord_client = cd_client_new();

        // Awful hack: If libcolord gets unloaded due to an unload of the Screen mex file,
        // then on a reload of Screen the cd_client_new() function above would try to reregister
        // some internal GType it forgot it already had registered. This would fail, would then
        // cause a failure of cd_client_new() and then a crash for ourselves. Yay!
        //
        // To prevent this we dlopen() libcolord again, for the sole purpose of getting
        // its reference count bumped by one, so a Screen() unload will not cause the
        // library to unload during the complete livetime of our hosting process, and
        // thereby prevent it from forgetting it already has registered that GType.
        if (colord_client && dladdr((void*) cd_client_new, &libcolord_info)) {
            // Yes, got the info.
            if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: libcolord.so path as detected by dladdr(): %s\n", libcolord_info.dli_fname);

            // Try to re-open libcolord.so. This will not really
            // reload it due to RTLD_NOLOAD, just increment its
            // refcount to prevent it from unloading when Screen
            // gets unloaded due to a "clear Screen" etc.:
            dlerror();
            if (dlopen(libcolord_info.dli_fname, RTLD_NOW | RTLD_NOLOAD | RTLD_GLOBAL)) {
                if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: Locked libcolord into memory.\n");
            }
            else {
                if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Failed to lock libcolord into memory [%s]! Expect color management failure after reloading Screen()!!\n", dlerror());
            }
            dlerror();
        }
        else {
            if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Failed to lock libcolord into memory! Expect color management failure when reloading Screen()!!\n");
        }

        error = NULL;
        if (!colord_client || !cd_client_connect_sync(colord_client, NULL, &error)) {
            if (colord_client) g_object_unref(colord_client);
            colord_client = NULL;
            if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-INFO: Couldn't connect to color daemon [%s]! Gamma table functions will be disabled.\n", error->message);
            if (error) g_error_free(error);
        }
        else {
            // Load Plug & Play device id database for proper matching of
            // wl_output's to colord managed display devices:
            load_pnp_ids();
            if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Connected to colord.\n");
        }
    }
    else {
        colord_client = NULL;
        if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: Disabled colord support via users setenv('PSYCH_DISABLE_COLORD', '1'). Gamma table functions disabled.\n");
    }

    // Get our own wl_registry, do the enumeration and binding:
    wl_registry = wl_display_get_registry(wl_display);
    if (!wl_registry) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: wl_display_get_registry failed\n");
        PsychErrorExitMsg(PsychError_system, "FATAL ERROR: Initialisation failed! Game over!");
    }

    if (wl_registry_add_listener(wl_registry, &wayland_registry_listener, wl_display) < 0) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: wl_registry_add_listener failed\n");
        PsychErrorExitMsg(PsychError_system, "FATAL ERROR: Initialisation failed! Game over!");
    }

    // Block until the Wayland server has processed all pending requests and
    // has sent out pending events on all event queues. This should ensure
    // that the registry listener has received announcement of the shell and
    // compositor.
    do {
        wayland_roundtrip_needed = FALSE;
        if (wl_display_roundtrip(wl_display) == -1) {
            if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: wl_display_roundtrip failed\n");
            PsychErrorExitMsg(PsychError_system, "FATAL ERROR: Initialisation failed! Game over!");
        }
        // Repeat until everything is enumerated.
    } while (wayland_roundtrip_needed);

    // Print info about all enumerated wl_seats:
    for (i = 0; i < numSeats; i++) {
        // Detailed info about enumerated output:
        if (PsychPrefStateGet_Verbosity() > 3) print_seat_info(waylandSeats[i]);
    }

    // Initialize screenId -> GPU headId mapping to identity mappings for numDisplays:
    PsychInitScreenToHeadMappings(numDisplays);

    // Setup the screenNumber --> Wayland display mappings:
    for (i = 0; i < numDisplays && i < kPsychMaxPossibleDisplays; i++) {
        // Store the wl_output* handle for retrieval by our client code:
        displayWaylandOutputs[i] = displayOutputs[i]->output;

        // Detailed info about enumerated output:
        if (PsychPrefStateGet_Verbosity() > 3) print_output_info(displayOutputs[i]);

        // Set reference crtc == our output info for primary output to always 0,
        // to create a screen == output mapping:
        PsychSetScreenToHead(i, 0, 0);
    }

    if (numDisplays > 1) printf("PTB-INFO: A total of %i Wayland display screens is available for use.\n", numDisplays);

    return;
}

void PsychCleanupDisplayGlue(void)
{
    int i;

    // Make sure the mmio mapping is shut down:
    PsychOSShutdownPsychtoolboxKernelDriverInterface();

    // Destroy libxkbcommon resources:
    if (xkb_context) {
        xkb_state_unref(xkb_state);
        xkb_keymap_unref(xkb_keymap);
        xkb_context_unref(xkb_context);
        xkb_context = NULL;
        xkb_state = NULL;
        xkb_keymap = NULL;
    }

    PsychLockDisplay();

    // Delete info about all enumerated wl_seats:
    for (i = 0; i < numSeats; i++) {
        if (waylandSeats[i]) destroy_seat_info(waylandSeats[i]);
        waylandSeats[i] = NULL;
    }
    numSeats = 0;

    // Clear all input device pointers for the now destroyed seats:
    for (i = 0; i< 3* kPsychMaxWaylandSeats; i++) waylandInputDevices[i] = NULL;
    numInputDevices = 0;

    // Delete cursor theme, if any:
    if (wayland_cursor_theme) {
        wl_cursor_theme_destroy(wayland_cursor_theme);
        wayland_cursor_theme = NULL;
    }

    // Go trough full screen list:
    for (i = 0; i < PsychGetNumDisplays(); i++) {
        if (displayOutputs[i]) {
            // Release all profiling inhibitions on colord managed displays:
            if (displayProfilingInhibit[i]) {
                PsychWaylandProfilingInhibit(i, FALSE);
            }
            destroy_output_info(displayOutputs[i]);
            displayOutputs[i] = NULL;
            displayWaylandOutputs[i] = NULL;
            displayCGIDs[i] = NULL;
        }
    }

    // Release wayland shared memory:
    if (wl_shm) wl_shm_destroy(wl_shm);
    wl_shm = NULL;

    // Reset our binding to presentation_feedback extension if this is our last onscreen window to close:
    wp_presentation_destroy(wayland_pres);
    wayland_pres = NULL;

    // Destroy our reference to the registry:
    wl_registry_destroy(wl_registry);
    wl_registry = NULL;

    // Release our shared waffle display connection:
    waffle_display_disconnect(waffle_display);
    waffle_display = NULL;

    // Close colord connection, delete device.
    if (colord_client) {
        g_object_unref(colord_client);
        colord_client = NULL;
        delete_pnp_ids();
    }

    PsychUnlockDisplay();

    // Destroy the display lock mutex, now that we're done with it for this Screen() session instance:
    PsychDestroyMutex(&displayLock);

    // All connections should be closed now. We can't NULL-out the display list, but
    // our scripting host environment will flush the Screen - Mexfile anyway...
    return;
}

void PsychGetScreenDepths(int screenNumber, PsychDepthType *depths)
{
    int* x11_depths = NULL;
    int  i, count = 0;

    if (screenNumber >= numDisplays || screenNumber < 0) PsychErrorExitMsg(PsychError_internal, "screenNumber is out of range");

    // Update out view of this screens configuration:
    PsychLockDisplay();
//      ProcessWaylandEvents(screenNumber);
//
//     if (displayCGIDs[screenNumber]) {
//         x11_depths = XListDepths(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber), &count);
//     }

    PsychUnlockDisplay();

    if (x11_depths && depths && count > 0) {
        // Query successful: Add all values to depth struct:
//         for(i=0; i<count; i++) PsychAddValueToDepthStruct(x11_depths[i], depths);
//         XFree(x11_depths);
    }
    else {
        // Query failed: Assume at least 32 bits is available.
        printf("PTB-WARNING: Couldn't query available display depths values! Returning a made up list...\n");
        PsychAddValueToDepthStruct(32, depths);
        PsychAddValueToDepthStruct(24, depths);
        PsychAddValueToDepthStruct(16, depths);
    }
}

double PsychOSVRefreshFromMode(XRRModeInfo *mode)
{
    // This routine is not really needed anymore. Only theoretically called from
    // SCREENResolutions.c, where it gets skipped anyway. Just define a no-op
    // implementation to avoid linker failure or ugly workarounds for the moment:
    return(0);
}

/*   PsychGetAllSupportedScreenSettings()
 *
 *     Queries the display system for a list of all supported display modes, ie. all valid combinations
 *     of resolution, pixeldepth and refresh rate. Allocates temporary arrays for storage of this list
 *     and returns it to the calling routine. This function is basically only used by Screen('Resolutions').
 */
int PsychGetAllSupportedScreenSettings(int screenNumber, int outputId, long** widths, long** heights, long** hz, long** bpp)
{
    int i, numPossibleModes;
    struct output_info* output_info = NULL;
    struct output_mode* mode = NULL;

    if (screenNumber >= numDisplays || screenNumber < 0) PsychErrorExit(PsychError_invalidScumber);

    PsychLockDisplay();
    ProcessWaylandEvents(screenNumber);
    PsychUnlockDisplay();

    if (outputId < 0) {
        // Count number of available modes:
        numPossibleModes = 0;
        output_info = displayOutputs[screenNumber];
        wl_list_for_each(mode, &output_info->modes, link) {
            numPossibleModes++;
        }

        // Allocate output arrays: These will get auto-released at exit from Screen():
        *widths  = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));
        *heights = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));
        *hz      = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));
        *bpp     = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));

        // Reiterate and fill all slots:
        i = 0;
        wl_list_for_each(mode, &output_info->modes, link) {
            if (i == numPossibleModes) break;

            (*widths)[i]  = (long) mode->width;
            (*heights)[i] = (long) mode->height;
            (*hz)[i]      = (long) ((float) mode->refresh / 1000.0 + 0.5);
            (*bpp)[i]     = (long) PsychGetScreenDepthValue(screenNumber);

            i++;
        }

        // Done:
        return(numPossibleModes);
    }

    if (outputId >= kPsychMaxPossibleCrtcs) PsychErrorExitMsg(PsychError_user, "Invalid output index provided! No such output for this screen!");

    // A bit dull: We only have one output per screen atm., as we go with the screen == output model for initial Wayland enablement:
    if (outputId >= 1 || outputId < 0) PsychErrorExitMsg(PsychError_user, "Invalid output index provided! No such output for this screen!");
    outputId = PsychScreenToHead(screenNumber, outputId);

    // Now have literally the same code as for the outputId < 0 case above,
    // Pretty dumb, but best we can do atm. given the screenid = outputid mapping.
    // TODO: Needs more thought...

    // Count number of available modes:
    numPossibleModes = 0;
    output_info = displayOutputs[screenNumber];
    wl_list_for_each(mode, &output_info->modes, link) {
        numPossibleModes++;
    }

    // Allocate output arrays: These will get auto-released at exit from Screen():
    *widths  = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));
    *heights = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));
    *hz      = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));
    *bpp     = (long*) PsychMallocTemp(numPossibleModes * sizeof(long));

    // Reiterate and fill all slots:
    i = 0;
    wl_list_for_each(mode, &output_info->modes, link) {
        if (i == numPossibleModes) break;

        (*widths)[i]  = (long) mode->width;
        (*heights)[i] = (long) mode->height;
        (*hz)[i]      = (long) ((float) mode->refresh / 1000.0 + 0.5);
        (*bpp)[i]     = (long) PsychGetScreenDepthValue(outputId);

        i++;
    }

    // Done:
    return(numPossibleModes);
}

/*
 * PsychGetCGModeFromVideoSetting()
 */
static psych_bool PsychGetCGModeFromVideoSetting(struct output_mode **cgMode, PsychScreenSettingsType *setting)
{
    struct output_info* output;
    struct output_mode *mode, *targetmode = NULL;

    // Extract parameters from setting struct:
    output = displayOutputs[setting->screenNumber];
    int width  = (int) PsychGetWidthFromRect(setting->rect);
    int height = (int) PsychGetHeightFromRect(setting->rect);
    int fps    = (int) (setting->nominalFrameRate + 0.5);

    if (PsychPrefStateGet_Verbosity() > 3) printf("PTB-INFO: Trying to validate/find closest video mode for requested spec: width = %i x height = %i, rate %i Hz.\n", width, height, fps);

    // Find matching mode:
    wl_list_for_each(mode, &output->modes, link) {
        if (PsychPrefStateGet_Verbosity() > 3) {
            printf("PTB-INFO: Testing against mode %p of resolution w x h = %i x %i with refresh rate: %i Hz.  ",
                   mode, mode->width, mode->height, (int) (((float) mode->refresh / 1000) + 0.5));
        }

        if ((width == mode->width) && (height == mode->height) && (fps == (int) (((float) mode->refresh / 1000) + 0.5))) {
            // Our requested mode is supported:
            if (PsychPrefStateGet_Verbosity() > 3) printf("--> Got it! Mode id %p. ", mode);
            targetmode = mode;
            break;
        }
        if (PsychPrefStateGet_Verbosity() > 3) printf("\n");
    }

    // Found valid settings?
    if (!targetmode) return(FALSE);

    *cgMode = targetmode;
    return(TRUE);
}

/*
 *    PsychCheckVideoSettings()
 *
 *    Check all available video display modes for the specified screen number and return true if the
 *    settings are valid and false otherwise.
 */
psych_bool PsychCheckVideoSettings(PsychScreenSettingsType *setting)
{
    struct output_mode *cgMode = NULL;
    return (PsychGetCGModeFromVideoSetting(&cgMode, setting));
}

/*
    PsychGetScreenDepth()

    The caller must allocate and initialize the depth struct.
*/
void PsychGetScreenDepth(int screenNumber, PsychDepthType *depth)
{
    if (screenNumber>=numDisplays || screenNumber < 0) PsychErrorExitMsg(PsychError_internal, "screenNumber is out of range"); //also checked within SCREENPixelSizes
// TODO: Make it work for real.
    PsychLockDisplay();
//     ProcessWaylandEvents(screenNumber);
//
//     if (displayCGIDs[screenNumber]) {
//         PsychAddValueToDepthStruct(DefaultDepth(displayCGIDs[screenNumber], PsychGetXScreenIdForScreen(screenNumber)), depth);
//     }
//     else {
        PsychAddValueToDepthStruct(24, depth);
//     }

    PsychUnlockDisplay();
}

float PsychGetNominalFramerate(int screenNumber)
{
    struct output_mode* mode = NULL;
    float vrefresh = 0;

    if (PsychPrefStateGet_ConserveVRAM() & kPsychIgnoreNominalFramerate) return(0);

    if (screenNumber >= numDisplays || screenNumber < 0)
        PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetNominalFramerate() is out of range");

    // Get current video mode for screenNumber:
    mode = PsychWaylandGetCurrentMode(screenNumber);
    if (mode == NULL) return(0);

    // Convert nominal refresh rate in milliHz to Hz and return it:
    return((float) mode->refresh / 1000.0);
}

float PsychSetNominalFramerate(int screenNumber, float requestedHz)
{
    // We start with a default vrefresh of zero, which means "couldn't query refresh from OS":
    float vrefresh = 0;

    if (screenNumber >= numDisplays || screenNumber < 0)
        PsychErrorExitMsg(PsychError_internal, "screenNumber is out of range");

    // Not available on non-X11:
    if (!displayOutputs[screenNumber]) return(0);

    // TODO: Could store target rate to set in a "pending new Hz" array,
    // and try to apply it via wl_shell fullscreen request...
    return(0);
}

/* Returns the physical display size in mm */
void PsychGetDisplaySize(int screenNumber, int *width, int *height)
{
    if (screenNumber >= numDisplays || screenNumber < 0)
        PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetDisplaySize() is out of range");

    if (!displayOutputs[screenNumber]) { *width = 0; *height = 0; return; }

    // Update XLib's view of this screens configuration:
    PsychLockDisplay();
    ProcessWaylandEvents(screenNumber);

    *width = (int) displayOutputs[screenNumber]->geometry.physical_width;
    *height = (int) displayOutputs[screenNumber]->geometry.physical_height;

    PsychUnlockDisplay();
}

// Get size of a video output in pixels - active scanout area in pixels:
void PsychGetScreenPixelSize(int screenNumber, long *width, long *height)
{
    struct output_mode* mode = NULL;

    if (screenNumber >= numDisplays || screenNumber < 0) PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetScreenPixelSize() is out of range");
    if (!displayOutputs[screenNumber]) PsychErrorExitMsg(PsychError_system, "Could not query screen size in PsychGetScreenPixelSize() for wanted screen");

    // Update XLib's view of this screens configuration:
    PsychLockDisplay();
    ProcessWaylandEvents(screenNumber);

    // Get size from current mode:
    // TODO: Make sure to take output scaling info provided via wl_output listener callback into account:
    mode = PsychWaylandGetCurrentMode(screenNumber);
    if (!mode) PsychErrorExitMsg(PsychError_system, "Could not query screen size in PsychGetScreenPixelSize() for wanted screen");

    *width = (int) mode->width;
    *height = (int) mode->height;

    PsychUnlockDisplay();
}

// Width and height of output in compositor space units (points?):
void PsychGetScreenSize(int screenNumber, long *width, long *height)
{
    struct output_mode* mode = NULL;

    if (screenNumber >= numDisplays || screenNumber < 0) PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychGetScreenSize() is out of range");
    if (!displayOutputs[screenNumber]) PsychErrorExitMsg(PsychError_system, "Could not query screen size in PsychGetScreenSize() for wanted screen");

    // Update XLib's view of this screens configuration:
    PsychLockDisplay();
    ProcessWaylandEvents(screenNumber);

    // Get size from current mode:
    // TODO: Make sure to take output scaling info provided via wl_output listener callback into account,
    // so we don't report wrong values on HiDPI / Retina style displays:
    mode = PsychWaylandGetCurrentMode(screenNumber);
    if (!mode) PsychErrorExitMsg(PsychError_system, "Could not query screen size in PsychGetScreenSize() for wanted screen");

    *width = mode->width;
    *height = mode->height;

    PsychUnlockDisplay();
}

// Global bounding rectangle of output in compositor space coordinates:
void PsychGetGlobalScreenRect(int screenNumber, double *rect)
{
    PsychGetScreenRect(screenNumber, rect);
    rect[kPsychLeft]   += (int) displayOutputs[screenNumber]->geometry.x;
    rect[kPsychRight]  += (int) displayOutputs[screenNumber]->geometry.x;
    rect[kPsychTop]    += (int) displayOutputs[screenNumber]->geometry.y;
    rect[kPsychBottom] += (int) displayOutputs[screenNumber]->geometry.y;;
}

// Bounding rectangle of output in compositor space units (points?):
void PsychGetScreenRect(int screenNumber, double *rect)
{
    long width, height;

    PsychGetScreenSize(screenNumber, &width, &height);
    rect[kPsychLeft] = 0;
    rect[kPsychTop] = 0;
    rect[kPsychRight] = (int) width;
    rect[kPsychBottom] = (int) height;
}

/*
    This is a place holder for a function which uncovers the number of dacbits.
    There seems to be no way to uncover the dacbits programatically.
    For now we just use pessimistic 0 bits to avoid false precision.
 */
int PsychGetDacBitsFromDisplay(int screenNumber)
{
    return(0);
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
//     XRRScreenResources *res = displayX11ScreenResources[screenNumber];

    // Need this later:
    PsychGetDisplaySize(screenNumber, &widthMM, &heightMM);

    if (has_xrandr_1_2 && (PsychScreenToHead(screenNumber, outputId) >= 0)) {
        PsychLockDisplay();
//         crtc_info = XRRGetCrtcInfo(dpy, res, res->crtcs[PsychScreenToHead(screenNumber, outputId)]);
        PsychUnlockDisplay();
    }
    else {
        // Failed!
        return(FALSE);
    }

//     // Disable auto-restore of screen settings - It would end badly...
//     displayOriginalCGSettingsValid[screenNumber] = FALSE;
//
//     // Find matching mode:
//     for (modeid = 0; modeid < res->nmode; modeid++) {
//         if (((int) res->modes[modeid].width == newWidth) && ((int) res->modes[modeid].height == newHeight) &&
//             (newHz == (int)(PsychOSVRefreshFromMode(&res->modes[modeid]) + 0.5))) {
//             break;
//         }
//     }
//
//     // Matching mode found for modesetting?
//     if (modeid < res->nmode) {
//         PsychLockDisplay();
//
//         // Assign default panning:
//         if (newX < 0) newX = crtc_info->x;
//         if (newY < 0) newY = crtc_info->y;
//
//         // Iterate over all outputs and compute the new screen bounding box:
//         maxw = maxh = 0;
//         for (output = 0; (PsychScreenToHead(screenNumber, output) >= 0) && (output < kPsychMaxPossibleCrtcs); output++) {
//             if (output == outputId) continue;
//
//             crtc_info2 = XRRGetCrtcInfo(dpy, res, res->crtcs[PsychScreenToHead(screenNumber, output)]);
//             if (crtc_info2->x + (int) crtc_info2->width > maxw) maxw = crtc_info2->x + (int) crtc_info2->width;
//             if (crtc_info2->y + (int) crtc_info2->height > maxh) maxh = crtc_info2->y + (int) crtc_info2->height;
//             XRRFreeCrtcInfo(crtc_info2);
//         }
//
//         // Incorporate our soon reconfigured crtc:
//         if (newX + newWidth  > maxw) maxw = newX + newWidth;
//         if (newY + newHeight > maxh) maxh = newY + newHeight;
//
//         // [0, 0, maxw, maxh] is the new bounding rectangle of the scanned out framebuffer. Set screen size accordingly:
//
//         // Prevent clients from getting confused by our config sequence:
//         // XGrabServer(dpy);
//
//         // Disable target crtc:
//         if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Disabling crtc %i.\n", outputId);
//         Status rc = XRRSetCrtcConfig(dpy, res, res->crtcs[PsychScreenToHead(screenNumber, outputId)], crtc_info->timestamp,
//                                     0, 0, None, RR_Rotate_0, NULL, 0);
//
//         // Resize screen: MK Don't! Skip this for now, use PsychSetScreenSettings() aka Screen('Resolution') to resize
//         // the screen without changing the crtc / output settings. More flexible...
//         // if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Resizing screen %i to %i x %i pixels.\n", screenNumber, maxw, maxh);
//         // XRRSetScreenSize(dpy, RootWindow(dpy, PsychGetXScreenIdForScreen(screenNumber)), maxw, maxh, widthMM, heightMM);
//
//         // Switch mode of target crtc and reenable it:
//         if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: Enabling crtc %i.\n", outputId);
//
//         crtc_info2 = XRRGetCrtcInfo(dpy, res, res->crtcs[PsychScreenToHead(screenNumber, outputId)]);
//         rc = XRRSetCrtcConfig(dpy, res, res->crtcs[PsychScreenToHead(screenNumber, outputId)], crtc_info2->timestamp,
//                                 newX, newY, res->modes[modeid].id, crtc_info->rotation,
//                                 crtc_info->outputs, crtc_info->noutput);
//         XRRFreeCrtcInfo(crtc_info);
//         XRRFreeCrtcInfo(crtc_info2);
//
//         // XUngrabServer(dpy);
//
//         // Make sure the screen change gets noticed by XLib:
//         ProcessWaylandEvents(screenNumber);
//
//         PsychUnlockDisplay();
//
//         return(TRUE);
//     } else {
//         XRRFreeCrtcInfo(crtc_info);
//         return(FALSE);
//     }
}

/*
    PsychSetScreenSettings()

    Accept a PsychScreenSettingsType structure holding a video mode and set the display mode accordingly.

    If we can not change the display settings because of a lock (set by open window or close window) then return false.

    SCREENOpenWindow should capture the display before it sets the video mode.  If it doesn't, then PsychSetScreenSettings will
    detect that and exit with an error.  SCREENClose should uncapture the display.

    The duties of SCREENOpenWindow are:
    -Lock the screen which serves the purpose of preventing changes in video setting with open Windows.
    -Capture the display which gives the application synchronous control of display parameters.

    TO DO: for 8-bit palletized mode there is probably more work to do.

*/

psych_bool PsychSetScreenSettings(psych_bool cacheSettings, PsychScreenSettingsType *settings)
{
    struct output_info* output;
    struct output_mode *mode;
    struct output_mode *cgMode = NULL;
    psych_bool isValid, isCaptured;

    if (settings->screenNumber >= numDisplays || settings->screenNumber < 0) PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychSetScreenSettings() is out of range");

    //Check to make sure that this display is captured, which OpenWindow should have done.  If it has not been done, then exit with an error.
    isCaptured = PsychIsScreenCaptured(settings->screenNumber);
    if (!isCaptured) PsychErrorExitMsg(PsychError_internal, "Attempt to change video settings without capturing the display");

    // Store the original display mode if this is the first time we have called this function.  The psychtoolbox will disregard changes in
    // the screen state made through the control panel after the Psychtoolbox was launched. That is, OpenWindow will by default continue to
    // open windows with finder settings which were in place at the first call of OpenWindow.  That's not intuitive, but not much of a problem
    // either.
    if (!displayOriginalCGSettingsValid[settings->screenNumber]) {
        PsychGetScreenSettings(settings->screenNumber, &displayOriginalCGSettings[settings->screenNumber]);
        displayOriginalCGSettingsValid[settings->screenNumber] = TRUE;
    }

    // Find video settings which correspond to settings as specified within by an abstracted psychsettings structure.
    isValid = PsychGetCGModeFromVideoSetting(&cgMode, settings);
    if (!isValid || !cgMode) {
        // This is an internal error because the caller is expected to check first.
        PsychErrorExitMsg(PsychError_user, "Attempt to set invalid video settings.");
    }

    // Change the display mode.
    PsychLockDisplay();

    // Extract parameters from setting struct:
    output = displayOutputs[settings->screenNumber];

    // Clear current mode flag on all modes:
    wl_list_for_each(mode, &output->modes, link) {
        mode->flags &= ~WL_OUTPUT_MODE_CURRENT;
    }

    // Then set current mode flag on our "new" mode:
    // This should guarantee that all our video mode reporting functions
    // will report the properties of our new current mode as current
    // settings, so our onscreen windows size will initialize to this
    // modes settings:
    cgMode->flags |= WL_OUTPUT_MODE_CURRENT;

    // Make sure the screen change gets noticed:
    ProcessWaylandEvents(settings->screenNumber);

    PsychUnlockDisplay();

    // Done:
    return(TRUE);
}

/*
    PsychRestoreScreenSettings()

    Restores video settings to the state set by the finder.  Returns TRUE if the settings can be restored or false if they
    can not be restored because a lock is in effect, which would mean that there are still open windows.

*/
psych_bool PsychRestoreScreenSettings(int screenNumber)
{
    struct output_info* output;
    struct output_mode *mode;
    struct output_mode *cgMode = NULL;
    psych_bool isValid, isCaptured;
    PsychScreenSettingsType *settings;

    if (screenNumber >= numDisplays || screenNumber < 0)
        PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychRestoreScreenSettings() is out of range"); //also checked within SCREENPixelSizes

    //Check to make sure that the original graphics settings were cached. If not,
    // it means that the settings were never changed, so we can just return true.
    if (!displayOriginalCGSettingsValid[screenNumber]) return(TRUE);

    //Check to make sure that this display is captured, which OpenWindow should have done.  If it has not been done, then exit with an error.
    isCaptured = PsychIsScreenCaptured(screenNumber);
    if (!isCaptured) PsychErrorExitMsg(PsychError_internal, "Attempt to change video settings without capturing the display");

    // Invalidate settings - we want a fresh game after restoring the resolution:
    displayOriginalCGSettingsValid[screenNumber] = FALSE;

    // End it here, for now. TODO: Need to first figure out how to do this
    // properly on Wayland. In principle would need to open a fullscreen
    // window, then close it again after applying restored settings below,
    // to force Wayland to restore original mode. But that seems cumbersome
    // and possibly not worth the extra trouble:
    return(TRUE);

    // Retrieve original screen settings which we should restore for this screen:
    settings = &displayOriginalCGSettings[screenNumber];

    //Find core graphics video settings which correspond to settings as specified withing by an abstracted psychsettings structure.
    isValid = PsychGetCGModeFromVideoSetting(&cgMode, settings);
    if (!isValid || !cgMode) {
        // This is an internal error because the caller is expected to check first.
        PsychErrorExitMsg(PsychError_user, "Attempt to restore invalid video settings.");
    }

    // Change the display mode.
    PsychLockDisplay();

    // Extract parameters from setting struct:
    output = displayOutputs[screenNumber];

    // Clear current mode flag on all modes:
    wl_list_for_each(mode, &output->modes, link) {
        mode->flags &= ~WL_OUTPUT_MODE_CURRENT;
    }

    // Then set current mode flag on our "new" mode:
    // This should guarantee that all our video mode reporting functions
    // will report the properties of our new current mode as current
    // settings, so our onscreen windows size will initialize to this
    // modes settings:
    cgMode->flags |= WL_OUTPUT_MODE_CURRENT;

    // Make sure the screen change gets noticed:
    ProcessWaylandEvents(screenNumber);

    PsychUnlockDisplay();

    // Done:
    return(TRUE);
}

void PsychOSDefineWaylandCursor(int screenNumber, int deviceId, const char* cursorName)
{
    struct seat_info *seat;

    // Check for valid screenNumber, although it will be ignored on Wayland:
    if ((screenNumber >= numDisplays) || (screenNumber < 0)) PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychOSDefineWaylandCursor() is out of range");

    // deviceId -1 means "auto-detected default pointer". Simply use the first found pointer device:
    if (deviceId < 0) {
        for (deviceId = 0; deviceId < numInputDevices; deviceId++) {
            if (waylandInputDevices[deviceId] && (waylandInputDevices[deviceId]->capabilities & WL_SEAT_CAPABILITY_POINTER))
                break;
        }
    }

    // Outside valid range for input devices?
    if (deviceId < 0 || deviceId >= numInputDevices) PsychErrorExitMsg(PsychError_user, "Invalid 'mouseIndex' provided. No such cursor pointer.");

    // No device under that id, or device isn't a pointer?
    if (!waylandInputDevices[deviceId] || !(waylandInputDevices[deviceId]->capabilities & WL_SEAT_CAPABILITY_POINTER)) {
        PsychErrorExitMsg(PsychError_user, "Invalid 'mouseIndex' provided. No such cursor pointer.");
    }

    seat = waylandInputDevices[deviceId];

    if (wayland_cursor_theme) {
        // Assign initial default cursor: The classic left-tilted arrow pointer:
        seat->current_cursor = wl_cursor_theme_get_cursor(wayland_cursor_theme, cursorName);
        if (!seat->current_cursor) {
            if (PsychPrefStateGet_Verbosity() > 0) {
                printf("PTB-ERROR: Loading Wayland cursor from theme for pointer device %i [seat %p] failed. Cursor support unavailable on this seat!\n", deviceId, seat);
            }
        }

        // Enforce update of cursor surface image:
        PsychShowCursor(screenNumber, deviceId);
    }

    return;
}

void PsychHideCursor(int screenNumber, int deviceIdx)
{
    struct seat_info *seat;

    // Check for valid screenNumber, although it will be ignored on Wayland:
    if ((screenNumber >= numDisplays) || (screenNumber < 0)) PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychHideCursor() is out of range");

    // deviceIdx -1 means "auto-detected default pointer". Simply use the first found pointer device:
    if (deviceIdx < 0) {
        for (deviceIdx = 0; deviceIdx < numInputDevices; deviceIdx++) {
            if (waylandInputDevices[deviceIdx] && (waylandInputDevices[deviceIdx]->capabilities & WL_SEAT_CAPABILITY_POINTER))
                break;
        }
    }

    // Outside valid range for input devices?
    if (deviceIdx < 0 || deviceIdx >= numInputDevices) PsychErrorExitMsg(PsychError_user, "Invalid 'mouseIndex' provided. No such cursor pointer.");

    // No device under that id, or device isn't a pointer?
    if (!waylandInputDevices[deviceIdx] || !(waylandInputDevices[deviceIdx]->capabilities & WL_SEAT_CAPABILITY_POINTER)) {
        PsychErrorExitMsg(PsychError_user, "Invalid 'mouseIndex' provided. No such cursor pointer.");
    }

    seat = waylandInputDevices[deviceIdx];
    seat->cursor_hidden = TRUE;

    if (seat->cursor_surface) {
        wl_surface_attach(seat->cursor_surface, NULL, 0, 0);
        wl_surface_commit(seat->cursor_surface);
    }

    ProcessWaylandEvents(screenNumber);

    return;
}

void PsychShowCursor(int screenNumber, int deviceIdx)
{
    struct wl_buffer *buffer;
    struct wl_cursor_image *image;
    struct seat_info *seat;

    // Check for valid screenNumber, although it will be ignored on Wayland:
    if ((screenNumber >= numDisplays) || (screenNumber < 0)) PsychErrorExitMsg(PsychError_internal, "screenNumber passed to PsychShowCursor() is out of range");

    // deviceIdx -1 means "auto-detected default pointer". Simply use the first found pointer device:
    if (deviceIdx < 0) {
        for (deviceIdx = 0; deviceIdx < numInputDevices; deviceIdx++) {
            if (waylandInputDevices[deviceIdx] && (waylandInputDevices[deviceIdx]->capabilities & WL_SEAT_CAPABILITY_POINTER))
                break;
        }
    }

    // Outside valid range for input devices?
    if (deviceIdx < 0 || deviceIdx >= numInputDevices) PsychErrorExitMsg(PsychError_user, "Invalid 'mouseIndex' provided. No such cursor pointer.");

    // No device under that id, or device isn't a pointer?
    if (!waylandInputDevices[deviceIdx] || !(waylandInputDevices[deviceIdx]->capabilities & WL_SEAT_CAPABILITY_POINTER)) {
        PsychErrorExitMsg(PsychError_user, "Invalid 'mouseIndex' provided. No such cursor pointer.");
    }

    seat = waylandInputDevices[deviceIdx];
    seat->cursor_hidden = FALSE;

    // Make cursor visible for this pointer: Any cursor assigned at the moment?
    if (seat->current_cursor && seat->cursor_surface) {
        // Yes, set it. Get first cursor image for current cursor:
        // There could be multiple, e.g., for animated cursors, but we only
        // care about one image - We don't do animated cursors:
        image = seat->current_cursor->images[0];

        // Get wl_buffer with cursor image. This buffer doesn't need to
        // be released, as it is memory managed by the wayland_cursor_theme:
        buffer = wl_cursor_image_get_buffer(image);

        // Assign cursor image to the cursor surface:
        wl_surface_attach(seat->cursor_surface, buffer, seat->hotspot_x - image->hotspot_x, seat->hotspot_y - image->hotspot_y);
        wl_surface_damage(seat->cursor_surface, 0, 0, image->width, image->height);
        wl_surface_commit(seat->cursor_surface);
        seat->hotspot_x = image->hotspot_x;
        seat->hotspot_y = image->hotspot_y;
    }

    ProcessWaylandEvents(screenNumber);

    return;
}

void PsychPositionCursor(int screenNumber, int x, int y, int deviceIdx)
{
    // Not available on Wayland, as of version 1.6 - January 2015:
    if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: SetCursor() request ignored, as this isn't currently supportable on Wayland.\n");
    return;
}

/*
    PsychReadNormalizedGammaTable()
*/
void PsychReadNormalizedGammaTable(int screenNumber, int outputId, int *numEntries, float **redTable, float **greenTable, float **blueTable)
{
    static float localRed[MAX_GAMMALUT_SIZE], localGreen[MAX_GAMMALUT_SIZE], localBlue[MAX_GAMMALUT_SIZE];
    int i, n;
    double maxval = 0.0;

    (void) outputId;

    // Initial assumption: Failed.
    *numEntries = 0;
    n = 0;

    if (!colord_client || !displayOutputs[screenNumber]->colord_device) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: ReadNormalizedGammatable: No color management for screen %i.\n", screenNumber);
        return;
    }

    CdDevice *device = displayOutputs[screenNumber]->colord_device;
    CdProfile *profile = cd_device_get_default_profile (device);
    if (!profile) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: ReadNormalizedGammatable: No color profile assigned for screen %i.\n", screenNumber);
        return;
    }

    if (!cd_profile_connect_sync(profile, NULL, NULL)) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: ReadNormalizedGammatable: Could not connect to color profile assigned for screen %i.\n", screenNumber);
        goto out_profile;
    }

    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: ReadNormalizedGammatable: Profile for screen %i is '%s'.\n", screenNumber, cd_profile_get_title(profile));
    if (!cd_profile_get_has_vcgt(profile)) {
        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: ReadNormalizedGammatable: Profile for screen %i does not have a gamma table. Synthesizing one...\n", screenNumber);
    }

    CdIcc *icc = cd_profile_load_icc(profile, CD_ICC_LOAD_FLAGS_ALL, NULL, NULL);
    if (!icc) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: ReadNormalizedGammatable: Could not load icc for screen %i.\n", screenNumber);
        goto out_profile;
    }

    GPtrArray *clut = cd_profile_get_has_vcgt(profile) ? cd_icc_get_vcgt(icc, 256, NULL) : cd_icc_get_response(icc, 256, NULL);
    if (!clut) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: ReadNormalizedGammatable: Could not get clut for screen %i.\n", screenNumber);
        goto out_icc;
    }

    n = 256;
    for (i = 0; i < n; i++) {
        CdColorRGB *rgb = g_ptr_array_index(clut, i);
        localRed[i]   = (float) rgb->R;
        localGreen[i] = (float) rgb->G;
        localBlue[i]  = (float) rgb->B;
        maxval = (maxval > localRed[i]) ? maxval : localRed[i];
        maxval = (maxval > localGreen[i]) ? maxval : localGreen[i];
        maxval = (maxval > localBlue[i]) ? maxval : localBlue[i];
    }

    g_ptr_array_free(clut, TRUE);

    if (!cd_profile_get_has_vcgt(profile) && (maxval > 0.0)) {
        // Synthesized pseudo gamma curves from cd_icc_get_response()
        // can have values > 1.0, so renormalize them to make their
        // maximum value map to 1.0:
        for (i = 0; i < n; i++) {
            localRed[i] /= maxval;
            localGreen[i] /= maxval;
            localBlue[i] /= maxval;
        }
    }

out_icc:
    g_object_unref(icc);

out_profile:
    g_object_unref(profile);

    // Assign output lut's:
    *redTable=localRed; *greenTable=localGreen; *blueTable=localBlue;

    // Assign size of LUT's::
    *numEntries = n;

    return;
}

unsigned int PsychLoadNormalizedGammaTable(int screenNumber, int outputId, int numEntries, float *redTable, float *greenTable, float *blueTable)
{
    int i, j, rc;
    CdColorRGB *rgb;
    GPtrArray *clut;
    CdIcc *icc;
    CdProfile *profile = NULL;
    CdProfile *oldprofile1 = NULL;
    CdProfile *oldprofile2 = NULL;
    GError *error = NULL;
    cmsContext context;
    cmsHPROFILE hProfile;
    char screenName[32];
    gchar *filename = NULL;
    GFile *file = NULL;
    int targetindex = 1;

    // Initial assumption: Failure.
    rc = 0;

    if (numEntries != 256)
        PsychErrorExitMsg(PsychError_user, "Provided gamma table has wrong number of slots!");

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

    // Color management supported on this screen?
    if (!colord_client || !displayOutputs[screenNumber]->colord_device) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: PsychLoadNormalizedGammaTable: No gamma table support on screen %i.\n", screenNumber);
        return(0);
    }

    // Does a profile with that name already exist? If so, remove it:
    sprintf(screenName, "icc-PTBScreen %i 1", screenNumber);
    oldprofile1 = cd_client_find_profile_by_property_sync(colord_client, "PTBName", screenName, NULL, NULL);
    if (oldprofile1 && !cd_profile_connect_sync(oldprofile1, NULL, NULL)) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychLoadNormalizedGammaTable: Could not connect to old color profile1 assigned for screen %i.\n", screenNumber);
        goto cleanup_profile;
    }

    sprintf(screenName, "icc-PTBScreen %i 2", screenNumber);
    oldprofile2 = cd_client_find_profile_by_property_sync(colord_client, "PTBName", screenName, NULL, NULL);
    if (oldprofile2 && !cd_profile_connect_sync(oldprofile2, NULL, NULL)) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychLoadNormalizedGammaTable: Could not connect to old color profile2 assigned for screen %i.\n", screenNumber);
        goto cleanup_profile;
    }

    profile = cd_device_get_default_profile(displayOutputs[screenNumber]->colord_device);
    if (profile && !cd_profile_connect_sync(profile, NULL, NULL)) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychLoadNormalizedGammaTable: Could not connect to current profile assigned for screen %i.\n", screenNumber);
        goto cleanup_profile;
    }

    // Current profile == oldprofile2?
    if (!profile || !oldprofile2 || !cd_profile_equal(profile, oldprofile2)) {
        // Yes. That means oldprofile2, if any, is unused and can be deleted now,
        // then recycled as next profile:
        targetindex = 2;
        sprintf(screenName, "icc-PTBScreen %i %i", screenNumber, targetindex);

        if (oldprofile2 && !cd_client_delete_profile_sync(colord_client, oldprofile2, NULL, &error)) {
            if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychLoadNormalizedGammaTable: Failed to remove old PTB profile 2 for screen %i: %s\n", screenNumber, error->message);
            g_error_free(error);
            goto cleanup_profile;
        }

        // Release our test profile:
        if (oldprofile2) g_object_unref(oldprofile2);
        oldprofile2 = NULL;
        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: PsychLoadNormalizedGammaTable: Old PTB profile 2 for screen %i removed.\n", screenNumber);
    }

    // Current profile == oldprofile1?
    if (!profile || !oldprofile1 || !cd_profile_equal(profile, oldprofile1)) {
        // Yes. That means oldprofile1, if any, is unused and can be deleted now,
        // then recycled as next profile:
        targetindex = 1;
        sprintf(screenName, "icc-PTBScreen %i %i", screenNumber, targetindex);

        if (oldprofile1 && !cd_client_delete_profile_sync(colord_client, oldprofile1, NULL, &error)) {
            if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychLoadNormalizedGammaTable: Failed to remove old PTB profile 1 for screen %i: %s\n", screenNumber, error->message);
            g_error_free(error);
            goto cleanup_profile;
        }

        // Release our test profile:
        if (oldprofile1) g_object_unref(oldprofile1);
        oldprofile1 = NULL;
        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: PsychLoadNormalizedGammaTable: Old PTB profile 1 for screen %i removed.\n", screenNumber);
    }

    g_object_unref(profile);
    profile = NULL;

    // Build clut from input arrays:
    clut = g_ptr_array_new_with_free_func((GDestroyNotify) cd_color_rgb_free);
    for (i = 0; i < numEntries; i++) {
        rgb = cd_color_rgb_new();
        cd_color_rgb_set(rgb, redTable[i], greenTable[i], blueTable[i]);
        g_ptr_array_add(clut, rgb);
    }

    // Create a new icc object with clut:
    icc = cd_icc_new();

    // Associate Little-CMS 2 color management system with it and create a sRGB profile
    // as a starting point:
    context = NULL; // In later versions of lcms2: cmsCreateContext(NULL, NULL);
    hProfile = cmsCreate_sRGBProfile();
    if (!hProfile || !cd_icc_load_handle(icc, hProfile, CD_ICC_LOAD_FLAGS_ALL,  &error)) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychLoadNormalizedGammaTable: Failed loading handle into icc object for screen %i. [hProfile=%p]: %s\n", screenNumber, hProfile, error->message);
        g_error_free(error);
        goto cleanup_iccset;
    }

    // Assign our clut as gamma table:
    if (!cd_icc_set_vcgt(icc, clut, &error)) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychLoadNormalizedGammaTable: Failed setting gamma table in icc object for screen %i: %s\n", screenNumber, error->message);
        g_error_free(error);
        goto cleanup_iccset;
    }

    // Add metainfo:
    cd_icc_set_colorspace(icc, CD_COLORSPACE_RGB);
    cd_icc_set_kind(icc, CD_PROFILE_KIND_DISPLAY_DEVICE);
    cd_icc_set_description(icc, NULL, screenName);
    cd_icc_add_metadata(icc, "PTBName", screenName);
    cd_icc_set_manufacturer(icc, NULL, displayOutputs[screenNumber]->geometry.make);
    cd_icc_set_model(icc, NULL, displayOutputs[screenNumber]->geometry.model);
    cd_icc_set_copyright(icc, NULL, "This profile is free of copyright, use as you wish.");

    // Save it to a ICC file:
    sprintf(screenName, "icc-PTBScreen_%i_%i.icc", screenNumber, targetindex);
    filename = g_build_filename(PsychRuntimeGetPsychtoolboxRoot(TRUE), screenName, NULL);
    file = g_file_new_for_path(filename);
    if (!cd_icc_save_file(icc, file, CD_ICC_SAVE_FLAGS_NONE, NULL,  &error)) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychLoadNormalizedGammaTable: Failed saving icc profile %i for screen %i to file '%s': %s\n", targetindex, screenNumber, filename, error->message);
        g_error_free(error);
        goto cleanup_file;
    }
    else {
        if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-DEBUG: PsychLoadNormalizedGammaTable: Saved icc profile %i for screen %i to file '%s'\n", targetindex, screenNumber, filename);
    }

    PsychYieldIntervalSeconds(0.1);

    // Import it into colord's database and the standard user ICC folder, get
    // a profile handle back:
    profile = cd_client_import_profile_sync(colord_client, file, NULL, &error);
    if (!profile) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychLoadNormalizedGammaTable: Failed import of icc profile %i for screen %i from file '%s': %s\n", targetindex, screenNumber, filename, error->message);
        g_error_free(error);
        goto cleanup_file;
    }

    // Now we need to add and set our profile as default profile for the device:
    PsychYieldIntervalSeconds(0.1);

    // First add the new instance of the profile to our target device:
    if (!cd_device_add_profile_sync(displayOutputs[screenNumber]->colord_device, CD_DEVICE_RELATION_HARD, profile, NULL, &error)) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychLoadNormalizedGammaTable: Failed to add icc profile %i to color device for screen %i: %s\n", targetindex, screenNumber, error->message);
        g_error_free(error);
        goto cleanup_file;
    }

    PsychYieldIntervalSeconds(0.1);

    if (!cd_device_make_profile_default_sync(displayOutputs[screenNumber]->colord_device, profile, NULL, &error)) {
        if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: PsychLoadNormalizedGammaTable: Failed to activate icc profile %i on color device for screen %i: %s\n", targetindex, screenNumber, error->message);
        g_error_free(error);
        goto cleanup_file;
    }

    PsychYieldIntervalSeconds(0.1);

    // Successfully done:
    if (PsychPrefStateGet_Verbosity() > 4) printf("PTB-INFO: PsychLoadNormalizedGammaTable: Gamma table updated via profile %i update for screen %i.\n", targetindex, screenNumber);
    rc = 1;

cleanup_file:
    g_object_unref(file);
    g_free(filename);

cleanup_iccset:
    g_object_unref(icc);
    g_ptr_array_free(clut, TRUE);

cleanup_profile:
    if (profile) g_object_unref(profile);
    if (oldprofile1) g_object_unref(oldprofile1);
    if (oldprofile2) g_object_unref(oldprofile2);

    // Return final status:
    return(rc);
}

/*
 * Enable or disable profiling inhibition for a screen on colord. Enabling this
 * will load an identity gamma table and prevent any further gamma table changes.
 * This is meant to make display calibration foolproof, but is also useful for us
 * to automatically enable pixel-perfect framebuffer passthrough for devices like
 * Bits+ or Datapixx etc. - Assuming it works, of course, and we don't run into the
 * problem of "non-identity" identity gamma tables. This needs to be verified on a
 * per driver + per GPU basis.
 */
psych_bool PsychWaylandProfilingInhibit(int screenNumber, psych_bool enableInhibit)
{
    GError *error = NULL;

    // Color management supported on this screen?
    if (!colord_client || !displayOutputs[screenNumber]->colord_device) {
        if (PsychPrefStateGet_Verbosity() > 1) printf("PTB-WARNING: PsychProfilingInhibit: No colord support on screen %i.\n", screenNumber);
        return(FALSE);
    }

    // Requested state already set? No op, if so:
    if (enableInhibit == displayProfilingInhibit[screenNumber])
        return(TRUE);

    if (enableInhibit) {
        if (!cd_device_profiling_inhibit_sync(displayOutputs[screenNumber]->colord_device, NULL, &error)) {
            if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Failed to enable profiling inhibition for screen %i: %s\n", screenNumber, error->message);
            g_error_free(error);
            return(FALSE);
        }
    }
    else {
        if (!cd_device_profiling_uninhibit_sync(displayOutputs[screenNumber]->colord_device, NULL, &error)) {
            if (PsychPrefStateGet_Verbosity() > 0) printf("PTB-ERROR: Failed to disable profiling inhibition for screen %i: %s\n", screenNumber, error->message);
            g_error_free(error);
            return(FALSE);
        }
    }

    // Done.
    displayProfilingInhibit[screenNumber] = enableInhibit;
    if (PsychPrefStateGet_Verbosity() > 2) printf("PTB-INFO: %s profiling inhibition (= identity pixel passthrough) for screen %i.\n", enableInhibit ? "Enabled" : "Disabled", screenNumber);

    return(TRUE);
}

// Return true (non-zero) if a desktop compositor is likely active on screen 'screenNumber':
int PsychOSIsDWMEnabled(int screenNumber)
{
    // This is Wayland. There is always a compositor active, but what this
    // function actually returns is if a compositor is active which is
    // incompatible with good visual presentation timing and timestamping.
    // Wayland should not have this problem, at least with the presentation
    // extension available and enabled, so return an optimistic "No" if the
    // extension is available and enabled:
    return((wayland_pres) ? 0 : 1);
}

// !PTB_USE_WAYLAND
#endif
