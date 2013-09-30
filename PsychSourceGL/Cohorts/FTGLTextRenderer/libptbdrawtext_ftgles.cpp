/*
 * libptbdrawtext_ftgles.cpp: Source of libptbdrawtext_ftgles dynamic shared library plugin
 * for Screen('DrawText') text rendering via external engine. This is the version for use
 * with OpenGL-ES 1.1 instead of desktop OpenGL.
 *
 * This works only on GNU/Linux.
 *
 * Features:
 *
 * - Texture mapped renderer, like on OS/X with ATSU Drawtext.
 * - Fast due to use of a texture atlas for fast glyph recycling.
 * - Good text layouting.
 * - Supports all Freetype-2 supported fonts, e.g., vectorgraphics TrueType fonts.
 * - Anti-Aliased drawing via Alpha-Blending.
 * - Full Unicode support.
 * - Text measuring support, e.g., bounding boxes.
 *
 * Requires:
 *
 * Modified FTGLES library forked from https://github.com/cdave1/ftgles, modified version
 * in https://github.com/kleinerm/ftgles in the linuxptb branch. The library is licensed
 * to us under a MIT style license.
 *
 * The library must be ./configured as follows to allow linking the plugin with the static lib:
 * ./configure --with-pic
 *
 * Requires the FTGL subfolder of that lib (ftgles/ftgles/src/FTGL) in our build directory, for
 * access to header definitions.
 *
 * Requires libfontconfig and libFreeType-2.
 *
 * Building the plugin for Linux:
 * Intel-64 Bit:
 * g++ -g -fPIC -I. -I/usr/include/ -I/usr/include/freetype2/ -I/usr/include/GLES -L/usr/lib libptbdrawtext_ftgles.cpp -lGLESv1_CM -lfontconfig -lfreetype -pie -shared -Wl,-Bsymbolic -Wl,-Bsymbolic-functions -Wl,--version-script=linuxexportlist.txt -o ../../../Psychtoolbox/PsychBasic/PsychPlugins/libptbdrawtext_ftgles64.so.1 libftgles_intel64.a
 *
 * ARM-32 Bit:
 * g++ -g -fPIC -I. -I/usr/include/ -I/usr/include/freetype2/ -I/usr/include/GLES -L/usr/lib libptbdrawtext_ftgles.cpp -lGLESv1_CM -lfontconfig -lfreetype -pie -shared -Wl,-Bsymbolic -Wl,-Bsymbolic-functions -Wl,--version-script=linuxexportlist.txt -o ../../../Psychtoolbox/PsychBasic/PsychPlugins/libptbdrawtext_ftgles_arm.so.1 libftgles_arm32.a
 *
 * libptbdrawtext_ftgles is copyright (c) 2013 by Mario Kleiner.
 * It is licensed to you as follows:
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

// Standard libs for debug output and unicode character string handling:
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <wchar.h>

// Include GLFTES and OpenGL-ES stuff, as well as GL on GLES emulation glue:
#include <gl.h>
#include <FTGL/ftgles.h>
#include <FTGL/ftglesGlue.h>

// Include fontconfig stuff:
#include "fontconfig/fontconfig.h"
#include "fontconfig/fcfreetype.h"


// Maximum number of cacheable font settings, combined over all onscreen windows
// font family names, sizes, styles and anti-aliasing settings. E.g., a setting of 40 means
// you could quickly switch between up to 40 different combinations of font name, size,
// style and anti-aliasing setting without speed penalty, due to caching. This is a tradeoff
// between memory consumption and switchability, so choose wisely.
// If usercode requests more than those number, caching for the given context (window)
// will use LRU replacement to recycle slots and switching between fonts in that context
// (aka window) will slow down by about a factor of 10x to 15x, unless the LRU will do a good
// job of finding a new steady state. It will still be reasonably fast though.
#define MAX_CACHE_SLOTS 40

// Guarantee correct text rendering on at least 10 different contexts (aka onscreen windows).
// This limit makes sure that as long as there aren't more than 10 contexts, none of them can
// starve of resources because another one is too greedy. When this limit is reached, the cache
// turns to agressive LRU replacement instead of allocating new slots for new font configurations.
// Note: 10 is a good number, because Screen() allows for maximum 10 screens, so therefore at most
// 10 fullscreen onscreen windows, so guaranteeing for 10 onscreen windows should be good enough.
#define MIN_GUARANTEED_CONTEXTS 10

unsigned int nowtime = 0;
unsigned int hitcount = 0;
unsigned int _verbosity = 2;
bool _firstCall = true;
bool _useOwnFontmapper = false;
unsigned _mapperFlags;
int _antiAliasing = 1;
double _vxs, _vys, _vw, _vh;
char _fontName[FILENAME_MAX] = { 0 };
unsigned int _fontStyle = 0;
double _fontSize = 0.0;
GLfloat _fgcolor[4];
GLfloat _bgcolor[4];

typedef struct fontCacheItem_t {
    int contextId;
    unsigned int timestamp;
    int antiAliasing;
    char fontName[FILENAME_MAX];
    char fontRealName[FILENAME_MAX];
    unsigned int fontStyle;
    double fontSize;
    FTFont *faceT;
    FTFont *faceM;
    FT_Face ft_face;
} fontCacheItem;
fontCacheItem cache[MAX_CACHE_SLOTS];

// Internal helper:
// Convert unicode text in 'text', with length 'textLen' into wchar_t text string and
// store it in internal global uniCodeText[] variable for use. Filter out unprintable
// characters:
#define MAX_TEXT_LENGTH 10000
wchar_t uniCodeText[MAX_TEXT_LENGTH + 1];
int UnicodeToWchar(int textInLen, double *text)
{
    int i, j = 0;

    // Synthesize Unicode wchar_t string from double vector:
    memset(&uniCodeText[0], 0, sizeof(uniCodeText));
    textInLen = (textInLen > MAX_TEXT_LENGTH) ? MAX_TEXT_LENGTH : textInLen;
    for (i = 0; i < textInLen; i++) {
        // Filter out unprintable characters like 10 (LF) and 13 (CR):
        if ((text[i] != 10) && (text[i] != 13)) uniCodeText[j++] = (wchar_t) text[i];
    }

    // Return new textLen after potential filtering of unprintable chars
    return j;
}

extern "C" {

    int PsychInitText(void);
    int PsychShutdownText(int context);
    int PsychRebuildFont(fontCacheItem* fi);
    int PsychSetTextFont(int context, const char* fontName);
    const char* PsychGetTextFont(int context);
    int PsychSetTextStyle(int context, unsigned int fontStyle);
    int PsychSetTextSize(int context, double fontSize);
    void PsychSetTextFGColor(int context, double* color);
    void PsychSetTextBGColor(int context, double* color);
    void PsychSetTextUseFontmapper(unsigned int useMapper, unsigned int mapperFlags);
    void PsychSetTextViewPort(int context, double xs, double ys, double w, double h);
    int PsychDrawText(int context, double xStart, double yStart, int textLen, double* text);
    int PsychMeasureText(int context, int textLen, double* text, float* xmin, float* ymin, float* xmax, float* ymax);
    void PsychSetTextVerbosity(unsigned int verbosity);
    void PsychSetTextAntiAliasing(int context, int antiAliasing);

    fontCacheItem* getForContext(int contextId)
    {
        int lruslotid = -1;
        unsigned int lruage = 0;
        int freeslot = -1;
        int freecount = 0;
        fontCacheItem* fi = NULL;

        // Update running "time" for LRU replacement:
        nowtime++;

        // Search for matching cached font object:
        for (int i = 0; i < MAX_CACHE_SLOTS; i++) {
            // Only look at slots for requested contextId:
            if (contextId == cache[i].contextId) {
                // This one is for our contextId.
                fi = &(cache[i]);

                //  LRU updating, in case we need to replace:
                if (nowtime - fi->timestamp > lruage) {
                    lruslotid = i;
                    lruage = nowtime - fi->timestamp;
                }

                // Matching attributes for current requested attributes?
                // We match requested fontName against both, the originally requested fontName for this cache slot,
                // and the real effective fontRealName that libfontconfig actually gave us. Otherwise, as fontRealName
                // is returned to Screen(), we could get into a funny loop which causes false cache misses if the loaded font
                // doesn't match exactly the required one.
                if ((fi->antiAliasing == _antiAliasing) && (fi->fontStyle == _fontStyle) && (fi->fontSize == _fontSize) &&
                    ((strcmp(fi->fontName, _fontName) == 0) || (strcmp(fi->fontRealName, _fontName) == 0))) {
                    // Match! We have cached OGLFT font objects for this font on this context. Return them:
                    hitcount++;

                    if (_verbosity > 3) fprintf(stderr, "libptbdrawtext_ftgles: Cache hit for contextId %i at slot %i. Hit ratio is %f%%\n", contextId, i, (double) hitcount / (double) nowtime * 100);

                    // Update last access timestamp for LRU:
                    fi->timestamp = nowtime;

                    return(fi);
                }
            }
            else if (cache[i].contextId == -1) {
                if (freeslot == -1) freeslot = i;
                freecount++;
            }
        }

        // No match. We need to (re-)create a matching object.

        // Free slots available?
        if ((freeslot >= 0) && ((lruslotid == -1) || (freecount > MIN_GUARANTEED_CONTEXTS))) {
            // Yes. Fill it with new font object of matching properties:
            fi = &(cache[freeslot]);
            if (_verbosity > 3) fprintf(stderr, "libptbdrawtext_ftgles: Nothing cached for contextId %i. Using new slot %i. %i free slots remaining.\n", contextId, freeslot, freecount);
        }
        else if (lruslotid >= 0) {
            // No. Overwrite least recently used font object for current contextId:
            fi = &(cache[lruslotid]);
            if (_verbosity > 3) fprintf(stderr, "libptbdrawtext_ftgles: Nothing cached for contextId %i but cache full. LRU replacing slot %i, age %i. %i free slots remaining.\n", contextId, lruslotid, lruage, freecount);
        }
        else {
            // Cache full, with no possibility to even LRU replace on this new context (aka window). Game over!
            if (_verbosity > 0) fprintf(stderr, "libptbdrawtext_ftgles: ERROR when trying to setup new drawtext context %i: Font cache full with no way to free up resources! Text drawing will fail!\n", contextId);
            return(NULL);
        }

        // Rebuild or create font objects for slot fi. Return NULL on failure:
        if (PsychRebuildFont(fi)) return(NULL);

        // Update tags:
        fi->contextId = contextId;
        fi->antiAliasing = _antiAliasing;
        fi->fontStyle = _fontStyle;
        fi->fontSize = _fontSize;
        strcpy(fi->fontName, _fontName);

        // Update last access timestamp for LRU:
        fi->timestamp = nowtime;

        // Return new font objects:
        return(fi);
    }

    void PsychSetTextVerbosity(unsigned int verbosity)
    {
        _verbosity = verbosity;
        return;
    }

    void PsychSetTextAntiAliasing(int context, int antiAliasing)
    {
        static bool firstAntiAliasCall = true;

        // FTGLES as of April 2013 can't support non-anti-aliased rendering. Give a one-time warning:
        if ((antiAliasing == 0) && (_verbosity > 0) && firstAntiAliasCall) {
            firstAntiAliasCall = false;
            fprintf(stderr, "libptbdrawtext_ftgles: WARNING: Asked to disable text anti-aliasing for Screen('DrawText'), but this is not yet supported on OpenGL-ES! Ignored.\n");
        }

        _antiAliasing = antiAliasing;
        return;
    }

    int PsychRebuildFont(fontCacheItem* fi)
    {
        int faceIndex = 0;
        char fontFileName[FILENAME_MAX] = { 0 };

        // Destroy old font object, if any:
        if (fi->faceT || fi->faceM) {
            // Delete OGLFT face object:
            if (fi->faceT) delete(fi->faceT);
            fi->faceT = NULL;

            if (fi->faceM) delete(fi->faceM);
            fi->faceM = NULL;

            if (_verbosity > 3) fprintf(stderr, "libptbdrawtext_ftgles: Destroying old font face...\n");

            // Delete underlying FreeType representation:
            FT_Done_Face(fi->ft_face);
            fi->ft_face = NULL;
        }

        if (_useOwnFontmapper) {
            FcResult result = FcResultMatch;	// Must init this due to weirdness in libfontconfig...
            FcPattern *target = NULL;

            if (_fontName[0] == '-') {
                // _fontName starts with a '-' dash: This is not a simple font family string but a special
                // fontspec string in FontConfig's special format. It contains many possible required font
                // properties encoded in the string. Parse it into a font matching pattern:
                target = FcNameParse((FcChar8 *) & (_fontName[1]));

                // Need to manually add the current _fontSize, otherwise inconsistent stuff may happen:
                FcPatternAddDouble(target, FC_PIXEL_SIZE, _fontSize);
            } else {
                // _fontName contains only font family name: Build matching pattern based on _fontSize and
                // the flags provided in _fontStyle, according to the conventions in Psychtoolbox Screen('TextStyle'):
                target = FcPatternBuild(0, FC_FAMILY, FcTypeString, _fontName, FC_PIXEL_SIZE, FcTypeDouble, _fontSize,
                                        FC_WEIGHT, FcTypeInteger, ((_fontStyle & 1) ? FC_WEIGHT_BOLD : FC_WEIGHT_NORMAL),
                                        FC_SLANT, FcTypeInteger, ((_fontStyle & 2) ? FC_SLANT_ITALIC : FC_SLANT_ROMAN),
                                        FC_OUTLINE, FcTypeBool, ((_fontStyle & 8) ? true : false),
                                        FC_WIDTH, FcTypeInteger, ((_fontStyle & 32) ? FC_WIDTH_CONDENSED : ((_fontStyle & 64) ? FC_WIDTH_EXPANDED : FC_WIDTH_NORMAL)),
                                        FC_DPI, FcTypeDouble, (double)72.0, FC_SCALABLE, FcTypeBool, true, FC_ANTIALIAS, FcTypeBool, ((_antiAliasing != 0) ? true : false), NULL);
            }

            // Set default settings for missing pattern properties:
            FcDefaultSubstitute(target);
            if (!FcConfigSubstitute(NULL, target, FcMatchPattern)) {
                // Failed!
                if (_verbosity > 1)
                    fprintf(stderr, "libptbdrawtext_ftgles: FontConfig failed to substitute default properties for family %s, size %f pts and style flags %i.\n", _fontName, (float)_fontSize, _fontStyle);
                FcPatternDestroy(target);
                return (1);
            }

            // Have a matching pattern:
            if (_verbosity > 3) {
                fprintf(stderr, "libptbdrawtext_ftgles: Trying to find font that closely matches following specification:\n");
                FcPatternPrint(target);
            }
            // Perform font matching for the font in the default configuration (0) that best matches the
            // specified target pattern:
            FcPattern *matched = FcFontMatch(NULL, target, &result);
            if ((matched == NULL) || (result == FcResultNoMatch)) {
                // Failed!
                if (_verbosity > 1)
                    fprintf(stderr, "libptbdrawtext_ftgles: FontConfig failed to find a matching font for family %s, size %f pts and style flags %i.\n", _fontName, (float)_fontSize, _fontStyle);
                FcPatternDestroy(target);
                return (1);
            }

            // Success: Extract relevant information for Freetype-2, the font filename and faceIndex:
            if (_verbosity > 3) {
                fprintf(stderr, "libptbdrawtext_ftgles: Best matching font which will be selected for drawing has following specs:\n");
                FcPatternPrint(matched);
            }

            // Retrieve font filename for matched font:
            FcChar8 *localfontFileName = NULL;
            if (FcPatternGetString(matched, FC_FILE, 0, (FcChar8 **) & localfontFileName) != FcResultMatch) {
                // Failed!
                if (_verbosity > 1)
                    fprintf(stderr, "libptbdrawtext_ftgles: FontConfig did not find filename for font with family %s, size %f pts and style flags %i.\n", _fontName, (float)_fontSize, _fontStyle);
                FcPatternDestroy(target);
                FcPatternDestroy(matched);
                return (1);
            }

            strcpy(fontFileName, (char *)localfontFileName);

            // Retrieve faceIndex within fontfile:
            if (FcPatternGetInteger(matched, FC_INDEX, 0, &faceIndex) != FcResultMatch) {
                // Failed!
                if (_verbosity > 1)
                    fprintf(stderr, "libptbdrawtext_ftgles: FontConfig did not find faceIndex for font file %s, family %s, size %f pts and style flags %i.\n", fontFileName, _fontName, (float)_fontSize, _fontStyle);
                FcPatternDestroy(target);
                FcPatternDestroy(matched);
                return (1);
            }

            // Retrieve font family name for matched font:
            if (FcPatternGetString(matched, FC_FAMILY, 0, (FcChar8**) &localfontFileName) != FcResultMatch) {
                // Failed!
                if (_verbosity > 1) fprintf(stderr, "libptbdrawtext_ftgles: FontConfig did not return actual font family name for font with requested family %s, size %f pts and style flags %i.\n", _fontName, (float) _fontSize, _fontStyle);
                FcPatternDestroy(target);
                FcPatternDestroy(matched);
                return(1);
            }

            // Store it as actual name in fi:
            strcpy(fi->fontRealName, (char*) localfontFileName);

            // Release target pattern and matched pattern objects:
            FcPatternDestroy(target);
            FcPatternDestroy(matched);
        } else {
            // Use "raw" values as passed by calling client code:
            strcpy(fontFileName, _fontName);
            strcpy(fi->fontRealName, fontFileName);
            faceIndex = (int)_fontStyle;
        }

        /*
        // MK TODO FIXME: THIS DOES NOT WORK WITH FTGLES: It does not support a simple method to pass in the ft_face and build
        // a suitable OpenGL-ES representation out of it. We need to leave control to the constructors below, solely
        // based on filename of the font file, and hope for the best for now.
        //
        // Load & Create new font and face object, based on current spec settings:
        // We directly use the Freetype library, so we can spec the faceIndex for selection of textstyle, which wouldn't be
        // possible with the higher-level OGLFT constructor...
        FT_Error error = FT_New_Face(FTLibrary::Instance(), fontFileName, faceIndex, &ft_face );
        if (error) {
        if (_verbosity > 1) fprintf(stderr, "libptbdrawtext_ftgles: Freetype did not load face with index %i from font file %s.\n", faceIndex, fontFileName);
        return(1);
        }
        else {
        if (_verbosity > 3) fprintf(stderr, "libptbdrawtext_ftgles: Freetype loaded face %p with index %i from font file %s.\n", ft_face, faceIndex, fontFileName);
        }
        */

        // Create FTGLES face from Freetype face with given size and a 72 DPI resolution, aka _fontSize == pixelsize:
        if (_antiAliasing != 0) {
            // Font glyphs stored in textures - or better, rectangular sybregions of a texture atlas
            // where one or few textures store all glyphs of a font face. Textures are GL_ALPHA textures,
            // and blitted with alpha blending hard-coded to src_alpha, one_minus_src_alpha, so alpha blending
            // takes care of anti-aliased rendering.
            fi->faceT = new FTTextureFont(fontFileName);
            // Test the created face to make sure it will work correctly:
            if (fi->faceT->Error()) {
                delete(fi->faceT);
                fi->faceT = NULL;
                if (_verbosity > 1)
                    fprintf(stderr, "libptbdrawtext_ftgles: FTGLES did not recognize %s as a font file.\n", _fontName);
                return (1);
            }
            // Set size of font in pixels for a DPI of 72.
            fi->faceT->FaceSize(_fontSize, 72);

            // Use Unicode character encoding:
            fi->faceT->CharMap(ft_encoding_unicode);
        } else {
            // The current FTGLES implementation does not support BitMap or PixMap fonts without anti-aliasing.
            // The polygon/geometry based variants are not useful to us, or they are dysfunctional.
            // Therefore we are left with anti-aliased, alpha-blended texture fonts, iow., there ain't no
            // easy way to disable anti-aliasing. For now we "solve" this problem by ignoring it and just
            // do the same anti-aliased rendering, effectively ignoring the anti-aliasing enable/disable switch.
            fi->faceM = new FTTextureFont(fontFileName);
            // Test the created face to make sure it will work correctly:
            if (fi->faceM->Error()) {
                delete(fi->faceM);
                fi->faceM = NULL;
                if (_verbosity > 1)
                    fprintf(stderr, "libptbdrawtext_ftgles: FTGLES did not recognize %s as a font file.\n", _fontName);
                return (1);
            }
            // Set size of font in pixels for a DPI of 72.
            fi->faceM->FaceSize(_fontSize, 72);

            // Use Unicode character encoding:
            fi->faceM->CharMap(ft_encoding_unicode);
        }

        return (0);
    }


    void PsychSetTextFGColor(int context, double* color)
    {
        _fgcolor[0] = color[0];
        _fgcolor[1] = color[1];
        _fgcolor[2] = color[2];
        _fgcolor[3] = color[3];

        return;
    }

    void PsychSetTextBGColor(int context, double* color)
    {
        _bgcolor[0] = color[0];
        _bgcolor[1] = color[1];
        _bgcolor[2] = color[2];
        _bgcolor[3] = color[3];

        return;
    }

    // Assign dimensions and position of target viewport for rendering:
    // Assumption is that origin is bottom-left, x-axis pointing right, y-axis pointint up.
    // The plugin needs to remap or setup proper OpenGL transformations if its underlying
    // renderbackend needs a different geometric setup:
    void PsychSetTextViewPort(int context, double xs, double ys, double w, double h)
    {
        _vxs = xs;
        _vys = ys;
        _vw  = w;
        _vh  = h;

        return;
    }

    // Should the plugin use its own fontMapper implementation, or can it rely at least
    // partially - as defined by mapperFlags - on the PTB internal mapper?
    void PsychSetTextUseFontmapper(unsigned int useMapper, unsigned int mapperFlags)
    {
        _useOwnFontmapper = (useMapper > 0) ? true : false;
        _mapperFlags      = mapperFlags;

        return;
    }

    // Select font Name:
    int PsychSetTextFont(int context, const char* fontName)
    {
        strcpy(_fontName, fontName);
        return(0);
    }

    // Select font Style:
    int PsychSetTextStyle(int context, unsigned int fontStyle)
    {
        _fontStyle = fontStyle;
        return(0);
    }

    // Select font Size:
    int PsychSetTextSize(int context, double fontSize)
    {
        _fontSize = fontSize;
        return(0);
    }

    const char* PsychGetTextFont(int context)
    {
        // Check if rebuild of font face needed due to parameter
        // change. Reload/Rebuild font face if so, check for errors:
        fontCacheItem *fi = getForContext(context);
        if (!fi) return(NULL);

        return(fi->fontRealName);
    }

    int PsychDrawText(int context, double xStart, double yStart, int textLen, double* text)
    {
        // Check if rebuild of font face needed due to parameter
        // change. Reload/Rebuild font face if so, check for errors:
        fontCacheItem *fi = getForContext(context);
        if (!fi) return(1);

        // Synthesize Unicode wchar_t string from double vector:
        UnicodeToWchar(textLen, text);

        //  Done by FTGLES: glEnable( GL_TEXTURE_2D );
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
        int oldBlendEnable = glIsEnabled(GL_BLEND);

        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glOrthof(_vxs, _vxs + _vw, _vys, _vys + _vh, -1, +1);

        glMatrixMode(GL_MODELVIEW);

        // Rendering of background quad requested? -- True if background alpha > 0.
        if (_bgcolor[3] > 0) {
            // Yes. Compute bounding box of "to be drawn" text and render a quad in background color:
            float xmin, ymin, xmax, ymax;
            PsychMeasureText(context, textLen, text, &xmin, &ymin, &xmax, &ymax);

            // Compensate for 1-pixel off problems with some fonts:
            xmin--;
            xmax++;
            ymin--;
            ymax++;

            // Use immediate mode emulation functions from glesGlue.h/c.
            // They implemented the same trick i use in the Screen() core
            // code to get immediate mode convenience back - Great minds
            // think alike ;-) -- Their implementation is slightly more
            // complete and refined though. While this is not needed in
            // Screen() core code, we might use their glesGlue in other
            // parts of ptb, e.g., our moglcore for OpenGL-ES support.
            ftglBegin(GL_QUADS);
            ftglColor4f(_bgcolor[0], _bgcolor[1], _bgcolor[2], _bgcolor[3]);
            ftglVertex2f(xmin + xStart, ymin + yStart);
            ftglVertex2f(xmin + xStart, ymax + yStart);
            ftglVertex2f(xmax + xStart, ymax + yStart);
            ftglVertex2f(xmax + xStart, ymin + yStart);
            ftglEnd();
        }
        // Set text color:
        glColor4f(_fgcolor[0], _fgcolor[1], _fgcolor[2], _fgcolor[3]);

        // Set drawing start position:
        glPushMatrix();
        glTranslatef(xStart, yStart, 0.0f);

        // Draw the text:
        if (fi->faceT) {
            fi->faceT->Render(uniCodeText);
        } else {
            fi->faceM->Render(uniCodeText);
        }

        // Reset drawing position:
        glPopMatrix();

        glMatrixMode(GL_PROJECTION);
        glPopMatrix();

        // Back to standard modelview matrix:
        glMatrixMode(GL_MODELVIEW);

        // Restore blending enable state:
        if (oldBlendEnable) {
            glEnable(GL_BLEND);
        } else {
            glDisable(GL_BLEND);
        }

        // Ready!
        return (0);
    }

    int PsychMeasureText(int context, int textLen, double* text, float* xmin, float* ymin, float* xmax, float* ymax)
    {
        // Check if rebuild of font face needed due to parameter
        // change. Reload/Rebuild font face if so, check for errors:
        fontCacheItem *fi = getForContext(context);
        if (!fi) return(1);

        // Convert text to unicode wchar_t, filtering for unprintable characters:
        textLen = UnicodeToWchar(textLen, text);

        // Compute its bounding box:
        FTBBox box = (fi->faceT) ? fi->faceT->BBox(uniCodeText, textLen) : fi->faceM->BBox(uniCodeText, textLen);

        *xmin = box.Lower().X();
        *ymin = box.Lower().Y();
        *xmax = box.Upper().X();
        *ymax = box.Upper().Y();

        return (0);
    }

    int PsychInitText(void)
    {
        _firstCall = true;

        // Try to initialize libfontconfig - our fontMapper library for font matching and selection:
        if (!FcInit()) {
            if (_verbosity > 0)
                fprintf(stderr, "libptbdrawtext_ftgles: FontMapper initialization failed!\n");
            return (1);
        }

        // Clear cache of all fonts instances:
        memset(&cache, 0, sizeof(cache));
        for (int i = 0; i < MAX_CACHE_SLOTS; i++) cache[i].contextId = -1;

        if (_verbosity > 2) {
            fprintf(stderr, "libptbdrawtext_ftgles: External 'DrawText' text rendering plugin initialized.\n");
            fprintf(stderr, "libptbdrawtext_ftgles: Maximum number of cacheable fonts is %i, minimum number of supported concurrent windows is %i.\n", MAX_CACHE_SLOTS, MIN_GUARANTEED_CONTEXTS);
            fprintf(stderr, "libptbdrawtext_ftgles: This plugin uses multiple excellent free software libraries to do its work:\n");
            fprintf(stderr, "libptbdrawtext_ftgles: FTGLES (https://github.com/cdave1/ftgles) the FreeType OpenGL-ES library.\n");
            fprintf(stderr, "libptbdrawtext_ftgles: The FreeType-2 (http://freetype.sourceforge.net/) library.\n");
            fprintf(stderr, "libptbdrawtext_ftgles: The FontConfig (http://www.fontconfig.org) library.\n");
            fprintf(stderr, "libptbdrawtext_ftgles: Thanks!\n\n");
        }

        return (0);
    }

    int PsychShutdownText(int context)
    {
        // Delete specific context?
        if (context >= 0) {
            // Yes. Delete all objects for this context:
            for (int i = 0; i < MAX_CACHE_SLOTS; i++) {
                // Is this slot to be destructed?
                if (cache[i].contextId == context) {
                    fontCacheItem *fi = &(cache[i]);
                    fi->contextId = -1;

                    if (fi->faceT || fi->faceM) {
                        if (_verbosity > 3) fprintf(stderr, "libptbdrawtext_ftgles: In shutdown for context %i, slot %i:  faceT = %p faceM = %p\n", context, i, fi->faceT, fi->faceM);

                        // Delete OGLFT face objects:
                        if (fi->faceT) delete(fi->faceT);
                        fi->faceT = NULL;

                        if (fi->faceM) delete(fi->faceM);
                        fi->faceM = NULL;

                        // Delete Freetype face object:
                        if (fi->ft_face) FT_Done_Face(fi->ft_face);
                        fi->ft_face = NULL;
                    }
                }
            }

            return(0);
        }

        // Complete shutdown for the plugin:
        if (_verbosity > 3) fprintf(stderr, "libptbdrawtext_ftgles: Shutting down. Overall cache hit ratio was %f%%\n", (double) hitcount / (double) nowtime * 100);
        _firstCall = false;

        // Shutdown fontmapper library:
        // Actually, don't! Some versions of octave also use fontconfig internally, and there is only
        // one shared library instance in the process. Calling FcFini() here will shutdown that instance
        // and wreak havoc if octave or other clients later try to access that shutdown instance after
        // we've closed down. Keeping it alive is ok, let octave/matlab/the os whatever do the cleanup.
        // Luckily we can always call FcInit() in our init path, as it turns into a no-op if fontconfig
        // has been already brought online by some other client, e.g., octave or matlab.
        // Should fix bug from forum msg #12560
        // FcFini();

        return (0);
    }
} // extern "C"
