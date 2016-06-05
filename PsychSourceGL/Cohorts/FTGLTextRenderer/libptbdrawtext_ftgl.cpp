/*
 * libptbdrawtext_ftgl.cpp: Source of libptbdrawtext_ftgl dynamic shared library plugin
 * for Screen('DrawText') text rendering via external engine.
 *
 * Features:
 *
 * - Texture mapped renderer.
 * - Fast due to texture object and display list caching for fast glyph recycling.
 * - Good text layouting.
 * - Supports all Freetype-2 supported fonts, e.g., vectorgraphics TrueType fonts.
 * - Anti-Aliased drawing via Alpha-Blending.
 * - Full Unicode support.
 * - Text measuring support, e.g., bounding boxes.
 *
 * New style:   Includes our own - slightly modified - version of OGLFT.h/cpp for libglft support, thereby avoiding runtime
 *              dependencies on a properly configured & built libglft (Many binary installs of this library are built without unicode support).
 *
 *              Includes qstringqcharemulation.cpp - our own super-minimalistic reimplementation of QT-4 QChar and QString and QRgb support. This
 *              reimplementation only contains the functionality and C++ functions that are needed for OGLFT's internal Unicode character and string
 *              processing and for communication of Unicode strings from libptbdrawtext_ftgl.cpp to OGLFT. The few exported interfaces are API compatible
 *              with QT-4, so one can choose to compile either against this wrapper or the real QT-4 toolkit, but the implementation is not based on or
 *              ABI compatible with QT-4, it just emulates its behaviour. Also almost all functionality is missing, except the few bits that we need.
 *
 *              Rationale:    Requiring users to install a multiple-hundred megabyte size toolkit and paying the cost of linking against a 15 MB library,
 *                            plus possible version conflicts if the runtime (Octave or Matlab) should use QT itself is just too much trouble and overhead
 *                            for simple Unicode string handling.
 *
 * The OS/X 64-Bit build assumes the default install locations for libfontconfig and libfreetype in /usr/X11/include and /usr/X11/lib:
 *
 * g++ -g -DHAVE_OPENGL_DIR -I.  -I/usr/X11/include/ -I/usr/X11/include/freetype2/ -L/usr/X11/lib/ -framework OpenGL -l fontconfig -l freetype -dynamiclib -o libptbdrawtext_ftgl64.dylib libptbdrawtext_ftgl.cpp qstringqcharemulation.cpp OGLFT.cpp
 *
 * Building for Linux:
 *
 * g++ -g -fPIC -I. -I/usr/include/ -I/usr/include/freetype2/ -L/usr/lib -pie -shared -Wl,-Bsymbolic -Wl,-Bsymbolic-functions -Wl,--version-script=linuxexportlist.txt -o libptbdrawtext_ftgl.so.1 libptbdrawtext_ftgl.cpp qstringqcharemulation.cpp OGLFT.cpp -lGL -lGLU -lfontconfig -lfreetype
 *
 * libptbdrawtext_ftgl is copyright (c) 2010-2016 by Mario Kleiner.
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

// Standard libs for debug output:
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Include all GLFT and QT stuff:
#include "OGLFT.h"

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
FT_Matrix _matrix;
FT_Vector _vector;
double _xp;
double _yp;

typedef struct fontCacheItem_t {
    int contextId;
    unsigned int timestamp;
    int antiAliasing;
    char fontName[FILENAME_MAX];
    char fontRealName[FILENAME_MAX];
    unsigned int fontStyle;
    double fontSize;
    FT_Matrix matrix;
    FT_Vector vector;
    OGLFT::TranslucentTexture    *faceT;
    OGLFT::MonochromeTexture    *faceM;
    FT_Face ft_face;
} fontCacheItem;
fontCacheItem cache[MAX_CACHE_SLOTS];

#ifdef _MSC_VER
#ifdef OGLFT_BUILD
#define OGLFT_API __declspec(dllexport)
#else
#define OGLFT_API __declspec(dllimport)
#endif
#else
#define OGLFT_API
#endif

extern "C" {

OGLFT_API int PsychInitText(void);
OGLFT_API int PsychShutdownText(int context);
OGLFT_API int PsychRebuildFont(fontCacheItem* fi);
OGLFT_API int PsychSetTextFont(int context, const char* fontName);
OGLFT_API const char* PsychGetTextFont(int context);
OGLFT_API int PsychSetTextStyle(int context, unsigned int fontStyle);
OGLFT_API int PsychSetTextSize(int context, double fontSize);
OGLFT_API void PsychSetTextFGColor(int context, double* color);
OGLFT_API void PsychSetTextBGColor(int context, double* color);
OGLFT_API void PsychSetTextUseFontmapper(unsigned int useMapper, unsigned int mapperFlags);
OGLFT_API void PsychSetTextViewPort(int context, double xs, double ys, double w, double h);
OGLFT_API int PsychDrawText(int context, double xStart, double yStart, int textLen, double* text);
OGLFT_API int PsychMeasureText(int context, int textLen, double* text, float* xmin, float* ymin, float* xmax, float* ymax, float* xadvance);
OGLFT_API void PsychSetTextVerbosity(unsigned int verbosity);
OGLFT_API void PsychSetTextAntiAliasing(int context, int antiAliasing);
OGLFT_API void PsychSetAffineTransformMatrix(int context, double matrix[2][3]);
OGLFT_API void PsychGetTextCursor(int context, double* xp, double* yp, double* height);

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
                ((strcmp(fi->fontName, _fontName) == 0) || (strcmp(fi->fontRealName, _fontName) == 0)) &&
                (fi->matrix.xx == _matrix.xx && fi->matrix.xy == _matrix.xy && fi->matrix.yx == _matrix.yx && fi->matrix.yy == _matrix.yy &&
                 fi->vector.x == _vector.x && fi->vector.y == _vector.y)) {
                // Match! We have cached OGLFT font objects for this font on this context. Return them:
                hitcount++;

                if (_verbosity > 15) fprintf(stdout, "libptbdrawtext_ftgl: Cache hit for contextId %i at slot %i. Hit ratio is %f%%\n", contextId, i, (double) hitcount / (double) nowtime * 100);

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
        if (_verbosity > 12) fprintf(stdout, "libptbdrawtext_ftgl: Nothing cached for contextId %i. Using new slot %i. %i free slots remaining.\n", contextId, freeslot, freecount);
    }
    else if (lruslotid >= 0) {
        // No. Overwrite least recently used font object for current contextId:
        fi = &(cache[lruslotid]);
        if (_verbosity > 12) fprintf(stdout, "libptbdrawtext_ftgl: Nothing cached for contextId %i but cache full. LRU replacing slot %i, age %i. %i free slots remaining.\n", contextId, lruslotid, lruage, freecount);
    }
    else {
        // Cache full, with no possibility to even LRU replace on this new context (aka window). Game over!
        if (_verbosity > 0) fprintf(stdout, "libptbdrawtext_ftgl: ERROR when trying to setup new drawtext context %i: Font cache full with no way to free up resources! Text drawing will fail!\n", contextId);
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
    fi->matrix = _matrix;
    fi->vector = _vector;

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

        if (_verbosity > 12) fprintf(stdout, "libptbdrawtext_ftgl: Destroying old font face...\n");

        // Delete underlying FreeType representation:
        FT_Done_Face(fi->ft_face);
        fi->ft_face = NULL;
    }

    if (_useOwnFontmapper) {
        FcResult result = FcResultMatch; // Must init this due to weirdness in libfontconfig...
        FcPattern* target = NULL;

        if (_fontName[0] == '-') {
            // _fontName starts with a '-' dash: This is not a simple font family string but a special
            // fontspec string in FontConfig's special format. It contains many possible required font
            // properties encoded in the string. Parse it into a font matching pattern:
            target = FcNameParse((FcChar8*) &(_fontName[1]));

            // Need to manually add the current _fontSize, otherwise inconsistent stuff may happen:
            FcPatternAddDouble(target, FC_PIXEL_SIZE, _fontSize);
        }
        else {
            // _fontName contains only font family name: Build matching pattern based on _fontSize and
            // the flags provided in _fontStyle, according to the conventions in Psychtoolbox Screen('TextStyle'):
            target = FcPatternBuild (0, FC_FAMILY, FcTypeString, _fontName, FC_PIXEL_SIZE, FcTypeDouble, _fontSize,
                                     FC_WEIGHT, FcTypeInteger, ((_fontStyle & 1) ? FC_WEIGHT_BOLD : FC_WEIGHT_NORMAL),
                                     FC_SLANT, FcTypeInteger, ((_fontStyle & 2) ?  FC_SLANT_ITALIC : FC_SLANT_ROMAN),
                                     FC_OUTLINE, FcTypeBool, ((_fontStyle & 8) ? true : false),
                                     FC_WIDTH, FcTypeInteger, ( (_fontStyle & 32) ?  FC_WIDTH_CONDENSED : ((_fontStyle & 64) ?  FC_WIDTH_EXPANDED : FC_WIDTH_NORMAL) ),
                                     FC_DPI, FcTypeDouble, (double) 72.0,
                                     FC_SCALABLE, FcTypeBool, true,
                                     FC_ANTIALIAS, FcTypeBool, ((_antiAliasing != 0) ? true : false),
                                     NULL);
        }

        // Set default settings for missing pattern properties:
        FcDefaultSubstitute(target);
        if (!FcConfigSubstitute(NULL, target, FcMatchPattern)) {
            // Failed!
            if (_verbosity > 1) fprintf(stdout, "libptbdrawtext_ftgl: FontConfig failed to substitute default properties for family %s, size %f pts and style flags %i.\n", _fontName, (float) _fontSize, _fontStyle);
            FcPatternDestroy(target);
            return(1);
        }

        // Have a matching pattern:
        if (_verbosity > 5) {
            fprintf(stdout, "libptbdrawtext_ftgl: Trying to find font that closely matches following specification:\n");
            FcPatternPrint(target);
        }

        // Perform font matching for the font in the default configuration (0) that best matches the
        // specified target pattern:
        FcPattern* matched = FcFontMatch(NULL, target, &result);
        if ((matched == NULL) || (result == FcResultNoMatch)) {
            // Failed!
            if (_verbosity > 1) fprintf(stdout, "libptbdrawtext_ftgl: FontConfig failed to find a matching font for family %s, size %f pts and style flags %i.\n", _fontName, (float) _fontSize, _fontStyle);
            FcPatternDestroy(target);
            return(1);
        }

        // Success: Extract relevant information for Freetype-2, the font filename and faceIndex:
        if (_verbosity > 5) {
            fprintf(stdout, "libptbdrawtext_ftgl: Best matching font which will be selected for drawing has following specs:\n");
            FcPatternPrint(matched);
        }

        // Retrieve font filename for matched font:
        FcChar8* localfontFileName = NULL;
        if (FcPatternGetString(matched, FC_FILE, 0, (FcChar8**) &localfontFileName) != FcResultMatch) {
            // Failed!
            if (_verbosity > 1) fprintf(stdout, "libptbdrawtext_ftgl: FontConfig did not find filename for font with family %s, size %f pts and style flags %i.\n", _fontName, (float) _fontSize, _fontStyle);
            FcPatternDestroy(target);
            FcPatternDestroy(matched);
            return(1);
        }

        strcpy(fontFileName, (char*) localfontFileName);

        // Retrieve faceIndex within fontfile:
        if (FcPatternGetInteger(matched, FC_INDEX, 0, &faceIndex) != FcResultMatch)  {
            // Failed!
            if (_verbosity > 1) fprintf(stdout, "libptbdrawtext_ftgl: FontConfig did not find faceIndex for font file %s, family %s, size %f pts and style flags %i.\n", fontFileName, _fontName, (float) _fontSize, _fontStyle);
            FcPatternDestroy(target);
            FcPatternDestroy(matched);
            return(1);
        }

        // Retrieve font family name for matched font:
        if (FcPatternGetString(matched, FC_FAMILY, 0, (FcChar8**) &localfontFileName) != FcResultMatch) {
            // Failed!
            if (_verbosity > 1) fprintf(stdout, "libptbdrawtext_ftgl: FontConfig did not return actual font family name for font with requested family %s, size %f pts and style flags %i.\n", _fontName, (float) _fontSize, _fontStyle);
            FcPatternDestroy(target);
            FcPatternDestroy(matched);
            return(1);
        }

        // Store it as actual name in fi:
        strcpy(fi->fontRealName, (char*) localfontFileName);

        // Release target pattern and matched pattern objects:
        FcPatternDestroy(target);
        FcPatternDestroy(matched);
    }
    else {
        // Use "raw" values as passed by calling client code:
        strcpy(fontFileName, _fontName);
        strcpy(fi->fontRealName, fontFileName);
        faceIndex = (int) _fontStyle;
    }

    // Load & Create new font and face object, based on current spec settings:
    // We directly use the Freetype library, so we can spec the faceIndex for selection of textstyle, which wouldn't be
    // possible with the higher-level OGLFT constructor...
    FT_Error error = FT_New_Face( OGLFT::Library::instance(), fontFileName, faceIndex, &fi->ft_face );
    if (error) {
        if (_verbosity > 1) fprintf(stdout, "libptbdrawtext_ftgl: Freetype did not load face with index %i from font file %s.\n", faceIndex, fontFileName);
        return(1);
    }
    else {
        if (_verbosity > 5) fprintf(stdout, "libptbdrawtext_ftgl: Freetype loaded face %p with index %i from font file %s.\n", fi->ft_face, faceIndex, fontFileName);
    }

    // Apply affine transformations, if any:
    FT_Set_Transform(fi->ft_face, &_matrix, &_vector);

    // Create FTGL face from Freetype face with given size and a 72 DPI resolution, aka _fontSize == pixelsize:
    if (_antiAliasing != 0) {
        fi->faceT = new OGLFT::TranslucentTexture(fi->ft_face, _fontSize, 72);
        // Test the created face to make sure it will work correctly:
        if (!fi->faceT->isValid()) {
            if (_verbosity > 1) fprintf(stdout, "libptbdrawtext_ftgl: Freetype did not recognize %s as a font file.\n", _fontName);
            delete(fi->faceT);
            fi->faceT = NULL;
            return(1);
        }

        // Do not reset GL_MODELVIEW matrix after text drawing, so
        // we can read out the post-draw text cursor position from it:
        fi->faceT->setAdvance(true);
    }
    else {
        fi->faceM = new OGLFT::MonochromeTexture(fi->ft_face, _fontSize, 72);
        // Test the created face to make sure it will work correctly:
        if (!fi->faceM->isValid()) {
            if (_verbosity > 1) fprintf(stdout, "libptbdrawtext_ftgl: Freetype did not recognize %s as a font file.\n", _fontName);
            delete(fi->faceM);
            fi->faceM = NULL;
            return(1);
        }

        // Do not reset GL_MODELVIEW matrix after text drawing, so
        // we can read out the post-draw text cursor position from it:
        fi->faceM->setAdvance(true);
    }

    // Ready!
    return(0);
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

// Assign affine 2D text transformation matrix:
void PsychSetAffineTransformMatrix(int context, double matrix[2][3])
{
    // Matrix in 16.16 format:
    _matrix.xx = (FT_Fixed) (matrix[0][0] * 65536);
    _matrix.xy = (FT_Fixed) (matrix[0][1] * 65536);
    _matrix.yx = (FT_Fixed) (matrix[1][0] * 65536);
    _matrix.yy = (FT_Fixed) (matrix[1][1] * 65536);
    // Translation vector in 26.6 format:
    _vector.x  = (FT_Pos)   (matrix[0][2] * 64);
    _vector.y  = (FT_Pos)   (matrix[1][2] * 64);
}

void PsychGetTextCursor(int context, double* xp, double* yp, double* height)
{
    fontCacheItem *fi = getForContext(context);
    *xp = _xp;
    *yp = _yp;

    // Try to get line height of font:
    if (fi) {
        if (fi->faceT)
            *height = fi->faceT->height();
        else
            *height = fi->faceM->height();
    }
    else {
        *height = 0;
    }
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
    int i;
    GLuint ti;
    QChar* myUniChars;
    GLdouble modelview[4][4];

    // On first invocation after init we need to generate a useless texture object.
    // This is a weird workaround for some weird bug somewhere in FTGL...
    if (_firstCall) {
        _firstCall = false;
        glGenTextures(1, &ti);
    }

    // Check if rebuild of font face needed due to parameter
    // change. Reload/Rebuild font face if so, check for errors:
    fontCacheItem *fi = getForContext(context);
    if (!fi) return(1);

    // Synthesize Unicode QString from double vector:
    myUniChars = new QChar[textLen];
    for(i = 0; i < textLen; i++) {
        myUniChars[i] = QChar((unsigned int) text[i]);
    }

    QString uniCodeText = QString(myUniChars, textLen);
    delete [] myUniChars;

    glPushClientAttrib(GL_CLIENT_ALL_ATTRIB_BITS);
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glPixelStorei( GL_UNPACK_ALIGNMENT, 1);
    glEnable( GL_TEXTURE_2D );
    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(_vxs, _vxs + _vw, _vys, _vys + _vh);
    glMatrixMode(GL_MODELVIEW);

    // Set text color: This will be filtered by OGLFT for redundant settings:
    if (fi->faceT) {
        fi->faceT->setForegroundColor( _fgcolor[0], _fgcolor[1], _fgcolor[2], _fgcolor[3]);
    }
    else {
        fi->faceM->setForegroundColor( _fgcolor[0], _fgcolor[1], _fgcolor[2], _fgcolor[3]);
    }

    // Rendering of background quad requested? -- True if background alpha > 0.
    if (_bgcolor[3] > 0) {
        // Yes. Compute bounding box of "to be drawn" text and render a quad in background color:
        float xmin, ymin, xmax, ymax, xadvance;
        PsychMeasureText(context, textLen, text, &xmin, &ymin, &xmax, &ymax, &xadvance);
        glColor4fv(&(_bgcolor[0]));
        glRectf(xmin + xStart, ymin + yStart, xmax + xStart, ymax + yStart);
    }

    // Enable alpha-test against an alpha-value greater zero during draw.
    // This way, non-text pixels (with alpha equal to zero) are discarded.
    glEnable(GL_ALPHA_TEST);
    glAlphaFunc(GL_GREATER, 0);

    // Draw the text at selected start location:
    glPushMatrix();
    if (fi->faceT) {
        fi->faceT->draw(xStart, yStart, uniCodeText);
    }
    else {
        fi->faceM->draw(xStart, yStart, uniCodeText);
    }

    // Extract final text cursor position from GL_MODELVIEW matrix:
    glGetDoublev(GL_MODELVIEW_MATRIX, &(modelview[0][0]));
    _xp = modelview[3][0];
    _yp = modelview[3][1];
    glPopMatrix();

    // Disable alpha test after blit:
    glDisable(GL_ALPHA_TEST);

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glDisable( GL_TEXTURE_2D );
    glPopAttrib();
    glPopClientAttrib();

    // Ready!
    return(0);
}

int PsychMeasureText(int context, int textLen, double* text, float* xmin, float* ymin, float* xmax, float* ymax, float* xadvance)
{
    int i;
    QChar* myUniChars;

    // Check if rebuild of font face needed due to parameter
    // change. Reload/Rebuild font face if so, check for errors:
    fontCacheItem *fi = getForContext(context);
    if (!fi) return(1);

    // Synthesize Unicode QString from double vector:
    myUniChars = new QChar[textLen];
    for(i = 0; i < textLen; i++) {
        myUniChars[i] = QChar((unsigned int) text[i]);
    }

    QString uniCodeText = QString(myUniChars, textLen);
    delete [] myUniChars;

    // Compute its bounding box:
    glPushMatrix();
    OGLFT::BBox box = (fi->faceT) ? fi->faceT->measure(uniCodeText) : fi->faceM->measure(uniCodeText);
    glPopMatrix();

    *xmin = box.x_min_;
    *ymin = box.y_min_;
    *xmax = box.x_max_;
    *ymax = box.y_max_;

    // xadvance = How far to horizontally advance the x text
    // drawing cursor position after drawing the text string:
    // newXPos = oldXPos + xadvance
    *xadvance = box.advance_.dx_;

    return(0);
}

int PsychInitText(void)
{
    _firstCall = true;

    // Try to initialize libfontconfig - our fontMapper library for font matching and selection:
    if (!FcInit()) {
        if (_verbosity > 0) fprintf(stdout, "libptbdrawtext_ftgl: FontMapper initialization failed!\n");
        return(1);
    }

    // Clear cache of all fonts instances:
    memset(&cache, 0, sizeof(cache));
    for (int i = 0; i < MAX_CACHE_SLOTS; i++) cache[i].contextId = -1;

    if (_verbosity > 2)    {
        fprintf(stdout, "libptbdrawtext_ftgl: External 'DrawText' text rendering plugin initialized.\n");
        fprintf(stdout, "libptbdrawtext_ftgl: Maximum number of cacheable fonts is %i, minimum number of supported concurrent windows is %i.\n", MAX_CACHE_SLOTS, MIN_GUARANTEED_CONTEXTS);
        fprintf(stdout, "libptbdrawtext_ftgl: This plugin uses multiple excellent free software libraries to do its work:\n");
        fprintf(stdout, "libptbdrawtext_ftgl: OGLFT (http://oglft.sourceforge.net/) the OpenGL-FreeType library.\n");
        fprintf(stdout, "libptbdrawtext_ftgl: The FreeType-2 (http://freetype.sourceforge.net/) library.\n");
        fprintf(stdout, "libptbdrawtext_ftgl: The FontConfig (http://www.fontconfig.org) library.\n");
        fprintf(stdout, "libptbdrawtext_ftgl: Thanks!\n\n");
    }

    return(0);
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
                    if (_verbosity > 5) fprintf(stdout, "libptbdrawtext_ftgl: In shutdown for context %i, slot %i:  faceT = %p faceM = %p\n", context, i, fi->faceT, fi->faceM);

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
    if (_verbosity > 5) fprintf(stdout, "libptbdrawtext_ftgl: Shutting down. Overall cache hit ratio was %f%%\n", (double) hitcount / (double) nowtime * 100);
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

    return(0);
}

} // extern "C"
