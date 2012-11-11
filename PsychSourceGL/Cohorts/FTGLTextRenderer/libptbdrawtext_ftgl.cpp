/*
 * libptbdrawtext_ftgl.cpp: Source of libptbdrawtext_ftgl dynamic shared library plugin
 * for Screen('DrawText') text rendering via external engine.
 *
 * This works on OS/X but is mostly intended for GNU/Linux.
 *
 * Features:
 *
 * - Texture mapped renderer, like on OS/X with ATSU Drawtext.
 * - Fast due to texture object and display list caching for fast glyph recycling.
 * - Good text layouting.
 * - Supports all Freetype-2 supported fonts, e.g., vectorgraphics TrueType fonts.
 * - Anti-Aliased drawing via Alpha-Blending.
 * - Full Unicode support.
 * - Text measuring support, e.g., bounding boxes.
 *
 * Building for OS/X:
 *
 * Old style: Needs runtimes and development headers of libfreetype2, libfontconfig, libglft build with Unicode+QT4 support, and the QT4 SDK installed:
 *
 * g++ -g -framework QtCore -DQT_CORE_LIB -DQT_GUI_LIB -DOGLFT_QT_VERSION=4 -I /usr/local/include/OGLFT/ -I/opt/local/include/ -I /usr/X11R6/include/freetype2/ -I /Library/Frameworks/QtGui.framework/Headers/ -I /Library/Frameworks/QtCore.framework/Headers/ -L/opt/local/lib/ -framework OpenGL -l fontconfig -l oglft -l freetype -dynamiclib -o libptbdrawtext_ftgl.dylib libptbdrawtext_ftgl.cpp 
 *
 * New style:	Includes our own - slightly modified - version of OGLFT.h/cpp for libglft support, thereby avoiding runtime
 *				dependencies on a properly configured & built libglft (Many binary installs of this library are built without unicode support).
 *
 *				Includes qstringqcharemulation.cpp - our own super-minimalistic reimplementation of QT-4 QChar and QString and QRgb support. This
 *				reimplementation only contains the functionality and C++ functions that are needed for OGLFT's internal Unicode character and string
 *				processing and for communication of Unicode strings from libptbdrawtext_ftgl.cpp to OGLFT. The few exported interfaces are API compatible
 *				with QT-4, so one can choose to compile either against this wrapper or the real QT-4 toolkit, but the implementation is not based on or
 *				ABI compatible with QT-4, it just emulates its behaviour. Also almost all functionality is missing, except the few bits that we need.
 *
 *				Rationale:	Requiring users to install a multiple-hundred megabyte size toolkit and paying the cost of linking against a 15 MB library,
 *							plus possible version conflicts if the runtime (Octave or Matlab) should use QT itself is just too much trouble and overhead
 *							for simple Unicode string handling.
 *
 *
 * The OS/X build assumes you have libfreetype2 and libfontconfig from the OS/X DarwinPorts project installed. See: http://darwinports.com/
 * 
 * g++ -g -DHAVE_OPENGL_DIR -DGLU_TESS_CALLBACK_TRIPLEDOT -I.  -I/opt/local/include/ -I/opt/local/include/freetype2/ -L/opt/local/lib/ -framework OpenGL -l fontconfig -l freetype -dynamiclib -o libptbdrawtext_ftgl.dylib libptbdrawtext_ftgl.cpp qstringqcharemulation.cpp OGLFT.cpp
 *
 * The alternative OS/X 64-Bit build on 10.7 Lion assumes the default install locations for libfontconfig and libfreetype in /usr/X11R6/include and /usr/X11R6/lib:
 *
 * g++ -g -DHAVE_OPENGL_DIR -I.  -I/usr/X11R6/include/ -I/usr/X11R6/include/freetype2/ -L/usr/X11R6/lib/ -framework OpenGL -l fontconfig -l freetype -dynamiclib -o libptbdrawtext_ftgl64.dylib libptbdrawtext_ftgl.cpp qstringqcharemulation.cpp OGLFT.cpp
 *
 * Building for Linux:
 * 
 * g++ -g -fPIC -I. -I/usr/include/ -I/usr/include/freetype2/ -L/usr/lib -l GL -l GLU -l fontconfig -l freetype -pie -shared -Wl,-Bsymbolic -Wl,-Bsymbolic-functions -Wl,--version-script=linuxexportlist.txt -o libptbdrawtext_ftgl.so.1 libptbdrawtext_ftgl.cpp qstringqcharemulation.cpp OGLFT.cpp
 *
 * libptbdrawtext_ftgl is copyright (c) 2010 by Mario Kleiner.
 * It is licensed to you under the LGPL license as follows:
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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

int _antiAliasing = 1;
unsigned int _verbosity = 2;
bool _firstCall = true;
bool _needsRebuild = true;
bool _useOwnFontmapper = false;
unsigned _mapperFlags;
double _vxs, _vys, _vw, _vh;
char _fontName[4096] = { 0 };
unsigned int _fontStyle = 0;
char _fontFileName[4096] = { 0 };
int _faceIndex = 0;
double _fontSize = 0.0;
GLfloat _fgcolor[4];
GLfloat _bgcolor[4];

OGLFT::TranslucentTexture	*faceT = NULL;
OGLFT::MonochromeTexture	*faceM = NULL;
static FT_Face ft_face = NULL;

extern "C" {

int PsychInitText(void);
int PsychShutdownText(void);
int PsychRebuiltFont(void);
int PsychSetTextFont(const char* fontName);
int PsychSetTextStyle(unsigned int fontStyle);
int PsychSetTextSize(double fontSize);
void PsychSetTextFGColor(double* color);
void PsychSetTextBGColor(double* color);
void PsychSetTextUseFontmapper(unsigned int useMapper, unsigned int mapperFlags);
void PsychSetTextViewPort(double xs, double ys, double w, double h);
int PsychDrawText(double xStart, double yStart, int textLen, double* text);
int PsychMeasureText(int textLen, double* text, float* xmin, float* ymin, float* xmax, float* ymax);
void PsychSetTextVerbosity(unsigned int verbosity);
void PsychSetTextAntiAliasing(int antiAliasing);

void PsychSetTextVerbosity(unsigned int verbosity)
{
	_verbosity = verbosity;
	return;
}

void PsychSetTextAntiAliasing(int antiAliasing)
{
	if (_antiAliasing != antiAliasing) _needsRebuild = true;
	
	_antiAliasing = antiAliasing;
	return;
}

int PsychRebuildFont(void)
{
	// Destroy old font object, if any:
	if (faceT || faceM) {
		// Delete OGLFT face object:
		if (faceT) delete(faceT);
		faceT = NULL;
		
		if (faceM) delete(faceM);
		faceM = NULL;

		if (_verbosity > 3) fprintf(stderr, "libptbdrawtext_ftgl: Destroying old font face...\n");
		
		// Delete underlying FreeType representation:
		FT_Done_Face(ft_face);
		ft_face = NULL;
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
			if (_verbosity > 1) fprintf(stderr, "libptbdrawtext_ftgl: FontConfig failed to substitute default properties for family %s, size %f pts and style flags %i.\n", _fontName, (float) _fontSize, _fontStyle);
			FcPatternDestroy(target);
			return(1);
		}
		
		// Have a matching pattern:
		if (_verbosity > 3) {
			fprintf(stderr, "libptbdrawtext_ftgl: Trying to find font that closely matches following specification:\n");
			FcPatternPrint(target);
		}

		// Perform font matching for the font in the default configuration (0) that best matches the
		// specified target pattern:
		FcPattern* matched = FcFontMatch(NULL, target, &result);
		if ((matched == NULL) || (result == FcResultNoMatch)) {
			// Failed!
			if (_verbosity > 1) fprintf(stderr, "libptbdrawtext_ftgl: FontConfig failed to find a matching font for family %s, size %f pts and style flags %i.\n", _fontName, (float) _fontSize, _fontStyle);
			FcPatternDestroy(target);
			return(1);
		}

		// Success: Extract relevant information for Freetype-2, the font filename and faceIndex:
		if (_verbosity > 3) {
			fprintf(stderr, "libptbdrawtext_ftgl: Best matching font which will be selected for drawing has following specs:\n");
			FcPatternPrint(matched);
		}

		// Retrieve font filename for matched font:
		FcChar8* localfontFileName = NULL;
		if (FcPatternGetString(matched, FC_FILE, 0, (FcChar8**) &localfontFileName) != FcResultMatch) {
			// Failed!
			if (_verbosity > 1) fprintf(stderr, "libptbdrawtext_ftgl: FontConfig did not find filename for font with family %s, size %f pts and style flags %i.\n", _fontName, (float) _fontSize, _fontStyle);
			FcPatternDestroy(target);
			FcPatternDestroy(matched);
			return(1);
		}

		strcpy(_fontFileName, (char*) localfontFileName);

		// Retrieve faceIndex within fontfile:
		if (FcPatternGetInteger(matched, FC_INDEX, 0, &_faceIndex) != FcResultMatch)  {
			// Failed!
			if (_verbosity > 1) fprintf(stderr, "libptbdrawtext_ftgl: FontConfig did not find faceIndex for font file %s, family %s, size %f pts and style flags %i.\n", _fontFileName, _fontName, (float) _fontSize, _fontStyle);
			FcPatternDestroy(target);
			FcPatternDestroy(matched);
			return(1);
		}

		// Release target pattern and matched pattern objects:
		FcPatternDestroy(target);
		FcPatternDestroy(matched);
	}
	else {
		// Use "raw" values as passed by calling client code:
		strcpy(_fontFileName, _fontName);
		_faceIndex = (int) _fontStyle;
	}
	
	// Load & Create new font and face object, based on current spec settings:
	// We directly use the Freetype library, so we can spec the faceIndex for selection of textstyle, which wouldn't be
	// possible with the higher-level OGLFT constructor...
	FT_Error error = FT_New_Face( OGLFT::Library::instance(), _fontFileName, _faceIndex, &ft_face );
	if (error) {
		if (_verbosity > 1) fprintf(stderr, "libptbdrawtext_ftgl: Freetype did not load face with index %i from font file %s.\n", _faceIndex, _fontFileName);
		return(1);
	}
	else {
		if (_verbosity > 3) fprintf(stderr, "libptbdrawtext_ftgl: Freetype loaded face %p with index %i from font file %s.\n", ft_face, _faceIndex, _fontFileName);
	}

	// Create FTGL face from Freetype face with given size and a 72 DPI resolution, aka _fontSize == pixelsize:
	if (_antiAliasing != 0) {
		faceT = new OGLFT::TranslucentTexture(ft_face, _fontSize, 72);
		// Test the created face to make sure it will work correctly:
		if (!faceT->isValid()) {
			if (_verbosity > 1) fprintf(stderr, "libptbdrawtext_ftgl: Freetype did not recognize %s as a font file.\n", _fontName);
			return(1);
		}
	}
	else {
		faceM = new OGLFT::MonochromeTexture(ft_face, _fontSize, 72);
		// Test the created face to make sure it will work correctly:
		if (!faceM->isValid()) {
			if (_verbosity > 1) fprintf(stderr, "libptbdrawtext_ftgl: Freetype did not recognize %s as a font file.\n", _fontName);
			return(1);
		}
	}

	// Ready!
	_needsRebuild = false;
	
	return(0);
}

void PsychSetTextFGColor(double* color)
{
	_fgcolor[0] = color[0];
	_fgcolor[1] = color[1];
	_fgcolor[2] = color[2];
	_fgcolor[3] = color[3];

	return;
}

void PsychSetTextBGColor(double* color)
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
void PsychSetTextViewPort(double xs, double ys, double w, double h)
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
int PsychSetTextFont(const char* fontName)
{
	if (strcmp(_fontName, fontName)) _needsRebuild = true;
	strcpy(_fontName, fontName);
	return(0);
}

// Select font Style:
int PsychSetTextStyle(unsigned int fontStyle)
{
	if (_fontStyle != fontStyle) _needsRebuild = true;
	_fontStyle = fontStyle;
	return(0);
}

// Select font Size:
int PsychSetTextSize(double fontSize)
{
	if (_fontSize != fontSize) _needsRebuild = true;
	_fontSize = fontSize;
	return(0);
}

int PsychDrawText(double xStart, double yStart, int textLen, double* text)
{
	int i;
	GLuint ti;
	QChar* myUniChars = new QChar[textLen];
	
	// On first invocation after init we need to generate a useless texture object.
	// This is a weird workaround for some weird bug somewhere in FTGL...
	if (_firstCall) {
		_firstCall = false;
		glGenTextures(1, &ti);
	}
	
	// Check if rebuild of font face needed due to parameter
	// change. Reload/Rebuild font face if so, check for errors:
	if (_needsRebuild && PsychRebuildFont()) return(1);

	// Synthesize Unicode QString from double vector:
	for(i = 0; i < textLen; i++) {
		myUniChars[i] = QChar((unsigned int) text[i]); 
	}	
	QString	uniCodeText = QString(myUniChars, textLen);  
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
	if (faceT) {
		faceT->setForegroundColor( _fgcolor[0], _fgcolor[1], _fgcolor[2], _fgcolor[3]);
	}
	else {
		faceM->setForegroundColor( _fgcolor[0], _fgcolor[1], _fgcolor[2], _fgcolor[3]);
	}

	// Rendering of background quad requested? -- True if background alpha > 0.
	if (_bgcolor[3] > 0) {
		// Yes. Compute bounding box of "to be drawn" text and render a quad in background color:
		float xmin, ymin, xmax, ymax;
		PsychMeasureText(textLen, text, &xmin, &ymin, &xmax, &ymax);
		glColor4fv(&(_bgcolor[0]));
		glRectf(xmin + xStart, ymin + yStart, xmax + xStart, ymax + yStart);
	}
	
	// Draw the text at selected start location:
	if (faceT) {
		faceT->draw(xStart, yStart, uniCodeText);
	}
	else {
		faceM->draw(xStart, yStart, uniCodeText);
	}
	
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glDisable( GL_TEXTURE_2D );
	glPopAttrib();
	glPopClientAttrib();
	
	// Ready!
	return(0);	
}

int PsychMeasureText(int textLen, double* text, float* xmin, float* ymin, float* xmax, float* ymax)
{
	int i;
	QChar* myUniChars = new QChar[textLen];
	
	// Check if rebuild of font face needed due to parameter
	// chage. Reload/Rebuild font face if so, check for errors:
	if (_needsRebuild && PsychRebuildFont()) return(1);

	// Synthesize Unicode QString from double vector:
	for(i = 0; i < textLen; i++) {
		myUniChars[i] = QChar((unsigned int) text[i]); 
	}	
	QString	uniCodeText = QString(myUniChars, textLen);  
	delete [] myUniChars;

	// Compute its bounding box:
	OGLFT::BBox box = (faceT) ? faceT->measure(uniCodeText) : faceM->measure(uniCodeText);
	
	*xmin = box.x_min_;
	*ymin = box.y_min_;
	*xmax = box.x_max_;
	*ymax = box.y_max_;

	return(0);
}

int PsychInitText(void)
{
	_firstCall = true;
	_needsRebuild = true;
	faceT = NULL;
	faceM = NULL;
	ft_face = NULL;

	// Try to initialize libfontconfig - our fontMapper library for font matching and selection:
	if (!FcInit()) {
		if (_verbosity > 0) fprintf(stderr, "libptbdrawtext_ftgl: FontMapper initialization failed!\n");
		return(1);
	}
	
	if (_verbosity > 2)	{
		fprintf(stderr, "libptbdrawtext_ftgl: External 'DrawText' text rendering plugin initialized.\n");
		fprintf(stderr, "libptbdrawtext_ftgl: This plugin uses multiple excellent free software libraries to do its work:\n");
		fprintf(stderr, "libptbdrawtext_ftgl: OGLFT (http://oglft.sourceforge.net/) the OpenGL-FreeType library.\n");
		fprintf(stderr, "libptbdrawtext_ftgl: The FreeType-2 (http://freetype.sourceforge.net/) library.\n");
		fprintf(stderr, "libptbdrawtext_ftgl: The FontConfig (http://www.fontconfig.org) library.\n");
		fprintf(stderr, "libptbdrawtext_ftgl: Thanks!\n\n");
	}

	return(0);
}

int PsychShutdownText(void)
{
	if (faceT || faceM) {
		if (_verbosity > 3) fprintf(stderr, "libptbdrawtext_ftgl: In shutdown: faceT = %p faceM = %p\n", faceT, faceM);
	
		// Delete OGLFT face objects:
		if (faceT) delete(faceT);
		faceT = NULL;
		
		if (faceM) delete(faceM);
		faceM = NULL;

		// Delete Freetype face object:
		if (ft_face) FT_Done_Face(ft_face);
		ft_face = NULL;
		if (_verbosity > 3) fprintf(stderr, "libptbdrawtext_ftgl: Shutting down.\n");
	}
	
	_needsRebuild = true;
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
