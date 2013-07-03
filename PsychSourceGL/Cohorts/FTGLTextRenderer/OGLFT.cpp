/*
 * OGLFT: A library for drawing text with OpenGL using the FreeType library
 * Copyright (C) 2002 lignum Computing, Inc. <oglft@lignumcomputing.com>
 * Copyright (C) 2008 Allen Barnett
 * $Id:$
 *
 * Contains small modifications by Mario Kleiner. Search for text "MK" for
 * comments on the minimal changes needed to provide Unicode support while
 * removing dependencies on the QT toolkit.
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

#include <iostream>
#include <iomanip>

#include "OGLFT.h"

#ifndef OGLFT_NO_QT
// MK changed:
#if OGLFT_QT_VERSION != 0
#include <qregexp.h>
#endif
#endif
#include <freetype/ftbitmap.h>

namespace OGLFT {

  // This is the static instance of the FreeType library wrapper ...

  Library Library::library;

  // ... and this is the FreeType library handle itself.

  FT_Library Library::library_;

  // The static instance above causes this constructor to be called
  // when the object module is loaded.

  Library::Library ( void )
  {

    FT_Error error = FT_Init_FreeType( &library_ );

    if ( error != 0 ) {
      std::cerr << "Could not initialize the FreeType library. Exiting." << std::endl;
      exit( 1 );
    }
  }

  Library::~Library ( void )
  {
    FT_Error error = FT_Done_FreeType( library_ );

    if ( error != 0 ) {
      std::cerr << "Could not terminate the FreeType library." << std::endl;
    }
  }

  // Return the only instance in the process

  FT_Library& Library::instance ( void )
  {
    return library_;
  }

  // Load a new face from file

  Face::Face ( const char* filename, float point_size, FT_UInt resolution )
    : point_size_( point_size ), resolution_( resolution )
  {
    valid_ = true; // Assume the best :-)

    FT_Face ft_face;

    FT_Error error = FT_New_Face( Library::instance(), filename, 0, &ft_face );

    if ( error != 0 ) {
      valid_ = false;
      return;
    }

    // As of FreeType 2.1: only a UNICODE charmap is automatically activated.
    // If no charmap is activated automatically, just use the first one.
    if ( ft_face->charmap == 0 && ft_face->num_charmaps > 0 )
      FT_Select_Charmap( ft_face, ft_face->charmaps[0]->encoding );

    faces_.push_back( FaceData( ft_face ) );

    init();
  }

  // Load a new face from memory

  Face::Face ( const FT_Byte* data_base, const FT_Long data_size, float point_size, FT_UInt resolution )
    : point_size_( point_size ), resolution_( resolution )
  {
    valid_ = true; // Assume the best :-)

    FT_Face ft_face;

    FT_Error error = FT_New_Memory_Face( Library::instance(), data_base, data_size, 0, &ft_face );

    if ( error != 0 ) {
      valid_ = false;
      return;
    }

    // As of FreeType 2.1: only a UNICODE charmap is automatically activated.
    // If no charmap is activated automatically, just use the first one.
    if ( ft_face->charmap == 0 && ft_face->num_charmaps > 0 )
      FT_Select_Charmap( ft_face, ft_face->charmaps[0]->encoding );

    faces_.push_back( FaceData( ft_face ) );

    init();
  }

  // Go with a face that the user has already opened.

  Face::Face ( FT_Face face, float point_size, FT_UInt resolution )
    : point_size_( point_size ), resolution_( resolution )
  {
    valid_ = true;

    // As of FreeType 2.1: only a UNICODE charmap is automatically activated.
    // If no charmap is activated automatically, just use the first one.
    if ( face->charmap == 0 && face->num_charmaps > 0 )
      FT_Select_Charmap( face, face->charmaps[0]->encoding );

    faces_.push_back( FaceData( face, false ) );

    init();
  }

  // Standard initialization behavior once the font file is opened.

  void Face::init ( void )
  {
    // By default, each glyph is compiled into a display list the first
    // time it is encountered

    compile_mode_ = COMPILE;

    // By default, all drawing is wrapped with push/pop matrix so that the
    // MODELVIEW matrix is not modified. If advance_ is set, then subsequent
    // drawings follow from the advance of the last glyph rendered.

    advance_ = false;

    // Initialize the default colors

    foreground_color_[R] = 0.;
    foreground_color_[G] = 0.;
    foreground_color_[B] = 0.;
    foreground_color_[A] = 1.;
    
    background_color_[R] = 1.;
    background_color_[G] = 1.;
    background_color_[B] = 1.;
    background_color_[A] = 0.;

    // The default positioning of the text is at the origin of the first glyph
    horizontal_justification_ = ORIGIN;
    vertical_justification_ = BASELINE;

    // By default, strings are rendered in their nominal direction
    string_rotation_ = 0;

    // setCharacterRotationReference calls the virtual function clearCaches()
    // so it is up to a subclass to set the real default
    rotation_reference_glyph_ = 0;
    rotation_reference_face_ = 0;
    rotation_offset_y_ = 0.;
  }

  Face::~Face ( void )
  {
    for ( unsigned int i = 0; i < faces_.size(); i++ )
      if ( faces_[i].free_on_exit_ )
	FT_Done_Face( faces_[i].face_ );
  }

  // Add another Face to select characters from file.

  bool Face::addAuxiliaryFace ( const char* filename )
  {
    FT_Face ft_face;

    FT_Error error = FT_New_Face( Library::instance(), filename, 0, &ft_face );

    if ( error != 0 )
      return false;

    faces_.push_back( FaceData( ft_face ) );

    setCharSize();

    return true;
  }

  // Add another Face to select characters from memory.

  bool Face::addAuxiliaryFace ( const FT_Byte* data_base, const FT_Long data_size )
  {
    FT_Face ft_face;

    FT_Error error = FT_New_Memory_Face( Library::instance(), data_base, data_size, 0, &ft_face );

    if ( error != 0 )
      return false;

    faces_.push_back( FaceData( ft_face ) );

    setCharSize();

    return true;
  }

  // Add another Face to select characters from (face)

  bool Face::addAuxiliaryFace ( FT_Face face )
  {
    faces_.push_back( FaceData( face, false ) );

    setCharSize();

    return true;
  }

  // Note: Changing the point size also clears the display list cache

  void Face::setPointSize ( float point_size )
  {
    if ( point_size != point_size_ ) {

      point_size_ = point_size;

      clearCaches();

      setCharSize();
    }
  }

  // Note: Changing the resolution also clears the display list cache

  void Face::setResolution ( FT_UInt resolution )
  {
    if ( resolution != resolution_ ) {

      resolution_ = resolution;

      clearCaches();

      setCharSize();
    }
  }

  // Note: Changing the background color also clears the display list cache.

  void Face::setBackgroundColor ( GLfloat red, GLfloat green, GLfloat blue,
				  GLfloat alpha )
  {
    if ( background_color_[R] != red ||
	 background_color_[G] != green ||
	 background_color_[B] != blue ||
	 background_color_[A] != alpha ) {

      background_color_[R] = red;
      background_color_[G] = green;
      background_color_[B] = blue;
      background_color_[A] = alpha;

      clearCaches();
    }
  }

  // Note: Changing the foreground color also clears the display list cache.

  void Face::setForegroundColor ( GLfloat red, GLfloat green, GLfloat blue,
				  GLfloat alpha )
  {
    if ( foreground_color_[R] != red ||
	 foreground_color_[G] != green ||
	 foreground_color_[B] != blue ||
	 foreground_color_[A] != alpha ) {

      foreground_color_[R] = red;
      foreground_color_[G] = green;
      foreground_color_[B] = blue;
      foreground_color_[A] = alpha;

      clearCaches();
    }
  }

  // Note: Changing the foreground color also clears the display list cache.

  void Face::setForegroundColor ( const GLfloat foreground_color[4] )
  {
    if ( foreground_color_[R] != foreground_color[R] ||
	 foreground_color_[G] != foreground_color[G] ||
	 foreground_color_[B] != foreground_color[B] ||
	 foreground_color_[A] != foreground_color[A] ) {

      foreground_color_[R] = foreground_color[R];
      foreground_color_[G] = foreground_color[G];
      foreground_color_[B] = foreground_color[B];
      foreground_color_[A] = foreground_color[A];

      clearCaches();
    }
  }

  // Note: Changing the background color also clears the display list cache.

  void Face::setBackgroundColor ( const GLfloat background_color[4] )
  {
    if ( background_color_[R] != background_color[R] ||
	 background_color_[G] != background_color[G] ||
	 background_color_[B] != background_color[B] ||
	 background_color_[A] != background_color[A] ) {

      background_color_[R] = background_color[R];
      background_color_[G] = background_color[G];
      background_color_[B] = background_color[B];
      background_color_[A] = background_color[A];

      clearCaches();
    }
  }
#ifndef OGLFT_NO_QT
  // Note: Changing the foreground color also clears the display list cache.

  void Face::setForegroundColor ( const QRgb foreground_rgba )
  {
    GLfloat foreground_color[4];
    foreground_color[R] = qRed( foreground_rgba ) / 255.;
    foreground_color[G] = qGreen( foreground_rgba ) / 255.;
    foreground_color[B] = qBlue( foreground_rgba ) / 255.;
    foreground_color[A] = qAlpha( foreground_rgba ) / 255.;

    if ( foreground_color_[R] != foreground_color[R] ||
	 foreground_color_[G] != foreground_color[G] ||
	 foreground_color_[B] != foreground_color[B] ||
	 foreground_color_[A] != foreground_color[A] ) {

      foreground_color_[R] = foreground_color[R];
      foreground_color_[G] = foreground_color[G];
      foreground_color_[B] = foreground_color[B];
      foreground_color_[A] = foreground_color[A];

      clearCaches();
    }
  }

  // Note: Changing the background color also clears the display list cache.

  void Face::setBackgroundColor ( const QRgb background_rgba )
  {
    GLfloat background_color[4];
    background_color[R] = qRed( background_rgba ) / 255.;
    background_color[G] = qGreen( background_rgba ) / 255.;
    background_color[B] = qBlue( background_rgba ) / 255.;
    background_color[A] = qAlpha( background_rgba ) / 255.;

    if ( background_color_[R] != background_color[R] ||
	 background_color_[G] != background_color[G] ||
	 background_color_[B] != background_color[B] ||
	 background_color_[A] != background_color[A] ) {

      background_color_[R] = background_color[R];
      background_color_[G] = background_color[G];
      background_color_[B] = background_color[B];
      background_color_[A] = background_color[A];

      clearCaches();
    }
  }
#endif /* OGLFT_NO_QT */
  // Note: Changing the string rotation angle clears the display list cache

  void Face::setStringRotation ( GLfloat string_rotation )
  {
    if ( string_rotation != string_rotation_ ) {
      string_rotation_ = string_rotation;

      clearCaches();

      // Note that this affects ALL glyphs accessed through
      // the Face, both the vector and the raster glyphs. Very nice!

      if ( string_rotation_ != 0. ) {
	float angle;
	if ( string_rotation_ < 0. ) {
	  angle = 360.f - fmod( fabs( string_rotation_ ), 360.f );
	}
	else {
	  angle = fmod( string_rotation_, 360.f );
	}

	FT_Matrix rotation_matrix;
	FT_Vector sinus;

	FT_Vector_Unit( &sinus, (FT_Angle)(angle * 0x10000L) );

	rotation_matrix.xx = sinus.x;
	rotation_matrix.xy = -sinus.y;
	rotation_matrix.yx = sinus.y;
	rotation_matrix.yy = sinus.x;

	for ( unsigned int i = 0; i < faces_.size(); i++ )
	  FT_Set_Transform( faces_[i].face_, &rotation_matrix, 0 );
      }
      else
	for ( unsigned int i = 0; i < faces_.size(); i++ )
	  FT_Set_Transform( faces_[i].face_, 0, 0 );
    }
  }

  // Note: Changing the rotation reference character clears the display list cache.

  void Face::setCharacterRotationReference ( unsigned char c )
  {
    unsigned int f;
    FT_UInt glyph_index = 0;

    for ( f = 0; f < faces_.size(); f++ ) {
      glyph_index = FT_Get_Char_Index( faces_[f].face_, c );
      if ( glyph_index != 0 ) break;
    }

    if ( f < faces_.size() && glyph_index != rotation_reference_glyph_ ) {

      FT_Error error = FT_Load_Glyph( faces_[f].face_, glyph_index,
				      FT_LOAD_DEFAULT );

      if ( error != 0 ) return;

      rotation_reference_glyph_ = glyph_index;

      rotation_reference_face_ = faces_[f].face_;

      setRotationOffset();

      clearCaches();
    }
  }

  BBox Face::measure ( const char* s )
  {
    BBox bbox;
    char c;

    if ( ( c = *s++ ) != 0 ) {

      bbox = measure( c );

      for ( c = *s; c != 0; c = *++s ) {

	BBox char_bbox = measure( c );

	bbox += char_bbox;
      }
    }

    return bbox;
  }

  BBox Face::measureRaw ( const char* s )
  {
    BBox bbox;

    for ( char c = *s; c != 0; c = *++s ) {
      BBox char_bbox;

      unsigned int f;
      FT_UInt glyph_index = 0;

      for ( f = 0; f < faces_.size(); f++ ) {
	glyph_index = FT_Get_Char_Index( faces_[f].face_, c );
	if ( glyph_index != 0 ) break;
      }

      if ( glyph_index == 0 ) continue;

      FT_Error error = FT_Load_Glyph( faces_[f].face_, glyph_index,
				      FT_LOAD_DEFAULT );
      if ( error != 0 ) continue;

      FT_Glyph glyph;
      error = FT_Get_Glyph( faces_[f].face_->glyph, &glyph );
      if ( error != 0 ) continue;

      FT_BBox ft_bbox;
      FT_Glyph_Get_CBox( glyph, ft_glyph_bbox_unscaled, &ft_bbox );

      FT_Done_Glyph( glyph );

      char_bbox = ft_bbox;
      char_bbox.advance_ = faces_[f].face_->glyph->advance;

      bbox += char_bbox;
    }

    return bbox;
  }

#ifndef OGLFT_NO_QT
  BBox Face::measure ( const QString& s )
  {
    BBox bbox;

    if ( s.length() > 0 ) {

      bbox = measure( s.at( 0 ) );

      for ( unsigned int i = 1; i < s.length(); i++ ) {

	BBox char_bbox = measure( s.at( i ) );

	bbox += char_bbox;
      }
    }

    return bbox;
  }

  BBox Face::measure ( const QString& format, double number )
  {
    return measure( format_number( format, number ) );
  }

  BBox Face::measureRaw ( const QString& s )
  {
    BBox bbox;

    for ( unsigned int i = 0; i < s.length(); i++ ) {
      BBox char_bbox;

      unsigned int f;
      FT_UInt glyph_index = 0;

      for ( f = 0; f < faces_.size(); f++ ) {
	glyph_index = FT_Get_Char_Index( faces_[f].face_, s.at( i ).unicode() );
	if ( glyph_index != 0 ) break;
      }

      if ( glyph_index == 0 ) {
	continue;
      }

      FT_Error error = FT_Load_Glyph( faces_[f].face_, glyph_index,
				      FT_LOAD_DEFAULT );
      if ( error != 0 ) continue;

      FT_Glyph glyph;
      error = FT_Get_Glyph( faces_[f].face_->glyph, &glyph );
      if ( error != 0 ) continue;

      FT_BBox ft_bbox;
      FT_Glyph_Get_CBox( glyph, ft_glyph_bbox_unscaled, &ft_bbox );

      FT_Done_Glyph( glyph );

      char_bbox = ft_bbox;
      char_bbox.advance_ = faces_[f].face_->glyph->advance;

      bbox += char_bbox;
    }

    return bbox;
  }
#endif /* OGLFT_NO_QT */

  // Measure the bounding box as if the (latin1) string were not rotated

  BBox Face::measure_nominal ( const char* s )
  {
    if ( string_rotation_ == 0. )
      return measure( s );
    
    for ( unsigned int f = 0; f < faces_.size(); f++ )
      FT_Set_Transform( faces_[f].face_, 0, 0 );

    BBox bbox = measure( s );

    float angle;
    if ( string_rotation_ < 0. ) {
      angle = 360.f - fmod( fabs( string_rotation_ ), 360.f );
    }
    else {
      angle = fmod( string_rotation_, 360.f );
    }

    FT_Matrix rotation_matrix;
    FT_Vector sinus;

    FT_Vector_Unit( &sinus, (FT_Angle)(angle * 0x10000L) );

    rotation_matrix.xx = sinus.x;
    rotation_matrix.xy = -sinus.y;
    rotation_matrix.yx = sinus.y;
    rotation_matrix.yy = sinus.x;

    for ( unsigned int f = 0; f < faces_.size(); f++ )
      FT_Set_Transform( faces_[f].face_, &rotation_matrix, 0 );

    return bbox;
  }

#ifndef OGLFT_NO_QT
  // Measure the bounding box as if the (UNICODE) string were not rotated

  BBox Face::measure_nominal ( const QString& s )
  {
    if ( string_rotation_ == 0. )
      return measure( s );
    
    for ( unsigned int f = 0; f < faces_.size(); f++ )
      FT_Set_Transform( faces_[f].face_, 0, 0 );

    BBox bbox = measure( s );

    float angle;
    if ( string_rotation_ < 0. ) {
      angle = 360. - fmod( fabs( string_rotation_ ), 360.f );
    }
    else {
      angle = fmod( string_rotation_, 360.f );
    }

    FT_Matrix rotation_matrix;
    FT_Vector sinus;

    FT_Vector_Unit( &sinus, (FT_Angle)(angle * 0x10000L) );

    rotation_matrix.xx = sinus.x;
    rotation_matrix.xy = -sinus.y;
    rotation_matrix.yx = sinus.y;
    rotation_matrix.yy = sinus.x;

    for ( unsigned int f = 0; f < faces_.size(); f++ )
      FT_Set_Transform( faces_[f].face_, &rotation_matrix, 0 );

    return bbox;
  }

// MK changed:
#if OGLFT_QT_VERSION != 0

  // Format the number per the given format. Mostly pointless
  // for the standard formats, e.g. %12e. You can use the regular
  // Qt functions to format such a string and avoid the parsing
  // which is done here.

  QString Face::format_number ( const QString& format, double number )
  {
    // This regexp says:
    // 1. optionally match any thing up to a format,
    // 2. the optional format (%...), and
    // 3. optionally anything after it.
    // Note that since everything is optional, the match always succeeds.
    QRegExp format_regexp("((?:[^%]|%%)*)(%[0-9]*\\.?[0-9]*[efgp])?((?:[^%]|%%)*)");
#if OGLFT_QT_VERSION == 3
    /*int pos = */ format_regexp.search( format );
#elif OGLFT_QT_VERSION == 4
    /*int pos = */ format_regexp.exactMatch( format );
#endif

    QStringList list = format_regexp.capturedTexts();

    QStringList::Iterator it = list.begin();

#if OGLFT_QT_VERSION == 3
    it = list.remove( it );	// Remove the "matched" string, leaving the pieces
#elif OGLFT_QT_VERSION == 4
    it = list.erase( it );	// Remove the "matched" string, leaving the pieces
#endif

    if ( it == list.end() ) return QString::null; // Probably an error

    // Extract each piece from the list

    QString prefix, value_format, postfix;
    char type = '\0';

    if ( !(*it).isEmpty() )
      prefix = *it;

    ++it;

    if ( it != list.end() ) {
      if ( !(*it).isEmpty() ) {
	// Reparse this to extract the details of the format
	QRegExp specifier_regexp( "([0-9]*)\\.?([0-9]*)([efgp])" );
#if OGLFT_QT_VERSION == 3
	(void)specifier_regexp.search( *it );
#elif OGLFT_QT_VERSION == 4
	(void)specifier_regexp.exactMatch( *it );
#endif
	QStringList specifier_list = specifier_regexp.capturedTexts();

	QStringList::Iterator sit = specifier_list.begin();
#if OGLFT_QT_VERSION == 3
	sit = specifier_list.remove( sit );
#elif OGLFT_QT_VERSION == 4
	sit = specifier_list.erase( sit );
#endif
	int width = (*sit).toInt();
	++sit;
	int precision = (*sit).toInt();
	++sit;
#if OGLFT_QT_VERSION == 3
	type = (*sit).at(0).latin1();
#elif OGLFT_QT_VERSION == 4
	type = (*sit).at(0).toLatin1();
#endif
	// The regular formats just use Qt's number formatting capability
	if ( type == 'e' || type == 'f' || type == 'g' )
	  value_format = QString( "%1" ).arg( number, width, type, precision );

	// For the fraction, though, we have to convert it the special
	// UNICODE encoding
	else if ( type == 'p' ) {
	  // Fixed for now...
	  if ( fabs( number ) < 1./256. )
	    value_format = "0";
	  else {
	    // Extract the integral part
	    int a = (int)number;

	    if ( a != 0 )
	      value_format = QString::number( a );

	    // Extract the fractional part: NOTE: THIS IS LIMITED TO
	    // REPRESENTING ALL FRACTIONS AS n/256
	    int b = (int)rint( 256. * fabs( number - a ) );

	    // If b is exactly 256, then the original number was
	    // essentially an integer (to within 1/256-th)
	    if ( b == 256 )
	      value_format = QString::number( rint( number ) );

	    else if ( b != 0 ) {
	      int c = 256;
	      // Remove common factors of two from the numerator and denominator
	      for ( ; ( b & 0x1 ) == 0; b >>= 1, c >>= 1 );

	      // Format the numerator and shift to 0xE000 sequence
	      QString numerator = QString::number( b );
	      for ( uint i = 0; i < numerator.length(); i++ ) {
		numerator[i] = QChar( numerator.at(i).unicode() -
				      QChar('0').unicode() +
				      0xE000 );
	      }
	      value_format += numerator;
	      value_format += QChar( 0xE00a ); // The '/'
	      // Format the denominator and shift to 0xE010 sequence
	      QString denominator = QString::number( c );
	      for ( uint i = 0; i < denominator.length(); i++ ) {
		denominator[i] = QChar( denominator.at(i).unicode() -
					QChar('0').unicode() +
					0xE010 );
	      }
	      value_format += denominator;
	    }
	  }
	}
      }

      ++it;
      
      if ( it != list.end() && !(*it).isEmpty() )
	postfix = *it;
    }

    return prefix + value_format + postfix;
  }

// MK changed:
#else
// Noop function to make linker happy:
QString Face::format_number ( const QString& format, double number ) { return( QString() ); }
#endif /* #if OGLFT_QT_VERSION != 0 */

#endif /* OGLFT_NO_QT */

  // Compile a (latin1) string into a display list

  GLuint Face::compile ( const char* s )
  {
    // First, make sure all the characters in the string are themselves
    // in display lists
    const char* s_tmp = s;

    for ( char c = *s_tmp; c != 0; c = *++s_tmp ) {
      compile( c );
    }
    
    GLuint dlist = glGenLists( 1 );
    glNewList( dlist, GL_COMPILE );

    glColor4f( foreground_color_[R], foreground_color_[G], foreground_color_[B],
	       foreground_color_[A] );
    if ( !advance_ )
      glPushMatrix();

    draw( s );

    if ( !advance_ )
      glPopMatrix();

    glEndList();

    return dlist;
  }
#ifndef OGLFT_NO_QT
  // Compile a (UNICODE) string into a display list

  GLuint Face::compile ( const QString& s )
  {
    // First, make sure all the characters in the string are themselves
    // in display lists
    for ( unsigned int i = 0; i < s.length(); i++ ) {
      compile( s.at( i ) );
    }
    
    GLuint dlist = glGenLists( 1 );
    glNewList( dlist, GL_COMPILE );

    glColor4f( foreground_color_[R], foreground_color_[G], foreground_color_[B],
	       foreground_color_[A] );
    if ( !advance_ )
      glPushMatrix();

    draw( s );

    if ( !advance_ )
      glPopMatrix();

    glEndList();

    return dlist;
  }
#endif /* OGLFT_NO_QT */
  // Compile a (latin1) character glyph into a display list and cache
  // it for later

  GLuint Face::compile ( unsigned char c )
  {
    // See if we've done it already

    GDLCI fgi = glyph_dlists_.find( c );

    if ( fgi != glyph_dlists_.end() )
      return fgi->second;

    unsigned int f;
    FT_UInt glyph_index = 0;

    for ( f = 0; f < faces_.size(); f++ ) {
      glyph_index = FT_Get_Char_Index( faces_[f].face_, c );
      if ( glyph_index != 0 ) break;
    }

    if ( glyph_index == 0 )
      return 0;

    GLuint dlist = compileGlyph( faces_[f].face_, glyph_index );

    glyph_dlists_[ c ] = dlist;

    return dlist;
  }

#ifndef OGLFT_NO_QT
  // Compile a (UNICODE) character glyph into a display list and cache
  // it for later

  GLuint Face::compile ( const QChar c )
  {
    // See if we've done it already

    GDLCI fgi = glyph_dlists_.find( c.unicode() );

    if ( fgi != glyph_dlists_.end() )
      return fgi->second;

    unsigned int f;
    FT_UInt glyph_index = 0;

    for ( f = 0; f < faces_.size(); f++ ) {
      glyph_index = FT_Get_Char_Index( faces_[f].face_, c.unicode() );
      if ( glyph_index != 0 ) break;
    }

    if ( glyph_index == 0 )
      return 0;

    GLuint dlist = compileGlyph( faces_[f].face_, glyph_index );

    glyph_dlists_[ c.unicode() ] = dlist;

    return dlist;
  }
#endif /* OGLFT_NO_QT */
  // Assume the MODELVIEW matrix is already set and draw the (latin1)
  // string.  Note: this routine now ignores almost all settings:
  // including the position (both modelview and raster), color,
  // justification and advance settings. Consider this to be the raw
  // drawing routine for which you are responsible for most of the
  // setup.

  void Face::draw ( const char* s )
  {
    DLCI character_display_list = character_display_lists_.begin();

    for ( char c = *s; c != 0; c = *++s ) {

      if ( character_display_list != character_display_lists_.end() ) {
	glCallList( *character_display_list );
	character_display_list++;
      }

      draw( c );
    }
  }
#ifndef OGLFT_NO_QT
  // Assume the MODELVIEW matrix is already set and draw the (UNICODE)
  // string.  Note: this routine now ignores almost all settings:
  // including the position (both modelview and raster), color,
  // justification and advance settings. Consider this to be the raw
  // drawing routine for which you are responsible for most of the
  // setup.

  void Face::draw ( const QString& s )
  {
    DLCI character_display_list = character_display_lists_.begin();

    for ( unsigned int i = 0; i < s.length(); i++ ) {

      if ( character_display_list != character_display_lists_.end() ) {
	glCallList( *character_display_list );
	character_display_list++;
      }

      draw( s.at( i ) );
    }
  }
#endif /* OGLFT_NO_QT */

  // Assume the MODELVIEW matrix is already setup and draw the
  // (latin1) character.

  void Face::draw ( unsigned char c )
  {
    // See if we've done it already

    GDLCI fgi = glyph_dlists_.find( c );

    if ( fgi != glyph_dlists_.end( ) ) {
      glCallList( fgi->second );
      return;
    }

    unsigned int f;
    FT_UInt glyph_index = 0;

    for ( f = 0; f < faces_.size(); f++ ) {
      glyph_index = FT_Get_Char_Index( faces_[f].face_, c );
      if ( glyph_index != 0 ) break;
    }

    if ( glyph_index == 0 )
      return;

    // Otherwise, either compile it (and call it) or ...

    else if ( compile_mode_ == COMPILE ) {
      GLuint dlist = compile( c );
      glCallList( dlist );
    }

    // ... render it immediately

    else {
      renderGlyph( faces_[f].face_, glyph_index );
    }
  }
#ifndef OGLFT_NO_QT
  // Assume the MODELVIEW matrix is already setup and draw the
  // (UNICODE) character.

  void Face::draw ( const QChar c )
  {
    // See if we've done it already

    GDLCI fgi = glyph_dlists_.find( c.unicode() );

    if ( fgi != glyph_dlists_.end( ) ) {
      glCallList( fgi->second );
      return;
    }

    unsigned int f;
    FT_UInt glyph_index = 0;

    for ( f = 0; f < faces_.size(); f++ ) {
      glyph_index = FT_Get_Char_Index( faces_[f].face_, c.unicode() );
      if ( glyph_index != 0 ) {
	break;
      }
    }

    if ( glyph_index == 0 )
      return;

    // Otherwise, either compile it (and call it) or ...

    if ( compile_mode_ == COMPILE ) {
      GLuint dlist = compile( c );
      glCallList( dlist );
    }

    // ... render it immediately

    else {
      renderGlyph( faces_[f].face_, glyph_index );
    }
  }
#endif /* OGLFT_NO_QT */
  // Draw the (latin1) character at the given position. The MODELVIEW
  // matrix is modified by the glyph advance.

  void Face::draw ( GLfloat x, GLfloat y, unsigned char c )
  {
    glTranslatef( x, y, 0. );

    glColor4f( foreground_color_[R], foreground_color_[G], foreground_color_[B],
	       foreground_color_[A] );

    glRasterPos2i( 0, 0 );

    draw( c );
  }
#ifndef OGLFT_NO_QT
  // Draw the (UNICODE) character at the given position. The MODELVIEW
  // matrix is modified by the glyph advance.

  void Face::draw ( GLfloat x, GLfloat y, QChar c )
  {
    glTranslatef( x, y, 0. );

    glColor4f( foreground_color_[R], foreground_color_[G], foreground_color_[B],
	       foreground_color_[A] );

    glRasterPos2i( 0, 0 );

    draw( c );
  }
#endif /* OGLFT_NO_QT */
  // Draw the (latin1) string at the given position.

  void Face::draw ( GLfloat x, GLfloat y, const char* s )
  {
    if ( !advance_ )
      glPushMatrix();

    if ( horizontal_justification_ != ORIGIN ||
	 vertical_justification_ != BASELINE ) {
      glPushMatrix();

      BBox bbox = measure_nominal( s );

      GLfloat dx = 0, dy = 0;

      switch ( horizontal_justification_ ) {
      case LEFT:
	dx = -bbox.x_min_; break;
      case CENTER:
	dx = -( bbox.x_min_ + bbox.x_max_ ) / 2.f; break;
      case RIGHT:
	dx = -bbox.x_max_; break;
      default:
	break;
      }
      switch ( vertical_justification_ ) {
      case BOTTOM:
	dy = -bbox.y_min_; break;
      case MIDDLE:
	dy = -( bbox.y_min_ + bbox.y_max_ ) / 2.f; break;
      case TOP:
	dy = -bbox.y_max_; break;
      default:
	break;
      }

      // There is probably a less expensive way to compute this

      glRotatef( string_rotation_, 0., 0., 1. );
      glTranslatef( dx, dy, 0 );
      glRotatef( -string_rotation_, 0., 0., 1. );
    }

    glTranslatef( x, y, 0. );

    glColor4f( foreground_color_[R], foreground_color_[G], foreground_color_[B],
	       foreground_color_[A] );

    glRasterPos2i( 0, 0 );

    draw( s );

    if ( horizontal_justification_ != ORIGIN ||
	 vertical_justification_ != BASELINE )
      glPopMatrix();

    if ( !advance_ )
      glPopMatrix();
  }

#ifndef OGLFT_NO_QT
  // Draw the (UNICODE) string at the given position.

  void Face::draw ( GLfloat x, GLfloat y, const QString& s )
  {
    if ( !advance_ )
      glPushMatrix();

    if ( horizontal_justification_ != ORIGIN ||
	 vertical_justification_ != BASELINE ) {
      glPushMatrix();

      BBox bbox = measure_nominal( s );

      GLfloat dx = 0, dy = 0;

      switch ( horizontal_justification_ ) {
      case LEFT:
	dx = -bbox.x_min_; break;
      case CENTER:
	dx = -( bbox.x_min_ + bbox.x_max_ ) / 2.; break;
      case RIGHT:
	dx = -bbox.x_max_; break;
      case ORIGIN:
	break;
      }
      switch ( vertical_justification_ ) {
      case BOTTOM:
	dy = -bbox.y_min_; break;
      case MIDDLE:
	dy = -( bbox.y_min_ + bbox.y_max_ ) / 2.; break;
      case TOP:
	dy = -bbox.y_max_; break;
      case BASELINE:
	break;
      }

      // There is probably a less expensive way to compute this

      glRotatef( string_rotation_, 0., 0., 1. );
      glTranslatef( dx, dy, 0 );
      glRotatef( -string_rotation_, 0., 0., 1. );
    }

    glTranslatef( x, y, 0. );

    glColor4f( foreground_color_[R], foreground_color_[G], foreground_color_[B],
	       foreground_color_[A] );

    glRasterPos2i( 0, 0 );

    draw( s );

    if ( horizontal_justification_ != ORIGIN ||
	 vertical_justification_ != BASELINE )
      glPopMatrix();

    if ( !advance_ )
      glPopMatrix();
  }

  // Draw the number at the given position per the given format.

  void Face::draw ( GLfloat x, GLfloat y, const QString& format, double number )
  {
    draw( x, y, format_number( format, number ) );
  }
#endif /* OGLFT_NO_QT */

  Raster::Raster ( const char* filename, float point_size, FT_UInt resolution )
    : Face( filename, point_size, resolution )
  {
    if ( !isValid() ) return;

    init();
  }

  Raster::Raster ( const FT_Byte* data_base, const FT_Long data_size,
		   float point_size, FT_UInt resolution )
    : Face( data_base, data_size, point_size, resolution )
  {
    if ( !isValid() ) return;

    init();
  }

  Raster::Raster ( FT_Face face, float point_size, FT_UInt resolution )
    : Face( face, point_size, resolution )
  {
    init();
  }

  void Raster::init ( void )
  {
    character_rotation_z_ = 0;

    setCharSize();

    setCharacterRotationReference( 'o' );
  }

  Raster::~Raster ( void )
  {
    clearCaches();
  }

  void Raster::setCharacterRotationZ ( GLfloat character_rotation_z )
  {
    if ( character_rotation_z != character_rotation_z_ ) {
      character_rotation_z_ = character_rotation_z;

      clearCaches();
    }
  }

  double Raster::height ( void ) const
  {
    if ( faces_[0].face_->height > 0 )
      return faces_[0].face_->height / 64.;
    else
      return faces_[0].face_->size->metrics.y_ppem;
  }

  BBox Raster::measure ( unsigned char c )
  {
    BBox bbox;
    // For starters, just get the unscaled glyph bounding box
    unsigned int f;
    FT_UInt glyph_index = 0;

    for ( f = 0; f < faces_.size(); f++ ) {
      glyph_index = FT_Get_Char_Index( faces_[f].face_, c );
      if ( glyph_index != 0 ) break;
    }

    if ( glyph_index == 0 )
      return bbox;

    FT_Error error = FT_Load_Glyph( faces_[f].face_, glyph_index,
				    FT_LOAD_DEFAULT );
    if ( error != 0 )
      return bbox;

    FT_Glyph glyph;
    error = FT_Get_Glyph( faces_[f].face_->glyph, &glyph );
    if ( error != 0 )
      return bbox;

    FT_BBox ft_bbox;
    FT_Glyph_Get_CBox( glyph, ft_glyph_bbox_unscaled, &ft_bbox );

    FT_Done_Glyph( glyph );

    bbox = ft_bbox;
    bbox.advance_ = faces_[f].face_->glyph->advance;

    // In order to be accurate regarding the placement of text not
    // aligned at the glyph's origin (CENTER/MIDDLE), the bounding box
    // of the raster format has to be projected back into the
    // view's coordinates

    GLint viewport[4];
    GLdouble modelview[16], projection[16];

    glGetIntegerv( GL_VIEWPORT, viewport );
    glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
    glGetDoublev( GL_PROJECTION_MATRIX, projection );

    // Well, first we have to get the Origin, since that is the basis
    // of the bounding box
    GLdouble x0, y0, z0;
    gluUnProject( 0., 0., 0., modelview, projection, viewport, &x0, &y0, &z0 );

    GLdouble x, y, z;
    gluUnProject( bbox.x_min_, bbox.y_min_, 0., modelview, projection, viewport,
		  &x, &y, &z );
    bbox.x_min_ = (float)( x - x0 );
    bbox.y_min_ = (float)( y - y0 );

    gluUnProject( bbox.x_max_, bbox.y_max_, 0., modelview, projection, viewport,
		  &x, &y, &z );
    bbox.x_max_ = (float)( x - x0 );
    bbox.y_max_ = (float)( y - y0 );

    gluUnProject( bbox.advance_.dx_, bbox.advance_.dy_, 0., modelview, projection,
		  viewport,
		  &x, &y, &z );
    bbox.advance_.dx_ = (float)( x - x0 );
    bbox.advance_.dy_ = (float)( y - y0 );

    return bbox;
  }

#ifndef OGLFT_NO_QT
  BBox Raster::measure ( const QChar c )
  {
    BBox bbox;
    // For starters, just get the unscaled glyph bounding box
    unsigned int f;
    FT_UInt glyph_index = 0;

    for ( f = 0; f < faces_.size(); f++ ) {
      glyph_index = FT_Get_Char_Index( faces_[f].face_, c.unicode() );
      if ( glyph_index != 0 ) break;
    }

    if ( glyph_index == 0 )
      return bbox;

    FT_Error error = FT_Load_Glyph( faces_[f].face_, glyph_index,
				    FT_LOAD_DEFAULT );
    if ( error != 0 )
      return bbox;

    FT_Glyph glyph;
    error = FT_Get_Glyph( faces_[f].face_->glyph, &glyph );
    if ( error != 0 )
      return bbox;

    FT_BBox ft_bbox;
    FT_Glyph_Get_CBox( glyph, ft_glyph_bbox_unscaled, &ft_bbox );

    FT_Done_Glyph( glyph );

    bbox = ft_bbox;
    bbox.advance_ = faces_[f].face_->glyph->advance;

    // In order to be accurate regarding the placement of text not
    // aligned at the glyph's origin (CENTER/MIDDLE), the bounding box
    // of the raster format has to be projected back into the
    // view's coordinates

    GLint viewport[4];
    GLdouble modelview[16], projection[16];

    glGetIntegerv( GL_VIEWPORT, viewport );
    glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
    glGetDoublev( GL_PROJECTION_MATRIX, projection );

    // Well, first we have to get the Origin, since that is the basis
    // of the bounding box
    GLdouble x0, y0, z0;
    gluUnProject( 0., 0., 0., modelview, projection, viewport, &x0, &y0, &z0 );

    GLdouble x, y, z;
    gluUnProject( bbox.x_min_, bbox.y_min_, 0., modelview, projection, viewport,
		  &x, &y, &z );
    bbox.x_min_ = x - x0;
    bbox.y_min_ = y - y0;

    gluUnProject( bbox.x_max_, bbox.y_max_, 0., modelview, projection, viewport,
		  &x, &y, &z );
    bbox.x_max_ = x - x0;
    bbox.y_max_ = y - y0;

    gluUnProject( bbox.advance_.dx_, bbox.advance_.dy_, 0., modelview, projection,
		  viewport,
		  &x, &y, &z );
    bbox.advance_.dx_ = x - x0;
    bbox.advance_.dy_ = y - y0;

    return bbox;
  }
  
  BBox Raster::measure ( const QString& format, double number )
  {
    return Face::measure( format, number );
  } 
#endif /* OGLFT_NO_QT */

  GLuint Raster::compileGlyph ( FT_Face face, FT_UInt glyph_index )
  {
    GLuint dlist = glGenLists( 1 );
    glNewList( dlist, GL_COMPILE );

    renderGlyph( face, glyph_index );

    glEndList( );

    return dlist;
  }

  void Raster::setCharSize ( void )
  {
    FT_Error error;
    for ( unsigned int i = 0; i < faces_.size(); i++ ) {
      error = FT_Set_Char_Size( faces_[i].face_,
				(FT_F26Dot6)( point_size_ * 64 ),
				(FT_F26Dot6)( point_size_ * 64 ),
				resolution_,
				resolution_ );
      if ( error != 0 ) return;
    }

    if ( rotation_reference_glyph_ != 0 )
      setRotationOffset();
  }

  void Raster::setRotationOffset ( void )
  {
    FT_Error error = FT_Load_Glyph( rotation_reference_face_,
				    rotation_reference_glyph_,
				    FT_LOAD_RENDER );

    if ( error != 0 )
      return;

    rotation_offset_y_ = rotation_reference_face_->glyph->bitmap.rows / 2.f;
  }

  void Raster::clearCaches ( void )
  {
    GDLI fgi = glyph_dlists_.begin();

    for ( ; fgi != glyph_dlists_.end(); ++fgi ) {
      glDeleteLists( fgi->second, 1 );
    }

    glyph_dlists_.clear();
  }

  Monochrome::Monochrome ( const char* filename, float point_size,
			   FT_UInt resolution )
    : Raster( filename, point_size, resolution )
  {}

  Monochrome::Monochrome ( const FT_Byte* data_base, const FT_Long data_size,
			   float point_size, FT_UInt resolution )
    : Raster( data_base, data_size, point_size, resolution )
  {}

  Monochrome::Monochrome ( FT_Face face, float point_size, FT_UInt resolution )
    : Raster( face, point_size, resolution )
  {}

  Monochrome::~Monochrome ( void )
  {}

  GLubyte* Monochrome::invertBitmap ( const FT_Bitmap& bitmap )
  {
    // In FreeType 2.0.9, the pitch of bitmaps was rounded up to an
    // even number. In general, this disagrees with what we had been
    // using for OpenGL.

    int width = bitmap.width / 8 + ( ( bitmap.width & 7 ) > 0 ? 1 : 0 );

    GLubyte* inverse = new GLubyte[ bitmap.rows * width ];
    GLubyte* inverse_ptr = inverse;

    for ( int r = 0; r < bitmap.rows; r++ ) {

      GLubyte* bitmap_ptr = &bitmap.buffer[bitmap.pitch * ( bitmap.rows - r - 1 )];

      for ( int p = 0; p < width; p++ )
	*inverse_ptr++ = *bitmap_ptr++;
    }

    return inverse;
  }

  void Monochrome::renderGlyph ( FT_Face face, FT_UInt glyph_index )
  {
    // Start by retrieving the glyph's data.

    FT_Error error = FT_Load_Glyph( face, glyph_index, FT_LOAD_DEFAULT );

    if ( error != 0 )
      return;

    FT_Glyph original_glyph;
    FT_Glyph glyph;

    error = FT_Get_Glyph( face->glyph, &original_glyph );

    if ( error != 0 )
      return;

    error = FT_Glyph_Copy( original_glyph, &glyph );

    FT_Done_Glyph( original_glyph );

    if ( error != 0 )
      return;

    // If the individual characters are rotated (as distinct from string
    // rotation), then apply that extra rotation here. This is equivalent
    // to the sequence
    // glTranslate(x_center,y_center);
    // glRotate(angle);
    // glTranslate(-x_center,-y_center);
    // which is used for the polygonal styles. The deal with the raster
    // styles is that you must retain the advance from the string rotation
    // so that the glyphs are laid out properly. So, we make a copy of
    // the string rotated glyph, and then rotate that and add back an
    // additional offset to (in effect) restore the proper origin and
    // advance of the glyph.

    if ( character_rotation_z_ != 0. ) {
      FT_Matrix rotation_matrix;
      FT_Vector sinus;

      FT_Vector_Unit( &sinus, (FT_Angle)(character_rotation_z_ * 0x10000L) );

      rotation_matrix.xx = sinus.x;
      rotation_matrix.xy = -sinus.y;
      rotation_matrix.yx = sinus.y;
      rotation_matrix.yy = sinus.x;

      FT_Vector original_offset, rotation_offset;

      original_offset.x = ( face->glyph->metrics.width / 2
	+ face->glyph->metrics.horiBearingX ) / 64 * 0x10000L;
      original_offset.y = (FT_Pos)(rotation_offset_y_ * 0x10000L);

      rotation_offset = original_offset;

      FT_Vector_Rotate( &rotation_offset,
			(FT_Angle)(character_rotation_z_ * 0x10000L) );

      rotation_offset.x = original_offset.x - rotation_offset.x;
      rotation_offset.y = original_offset.y - rotation_offset.y;

      rotation_offset.x /= 1024;
      rotation_offset.y /= 1024;

      error = FT_Glyph_Transform( glyph, &rotation_matrix, &rotation_offset );
    }

    error = FT_Glyph_To_Bitmap( &glyph, ft_render_mode_mono, 0, 1 );

    if ( error != 0 ) {
      FT_Done_Glyph( glyph );
      return;
    }

    FT_BitmapGlyph bitmap_glyph = (FT_BitmapGlyph)glyph;

    // Evidently, in FreeType2, you can only get "upside-down" bitmaps and
    // OpenGL won't invert a bitmap with PixelZoom, so we have to invert the
    // glyph's bitmap ourselves.

    GLubyte* inverted_bitmap = invertBitmap( bitmap_glyph->bitmap );

    glBitmap( bitmap_glyph->bitmap.width, bitmap_glyph->bitmap.rows,
	      (GLfloat)-bitmap_glyph->left,
	      (GLfloat)( bitmap_glyph->bitmap.rows - bitmap_glyph->top ),
	      face->glyph->advance.x / 64.f,
	      face->glyph->advance.y / 64.f,
	      inverted_bitmap );

    FT_Done_Glyph( glyph );

    delete[] inverted_bitmap;
  }

  Grayscale::Grayscale ( const char* filename, float point_size,
			 FT_UInt resolution )
    : Raster( filename, point_size, resolution )
  {}

  Grayscale::Grayscale ( const FT_Byte* data_base, const FT_Long data_size,
			 float point_size, FT_UInt resolution )
    : Raster( data_base, data_size, point_size, resolution )
  {}

  Grayscale::Grayscale ( FT_Face face, float point_size, FT_UInt resolution )
    : Raster( face, point_size, resolution )
  {}

  Grayscale::~Grayscale ( void )
  {}

  GLubyte* Grayscale::invertPixmap ( const FT_Bitmap& bitmap )
  {
    GLubyte* inverse = new GLubyte[ bitmap.rows * bitmap.pitch ];
    GLubyte* inverse_ptr = inverse;

    for ( int r = 0; r < bitmap.rows; r++ ) {

      GLubyte* bitmap_ptr = &bitmap.buffer[bitmap.pitch * ( bitmap.rows - r - 1 )];

      for ( int p = 0; p < bitmap.pitch; p++ ) {
	*inverse_ptr++ = *bitmap_ptr++;
      }
    }

    return inverse;
  }

  void Grayscale::renderGlyph ( FT_Face face, FT_UInt glyph_index )
  {
    FT_Error error = FT_Load_Glyph( face, glyph_index, FT_LOAD_DEFAULT );

    if ( error != 0 )
      return;

    FT_Glyph original_glyph;
    FT_Glyph glyph;

    error = FT_Get_Glyph( face->glyph, &original_glyph );

    if ( error != 0 ) return;

    error = FT_Glyph_Copy( original_glyph, &glyph );

    FT_Done_Glyph( original_glyph );

    if ( error != 0 ) return;

    if ( character_rotation_z_ != 0. ) {
      FT_Matrix rotation_matrix;
      FT_Vector sinus;

      FT_Vector_Unit( &sinus, (FT_Angle)(character_rotation_z_ * 0x10000L) );

      rotation_matrix.xx = sinus.x;
      rotation_matrix.xy = -sinus.y;
      rotation_matrix.yx = sinus.y;
      rotation_matrix.yy = sinus.x;

      FT_Vector original_offset, rotation_offset;

      original_offset.x = ( face->glyph->metrics.width / 2
	+ face->glyph->metrics.horiBearingX ) / 64 * 0x10000L;
      original_offset.y = (FT_Pos)(rotation_offset_y_ * 0x10000L);

      rotation_offset = original_offset;

      FT_Vector_Rotate( &rotation_offset,
			(FT_Angle)(character_rotation_z_ * 0x10000L) );

      rotation_offset.x = original_offset.x - rotation_offset.x;
      rotation_offset.y = original_offset.y - rotation_offset.y;

      rotation_offset.x /= 1024;
      rotation_offset.y /= 1024;

      error = FT_Glyph_Transform( glyph, &rotation_matrix, &rotation_offset );
    }

    error = FT_Glyph_To_Bitmap( &glyph, ft_render_mode_normal, 0, 1 );

    if ( error != 0 ) {
      FT_Done_Glyph( glyph );
      return;
    }

    FT_BitmapGlyph bitmap_glyph = (FT_BitmapGlyph)glyph;

    // Evidently, in FreeType2, you can only get "upside-down" bitmaps
    // (this could be cured with PixelZoom, but that an additional function)

    GLubyte* inverted_pixmap = invertPixmap( bitmap_glyph->bitmap );

    // :-( If this is compiled in a display list, it may or not be in effect
    // later when the list is actually called. So, the client should be alerted
    // to this fact: unpack alignment must be 1

    glPushAttrib( GL_PIXEL_MODE_BIT );
    glPixelTransferf( GL_RED_SCALE, foreground_color_[R] - background_color_[R] );
    glPixelTransferf( GL_GREEN_SCALE, foreground_color_[G] - background_color_[G] );
    glPixelTransferf( GL_BLUE_SCALE, foreground_color_[B] - background_color_[B] );
    glPixelTransferf( GL_ALPHA_SCALE, foreground_color_[A] );
    glPixelTransferf( GL_RED_BIAS, background_color_[R] );
    glPixelTransferf( GL_GREEN_BIAS, background_color_[G] );
    glPixelTransferf( GL_BLUE_BIAS, background_color_[B] );
    glPixelTransferf( GL_ALPHA_BIAS, background_color_[A] );

    glBitmap( 0, 0, 0, 0,
	      (GLfloat)bitmap_glyph->left,
	      (GLfloat)( bitmap_glyph->top - bitmap_glyph->bitmap.rows ),
	      0 );

    glDrawPixels( bitmap_glyph->bitmap.width, bitmap_glyph->bitmap.rows,
		  GL_LUMINANCE, GL_UNSIGNED_BYTE,
		  inverted_pixmap );

    // This is how you advance the raster position when drawing PIXMAPS
    // (without querying the state)

    glBitmap( 0, 0, 0, 0,
	      (GLfloat)( -bitmap_glyph->left + face->glyph->advance.x / 64.f ),
	      (GLfloat)( bitmap_glyph->bitmap.rows - bitmap_glyph->top +
	      face->glyph->advance.y / 64. ),
	      0 );

    FT_Done_Glyph( glyph );

    glPopAttrib();

    delete[] inverted_pixmap;
  }

  Translucent::Translucent ( const char* filename, float point_size,
			     FT_UInt resolution )
    : Raster( filename, point_size, resolution )
  {}

  Translucent::Translucent ( const FT_Byte* data_base, const FT_Long data_size,
			     float point_size, FT_UInt resolution )
    : Raster( data_base, data_size, point_size, resolution )
  {}

  Translucent::Translucent ( FT_Face face, float point_size, FT_UInt resolution )
    : Raster( face, point_size, resolution )
  {}

  Translucent::~Translucent ( void )
  {}

  // The simplest format which glDrawPixels can render with (varying) transparency
  // is GL_LUMINANCE_ALPHA; so, we take the grayscale bitmap from FreeType
  // and treat all non-zero values as full luminance (basically the mask for
  // rendering) and duplicate the grayscale values as alpha values
  // (as well as turn it upside-down).

  GLubyte* Translucent::invertPixmapWithAlpha ( const FT_Bitmap& bitmap )
  {
    GLubyte* inverse = new GLubyte[ 2 * bitmap.rows * bitmap.pitch ];
    GLubyte* inverse_ptr = inverse;

    for ( int r = 0; r < bitmap.rows; r++ ) {

      GLubyte* bitmap_ptr = &bitmap.buffer[bitmap.pitch * ( bitmap.rows - r - 1 )];

      for ( int p = 0; p < bitmap.pitch; p++ ) {
	*inverse_ptr++ = *bitmap_ptr ? 255 : 0;
	*inverse_ptr++ = *bitmap_ptr++;
      }
    }

    return inverse;
  }

  void Translucent::renderGlyph ( FT_Face face, FT_UInt glyph_index )
  {
    FT_Error error = FT_Load_Glyph( face, glyph_index, FT_LOAD_DEFAULT );

    if ( error != 0 )
      return;

    FT_Glyph original_glyph;
    FT_Glyph glyph;

    error = FT_Get_Glyph( face->glyph, &original_glyph );

    if ( error != 0 ) return;

    error = FT_Glyph_Copy( original_glyph, &glyph );

    FT_Done_Glyph( original_glyph );

    if ( error != 0 ) return;

    if ( character_rotation_z_ != 0. ) {
      FT_Matrix rotation_matrix;
      FT_Vector sinus;

      FT_Vector_Unit( &sinus, (FT_Angle)(character_rotation_z_ * 0x10000L) );

      rotation_matrix.xx = sinus.x;
      rotation_matrix.xy = -sinus.y;
      rotation_matrix.yx = sinus.y;
      rotation_matrix.yy = sinus.x;

      FT_Vector original_offset, rotation_offset;

      original_offset.x = ( face->glyph->metrics.width / 2
	+ face->glyph->metrics.horiBearingX ) / 64 * 0x10000L;
      original_offset.y = (FT_Pos)(rotation_offset_y_ * 0x10000L);

      rotation_offset = original_offset;

      FT_Vector_Rotate( &rotation_offset,
			(FT_Angle)(character_rotation_z_ * 0x10000L) );

      rotation_offset.x = original_offset.x - rotation_offset.x;
      rotation_offset.y = original_offset.y - rotation_offset.y;

      rotation_offset.x /= 1024;
      rotation_offset.y /= 1024;

      error = FT_Glyph_Transform( glyph, &rotation_matrix, &rotation_offset );
    }

    error = FT_Glyph_To_Bitmap( &glyph, ft_render_mode_normal, 0, 1 );

    if ( error != 0 ) {
      FT_Done_Glyph( glyph );
      return;
    }

    FT_BitmapGlyph bitmap_glyph = (FT_BitmapGlyph)glyph;

    // Evidently, in FreeType2, you can only get "upside-down" bitmaps. For
    // translucency, the grayscale bitmap generated by FreeType is expanded
    // to include an alpha value (and the non-zero values of the
    // grayscale bitmap are saturated to provide a "mask" of the glyph).

    GLubyte* inverted_pixmap = invertPixmapWithAlpha( bitmap_glyph->bitmap );

    glPushAttrib( GL_PIXEL_MODE_BIT );
    glPixelTransferf( GL_RED_SCALE, foreground_color_[R] - background_color_[R] );
    glPixelTransferf( GL_GREEN_SCALE, foreground_color_[G] -background_color_[G] );
    glPixelTransferf( GL_BLUE_SCALE, foreground_color_[B] - background_color_[B] );
    glPixelTransferf( GL_ALPHA_SCALE, foreground_color_[A] );
    glPixelTransferf( GL_RED_BIAS, background_color_[R] );
    glPixelTransferf( GL_GREEN_BIAS, background_color_[G] );
    glPixelTransferf( GL_BLUE_BIAS, background_color_[B] );
    glPixelTransferf( GL_ALPHA_BIAS, background_color_[A] );

    // Set the proper raster position for rendering this glyph (why doesn't
    // OpenGL have a similar function for pixmaps?)

    glBitmap( 0, 0, 0, 0,
	      (GLfloat)bitmap_glyph->left,
	      (GLfloat)( bitmap_glyph->top - bitmap_glyph->bitmap.rows ),
	      0 );

    glDrawPixels( bitmap_glyph->bitmap.width, bitmap_glyph->bitmap.rows,
		  GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE,
		  inverted_pixmap );

    // This is how you advance the raster position when drawing PIXMAPS
    // (without querying the state)

    glBitmap( 0, 0, 0, 0,
	      -bitmap_glyph->left + face->glyph->advance.x / 64.f,
	      bitmap_glyph->bitmap.rows - bitmap_glyph->top +
	      face->glyph->advance.y / 64.f,
	      0 );

    FT_Done_Glyph( glyph );

    glPopAttrib();

    delete[] inverted_pixmap;
  }

  Polygonal::Polygonal ( const char* filename, float point_size, FT_UInt resolution )
    : Face( filename, point_size, resolution )
  {
    if ( !isValid() ) return;

    init();
  }

  Polygonal::Polygonal ( const FT_Byte* data_base, const FT_Long data_size, float point_size, FT_UInt resolution)
    : Face( data_base, data_size, point_size, resolution )
  {
    if ( !isValid() ) return;

    init();
  }

  Polygonal::Polygonal ( FT_Face face, float point_size, FT_UInt resolution )
    : Face( face, point_size, resolution )
  {
    init();
  }

  void Polygonal::init ( void )
  {
    character_rotation_.active_ = false;
    character_rotation_.x_ = 0;
    character_rotation_.y_ = 0;
    character_rotation_.z_ = 0;

    tessellation_steps_ =  DEFAULT_TESSELLATION_STEPS;

    delta_ = 1. / (double)tessellation_steps_;
    delta2_ = delta_ * delta_;
    delta3_ = delta2_ * delta_;

    // For vector rendition modes, FreeType is allowed to generate the
    // lines and arcs at the original face definition resolution. To
    // get to the proper glyph size, the vertices are scaled before
    // they're passed to the GLU tessellation routines.

    if ( resolution_ != 0 )
      vector_scale_ = ( point_size_ * resolution_ ) /
	(float)( faces_.front().face_->units_per_EM * 72 );
    else // According to the FreeType documentation, resolution == 0 -> 72 DPI
      vector_scale_ = ( point_size_ ) /
	(float)( faces_.front().face_->units_per_EM );

    color_tess_ = 0;
    texture_tess_ = 0;

    setCharSize();

    // Can't call this until a valid character size is set!

    setCharacterRotationReference( 'o' );
  }

  Polygonal::~Polygonal ( void )
  {
    clearCaches();
  }

  // Note: Changing the color tessellation object also clears the
  // display list cache

  void Polygonal::setColorTess ( ColorTess* color_tess )
  {
    color_tess_ = color_tess;

    clearCaches();
  }

  // Note: Changing the texture coordinate tessellation object also
  // clears the display list cache

  void Polygonal::setTextureTess ( TextureTess* texture_tess )
  {
    texture_tess_ = texture_tess;

    clearCaches();
  }

  // Note: Changing the appoximation steps also clears the display list cache

  void Polygonal::setTessellationSteps ( unsigned int tessellation_steps )
  {
    if ( tessellation_steps != tessellation_steps_ ) {

      tessellation_steps_ = tessellation_steps;

      delta_ = 1. / (double)tessellation_steps_;
      delta2_ = delta_ * delta_;
      delta3_ = delta2_ * delta_;

      clearCaches();
    }
  }

  // Note: Changing the character rotation also clears the display list cache.

  void Polygonal::setCharacterRotationX ( GLfloat character_rotation_x )
  {
    if ( character_rotation_x != character_rotation_.x_ ) {
      character_rotation_.x_ = character_rotation_x;

      if ( character_rotation_.x_ != 0. || character_rotation_.y_ != 0. ||
	   character_rotation_.z_ != 0. )
	character_rotation_.active_ = true;
      else
	character_rotation_.active_ = false;

      clearCaches();
    }
  }

  void Polygonal::setCharacterRotationY ( GLfloat character_rotation_y )
  {
    if ( character_rotation_y != character_rotation_.y_ ) {
      character_rotation_.y_ = character_rotation_y;

      if ( character_rotation_.x_ != 0. || character_rotation_.y_ != 0. ||
	   character_rotation_.z_ != 0. )
	character_rotation_.active_ = true;
      else
	character_rotation_.active_ = false;

      clearCaches();
    }
  }

  void Polygonal::setCharacterRotationZ ( GLfloat character_rotation_z )
  {
    if ( character_rotation_z != character_rotation_.z_ ) {
      character_rotation_.z_ = character_rotation_z;

      if ( character_rotation_.x_ != 0. || character_rotation_.y_ != 0. ||
	   character_rotation_.z_ != 0. )
	character_rotation_.active_ = true;
      else
	character_rotation_.active_ = false;

      clearCaches();
    }
  }

  void Polygonal::setCharSize ( void )
  {
    for ( unsigned int i = 0; i < faces_.size(); i++ ) {
      FT_Error error = FT_Set_Char_Size( faces_[i].face_,
					 0,
					 faces_[i].face_->units_per_EM * 64,
					 0,
					 0 );
      if ( error != 0 ) return;
    }

    if ( rotation_reference_glyph_ != 0 )
      setRotationOffset();
  }

  void Polygonal::setRotationOffset ( void )
  {
    FT_Error error = FT_Load_Glyph( rotation_reference_face_,
				    rotation_reference_glyph_,
				    FT_LOAD_RENDER );

    if ( error != 0 )
      return;

    vector_scale_ = ( point_size_ * resolution_ ) /
      ( 72.f * rotation_reference_face_->units_per_EM );

    rotation_offset_y_ =
      ( rotation_reference_face_->glyph->metrics.horiBearingY / 2.f ) / 64.f
      * vector_scale_;
  }

  double Polygonal::height ( void ) const
  {
    if ( faces_[0].face_->height > 0 )
      return ( faces_[0].face_->height * point_size_ * resolution_ ) /
	( 72. * faces_[0].face_->units_per_EM );
    else
      return ( faces_[0].face_->size->metrics.y_ppem * point_size_ * resolution_ ) /
	( 72. * faces_[0].face_->units_per_EM );
  }

  BBox Polygonal::measure ( unsigned char c )
  {
    BBox bbox;
    // For starters, just get the unscaled glyph bounding box
    unsigned int f;
    FT_UInt glyph_index = 0;

    for ( f = 0; f < faces_.size(); f++ ) {
      glyph_index = FT_Get_Char_Index( faces_[f].face_, c );
      if ( glyph_index != 0 ) break;
    }

    if ( glyph_index == 0 )
      return bbox;

    FT_Error error = FT_Load_Glyph( faces_[f].face_, glyph_index,
				    FT_LOAD_DEFAULT );
    if ( error != 0 )
      return bbox;

    FT_Glyph glyph;
    error = FT_Get_Glyph( faces_[f].face_->glyph, &glyph );
    if ( error != 0 )
      return bbox;

    FT_BBox ft_bbox;
    FT_Glyph_Get_CBox( glyph, ft_glyph_bbox_unscaled, &ft_bbox );

    FT_Done_Glyph( glyph );

    bbox = ft_bbox;
    bbox.advance_ = faces_[f].face_->glyph->advance;

    bbox *=
      ( point_size_ * resolution_ ) / ( 72.f * faces_[f].face_->units_per_EM );

    return bbox;
  }
#ifndef OGLFT_NO_QT
  BBox Polygonal::measure ( const QChar c )
  {
    BBox bbox;
    // For starters, just get the unscaled glyph bounding box
    unsigned int f;
    FT_UInt glyph_index = 0;

    for ( f = 0; f < faces_.size(); f++ ) {
      glyph_index = FT_Get_Char_Index( faces_[f].face_, c.unicode() );
      if ( glyph_index != 0 ) break;
    }

    if ( glyph_index == 0 )
      return bbox;

    FT_Error error = FT_Load_Glyph( faces_[f].face_, glyph_index,
				    FT_LOAD_DEFAULT );
    if ( error != 0 )
      return bbox;

    FT_Glyph glyph;
    error = FT_Get_Glyph( faces_[f].face_->glyph, &glyph );
    if ( error != 0 )
      return bbox;

    FT_BBox ft_bbox;
    FT_Glyph_Get_CBox( glyph, ft_glyph_bbox_unscaled, &ft_bbox );

    FT_Done_Glyph( glyph );

    bbox = ft_bbox;
    bbox.advance_ = faces_[f].face_->glyph->advance;

    bbox *= ( point_size_ * resolution_ ) / ( 72. * faces_[f].face_->units_per_EM );

    return bbox;
  }
#endif /* OGLFT_NO_QT */

  GLuint Polygonal::compileGlyph ( FT_Face face, FT_UInt glyph_index )
  {
    GLuint dlist = glGenLists( 1 );

    glNewList( dlist, GL_COMPILE );

    renderGlyph( face, glyph_index );

    glEndList( );

    return dlist;
  }

  void Polygonal::clearCaches ( void )
  {
    GDLI fgi = glyph_dlists_.begin();

    for ( ; fgi != glyph_dlists_.end(); ++fgi ) {
      glDeleteLists( fgi->second, 1 );
    }

    glyph_dlists_.clear();
  }

  Outline::Outline ( const char* filename, float point_size, FT_UInt resolution )
    : Polygonal( filename, point_size, resolution )
  {
    if ( !isValid() ) return;

    init();
  }

  Outline::Outline ( const FT_Byte* data_base, const FT_Long data_size, float point_size, FT_UInt resolution)
    : Polygonal( data_base, data_size, point_size, resolution )
  {
    if ( !isValid() ) return;

    init();
  }

  Outline::Outline ( FT_Face face, float point_size, FT_UInt resolution )
    : Polygonal( face, point_size, resolution )
  {
    init();
  }
  
  void Outline::init ( void )
  {
    interface_.move_to = (FT_Outline_MoveTo_Func)moveToCallback;
    interface_.line_to = (FT_Outline_LineTo_Func)lineToCallback;
    interface_.conic_to = (FT_Outline_ConicTo_Func)conicToCallback;
    interface_.cubic_to = (FT_Outline_CubicTo_Func)cubicToCallback;
    interface_.shift = 0;
    interface_.delta = 0;
  }

  Outline::~Outline ( void )
  {}

  void Outline::renderGlyph ( FT_Face face, FT_UInt glyph_index )
  {
    FT_Error error = FT_Load_Glyph( face, glyph_index, FT_LOAD_DEFAULT );

    if ( error != 0 )
      return;

    FT_Glyph g;

    error = FT_Get_Glyph( face->glyph, &g );

    if ( error != 0 || g->format != FT_GLYPH_FORMAT_OUTLINE )
      return;

    vector_scale_ = ( point_size_ * resolution_ ) /
      ( 72.f * face->units_per_EM );

    if ( character_rotation_.active_ ) {
      glPushMatrix();
      glTranslatef( ( face->glyph->metrics.width / 2.f +
		      face->glyph->metrics.horiBearingX ) / 64.f
		    * vector_scale_,
		    rotation_offset_y_,
		    0.f );

      if ( character_rotation_.x_ != 0. )
	glRotatef( character_rotation_.x_, 1., 0., 0. );

      if ( character_rotation_.y_ != 0. )
	glRotatef( character_rotation_.y_, 0., 1., 0. );

      if ( character_rotation_.z_ != 0. )
	glRotatef( character_rotation_.z_, 0., 0., 1. );

      glTranslatef( -( face->glyph->metrics.width / 2.f +
		      face->glyph->metrics.horiBearingX ) / 64.f
		    * vector_scale_,
		    -rotation_offset_y_,
		    0.f );
    }

    contour_open_ = false;

    // The Big Kahuna: the FreeType glyph decomposition routine traverses
    // the outlines of the font by calling the various routines stored in
    // outline_interface_. These routines in turn call the GL vertex routines.

    error = FT_Outline_Decompose( &((FT_OutlineGlyph)g)->outline,
				  &interface_, this );

    FT_Done_Glyph( (FT_Glyph)g );

    // Some glyphs may be empty (the 'blank' for instance!)

    if ( contour_open_ )
      glEnd( );


    if ( character_rotation_.active_ ) {
      glPopMatrix();
    }

    // Drawing a character always advances the MODELVIEW.

    glTranslatef( face->glyph->advance.x / 64.f * vector_scale_,
		  face->glyph->advance.y / 64.f * vector_scale_,
		  0.f );

    for ( VILI vili = vertices_.begin(); vili != vertices_.end(); vili++ )
      delete *vili;

    vertices_.clear();
  }

  int Outline::moveToCallback ( FT_Vector* to, Outline* outline )
  {
    if ( outline->contour_open_ ) {
      glEnd();
    }

    outline->last_vertex_ = VertexInfo( to,
					outline->colorTess(),
					outline->textureTess() );

    glBegin( GL_LINE_LOOP );

    outline->contour_open_ = true;

    return 0;
  }

  int Outline::lineToCallback ( FT_Vector* to, Outline* outline )
  {
    outline->last_vertex_ = VertexInfo( to,
					outline->colorTess(),
					outline->textureTess() );
    GLdouble g[2];

    g[X] = outline->last_vertex_.v_[X] * outline->vector_scale_;
    g[Y] = outline->last_vertex_.v_[Y] * outline->vector_scale_;

    glVertex2dv( g );

    return 0;
  }

  int Outline::conicToCallback ( FT_Vector* control, FT_Vector* to, Outline* outline )
  {
    // This is crude: Step off conics with a fixed number of increments

    VertexInfo to_vertex( to, outline->colorTess(), outline->textureTess() );
    VertexInfo control_vertex( control, outline->colorTess(), outline->textureTess() );

    double b[2], c[2], d[2], f[2], df[2], d2f[2];
    GLdouble g[3];

    g[Z] = 0.;

    b[X] = outline->last_vertex_.v_[X] - 2 * control_vertex.v_[X] +
      to_vertex.v_[X];
    b[Y] = outline->last_vertex_.v_[Y] - 2 * control_vertex.v_[Y] +
      to_vertex.v_[Y];

    c[X] = -2 * outline->last_vertex_.v_[X] + 2 * control_vertex.v_[X];
    c[Y] = -2 * outline->last_vertex_.v_[Y] + 2 * control_vertex.v_[Y];

    d[X] = outline->last_vertex_.v_[X];
    d[Y] = outline->last_vertex_.v_[Y];

    f[X] = d[X];
    f[Y] = d[Y];
    df[X] = c[X] * outline->delta_ + b[X] * outline->delta2_;
    df[Y] = c[Y] * outline->delta_ + b[Y] * outline->delta2_;
    d2f[X] = 2 * b[X] * outline->delta2_;
    d2f[Y] = 2 * b[Y] * outline->delta2_;

    for ( unsigned int i = 0; i < outline->tessellation_steps_-1; i++ ) {

      f[X] += df[X];
      f[Y] += df[Y];

      g[X] = f[X] * outline->vector_scale_;
      g[Y] = f[Y] * outline->vector_scale_;

      if ( outline->colorTess() )
	glColor4fv( outline->colorTess()->color( g ) );

      glVertex2dv( g );

      df[X] += d2f[X];
      df[Y] += d2f[Y];
    }

    g[X] = to_vertex.v_[X] * outline->vector_scale_;
    g[Y] = to_vertex.v_[Y] * outline->vector_scale_;

    if ( outline->colorTess() )
      glColor4fv( outline->colorTess()->color( g ) );

    glVertex2dv( g );

    outline->last_vertex_ = to_vertex;

    return 0;
  }

  int Outline::cubicToCallback ( FT_Vector* control1, FT_Vector* control2,
				 FT_Vector* to, Outline* outline )
  {
    // This is crude: Step off cubics with a fixed number of increments

    VertexInfo to_vertex( to, outline->colorTess(), outline->textureTess() );
    VertexInfo control1_vertex( control1, outline->colorTess(), outline->textureTess() );
    VertexInfo control2_vertex( control2, outline->colorTess(), outline->textureTess() );

    double a[2], b[2], c[2], d[2], f[2], df[2], d2f[2], d3f[2];
    GLdouble g[3];

    g[Z] = 0.;

    a[X] = -outline->last_vertex_.v_[X] + 3 * control1_vertex.v_[X]
      -3 * control2_vertex.v_[X] + to_vertex.v_[X];
    a[Y] = -outline->last_vertex_.v_[Y] + 3 * control1_vertex.v_[Y]
      -3 * control2_vertex.v_[Y] + to_vertex.v_[Y];

    b[X] = 3 * outline->last_vertex_.v_[X] - 6 * control1_vertex.v_[X] +
      3 * control2_vertex.v_[X];
    b[Y] = 3 * outline->last_vertex_.v_[Y] - 6 * control1_vertex.v_[Y] +
      3 * control2_vertex.v_[Y];

    c[X] = -3 * outline->last_vertex_.v_[X] + 3 * control1_vertex.v_[X];
    c[Y] = -3 * outline->last_vertex_.v_[Y] + 3 * control1_vertex.v_[Y];

    d[X] = outline->last_vertex_.v_[X];
    d[Y] = outline->last_vertex_.v_[Y];

    f[X] = d[X];
    f[Y] = d[Y];
    df[X] = c[X] * outline->delta_ + b[X] * outline->delta2_
      + a[X] * outline->delta3_;
    df[Y] = c[Y] * outline->delta_ + b[Y] * outline->delta2_
      + a[Y] * outline->delta3_;
    d2f[X] = 2 * b[X] * outline->delta2_ + 6 * a[X] * outline->delta3_;
    d2f[Y] = 2 * b[Y] * outline->delta2_ + 6 * a[Y] * outline->delta3_;
    d3f[X] = 6 * a[X] * outline->delta3_;
    d3f[Y] = 6 * a[Y] * outline->delta3_;

    for ( unsigned int i = 0; i < outline->tessellation_steps_-1; i++ ) {

      f[X] += df[X];
      f[Y] += df[Y];

      g[X] = f[X] * outline->vector_scale_;
      g[Y] = f[Y] * outline->vector_scale_;

      if ( outline->colorTess() )
	glColor4fv( outline->colorTess()->color( g ) );

      glVertex2dv( g );

      df[X] += d2f[X];
      df[Y] += d2f[Y];
      d2f[X] += d3f[X];
      d2f[Y] += d3f[Y];
    }

    g[X] = to_vertex.v_[X] * outline->vector_scale_;
    g[Y] = to_vertex.v_[Y] * outline->vector_scale_;

    if ( outline->colorTess() )
      glColor4fv( outline->colorTess()->color( g ) );

    glVertex2dv( g );

    outline->last_vertex_ = to_vertex;

    return 0;
  }

  Filled::Filled ( const char* filename, float point_size, FT_UInt resolution )
    : Polygonal( filename, point_size, resolution )
  {
    if ( !isValid() ) return;

    init();
  }

  Filled::Filled ( const FT_Byte* data_base, const FT_Long data_size,
		   float point_size, FT_UInt resolution )
    : Polygonal( data_base, data_size, point_size, resolution )
  {
    if ( !isValid() ) return;

    init();
  }

  Filled::Filled ( FT_Face face, float point_size, FT_UInt resolution )
    : Polygonal( face, point_size, resolution )
  {
    init();
  }

  void Filled::init ( void )
  {
    depth_offset_ = 0;

    interface_.move_to = (FT_Outline_MoveTo_Func)moveToCallback;
    interface_.line_to = (FT_Outline_LineTo_Func)lineToCallback;
    interface_.conic_to = (FT_Outline_ConicTo_Func)conicToCallback;
    interface_.cubic_to = (FT_Outline_CubicTo_Func)cubicToCallback;
    interface_.shift = 0;
    interface_.delta = 0;

    tess_obj_ = gluNewTess();

#if defined(WIN32)
    typedef void (CALLBACK*(CB))();
#else
    typedef GLUTessCallback CB;
#endif
    gluTessCallback( tess_obj_, GLU_TESS_VERTEX, CB(vertexCallback) );
    gluTessCallback( tess_obj_, GLU_TESS_BEGIN, CB(beginCallback) );
    gluTessCallback( tess_obj_, GLU_TESS_END, CB(endCallback) );
    gluTessCallback( tess_obj_, GLU_TESS_COMBINE_DATA, CB(combineCallback) );
    gluTessCallback( tess_obj_, GLU_TESS_ERROR, CB(errorCallback) );
  }

  Filled::~Filled ( void )
  {
    gluDeleteTess( tess_obj_ );
  }

  void Filled::renderGlyph ( FT_Face face, FT_UInt glyph_index )
  {
    FT_Error error = FT_Load_Glyph( face, glyph_index, FT_LOAD_DEFAULT );

    if ( error != 0 )
      return;

    FT_Glyph g;

    error = FT_Get_Glyph( face->glyph, &g );

    if ( error != 0 || g->format != FT_GLYPH_FORMAT_OUTLINE )
      return;

    vector_scale_ = ( point_size_ * resolution_ ) /
      ( 72.f * face->units_per_EM );

    if ( character_rotation_.active_ ) {
      glPushMatrix();
      glTranslatef( ( face->glyph->metrics.width / 2.f +
		      face->glyph->metrics.horiBearingX ) / 64.f
		    * vector_scale_,
		    rotation_offset_y_,
		    0. );

      if ( character_rotation_.x_ != 0. )
	glRotatef( character_rotation_.x_, 1., 0., 0. );

      if ( character_rotation_.y_ != 0. )
	glRotatef( character_rotation_.y_, 0., 1., 0. );

      if ( character_rotation_.z_ != 0. )
	glRotatef( character_rotation_.z_, 0., 0., 1. );

      glTranslatef( -( face->glyph->metrics.width / 2.f +
		      face->glyph->metrics.horiBearingX ) / 64.f
		    * vector_scale_,
		    -rotation_offset_y_,
		    0.f );
    }

    if ( depth_offset_ != 0. ) {
      glPushMatrix();
      glTranslatef( 0., 0., depth_offset_ );
      glNormal3f( 0., 0., 1. );
    }
    else {
      glNormal3f( 0., 0., -1. );
    }

    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

    contour_open_ = false;

    gluTessBeginPolygon( tess_obj_, this );

    // The Big Kahuna: the FreeType glyph decomposition routine traverses
    // the outlines of the font by calling the various routines stored in
    // interface_. These routines in turn call the GLU tessellation routines
    // to create OGL polygons.

    error = FT_Outline_Decompose( &((FT_OutlineGlyph)g)->outline,
				  &interface_, this );

    FT_Done_Glyph( (FT_Glyph)g );

    // Some glyphs may be empty (the 'blank' for instance!)

    if ( contour_open_ )
      gluTessEndContour( tess_obj_ );

    gluTessEndPolygon( tess_obj_ );

    if ( depth_offset_ != 0. ) {
      glPopMatrix();
    }
    if ( character_rotation_.active_ ) {
      glPopMatrix();
    }

    // Drawing a character always advances the MODELVIEW.

    glTranslatef( face->glyph->advance.x / 64 * vector_scale_,
		  face->glyph->advance.y / 64 * vector_scale_,
		  0. );

    for ( VILI vili = extra_vertices_.begin(); vili != extra_vertices_.end(); vili++ )
      delete *vili;

    extra_vertices_.clear();

    for ( VILI vili = vertices_.begin(); vili != vertices_.end(); vili++ )
      delete *vili;

    vertices_.clear();
  }

  int Filled::moveToCallback ( FT_Vector* to, Filled* filled )
  {
    if ( filled->contour_open_ ) {
      gluTessEndContour( filled->tess_obj_ );
    }

    filled->last_vertex_ = VertexInfo( to, filled->colorTess(), filled->textureTess() );

    gluTessBeginContour( filled->tess_obj_ );

    filled->contour_open_ = true;

    return 0;
  }

  int Filled::lineToCallback ( FT_Vector* to, Filled* filled )
  {
    filled->last_vertex_ = VertexInfo( to, filled->colorTess(), filled->textureTess() );

    VertexInfo* vertex = new VertexInfo( to, filled->colorTess(), filled->textureTess() );

    vertex->v_[X] *= filled->vector_scale_;
    vertex->v_[Y] *= filled->vector_scale_;

    gluTessVertex( filled->tess_obj_, vertex->v_, vertex );

    filled->vertices_.push_back( vertex );

    return 0;
  }

  int Filled::conicToCallback ( FT_Vector* control, FT_Vector* to, Filled* filled )
  {
    // This is crude: Step off conics with a fixed number of increments

    VertexInfo to_vertex( to, filled->colorTess(), filled->textureTess() );
    VertexInfo control_vertex( control, filled->colorTess(), filled->textureTess() );

    double b[2], c[2], d[2], f[2], df[2], d2f[2];

    b[X] = filled->last_vertex_.v_[X] - 2 * control_vertex.v_[X] +
      to_vertex.v_[X];
    b[Y] = filled->last_vertex_.v_[Y] - 2 * control_vertex.v_[Y] +
      to_vertex.v_[Y];

    c[X] = -2 * filled->last_vertex_.v_[X] + 2 * control_vertex.v_[X];
    c[Y] = -2 * filled->last_vertex_.v_[Y] + 2 * control_vertex.v_[Y];

    d[X] = filled->last_vertex_.v_[X];
    d[Y] = filled->last_vertex_.v_[Y];

    f[X] = d[X];
    f[Y] = d[Y];
    df[X] = c[X] * filled->delta_ + b[X] * filled->delta2_;
    df[Y] = c[Y] * filled->delta_ + b[Y] * filled->delta2_;
    d2f[X] = 2 * b[X] * filled->delta2_;
    d2f[Y] = 2 * b[Y] * filled->delta2_;

    for ( unsigned int i = 0; i < filled->tessellation_steps_-1; i++ ) {

      f[X] += df[X];
      f[Y] += df[Y];

      VertexInfo* vertex = new VertexInfo( f, filled->colorTess(), filled->textureTess() );

      vertex->v_[X] *= filled->vector_scale_;
      vertex->v_[Y] *= filled->vector_scale_;

      filled->vertices_.push_back( vertex );

      gluTessVertex( filled->tess_obj_, vertex->v_, vertex );

      df[X] += d2f[X];
      df[Y] += d2f[Y];
    }

    VertexInfo* vertex = new VertexInfo( to, filled->colorTess(), filled->textureTess() );
  
    vertex->v_[X] *= filled->vector_scale_;
    vertex->v_[Y] *= filled->vector_scale_;

    filled->vertices_.push_back( vertex );

    gluTessVertex( filled->tess_obj_, vertex->v_, vertex );

    filled->last_vertex_ = to_vertex;

    return 0;
  }

  int Filled::cubicToCallback ( FT_Vector* control1, FT_Vector* control2,
				FT_Vector* to, Filled* filled )
  {
    // This is crude: Step off cubics with a fixed number of increments

    VertexInfo to_vertex( to, filled->colorTess(), filled->textureTess() );
    VertexInfo control1_vertex( control1, filled->colorTess(), filled->textureTess() );
    VertexInfo control2_vertex( control2, filled->colorTess(), filled->textureTess() );

    double a[2], b[2], c[2], d[2], f[2], df[2], d2f[2], d3f[2];

    a[X] = -filled->last_vertex_.v_[X] + 3 * control1_vertex.v_[X]
      -3 * control2_vertex.v_[X] + to_vertex.v_[X];
    a[Y] = -filled->last_vertex_.v_[Y] + 3 * control1_vertex.v_[Y]
      -3 * control2_vertex.v_[Y] + to_vertex.v_[Y];

    b[X] = 3 * filled->last_vertex_.v_[X] - 6 * control1_vertex.v_[X] +
      3 * control2_vertex.v_[X];
    b[Y] = 3 * filled->last_vertex_.v_[Y] - 6 * control1_vertex.v_[Y] +
      3 * control2_vertex.v_[Y];

    c[X] = -3 * filled->last_vertex_.v_[X] + 3 * control1_vertex.v_[X];
    c[Y] = -3 * filled->last_vertex_.v_[Y] + 3 * control1_vertex.v_[Y];

    d[X] = filled->last_vertex_.v_[X];
    d[Y] = filled->last_vertex_.v_[Y];

    f[X] = d[X];
    f[Y] = d[Y];
    df[X] = c[X] * filled->delta_ + b[X] * filled->delta2_
      + a[X] * filled->delta3_;
    df[Y] = c[Y] * filled->delta_ + b[Y] * filled->delta2_
      + a[Y] * filled->delta3_;
    d2f[X] = 2 * b[X] * filled->delta2_ + 6 * a[X] * filled->delta3_;
    d2f[Y] = 2 * b[Y] * filled->delta2_ + 6 * a[Y] * filled->delta3_;
    d3f[X] = 6 * a[X] * filled->delta3_;
    d3f[Y] = 6 * a[Y] * filled->delta3_;

    for ( unsigned int i = 0; i < filled->tessellation_steps_-1; i++ ) {

      f[X] += df[X];
      f[Y] += df[Y];

      VertexInfo* vertex = new VertexInfo( f, filled->colorTess(), filled->textureTess() );

      vertex->v_[X] *= filled->vector_scale_;
      vertex->v_[Y] *= filled->vector_scale_;

      filled->vertices_.push_back( vertex );

      gluTessVertex( filled->tess_obj_, vertex->v_, vertex );

      df[X] += d2f[X];
      df[Y] += d2f[Y];
      d2f[X] += d3f[X];
      d2f[Y] += d3f[Y];
    }

    VertexInfo* vertex = new VertexInfo( to, filled->colorTess(), filled->textureTess() );
  
    vertex->v_[X] *= filled->vector_scale_;
    vertex->v_[Y] *= filled->vector_scale_;
  
    filled->vertices_.push_back( vertex );

    gluTessVertex( filled->tess_obj_, vertex->v_, vertex );

    filled->last_vertex_ = to_vertex;

    return 0;
  }

  void Filled::vertexCallback ( VertexInfo* vertex )
  {
    if ( vertex->color_tess_ != 0 )
      glColor4fv( vertex->color_tess_->color( vertex->v_ ) );

    if ( vertex->texture_tess_ != 0 )
      glTexCoord2fv( vertex->texture_tess_->texCoord( vertex->v_ ) );

    glVertex3dv( vertex->v_ );
  }

  void Filled::beginCallback ( GLenum which )
  {
    glBegin( which );
  }

  void Filled::endCallback ( void )
  {
    glEnd();
  }

  void Filled::combineCallback ( GLdouble coords[3], void* vertex_data[4],
				 GLfloat weight[4], void** out_data,
				 Filled* filled )
  {
    (void)vertex_data;
    (void)weight;
    //    std::cerr << "called combine" << std::endl;
    VertexInfo* vertex = new VertexInfo( coords );
    *out_data = vertex;
    filled->extraVertices().push_back( vertex );
  }

  void Filled::errorCallback ( GLenum error_code )
  {
    std::cerr << "hmm. error during tessellation?:" << gluErrorString( error_code ) << std::endl;
  }

#ifndef OGLFT_NO_SOLID
  Solid::Solid ( const char* filename, float point_size, FT_UInt resolution )
    : Filled( filename, point_size, resolution )
  {
    if ( !isValid() ) return;

    init();
  }

  Solid::Solid ( const FT_Byte* data_base, const FT_Long data_size,
		 float point_size, FT_UInt resolution )
    : Filled( data_base, data_size, point_size, resolution )
  {
    if ( !isValid() ) return;

    init();
  }

  Solid::Solid ( FT_Face face, float point_size, FT_UInt resolution )
    : Filled( face, point_size, resolution )
  {
    init();
  }

  void Solid::init ( void )
  {
    interface_.move_to = (FT_Outline_MoveTo_Func)moveToCallback;
    interface_.line_to = (FT_Outline_LineTo_Func)lineToCallback;
    interface_.conic_to = (FT_Outline_ConicTo_Func)conicToCallback;
    interface_.cubic_to = (FT_Outline_CubicTo_Func)cubicToCallback;
    interface_.shift = 0;
    interface_.delta = 0;

    // Set up for extrusion. Default depth is 1 (units of what?)
    extrusion_.depth_ = 1.;
    extrusion_.up_[X] = 0.;
    extrusion_.up_[Y] = 1.;
    extrusion_.up_[Z] = 0.;
    extrusion_.n_polyline_pts_ = N_POLYLINE_PTS;

    assign( extrusion_.point_array_[0], 0., 0., extrusion_.depth_ + 1. );
    assign( extrusion_.point_array_[1], 0., 0., extrusion_.depth_ );
    assign( extrusion_.point_array_[2], 0., 0., 0. );
    assign( extrusion_.point_array_[3], 0., 0., -1. );

    // Turn on closed contours and smooth vertices; turn off end capping

    gleSetJoinStyle( TUBE_JN_RAW | TUBE_CONTOUR_CLOSED | TUBE_NORM_EDGE );
  }

  Solid::~Solid ( void )
  {}

  // Note: as usual, setting this clears the caches

  void Solid::setDepth ( double depth )
  {
    if ( depth > 0. && depth != extrusion_.depth_ ) {
      extrusion_.depth_ = depth;

      assign( extrusion_.point_array_[0], 0., 0., extrusion_.depth_ + 1. );
      assign( extrusion_.point_array_[1], 0., 0., extrusion_.depth_ );

      clearCaches();
    }
  }

  void Solid::renderGlyph ( FT_Face face, FT_UInt glyph_index )
  {
    FT_Error error = FT_Load_Glyph( face, glyph_index, FT_LOAD_DEFAULT );

    if ( error != 0 )
      return;

    FT_OutlineGlyph g;

    error = FT_Get_Glyph( face->glyph, (FT_Glyph*)&g );

    if ( error != 0 )
      return;

    vector_scale_ = ( point_size_ * resolution_ ) / ( 72. * face->units_per_EM );

    if ( character_rotation_.active_ ) {
      glPushMatrix();

      glTranslatef( ( face->glyph->metrics.width / 2. +
		      face->glyph->metrics.horiBearingX ) / 64.
		    * vector_scale_,
		    rotation_offset_y_,
		    0. );

      if ( character_rotation_.x_ != 0. )
	glRotatef( character_rotation_.x_, 1., 0., 0. );

      if ( character_rotation_.y_ != 0. )
	glRotatef( character_rotation_.y_, 0., 1., 0. );

      if ( character_rotation_.z_ != 0. )
	glRotatef( character_rotation_.z_, 0., 0., 1. );

      glTranslatef( -( face->glyph->metrics.width / 2. +
		      face->glyph->metrics.horiBearingX ) / 64.
		    * vector_scale_,
		    -rotation_offset_y_,
		    0. );
    }

    contour_open_ = false;

    // In theory, TrueType contours are defined clockwise and Type1 contours
    // are defined counter-clockwise. Trust the flag set by FreeType to
    // indicate this since it is critical to getting the orientation of the
    // surface normals correct.
    if ( g->outline.flags & ft_outline_reverse_fill ) {
      extrusion_.normal_sign_.x_ = -1;
      extrusion_.normal_sign_.y_ = 1;
    }
    else {
      extrusion_.normal_sign_.x_ = 1;
      extrusion_.normal_sign_.y_ = -1;
    }
    // The Big Kahuna: the FreeType glyph decomposition routine traverses
    // the outlines of the font by calling the various routines stored in
    // extrude_interface_. These in turn call the gleExtrusion routine.

    error = FT_Outline_Decompose( &g->outline, &interface_, this );

    FT_Done_Glyph( (FT_Glyph)g );

    // Some glyphs may be empty (the 'blank' for instance!)

    if ( contour_open_ ) {
      extrusion_.contour_normals_.push_back( extrusion_.contour_normals_.front() );

      gleExtrusion( extrusion_.contour_.size(),
		    &extrusion_.contour_.begin()->p_,
		    &extrusion_.contour_normals_[1].p_,
		    extrusion_.up_,
		    extrusion_.n_polyline_pts_,
		    extrusion_.point_array_,
		    0 );

      extrusion_.contour_.clear();
      extrusion_.contour_normals_.clear();
    }

    if ( character_rotation_.active_ ) {
      glPopMatrix();
    }

    // Apply the front and back faces of the solid character (recall that
    // drawing a character advances the MODELVIEW, so defend against that
    // with the stack operations)

    glPushMatrix();
    depth_offset_ = 0.;
    Filled::renderGlyph( face, glyph_index );
    glPopMatrix();

    glPushMatrix();
    depth_offset_ = extrusion_.depth_;
    Filled::renderGlyph( face, glyph_index );
    glPopMatrix();

    // Drawing a character always advances the MODELVIEW.

    glTranslatef( face->glyph->advance.x / 64. * vector_scale_,
		  face->glyph->advance.y / 64. * vector_scale_,
		  0. );

    for ( VILI vili = vertices_.begin(); vili != vertices_.end(); vili++ )
      delete *vili;

    vertices_.clear();
  }

  int Solid::moveToCallback ( FT_Vector* to, Solid* solid )
  {
    if ( solid->contour_open_ ) {

      // A word of explanation: since you can't predict when the
      // contour is going to end (its end is signaled by calling this
      // routine, i.e., the contour ends when another is started
      // abruptly), only the lineTo and arcTo functions generate contour
      // points. The upshot is that the normals, which are computed for the
      // current segment, are one behind the segment described in the
      // the contour array. To make things match up at the end, the first
      // normal is copied to the end of the normal array and the extrusion
      // routine is passed the list of normals starting at the second entry.
      
      solid->extrusion_.contour_normals_.
	push_back( solid->extrusion_.contour_normals_.front() );
#if 1
      gleExtrusion( solid->extrusion_.contour_.size(),
		    &solid->extrusion_.contour_.begin()->p_,
		    &solid->extrusion_.contour_normals_[1].p_,
		    solid->extrusion_.up_,
		    solid->extrusion_.n_polyline_pts_,
		    solid->extrusion_.point_array_,
		    0 );
#endif
      solid->extrusion_.contour_.clear();
      solid->extrusion_.contour_normals_.clear();
    }

    solid->last_vertex_ = VertexInfo( to, solid->colorTess(), solid->textureTess() );

    solid->contour_open_ = true;

    return 0;
  }

  int Solid::lineToCallback ( FT_Vector* to, Solid* solid )
  {
    VertexInfo vertex( to, solid->colorTess(), solid->textureTess() );

    VertexInfo normal( solid->extrusion_.normal_sign_.y_ *
		       ( vertex.v_[Y] - solid->last_vertex_.v_[Y] ),
		       solid->extrusion_.normal_sign_.x_ *
		       ( vertex.v_[X] - solid->last_vertex_.v_[X] ) );

    solid->last_vertex_ = vertex;

    vertex.v_[X] *= solid->vector_scale_;
    vertex.v_[Y] *= solid->vector_scale_;

    normal.normalize();

    solid->extrusion_.contour_.push_back( vertex );
    solid->extrusion_.contour_normals_.push_back( normal );

    return 0;
  }

  int Solid::conicToCallback ( FT_Vector* control, FT_Vector* to, Solid* solid )
  {
    // This is crude: Step off conics with a fixed number of increments

    VertexInfo to_vertex( to, solid->colorTess(), solid->textureTess() );
    VertexInfo control_vertex( control, solid->colorTess(), solid->textureTess() );

    double b[2], c[2], d[2], f[2], df[2], d2f[2];

    b[X] = solid->last_vertex_.v_[X] - 2 * control_vertex.v_[X] +
      to_vertex.v_[X];
    b[Y] = solid->last_vertex_.v_[Y] - 2 * control_vertex.v_[Y] +
      to_vertex.v_[Y];

    c[X] = -2 * solid->last_vertex_.v_[X] + 2 * control_vertex.v_[X];
    c[Y] = -2 * solid->last_vertex_.v_[Y] + 2 * control_vertex.v_[Y];

    d[X] = solid->last_vertex_.v_[X];
    d[Y] = solid->last_vertex_.v_[Y];

    f[X] = d[X];
    f[Y] = d[Y];
    df[X] = c[X] * solid->delta_ + b[X] * solid->delta2_;
    df[Y] = c[Y] * solid->delta_ + b[Y] * solid->delta2_;
    d2f[X] = 2 * b[X] * solid->delta2_;
    d2f[Y] = 2 * b[Y] * solid->delta2_;

    for ( unsigned int i = 0; i < solid->tessellation_steps_-1; i++ ) {

      f[X] += df[X];
      f[Y] += df[Y];

      VertexInfo vertex( f, solid->colorTess(), solid->textureTess() );

      VertexInfo normal( solid->extrusion_.normal_sign_.y_ * df[Y],
			 solid->extrusion_.normal_sign_.x_ * df[X] );

      vertex.v_[X] *= solid->vector_scale_;
      vertex.v_[Y] *= solid->vector_scale_;

      normal.normalize();

      solid->extrusion_.contour_.push_back( vertex );
      solid->extrusion_.contour_normals_.push_back( normal );

      df[X] += d2f[X];
      df[Y] += d2f[Y];
    }

    VertexInfo vertex( to, solid->colorTess(), solid->textureTess() );

    VertexInfo normal( solid->extrusion_.normal_sign_.y_ * df[Y],
		       solid->extrusion_.normal_sign_.x_ * df[X] );
  
    vertex.v_[X] *= solid->vector_scale_;
    vertex.v_[Y] *= solid->vector_scale_;

    normal.normalize();

    solid->extrusion_.contour_.push_back( vertex );
    solid->extrusion_.contour_normals_.push_back( normal );

    solid->last_vertex_ = to_vertex;

    return 0;
  }

  int Solid::cubicToCallback ( FT_Vector* control1, FT_Vector* control2,
			       FT_Vector* to, Solid* solid )
  {
    // This is crude: Step off cubics with a fixed number of increments

    VertexInfo to_vertex( to, solid->colorTess(), solid->textureTess() );
    VertexInfo control1_vertex( control1, solid->colorTess(), solid->textureTess() );
    VertexInfo control2_vertex( control2, solid->colorTess(), solid->textureTess() );

    double a[2], b[2], c[2], d[2], f[2], df[2], d2f[2], d3f[2];

    a[X] = -solid->last_vertex_.v_[X] + 3 * control1_vertex.v_[X]
      -3 * control2_vertex.v_[X] + to_vertex.v_[X];
    a[Y] = -solid->last_vertex_.v_[Y] + 3 * control1_vertex.v_[Y]
      -3 * control2_vertex.v_[Y] + to_vertex.v_[Y];

    b[X] = 3 * solid->last_vertex_.v_[X] - 6 * control1_vertex.v_[X] +
      3 * control2_vertex.v_[X];
    b[Y] = 3 * solid->last_vertex_.v_[Y] - 6 * control1_vertex.v_[Y] +
      3 * control2_vertex.v_[Y];

    c[X] = -3 * solid->last_vertex_.v_[X] + 3 * control1_vertex.v_[X];
    c[Y] = -3 * solid->last_vertex_.v_[Y] + 3 * control1_vertex.v_[Y];

    d[X] = solid->last_vertex_.v_[X];
    d[Y] = solid->last_vertex_.v_[Y];

    f[X] = d[X];
    f[Y] = d[Y];
    df[X] = c[X] * solid->delta_ + b[X] * solid->delta2_
      + a[X] * solid->delta3_;
    df[Y] = c[Y] * solid->delta_ + b[Y] * solid->delta2_
      + a[Y] * solid->delta3_;
    d2f[X] = 2 * b[X] * solid->delta2_ + 6 * a[X] * solid->delta3_;
    d2f[Y] = 2 * b[Y] * solid->delta2_ + 6 * a[Y] * solid->delta3_;
    d3f[X] = 6 * a[X] * solid->delta3_;
    d3f[Y] = 6 * a[Y] * solid->delta3_;

    for ( unsigned int i = 0; i < solid->tessellation_steps_-1; i++ ) {

      f[X] += df[X];
      f[Y] += df[Y];

      VertexInfo vertex( f, solid->colorTess(), solid->textureTess() );

      VertexInfo normal( solid->extrusion_.normal_sign_.y_ * df[Y],
			 solid->extrusion_.normal_sign_.x_ * df[X] );

      vertex.v_[X] *= solid->vector_scale_;
      vertex.v_[Y] *= solid->vector_scale_;

      normal.normalize();

      solid->extrusion_.contour_.push_back( vertex );
      solid->extrusion_.contour_normals_.push_back( normal );

      df[X] += d2f[X];
      df[Y] += d2f[Y];
      d2f[X] += d3f[X];
      d2f[Y] += d3f[Y];
    }

    VertexInfo vertex( to, solid->colorTess(), solid->textureTess() );

    VertexInfo normal( solid->extrusion_.normal_sign_.y_ * df[Y],
		       solid->extrusion_.normal_sign_.x_ * df[X] );
  
    vertex.v_[X] *= solid->vector_scale_;
    vertex.v_[Y] *= solid->vector_scale_;

    normal.normalize();
  
    solid->extrusion_.contour_.push_back( vertex );
    solid->extrusion_.contour_normals_.push_back( normal );

    solid->last_vertex_ = to_vertex;

    return 0;
  }
#endif // OGLFT_NO_SOLID

  Texture::Texture ( const char* filename, float point_size, FT_UInt resolution )
    : Face( filename, point_size, resolution )
  {
    if ( !isValid() ) return;

    init();
  }

  Texture::Texture ( const FT_Byte* data_base, const FT_Long data_size,
		     float point_size, FT_UInt resolution )
    : Face( data_base, data_size, point_size, resolution )
  {
    if ( !isValid() ) return;

    init();
  }

  Texture::Texture ( FT_Face face, float point_size, FT_UInt resolution )
    : Face( face, point_size, resolution )
  {
    init();
  }

  void Texture::init ( void )
  {
    character_rotation_.active_ = false;
    character_rotation_.x_ = 0;
    character_rotation_.y_ = 0;
    character_rotation_.z_ = 0;

    setCharSize();

    setCharacterRotationReference( 'o' );
  }

  Texture::~Texture ( void )
  {
    clearCaches();
  }

  // Note: Changing the character rotation also clears the display list cache.

  void Texture::setCharacterRotationX ( GLfloat character_rotation_x )
  {
    if ( character_rotation_x != character_rotation_.x_ ) {
      character_rotation_.x_ = character_rotation_x;

      if ( character_rotation_.x_ != 0. || character_rotation_.y_ != 0. ||
	   character_rotation_.z_ != 0. )
	character_rotation_.active_ = true;
      else
	character_rotation_.active_ = false;

      clearCaches();
    }
  }

  void Texture::setCharacterRotationY ( GLfloat character_rotation_y )
  {
    if ( character_rotation_y != character_rotation_.y_ ) {
      character_rotation_.y_ = character_rotation_y;

      if ( character_rotation_.x_ != 0. || character_rotation_.y_ != 0. ||
	   character_rotation_.z_ != 0. )
	character_rotation_.active_ = true;
      else
	character_rotation_.active_ = false;

      clearCaches();
    }
  }

  void Texture::setCharacterRotationZ ( GLfloat character_rotation_z )
  {
    if ( character_rotation_z != character_rotation_.z_ ) {
      character_rotation_.z_ = character_rotation_z;

      if ( character_rotation_.x_ != 0. || character_rotation_.y_ != 0. ||
	   character_rotation_.z_ != 0. )
	character_rotation_.active_ = true;
      else
	character_rotation_.active_ = false;

      clearCaches();
    }
  }

  void Texture::setCharSize ( void )
  {
    for ( unsigned int f = 0; f < faces_.size(); f++ ) {
      FT_Error error = FT_Set_Char_Size( faces_[f].face_,
					 (FT_F26Dot6)( point_size_ * 64 ),
					 (FT_F26Dot6)( point_size_ * 64 ),
					 resolution_,
					 resolution_ );
      if ( error != 0 )
	return;
    }

    if ( rotation_reference_glyph_ != 0 )
      setRotationOffset();
  }

  void Texture::setRotationOffset ( void )
  {
    FT_Error error = FT_Load_Glyph( rotation_reference_face_,
				    rotation_reference_glyph_,
				    FT_LOAD_RENDER );

    if ( error != 0 )
      return;

    rotation_offset_y_ = rotation_reference_face_->glyph->bitmap.rows / 2.f;
  }

  BBox Texture::measure ( unsigned char c )
  {
    BBox bbox;
    // For starters, just get the unscaled glyph bounding box
    unsigned int f;
    FT_UInt glyph_index = 0;

    for ( f = 0; f < faces_.size(); f++ ) {
      glyph_index = FT_Get_Char_Index( faces_[f].face_, c );
      if ( glyph_index != 0 ) break;
    }

    if ( glyph_index == 0 )
      return bbox;

    FT_Error error = FT_Load_Glyph( faces_[f].face_, glyph_index,
				    FT_LOAD_DEFAULT );
    if ( error != 0 )
      return bbox;

    FT_Glyph glyph;
    error = FT_Get_Glyph( faces_[f].face_->glyph, &glyph );
    if ( error != 0 )
      return bbox;

    FT_BBox ft_bbox;
    FT_Glyph_Get_CBox( glyph, ft_glyph_bbox_unscaled, &ft_bbox );

    FT_Done_Glyph( glyph );

    bbox = ft_bbox;
    bbox.advance_ = faces_[f].face_->glyph->advance;

    return bbox;
  }

  double Texture::height ( void ) const
  {
    if ( faces_[0].face_->height > 0 )
      return faces_[0].face_->height / 64.;
    else
      return faces_[0].face_->size->metrics.y_ppem;
  }

#ifndef OGLFT_NO_QT

  BBox Texture::measure ( const QChar c )
  {
    BBox bbox;
    // For starters, just get the unscaled glyph bounding box
    unsigned int f;
    FT_UInt glyph_index = 0;

    for ( f = 0; f < faces_.size(); f++ ) {
      glyph_index = FT_Get_Char_Index( faces_[f].face_, c.unicode() );
      if ( glyph_index != 0 ) break;
    }

    if ( glyph_index == 0 )
      return bbox;

    FT_Error error = FT_Load_Glyph( faces_[f].face_, glyph_index,
				    FT_LOAD_DEFAULT );
    if ( error != 0 )
      return bbox;

    FT_Glyph glyph;
    error = FT_Get_Glyph( faces_[f].face_->glyph, &glyph );
    if ( error != 0 )
      return bbox;

    FT_BBox ft_bbox;
    FT_Glyph_Get_CBox( glyph, ft_glyph_bbox_unscaled, &ft_bbox );

    FT_Done_Glyph( glyph );

    bbox = ft_bbox;
    bbox.advance_ = faces_[f].face_->glyph->advance;

    return bbox;
  }
#endif /* OGLFT_NO_QT */
  GLuint Texture::compileGlyph ( FT_Face face, FT_UInt glyph_index )
  {
    bindTexture( face, glyph_index );

    GLuint dlist = glGenLists( 1 );
    glNewList( dlist, GL_COMPILE );

    renderGlyph( face, glyph_index );

    glEndList( );

    return dlist;
  }

  void Texture::renderGlyph ( FT_Face face, FT_UInt glyph_index )
  {
    FT_Error error = FT_Load_Glyph( face, glyph_index, FT_LOAD_DEFAULT );

    if ( error != 0 )
      return;

    TextureInfo texture_info;

    GTOCI texture_object = glyph_texobjs_.find( glyph_index );

    if ( texture_object == glyph_texobjs_.end() ) {

      bindTexture( face, glyph_index );

      texture_object = glyph_texobjs_.find( glyph_index );

      if ( texture_object == glyph_texobjs_.end() )
	return;
    }

    texture_info = texture_object->second;

    glBindTexture( GL_TEXTURE_2D, texture_info.texture_name_ );

    if ( character_rotation_.active_ ) {
      glPushMatrix();
      glTranslatef( ( texture_info.width_ / 2.f +
		      texture_info.left_bearing_ ),
		    rotation_offset_y_, 0.f );

      if ( character_rotation_.x_ != 0. )
	glRotatef( character_rotation_.x_, 1.f, 0.f, 0.f );

      if ( character_rotation_.y_ != 0. )
	glRotatef( character_rotation_.y_, 0.f, 1.f, 0.f );

      if ( character_rotation_.z_ != 0. )
	glRotatef( character_rotation_.z_, 0.f, 0.f, 1.f );

      glTranslatef( -( texture_info.width_ / 2.f +
		      texture_info.left_bearing_ ),
		    -rotation_offset_y_, 0.f );
    }

    glBegin( GL_QUADS );

    glTexCoord2i( 0, 0 );
    glVertex2i( texture_info.left_bearing_, texture_info.bottom_bearing_ );

    glTexCoord2f( texture_info.texture_s_, 0.f );
    glVertex2i( texture_info.left_bearing_ + texture_info.width_,
		texture_info.bottom_bearing_ );

    glTexCoord2f( texture_info.texture_s_, texture_info.texture_t_ );
    glVertex2i( texture_info.left_bearing_ + texture_info.width_,
		texture_info.bottom_bearing_ + texture_info.height_ );

    glTexCoord2f( 0.f, texture_info.texture_t_ );
    glVertex2i( texture_info.left_bearing_,
		texture_info.bottom_bearing_ + texture_info.height_ );
    
    glEnd();

    if ( character_rotation_.active_ ) {
      glPopMatrix();
    }

    // Drawing a character always advances the MODELVIEW.
    glTranslatef( texture_info.advance_.x / 64.f,
		  texture_info.advance_.y / 64.f,
		  0. );
  }

  void Texture::clearCaches ( void )
  {
    GDLI fgi = glyph_dlists_.begin();

    for ( ; fgi != glyph_dlists_.end(); ++fgi ) {
      glDeleteLists( fgi->second, 1 );
    }

    glyph_dlists_.clear();

    GTOI fti = glyph_texobjs_.begin();

    for ( ; fti != glyph_texobjs_.end(); ++fti ) {
      glDeleteTextures( 1, &fti->second.texture_name_ );
    }

    glyph_texobjs_.clear();
  }

  unsigned int Texture::nearestPowerCeil ( unsigned int a )
  {
    unsigned int b = a;
    unsigned int c = 1;

    if ( a == 0 ) return 1;

    // Take the log-2 of a
    for ( ; ; ) {
      if ( b == 1 )
	break;
      
      else if ( b == 3 ) {
	c *= 4;
	break;
      }

      b >>= 1;
      c *= 2;
    }
    // If it's too small, raise it another power
    if ( c < a ) c *= 2;

    return c;
  }

  MonochromeTexture::MonochromeTexture ( const char* filename, float point_size,
					 FT_UInt resolution )
    : Texture( filename, point_size, resolution )
  {}

  MonochromeTexture::MonochromeTexture ( const FT_Byte* data_base, const FT_Long data_size,
					 float point_size, FT_UInt resolution )
    : Texture( data_base, data_size, point_size, resolution )
  {}

  MonochromeTexture::MonochromeTexture ( FT_Face face, float point_size,
					 FT_UInt resolution )
    : Texture( face, point_size, resolution )
  {}

  MonochromeTexture::~MonochromeTexture ( void )
  {}

  // Round up the size of the image to a power of two, but otherwise
  // use the bitmap as is (i.e., don't expand it into separate
  // luminance and alpha components)

  GLubyte* MonochromeTexture::invertBitmap ( const FT_Bitmap& bitmap,
					     int* width, int* height )
  {
    *width = nearestPowerCeil( bitmap.width );
    *height = nearestPowerCeil( bitmap.rows );

    // MK: Alloc at least 128 Bytes to workaround potential gfx-driver out-of-bounds
    //     memory access bugs for very small 1x1 pixel textures (NVidia?)
    GLubyte* inverse = new GLubyte[ (( *width + 7) / 8 * *height) + 128 ];
    GLubyte* inverse_ptr = inverse;

    memset( inverse, 0, sizeof( GLubyte ) * ((( *width + 7 ) / 8 * *height) + 128) );

    for ( int r = 0; r < bitmap.rows; r++ ) {

      GLubyte* bitmap_ptr = &bitmap.buffer[bitmap.pitch * ( bitmap.rows - r - 1 )];

      for ( int p = 0; p < bitmap.pitch; p++ ) {

          *inverse_ptr++ = *bitmap_ptr++;
      }

      inverse_ptr += ( ( *width + 7 ) / 8 - bitmap.pitch );
    }

    return inverse;
  }

  // Hmm. This is the only routine which is different between the different
  // styles.

  void MonochromeTexture::bindTexture ( FT_Face face, FT_UInt glyph_index )
  {
    GTOCI texobj = glyph_texobjs_.find( glyph_index );

    if ( texobj != glyph_texobjs_.end() )
      return;

    // Retrieve the glyph's data.

    FT_Error error = FT_Load_Glyph( face, glyph_index, FT_LOAD_DEFAULT );

    if ( error != 0 )
      return;

    error = FT_Render_Glyph( face->glyph, ft_render_mode_mono );
	
    if ( error != 0 )
      return;

    // Double-Check that glyph bitmaps are 1bpp format, otherwise we'd crash. No-Op if this requirement is violated:
    if (face->glyph->bitmap.pixel_mode != FT_PIXEL_MODE_MONO) {
        printf("libptbdrawtext_ftgl: ERROR: Tried to draw text character without anti-aliasing, but source font does not\n");
        printf("libptbdrawtext_ftgl: ERROR: contain glyph bitmaps in 1 bpp mono format as required! Operation aborted.\n");
        return;
    }

    TextureInfo texture_info;

    glGenTextures( 1, &texture_info.texture_name_ );
    glBindTexture( GL_TEXTURE_2D, texture_info.texture_name_ );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );

    // Texture maps have be a power of 2 in size (is 1 a power of 2?), so
    // pad it out while flipping it over
    int width, height;
    GLubyte* inverted_pixmap =
      invertBitmap( face->glyph->bitmap, &width, &height );

    GLfloat red_map[2] = { background_color_[R], foreground_color_[R] };
    GLfloat green_map[2] = { background_color_[G], foreground_color_[G] };
    GLfloat blue_map[2] = { background_color_[B], foreground_color_[B] };
    GLfloat alpha_map[2] = { background_color_[A], foreground_color_[A] };

    glPixelMapfv( GL_PIXEL_MAP_I_TO_R, 2, red_map );
    glPixelMapfv( GL_PIXEL_MAP_I_TO_G, 2, green_map );
    glPixelMapfv( GL_PIXEL_MAP_I_TO_B, 2, blue_map );
    glPixelMapfv( GL_PIXEL_MAP_I_TO_A, 2, alpha_map );

    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, width, height,
		  0, GL_COLOR_INDEX, GL_BITMAP, inverted_pixmap );

    // Save a good bit of the data about this glyph
    texture_info.left_bearing_ = face->glyph->bitmap_left;
    texture_info.bottom_bearing_ = -( face->glyph->bitmap.rows
					- face->glyph->bitmap_top );
    texture_info.width_ = face->glyph->bitmap.width;
    texture_info.height_ = face->glyph->bitmap.rows;
    texture_info.texture_s_ = (GLfloat)texture_info.width_ / width;
    texture_info.texture_t_ = (GLfloat)texture_info.height_ / height;
    texture_info.advance_ = face->glyph->advance;

    glyph_texobjs_[ glyph_index ] = texture_info;

    delete[] inverted_pixmap;
  }

  GrayscaleTexture::GrayscaleTexture ( const char* filename, float point_size,
				       FT_UInt resolution )
    : Texture( filename, point_size, resolution )
  {}

  GrayscaleTexture::GrayscaleTexture ( const FT_Byte* data_base, const FT_Long data_size,
				       float point_size, FT_UInt resolution )
    : Texture( data_base, data_size, point_size, resolution )
  {}

  GrayscaleTexture::GrayscaleTexture ( FT_Face face, float point_size,
				       FT_UInt resolution )
    : Texture( face, point_size, resolution )
  {}

  GrayscaleTexture::~GrayscaleTexture ( void )
  {}

  // For the grayscale style, the luminance is the grayscale FreeType value,
  // so this just rounds up to a power of two and inverts the pixmap

  GLubyte* GrayscaleTexture::invertPixmap ( const FT_Bitmap& bitmap,
					    int* width, int* height )
  {
    *width = nearestPowerCeil( bitmap.width );
    *height = nearestPowerCeil( bitmap.rows );

    GLubyte* inverse = new GLubyte[ *width * *height ];
    GLubyte* inverse_ptr = inverse;

    memset(inverse, 0, *width * *height);

    for ( int r = 0; r < bitmap.rows; r++ ) {

      GLubyte* bitmap_ptr = &bitmap.buffer[bitmap.pitch * ( bitmap.rows - r - 1 )];

      for ( int p = 0; p < bitmap.width; p++ ) {
	*inverse_ptr++ = *bitmap_ptr++;
      }

      inverse_ptr += ( *width - bitmap.pitch );
    }
    return inverse;
  }

  // Hmm. This is the only routine which is different between the different
  // styles.

  void GrayscaleTexture::bindTexture ( FT_Face face, FT_UInt glyph_index )
  {
    GTOCI texobj = glyph_texobjs_.find( glyph_index );

    if ( texobj != glyph_texobjs_.end() )
      return;

    // Retrieve the glyph's data.

    FT_Error error = FT_Load_Glyph( face, glyph_index, FT_LOAD_DEFAULT );

    if ( error != 0 )
      return;

    error = FT_Render_Glyph( face->glyph, ft_render_mode_normal );
	
    if ( error != 0 )
      return;

    TextureInfo texture_info;

    glGenTextures( 1, &texture_info.texture_name_ );
    glBindTexture( GL_TEXTURE_2D, texture_info.texture_name_ );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );

    // Texture maps have be a power of 2 in size (is 1 a power of 2?), so
    // pad it out while flipping it over
    int width, height;
    GLubyte* inverted_pixmap =
      invertPixmap( face->glyph->bitmap, &width, &height );

    glPushAttrib( GL_PIXEL_MODE_BIT );
    glPixelTransferf( GL_RED_SCALE, foreground_color_[R] - background_color_[R] );
    glPixelTransferf( GL_GREEN_SCALE, foreground_color_[G]-background_color_[G] );
    glPixelTransferf( GL_BLUE_SCALE, foreground_color_[B]-background_color_[B] );
    glPixelTransferf( GL_ALPHA_SCALE, foreground_color_[A]-background_color_[A] );
    glPixelTransferf( GL_RED_BIAS, background_color_[R] );
    glPixelTransferf( GL_GREEN_BIAS, background_color_[G] );
    glPixelTransferf( GL_BLUE_BIAS, background_color_[B] );
    glPixelTransferf( GL_ALPHA_BIAS, background_color_[A] );

    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, width, height,
		  0, GL_LUMINANCE, GL_UNSIGNED_BYTE, inverted_pixmap );

    glPopAttrib();
    // Save a good bit of the data about this glyph
    texture_info.left_bearing_ = face->glyph->bitmap_left;
    texture_info.bottom_bearing_ = -( face->glyph->bitmap.rows
				      - face->glyph->bitmap_top );
    texture_info.width_ = face->glyph->bitmap.width;
    texture_info.height_ = face->glyph->bitmap.rows;
    texture_info.texture_s_ = (GLfloat)texture_info.width_ / width;
    texture_info.texture_t_ = (GLfloat)texture_info.height_ / height;
    texture_info.advance_ = face->glyph->advance;

    glyph_texobjs_[ glyph_index ] = texture_info;

    delete[] inverted_pixmap;
  }

  TranslucentTexture::TranslucentTexture ( const char* filename, float point_size,
					   FT_UInt resolution )
    : Texture( filename, point_size, resolution )
  {}

  TranslucentTexture::TranslucentTexture ( const FT_Byte* data_base, const FT_Long data_size,
					   float point_size, FT_UInt resolution )
    : Texture( data_base, data_size, point_size, resolution )
  {}

  TranslucentTexture::TranslucentTexture ( FT_Face face, float point_size,
					   FT_UInt resolution )
    : Texture( face, point_size, resolution )
  {}

  TranslucentTexture::~TranslucentTexture ( void )
  {}

  // For the translucent style, the luminance is saturated and alpha value
  // is the translucent FreeType value

  GLubyte* TranslucentTexture::invertPixmap ( const FT_Bitmap& bitmap,
					    int* width, int* height )
  {
    *width = nearestPowerCeil( bitmap.width );
    *height = nearestPowerCeil( bitmap.rows );

    // MK: Alloc at least 128 Bytes to workaround potential gfx-driver out-of-bounds
    //     memory access bugs for very small 1x1 pixel textures (NVidia?)
    GLubyte* inverse = new GLubyte[ (2 * *width * *height)  + 128];
    GLubyte* inverse_ptr = inverse;

    memset(inverse, 0, (size_t) (2 * *width * *height) + 128);

    for ( int r = 0; r < bitmap.rows; r++ ) {

      GLubyte* bitmap_ptr = &bitmap.buffer[bitmap.pitch * ( bitmap.rows - r - 1 )];

      for ( int p = 0; p < bitmap.width; p++ ) {
          *inverse_ptr++ = 0xff;
          *inverse_ptr++ = *bitmap_ptr++;
      }

      // MK: use bitmap.width instead of bitmap.pitch, as i think .pitch is a bug
      // and just worked by accident due to usually .pitch == .width.
      inverse_ptr += 2 * ( *width - bitmap.width );
    }

    return inverse;
  }

  // Hmm. This is the only routine which is different between the different
  // styles.

  void TranslucentTexture::bindTexture ( FT_Face face, FT_UInt glyph_index )
  {
    GTOCI texobj = glyph_texobjs_.find( glyph_index );

    if ( texobj != glyph_texobjs_.end() )
      return;

    // Retrieve the glyph's data.

    FT_Error error = FT_Load_Glyph( face, glyph_index, FT_LOAD_DEFAULT );

    if ( error != 0 )
      return;

    error = FT_Render_Glyph( face->glyph, ft_render_mode_normal );
	
    if ( error != 0 )
      return;

    TextureInfo texture_info;

    glGenTextures( 1, &texture_info.texture_name_ );
    glBindTexture( GL_TEXTURE_2D, texture_info.texture_name_ );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );

    // Texture maps have be a power of 2 in size (is 1 a power of 2?), so
    // pad it out while flipping it over
    int width, height;
    GLubyte* inverted_pixmap;
    
    // MK: We can only handle 8bpp grayscale bitmaps. If this bitmap isn't 8bpp,
    // convert it on the fly to keep going:
    if (face->glyph->bitmap.pixel_mode != FT_PIXEL_MODE_GRAY) {
        // Input bitmap is not 8bpp = 1 Byte per pixel anti-aliased grayscale image
        // as we expect. We need to convert it into 8bpp bitmaps without (aka 1 Byte) alignment:
        FT_Bitmap bitmap;
        
        // Initialize target temporary bitmap:
        FT_Bitmap_New( &bitmap );
        
        // Convert into temporary 8 bpp 'bitmap', with a 1 Byte alignment, ie. no alignment:
        if (0 != FT_Bitmap_Convert( OGLFT::Library::instance(), &(face->glyph->bitmap), &bitmap, 1 )) return;
        
        // Pass it to our inverter:
        inverted_pixmap = invertPixmap( bitmap, &width, &height );
        
        // Release temporary bitmap:
        FT_Bitmap_Done( OGLFT::Library::instance(), &bitmap );
    }
    else {
        // Input bitmap is 8bpp = 1 Byte per pixel anti-aliased grayscale image.
        // This is what we expect, so go for it:
        inverted_pixmap = invertPixmap( face->glyph->bitmap, &width, &height );
    }

    glPushAttrib( GL_PIXEL_MODE_BIT );
    glPixelTransferf( GL_RED_SCALE, foreground_color_[R] - background_color_[R] );
    glPixelTransferf( GL_GREEN_SCALE, foreground_color_[G]-background_color_[G] );
    glPixelTransferf( GL_BLUE_SCALE, foreground_color_[B]-background_color_[B] );
    glPixelTransferf( GL_ALPHA_SCALE, foreground_color_[A]-background_color_[A] );
    glPixelTransferf( GL_RED_BIAS, background_color_[R] );
    glPixelTransferf( GL_GREEN_BIAS, background_color_[G] );
    glPixelTransferf( GL_BLUE_BIAS, background_color_[B] );
    glPixelTransferf( GL_ALPHA_BIAS, background_color_[A] );

    glTexImage2D( GL_TEXTURE_2D, 0, GL_RGBA, width, height,
		  0, GL_LUMINANCE_ALPHA, GL_UNSIGNED_BYTE, inverted_pixmap );

    glPopAttrib();

    // Save a good bit of the data about this glyph
    texture_info.left_bearing_ = face->glyph->bitmap_left;
    texture_info.bottom_bearing_ = -( face->glyph->bitmap.rows
					- face->glyph->bitmap_top );
    texture_info.width_ = face->glyph->bitmap.width;
    texture_info.height_ = face->glyph->bitmap.rows;
    texture_info.texture_s_ = (GLfloat)texture_info.width_ / width;
    texture_info.texture_t_ = (GLfloat)texture_info.height_ / height;
    texture_info.advance_ = face->glyph->advance;

    glyph_texobjs_[ glyph_index ] = texture_info;

    delete[] inverted_pixmap;
  }

} // close OGLFT namespace
