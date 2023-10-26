/*
 *  bitmap.h
 *  PsychToolbox
 *
 *  Created by Nuha Jabakhanji on 10-10-01.
 *
 */

/*
 * Windows BMP file definitions for OpenGL.
 *
 * Written by Michael Sweet.
 */

#ifndef _BITMAP_H_
#define _BITMAP_H_

#  ifndef WIN32
	typedef struct                       /**** BMP file header structure ****/
    {
		unsigned short bfType;           /* Magic number for file */
		unsigned int   bfSize;           /* Size of file */
		unsigned short bfReserved1;      /* Reserved */
		unsigned short bfReserved2;      /* ... */
		unsigned int   bfOffBits;        /* Offset to bitmap data */
    } BITMAPFILEHEADER;
	
#  define BF_TYPE 0x4D42            
	
	typedef struct                       /**** BMP file info structure ****/
    {
		unsigned int   biSize;           /* Size of info header */
		int            biWidth;          /* Width of image */
		int            biHeight;         /* Height of image */
		unsigned short biPlanes;         /* Number of color planes */
		unsigned short biBitCount;       /* Number of bits per pixel */
		unsigned int   biCompression;    /* Type of compression to use */
		unsigned int   biSizeImage;      /* Size of image data */
		int            biXPelsPerMeter;  /* X pixels per meter */
		int            biYPelsPerMeter;  /* Y pixels per meter */
		unsigned int   biClrUsed;        /* Number of colors used */
		unsigned int   biClrImportant;   /* Number of important colors */
    } BITMAPINFOHEADER;
		
	typedef struct                       /**** Bitmap information structure ****/
    {
		BITMAPINFOHEADER bmiHeader;      /* Image header */
    } BITMAPINFO;
#  endif /* !WIN32 */

extern psych_uint8* LoadDIBitmap(const char *filename, BITMAPINFO **info);

#endif
