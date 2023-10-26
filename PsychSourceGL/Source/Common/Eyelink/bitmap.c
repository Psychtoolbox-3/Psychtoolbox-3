/*

 bitmap.c		
 
 AUTHORS:
 
 nuha@sr-research.com nj
 based on code by Michael Sweet 
 
 PLATFORMS:
 
 All.

*/

#include "PsychEyelink.h"

#include "bitmap.h"
#include <errno.h>

#ifdef WIN32

/*
 * 'LoadDIBitmap()' - Load a DIB/BMP file from disk.
 *
 * Returns a pointer to the bitmap if successful, NULL otherwise...
 */

psych_uint8* LoadDIBitmap(const char *filename, BITMAPINFO **info)    /* O - Bitmap information */
{
    FILE				*fp;          // Open file pointer
    psych_uint8			*bits;        // Bitmap pixel bits 
    int					bitsize;      // Size of bitmap 
    int					infosize;     // Size of header information     
	BITMAPFILEHEADER	header;       // File header 

	if ((fp = fopen(filename, "rb")) == NULL){
		printf("Bitmap Loader: can't open file\n");
		return (NULL);
	}
		
    if (fread(&header, sizeof(BITMAPFILEHEADER), 1, fp) < 1)
	{
		printf("Bitmap Loader: No header information\n");
		fclose(fp);
        return (NULL);
	}
	
    if (header.bfType != 'MB')	// Check for BM reversed... 
	{
		printf("Bitmap Loader: not a bitmap file\n");
        fclose(fp);
        return (NULL);
	}

    infosize = header.bfOffBits - sizeof(BITMAPFILEHEADER);
    if ((*info = (BITMAPINFO *)malloc(infosize)) == NULL)
	{
        printf("Bitmap Loader: bitmapinfo memory allocation failed\n");
        fclose(fp);
        return (NULL);
	}
		
    if (fread(*info, 1, infosize, fp) < infosize)
	{
        printf("Bitmap Loader: couldn't read memory\n");
        free(*info);
        fclose(fp);
        return (NULL);
	}

	if ((*info)->bmiHeader.biBitCount != 32 && (*info)->bmiHeader.biBitCount != 24){
		printf("Bitmap loading error: Bitmap must be 24 bits or 32 bits\n");
        free(*info);
        fclose(fp);
		return (NULL);
	}
    
	if ((*info)->bmiHeader.biCompression != 0){
		printf("Bitmap Loading Error: Bitmap must not be compressed\n");
        free(*info);
        fclose(fp);
		return (NULL);		
	}

    // Now that we have all the header info read in, allocate memory for 
	// the bitmap and read *it* in...                                    
    if ((bitsize = (*info)->bmiHeader.biSizeImage) == 0)
        bitsize = ((*info)->bmiHeader.biWidth * (*info)->bmiHeader.biBitCount + 7) / 8 * abs((*info)->bmiHeader.biHeight);
	else
		bitsize = (*info)->bmiHeader.biSizeImage;

    if ((bits = malloc(bitsize)) == NULL) {
        // Couldn't allocate memory - return NULL! 
		printf("Bitmap Loader:  Couldn't allocate memory for bitmap\n");
		
        free(*info);
        fclose(fp);
        return (NULL);
	}
	 
    if (fread(bits, 1, bitsize, fp) < bitsize) {
        // Couldn't read bitmap - free memory and return NULL!
		printf("Bitmap Loader:  Couldn't read bitmap\n");
		
        free(*info);
        free(bits);
        fclose(fp);
        return (NULL);
	}	

	// OK, everything went fine - return the allocated bitmap...
    fclose(fp);

	if (bits == NULL){
		printf("Bitmap Loader: output is NULL\n");
        free(*info);		
		return (NULL);
	}

    return (bits);
}

#else /* !WIN32 */

/*
 * Functions for reading and writing 16- and 32-bit little-endian integers.
 */

static unsigned short read_word(FILE *fp);
static unsigned int   read_dword(FILE *fp);
static int            read_long(FILE *fp);

psych_uint8* LoadDIBitmap(const char *filename, BITMAPINFO **info)
{	
	int              bitsize;      // Size of bitmap 
    int              infosize;     // Size of header information
    int				 index;
	int				 newIndex = 0;		
	int				 r; 
	int				 c;
	int				 pitch;	
	FILE             *fp;          // Open file pointer 
    psych_uint8      *bits;        // Bitmap pixel bits 
    BITMAPFILEHEADER header;       // File header 	

    if ((fp = fopen(filename, "rb")) == NULL) {        
		printf("Bitmap Loader: can't open file\n");
		return (NULL);
	}
	
    // Read the file header and any following bitmap information... 
    header.bfType      = read_word(fp);
    header.bfSize      = read_dword(fp);
    header.bfReserved1 = read_word(fp);
    header.bfReserved2 = read_word(fp);
    header.bfOffBits   = read_dword(fp);
	
    if (header.bfType != BF_TYPE) // Check for BM reversed... 
	{
        // Not a bitmap file - return NULL...
        fclose(fp);
		printf("Bitmap Loader: not a bitmap file\n");		
        return (NULL);
	}
	
    infosize = header.bfOffBits - 18;
    if ((*info = (BITMAPINFO *)malloc(sizeof (BITMAPINFO))) == NULL)
	{
        // Couldn't allocate memory for bitmap info - return NULL... 
        fclose(fp);
		printf("Bitmap Loader: bitmapinfo memory allocation failed\n");
		
        return (NULL);
	}
	
    (*info)->bmiHeader.biSize          = read_dword(fp);
    (*info)->bmiHeader.biWidth         = read_long(fp);
    (*info)->bmiHeader.biHeight        = read_long(fp);
    (*info)->bmiHeader.biPlanes        = read_word(fp);
    (*info)->bmiHeader.biBitCount      = read_word(fp);
    (*info)->bmiHeader.biCompression   = read_dword(fp);
    (*info)->bmiHeader.biSizeImage     = read_dword(fp);
    (*info)->bmiHeader.biXPelsPerMeter = read_long(fp);
    (*info)->bmiHeader.biYPelsPerMeter = read_long(fp);
    (*info)->bmiHeader.biClrUsed       = read_dword(fp);
    (*info)->bmiHeader.biClrImportant  = read_dword(fp);

	if ((*info)->bmiHeader.biBitCount != 32 && (*info)->bmiHeader.biBitCount != 24){
		printf("Bitmap loading error: Bitmap must be either 24 or 32 bits");
		free(*info);
		fclose(fp);
		return (NULL);
	}
	
	if ((*info)->bmiHeader.biCompression != 0){
		printf("Bitmap Loading Error: Bitmap must not be compressed");
		free(*info);
		fclose(fp);
		return (NULL);		
	}
	
    if ((bitsize = (*info)->bmiHeader.biSizeImage) == 0)
        bitsize = ((*info)->bmiHeader.biWidth * (*info)->bmiHeader.biBitCount + 7) / 8 * abs((*info)->bmiHeader.biHeight);
	else
		bitsize = (*info)->bmiHeader.biSizeImage;

	if ((bits = malloc(bitsize)) == NULL)
	{
		free(*info);
		fclose(fp);
		printf("Bitmap Loader:  memory allocation failed\n");		
		return (NULL);
	}

	if (fread(bits, 1, bitsize, fp) < bitsize)
	{
		free(*info);
		free(bits);
		fclose(fp);
		printf("Bitmap Loader: couldn't read memory\n");
		return (NULL);
	}

    fclose(fp);
    	
	return (bits);
}


/*
 * 'read_word()' - Read a 16-bit unsigned integer.
 */

static unsigned short read_word(FILE *fp)       /* I - File to read from */
{
    unsigned char b0, b1; 	
    b0 = getc(fp);
    b1 = getc(fp);
	
    return ((b1 << 8) | b0);
}


/*
 * 'read_dword()' - Read a 32-bit unsigned integer.
 */

static unsigned int  read_dword(FILE *fp)    
{
    unsigned char b0, b1, b2, b3; 
	
    b0 = getc(fp);
    b1 = getc(fp);
    b2 = getc(fp);
    b3 = getc(fp);
	
    return ((((((b3 << 8) | b2) << 8) | b1) << 8) | b0);
}


/*
 * 'read_long()' - Read a 32-bit signed integer.
 */
static int read_long(FILE *fp)       {
    unsigned char b0, b1, b2, b3; 
	
    b0 = getc(fp);
    b1 = getc(fp);
    b2 = getc(fp);
    b3 = getc(fp);
	
    return ((int)(((((b3 << 8) | b2) << 8) | b1) << 8) | b0);
}

#endif /* WIN32 */
