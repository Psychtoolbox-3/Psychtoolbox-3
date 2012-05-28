/* Memorybuffer2Texture.c -- Minimal example of how to
 * pass a pointer to a memory buffer to Psychtoolbox,
 * so that one can create a texture out of its content.
 *
 * This simple example creates a 640 x 480 pixels color
 * image in a memory buffer and returns a pointer to it
 * to Matlab, encoded as a double value: The way we encode
 * the pointer as double may look clumsy, but it guarantees
 * that it works even on 64-Bit systems.
 *
 * Compile with: mex Memorybuffer2Texture.c
 * Then run Memorybuffer2TextureDemo.m
 *
 * Written by Mario Kleiner.
 */

#include "mex.h"

/* Byte array with 640 x 480 RGB  pixels */
unsigned char rawimage[640 * 480 * 3];

/* Definition of a 64-bit number: Operating system dependend: */
#ifndef __WIN32__
/* Unix style: */
typedef unsigned long long  psych_uint64;
#else
/* Microsoft Windows style: */
typedef LONG  psych_uint64;
#endif


/* Convert a memory address pointer into a double value: */
double PtrToDouble(void* ptr)
{
  volatile psych_uint64 ival = (psych_uint64) ptr;
  volatile double* dptr = (double*) &ival;
  volatile double outval = *dptr;
  return(outval);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
  int x,y, pos;
  double dblimgptr;
  
  /* Fill our rawimage image array with content: */
  pos = 0;
  for(y=0; y<480; y++)
    for (x=0; x<640; x++) {
      // Red channel:
      rawimage[pos++] = (unsigned char) x;
      // Green channel:
      rawimage[pos++] = (unsigned char) y;
      // Blue channel:
      rawimage[pos++] = (unsigned char) 0;
    }

  /* Convert our memory pointer to something Psychtoolbox
     understands, i.e., convert it into a double value:
  */
  dblimgptr = PtrToDouble((void*) &rawimage);

  /* Return double-encoded pointer to Matlab: */
  plhs[0]=mxCreateDoubleMatrix(1, 1, mxREAL);
  *mxGetPr(plhs[0]) = dblimgptr;
}
