/*
 * memcpyCudaOpenGL.cpp -- Fast memory copy between CUDA and OpenGL resources.
 *
 * This function, compiled into a MEX file for Mathworks Matlab or GNU/Octave
 * allows to exchange data between NVidia CUDA memory buffers and other resources
 * and OpenGL resources, such as textures, renderbuffers, framebuffer objects,
 * and pixelbuffer objects. Data exchange can happen via copying data between
 * OpenGL memory and CUDA memory, e.g., for use with high-level GPGPU toolkits,
 * or via shared access to the same underlying OpenGL resource, e.g., for even
 * more efficient use with low-level GPGPU toolkits.
 *
 * To avoid redundant calls to expensive high-overhead CUDA functions like
 * CUDA-OpenGL resource (un)registration and resource (un)mapping, it uses a
 * LRU cache with a default size of 8 slots, configurable to up to 1024 slots.
 * The cache stores registration and mapping info for all OpenGL resources which
 * are shared between CUDA and OpenGL and tries to keep resources registered as
 * long as possible and mapped as long as possible. If the cache reaches its
 * capacity, it uses a LRU strategy to decide which items to evict and thereby
 * unmap and unregister.
 *
 * Efficient use of the cache to minimize such operations is important, as they
 * can be very expensive. Example: CUDA-5.0 + GPUmat toolkit on 64-Bit Matlab under
 * 64-Bit Ubuntu Linux 12.04.1 LTS on a 8-core MacPro with a NVidia GeForce-8800:
 *
 * Time for 1024x768 pixel RGBA 32 bpc floating point data copies between CUDA
 * (GPUmat GPUsingle() matrix) and OpenGL texture (PTB texture of equivalent size):
 * Pure copy: < 0.5 msecs. Corresponds to about 97% of GPU peak bandwidth.
 * With mapping and unmapping: ~ 1.0 msecs.
 * With registering and unregistering of resources: ~ 4.9 msecs!
 * Non-Cached OpenGL interop can cause a 10x performance hit just for communication,
 * so effective caching is important.
 *
 *
 * See accompanying M-File for building the mex file.
 * See M-Wrapper file for actually using this low-level function.
 *
 * This file was originally developed for Psychtoolbox-3 (www.psychtoolbox.org)
 * but can be used for interop between any Matlab or Octave based CUDA and
 * OpenGL toolkit.
 *
 *
 * Copyright (c) 2013 Mario Kleiner. Licensed to you under the MIT license:
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to permit
 * persons to whom the Software is furnished to do so, subject to the
 * following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
 * NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/* MATLAB / OCTAVE includes: */
#include "mex.h"

#include "string.h"

#ifdef _WIN32
#include <windows.h>
#endif

/* CUDA includes: */
#include "cuda_runtime_api.h"
#include "cuda_gl_interop.h"

/* Global verbosity level and firsttime init flag: */
int verbose = 0;
static int firsttime = 1;

/* Global CUDA error state: */
cudaError_t cudastatus;

/* LRU resource registration/mapping cache: */
typedef struct cudares {
  unsigned int glhandle;    /* OpenGL object handle: Used as search key in the cache. 0 = Free slot. */
  unsigned int gltarget;    /* OpenGL object type: Texture, renderbuffer, pixelbuffer etc. */
  cudaGraphicsResource_t resource;  /* CUDA resource descriptor if OpenGL object registered with CUDA, NULL otherwise. */
  unsigned int ismapped;    /* 0 = Resource not mapped, 1 = Resource mapped. */
  unsigned int mapflags;    /* Resource mapping flags used. */
  unsigned long long lastaccess;  /* Time in "cache-ticks" (0 = Never) when this resource was last accessed: For LRU replacement. */
} cudares;

/* The resource mapping cache: Has MAX_CACHE_SLOTS slots maximum capacity. */
#define MAX_CACHE_SLOTS 1024
struct cudares resourceCache[MAX_CACHE_SLOTS];

/* Effective cache capacity: */
int cachesize;

/* First time cache full warning already emitted? */
int firstLRUCycle;

/* Current cache reference time for LRU replacement: */
unsigned long long cacheclock;

void doCudaUnmap(int i)
{
  /* Unmap if it is actually mapped */
  if (resourceCache[i].resource && resourceCache[i].ismapped) {
    /* Mark as unmapped */
    resourceCache[i].ismapped = 0;
    
	  cudastatus = cudaGraphicsUnmapResources(1, &(resourceCache[i].resource), 0);
	  if (cudastatus != cudaSuccess) {
	    mexPrintf("\nmemcpyCudaOpenGL: ERROR in cudaGraphicsUnmapResources(): %s\n", cudaGetErrorString(cudastatus));
	  }
	  
    if (verbose) mexPrintf("\n%s: doCudaUnmap(%i).\n", mexFunctionName(), i);	  
  }
}

void doCudaUnregister(int i)
{
  /* Unmap if it is mappped: */
  doCudaUnmap(i);
  
  /* Unregister if it is registered: */
  if (resourceCache[i].resource) {
    cudastatus = cudaGraphicsUnregisterResource(resourceCache[i].resource);
    if (cudastatus != cudaSuccess) {
      mexPrintf("\nmemcpyCudaOpenGL: ERROR in cudaGraphicsUnregisterResource(): %s\n", cudaGetErrorString(cudastatus));
    }
    
    if (verbose) mexPrintf("\n%s: doCudaUnregister(%i).\n", mexFunctionName(), i);    
  }
  
  /* Clear out and thereby free the cache slot: */
  memset(&(resourceCache[i]), 0, sizeof(resourceCache[i]));
}

/* Search for an item in the cache: */
int cacheLookup(unsigned int glhhandle, unsigned int gltarget)
{
  int i;

  /* Iterate over all cache slots and search for (glhandle,gltarget) key: */
  for (i = 1; i < cachesize; i++) {
    /* In use? */
    if ((resourceCache[i].glhandle == glhhandle) && (resourceCache[i].gltarget == gltarget)) {
      /* Cache hit: Update access timestamp and return. */
      resourceCache[i].lastaccess = cacheclock;
      return(i);
    }
  }

  /* Cache miss: */
  return(-1);
}

/* Enter an item into the cache, unless it is already there, in which case its index is returned: */
int cacheInsert(unsigned int glhhandle, unsigned int gltarget, unsigned int mapflags)
{
  int i, lruvictim;
  unsigned long long curage;
  unsigned long long maxage = 0;

  /* Iterate over all cache slots and search for (glhandle,gltarget) key:
     We omit slot 0 as it is reserved. Otherwise we couldn't signal caching state
     to caller via the sign of the return value, as zero is ambiguous.
   */
  for (i = 1; i < cachesize; i++) {
    /* Found? */
    if ((resourceCache[i].glhandle == glhhandle) && (resourceCache[i].gltarget == gltarget)) {
      /* Cache hit: Update access timestamp and return slot. */
      resourceCache[i].lastaccess = cacheclock;
      
      /* Well, almost. We need to make sure that this resource is also registered
         with compatible 'mapflags', or bad things would happen.
       */
      if (mapflags != resourceCache[i].mapflags) {
        /* Incompatible flags. Need to unregister the resource and then return
         * its former host slot 'i' as "free" slot, so calling code thinks the
         * resource is not in cache and reregisters it with proper mapflags.
         */
         if (verbose) mexPrintf("\n%s: cacheInsert(%i): Incompatible cached map flags [%i vs. %i]. Triggering reregister.\n", mexFunctionName(), i, mapflags, resourceCache[i].mapflags);
         doCudaUnregister(i);
         return(-i);
      }
     
      /* Registered in a directly useable way: */ 
      return(i);
    }

    /* Nope. Keep track of first free slot, or slot which hasn't been accessed for
       the longest time, ie., the victim for lru replacement if neccessary. */
    curage = cacheclock - resourceCache[i].lastaccess;
    if (curage > maxage) {
      maxage = curage;
      lruvictim = i;
    }
  }

  /* If we reach this point, then the object isn't in the cache yet and needs
     to be entered into a free cache slot, if any. lruvictim is our candidate
     slot. If its lastaccess timestamp is zero, then we got a free slot for use.
     Otherwise our cache is full and lruvictim is the slot that wasn't used for
     the longest time which makes it our victim for expelling it from the cache.
   */
  if (resourceCache[lruvictim].lastaccess > 0) {
    /* Occupied slot. Free it up by unmapping and unregistering the cached
       OpenGL resource.
     */
    doCudaUnregister(lruvictim);
    if (firstLRUCycle) {
      firstLRUCycle = 0;
      mexPrintf("%s: Warning: CUDA<->OpenGL interop cache soft-limit of %i slots reached. Cache full!\n", mexFunctionName(), cachesize - 1);
      mexPrintf("%s: Warning: Will switch to LRU cache replacement mode now. Performance may be degraded!\n", mexFunctionName());
      mexPrintf("%s: Warning: Increase the cache soft-limit, or check your code for unintended resource leaks.\n", mexFunctionName());
    }

    if (verbose) mexPrintf("\n%s: cacheInsert(%i): LRU-Replacement done (Victims age %i).\n", mexFunctionName(), lruvictim, maxage);
  }
  else if (verbose) mexPrintf("\n%s: cacheInsert(%i): Inserted into free slot.\n", mexFunctionName(), lruvictim);

  /* Cache miss: Return the (now) free cache slot 'lruvictim' to be filled. Return
     a negative index to indicate this was a cache-miss and the cache needs to be filled.
   */
  return(-lruvictim);
}

/* Flush the resourceCache - Unmap and unregister all cached resources: */
void cacheFlush(void)
{
  int i;

  if (verbose) mexPrintf("\n%s: cacheFlush.\n", mexFunctionName());

  /* Iterate over all cache slots and unmap/unregister their resources: */
  for (i = 1; i < cachesize; i++) {
    /* In use? */
    if (resourceCache[i].glhandle) {
      doCudaUnregister(i);
    }
  }

  /* Nothing valid in cache anymore. Clear and reset it. */
  cacheclock = 0;
  memset(resourceCache, 0, sizeof(resourceCache[0]) * MAX_CACHE_SLOTS);

  /* Cache is in a pristine "as good as new" condition. */
  return;
}

/* Increment cache time and thereby age of all cached items for LRU ops: */
void ageCache(void)
{
  int i;
  unsigned long long lastaccess, mint = -1, maxt = 0;

  /* Increase age of items in cache by progressing time: */
  cacheclock++;

  /* Wraparound? */
  if (cacheclock == 0) {
    /* Clock wrapped around to zero. We need to recalculate a proper
       age and cacheclock value to keep LRU statistics intact:
       Note: This is mostly an exercise in correctness. In typical
       usage scenarios for this code, e.g., for visual stimulus presentation
       with Psychtoolbox or other OpenGL rendering toolkits, it would take
       anywhere between days and months of uninterruped operation to hit
       the wraparound. Few vision science sessions last more than a few hours.
     */
    for (i = 1; i < cachesize; i++) {
      lastaccess = resourceCache[i].lastaccess;
      if (lastaccess > 0) {
        mint = (lastaccess < mint) ? lastaccess : mint;
        maxt = (lastaccess > maxt) ? lastaccess : maxt;
      }
    }

    if (maxt > 0) {
      /* Shift access times of all valid items by mint - 1, so that
       * the oldest item has access time 1, and the newest item has
       * a much lower access time:
       */
      mint = mint - 1;
      maxt -= mint;
      for (i = 1; i < cachesize; i++) {
        if (resourceCache[i].lastaccess > 0) {
      	  resourceCache[i].lastaccess -= mint;
        }
      }

      /* Adjust cacheclock to maxt + 1. This way, the scale is compressed,
         but the relative order is preserved, ie., distances between cacheclock
         and lastaccess times of all active items in the cache are as before the
         wraparound and LRU semantics is preserved across the wraparound.
       */
      cacheclock = maxt + 1;
    }
    else {
      /* All items were inactive the whole time. Just set clock to 1
         and all is good, as all items have lastaccess time zero.
       */
      cacheclock = 1;
    }
  }
}

void unregisterResource(unsigned int glhandle, unsigned int gltarget)
{
  /* Look it up in cache. */
  int i = cacheLookup(glhandle, gltarget);
  if (i > 0) {
    /* There it is! Unmap and unregister it, release its cache slot: */
    doCudaUnregister(i);
  }
}

void unmapResource(unsigned int glhandle, unsigned int gltarget)
{
  /* Look it up in cache. */
  int i = cacheLookup(glhandle, gltarget);
  if (i > 0) {
    /* There it is! Unmap it: */
    doCudaUnmap(i);
  }
}

void usageExit(int errorAbort)
{
  mexPrintf("\nUsage:\n");
  mexPrintf("======\n\n");
  mexPrintf("%s(0); - Unmap and unregister all resources. Automatically called at shutdown.\n\n", mexFunctionName());
  mexPrintf("%s(1, glObject, glTarget); - Unmap resource glObject of type glTarget\n\n", mexFunctionName());
  mexPrintf("%s(2, glObject, glTarget); - Unregister (and unmap if mapped) resource glObject of type glTarget\n\n", mexFunctionName());
  mexPrintf("%s(3, glObject, glTarget, gpuptr, nrbytes, direction, keepmapped, mapflags);\n", mexFunctionName());
  mexPrintf("- Copy data in VRAM between a CUDA memory buffer and a OpenGL resource. The resource gets\n");
  mexPrintf("  registered for use with CUDA if it isn't registered already, and mapped for use if it\n");
  mexPrintf("  isn't mapped already. 'glObject' and 'glTarget' are OpenGL object and target type handles\n");
  mexPrintf("  of the OpenGL resource. 'gpuptr' is a CUDA memory pointer into the CUDA memory buffer.\n");
  mexPrintf("  'nrbytes' = Number of bytes to copy. 'direction' 0 = OpenGL->CUDA, 1 = CUDA->OpenGL.\n");
  mexPrintf("  'keepmapped' 0 = Unmap immediately after copy [default], 1 = Keep mapped for future use.\n");
  mexPrintf("  Mapping and unmapping comes at some overhead, so should be avoided if possible. However,\n");
  mexPrintf("  OpenGL can't operate on mapped buffers, so one must unmap before OpenGL object is used again\n");
  mexPrintf("  by OpenGL. For a processing flow where a specific object is used as input and output for a\n");
  mexPrintf("  CUDA computation and not used by OpenGL during the CUDA computation, it can make sense to\n");
  mexPrintf("  keep the object mapped inbetween, avoiding one extra map/unmap cycle.\n");
  mexPrintf("  'mapflags': Optional, auto-selected. If set to a >= 0 value, use them as flags for CUDA resource registration.\n\n");
  mexPrintf("%s(4, glObject, glTarget, gpuptr, nrbytes, direction, keepmapped, mapflags);\n\n", mexFunctionName());
  mexPrintf("- Like 3, but avoid data copy. Instead return a uint64 64-Bit memory pointer to the actual mapped resource.\n");
  mexPrintf("  This implies 'keepmapped' == 1, regardless of what you specify.\n");
  mexPrintf("  This allows C-Level plugins to access the mapped resource directly via CUDA pointers.\n\n");
  mexPrintf("%s(5, newlimit); - Increase soft-limit of LRU cache to 'newlimit'.\n\n", mexFunctionName());
  mexPrintf("%s(6, verbose); - Set 'verbose' flag. Defaults to zero. 1 = Verbose status output.\n\n", mexFunctionName());

  /* mexPrintf("%s(5) \n", mexFunctionName()); */

  if (errorAbort) mexErrMsgTxt("Wrong number of arguments. See usage help above.");

  return;
}

void mexExit(void)
{
  cacheFlush();
  firsttime = 1;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  unsigned int cmd, glhandle, gltarget, direction, flags, keepmapped, ispbo;
  unsigned long nrbytes;
  void *gpuptr;
  int slot = 0;
  cudaGraphicsResource_t resource = NULL;
  struct cudaArray *mappedArray = NULL;
  void* mappedPtr = NULL;
  size_t mappedSize = 0;

  /* Be optimistic, assume success unless told otherwise: */
  cudastatus = cudaSuccess;
  
  if (firsttime) {
      firsttime = 0;

      mexPrintf("\n%s: A simple CUDA <=> OpenGL interoperation interface.\n", mexFunctionName());
      mexPrintf("(c) 2013 by Mario Kleiner. Licensed to you under the MIT license.\n\n");

      /* Reset cache clock to zero and clear the cache: */
      cacheclock = 0;
      memset(resourceCache, 0, sizeof(resourceCache[0]) * MAX_CACHE_SLOTS);

      /* Start off with an effective cache capacity of 8 slots (1 slot is blocked from use): */
      cachesize = 8 + 1;
      firstLRUCycle = 1;

      /* Make sure the cache is flushed at mex file shutdown time: */
      mexAtExit(mexExit);
  }
  
  /* Retrieve command code: Give usage info if none given. */
  if (nrhs < 1) { usageExit(0); return; }

  cmd = (unsigned int) mxGetScalar(prhs[0]);

  /* Change of verbosity? */
  if (cmd == 6) {
    if (nrhs < 2) usageExit(1);  
    verbose = (unsigned int) mxGetScalar(prhs[1]);
    if (verbose) mexPrintf("\n%s: Verbose tracing of operations enabled.\n", mexFunctionName());
    return;
  }

  /* Resizing the LRU cache requested? */
  if (cmd == 5) {
    if (nrhs < 2) usageExit(1);

    /* Reset LRU cache full warning: */
    firstLRUCycle = 1;
    
    slot = (unsigned int) mxGetScalar(prhs[1]);
    
    /* Increment request by 1 to compensate for the "lost" slot 0: */
    slot = slot + 1;
    
    /* Child protections: */
    if (slot > MAX_CACHE_SLOTS) {
      mexPrintf("%s: Requested new softlimit %i for cache exceeds compiled in maximum %i. Will clamp to maximum.\n", mexFunctionName(), slot - 1, MAX_CACHE_SLOTS - 1);
      cachesize = MAX_CACHE_SLOTS;
      return;
    }
    
    if (slot < cachesize) {
      /* Shrinking the cache requested. This implies a full cache flush: */
      mexPrintf("%s: Requested new softlimit %i for cache is smaller than old softlimit %i. Will flush the cache before shrinking it.\n", mexFunctionName(), slot - 1, cachesize - 1);
      cacheFlush();
    }
    
    /* Set new softlimit: */
    cachesize = slot;
    mexPrintf("%s: New softlimit for LRU cache set to %i slots.\n", mexFunctionName(), cachesize - 1);
    
    return;
  }

  if (cmd == 0) {
    /* Cache flush requested: */
    cacheFlush();
    return;
  }

  /* Following ops require at least object handle and target type: */
  if (nrhs < 3) usageExit(1);

  /* Time to increment the age of our cached items by a clock tick: */
  ageCache();

  /* Retrieve OpenGL object handle to our image buffer: */
  glhandle = (unsigned int) mxGetScalar(prhs[1]);

  /* Get GLEnum target: */
  gltarget = (unsigned int) mxGetScalar(prhs[2]);
  
  if (cmd == 1) {
    /* Unmap resource if it is mapped: */
    unmapResource(glhandle, gltarget);
    return;
  }

  if (cmd == 2) {
    /* Unmap and unregister resource if it is mapped and/or registered: */
    unregisterResource(glhandle, gltarget);
    return;
  }

  if (nrhs < 6) usageExit(1);

  /* Retrieve CUDA memory pointer to source/destination CUDA memory buffer: */
  gpuptr = (void*) (unsigned long) mxGetScalar(prhs[3]);
  
  /* Retrieve number of bytes to copy: */
  nrbytes = (unsigned long) mxGetScalar(prhs[4]);
  
  /* Retrieve direction: 0 = OpenGL -> CUDA, 1 = CUDA -> OpenGL : */
  direction = (unsigned int) mxGetScalar(prhs[5]);

  /* Retrieve optional 'keepmapped' flag. */
  keepmapped = 0;
  if (nrhs >= 7) keepmapped = (unsigned int) mxGetScalar(prhs[6]);

  /* Define CUDA optimization flags, depending if this is a OpenGL->CUDA or
   * CUDA->OpenGL copy operation.
   */
  if ((nrhs >= 8) && (mxGetScalar(prhs[7]) >= 0)) {
    /* Override map flags provided. Use them: */
    flags = (unsigned int) mxGetScalar(prhs[7]);
  }
  else {
    /* Use auto-selected map flags: */
    flags = (direction) ? cudaGraphicsRegisterFlagsWriteDiscard : cudaGraphicsRegisterFlagsReadOnly;
  }

  /* Is gltarget a OpenGL pixelbuffer object? Check for gltarget == GL_PACK_BUFFER or GL_UNPACK_BUFFER. */
  ispbo = (gltarget == 35051 || gltarget == 35052) ? 1 : 0;

  /* Copy of data or mapped resource access pointer requested? */
  if (cmd == 3 || cmd == 4) {
    /* Register OpenGL object with CUDA as 'resource': */
    
    /* Already in cache? This would mean it is registered already with compatible mapping flags: */
    slot = cacheInsert(glhandle, gltarget, flags);
    if (slot < 0) {
      /* Not yet in cache. This means it is not registered at this time, either because it
         wasn't registered at all, or because it was registered with incompatible 'flags',
         so it just got unregistered and expelled from the cache. In any case, we need to
         insert it into the cache and register it. -slot is the free target slot for this
         purpose.
       */
       
      /* Turn slot into something useful: */
      slot = -slot;
    
      if (ispbo) {
        /* OpenGL Pixelbuffer object (GL_PACK_BUFFER or GL_UNPACK_BUFFER): */
        cudastatus = cudaGraphicsGLRegisterBuffer(&(resourceCache[slot].resource), glhandle, flags);
      }
      else {
        /* OpenGL texture or renderbuffer object: */
        cudastatus = cudaGraphicsGLRegisterImage(&(resourceCache[slot].resource), glhandle, gltarget, flags);
      }
      
      if (cudastatus != cudaSuccess) {
        mexPrintf("\nmemcpyCudaOpenGL: ERROR in %s(): %s\n", (ispbo) ? "cudaGraphicsGLRegisterBuffer" : "cudaGraphicsGLRegisterImage", cudaGetErrorString(cudastatus));
        resourceCache[slot].resource = NULL;
        goto err_final;
      }

      if (verbose) mexPrintf("\n%s: cacheInsert(%i): CUDA resource registered (globject %i, gltarget %i, flags %i).\n", mexFunctionName(), slot, glhandle, gltarget, flags);
      
      /* Fill cache slot: */
      resourceCache[slot].glhandle = glhandle;
      resourceCache[slot].gltarget = gltarget;
      resourceCache[slot].mapflags = flags;
      resourceCache[slot].lastaccess = cacheclock;
      resourceCache[slot].ismapped = 0;
    }
      
    /* At this point, the resource is stored in slot 'slot' of the cache and registered in a compatible way: */
    
    /* Map the 'resource', unless it is already mapped: */
    if (!resourceCache[slot].ismapped) {
      /* Map it: */
      cudastatus = cudaGraphicsMapResources(1, &(resourceCache[slot].resource), 0);
      if (cudastatus != cudaSuccess) {
        mexPrintf("\nmemcpyCudaOpenGL: ERROR in cudaGraphicsMapResources(): %s\n", cudaGetErrorString(cudastatus));
        goto err_unregister;
      }

      if (verbose) mexPrintf("\n%s: CUDA resource %i mapped (globject %i, gltarget %i, flags %i).\n", mexFunctionName(), slot, glhandle, gltarget, flags);
      
      /* Successfully mapped: */
      resourceCache[slot].ismapped = 1;
    }
    
    /* Get simpler handle: */
    resource = resourceCache[slot].resource;
    
    /* Get mapped resource image array handle or PBO pointer: */
    if (ispbo) {
      cudastatus = cudaGraphicsResourceGetMappedPointer(&mappedPtr, &mappedSize, resource);
    }
    else {
      cudastatus = cudaGraphicsSubResourceGetMappedArray(&mappedArray, resource, 0, 0);
    }
    
    if (cudastatus != cudaSuccess) {
      mexPrintf("\nmemcpyCudaOpenGL: ERROR in %s(): %s\n", (ispbo) ? "cudaGraphicsResourceGetMappedPointer" : "cudaGraphicsSubResourceGetMappedArray", cudaGetErrorString(cudastatus));
      goto err_unmap;
    }
  }
  
  /* Copy of PBO data between CUDA and OpenGL requested? */
  if (cmd == 3 && ispbo) {
    /* Copy from OpenGL PBO to CUDA buffer? */
    if (direction == 0) {
        /* OpenGL -> CUDA copy: */
        cudastatus = cudaMemcpyAsync(gpuptr, (const void*) mappedPtr, (size_t) nrbytes, cudaMemcpyDeviceToDevice, 0);
        if (cudastatus != cudaSuccess) {
            mexPrintf("\nmemcpyCudaOpenGL: ERROR in cudaMemcpyAsync(): %s\n", cudaGetErrorString(cudastatus));
            goto err_unmap;
        }
    }

    if (direction == 1) {
        /* CUDA -> OpenGL copy: */
        cudastatus = cudaMemcpyAsync(mappedPtr, (const void*) gpuptr, (size_t) nrbytes, cudaMemcpyDeviceToDevice, 0);
        if (cudastatus != cudaSuccess) {
            mexPrintf("\nmemcpyCudaOpenGL: ERROR in cudaMemcpyAsync(): %s\n", cudaGetErrorString(cudastatus));
            goto err_unmap;
        }
    }
  }

  /* Copy of texture or renderbuffer data between CUDA and OpenGL requested? */
  if (cmd == 3 && !ispbo) {
    /* Copy from OpenGL object to CUDA buffer? */
    if (direction == 0) {
        /* OpenGL -> CUDA copy: */
        cudastatus = cudaMemcpyFromArrayAsync(gpuptr, (const struct cudaArray*) mappedArray, 0, 0, (size_t) nrbytes, cudaMemcpyDeviceToDevice, 0);
        if (cudastatus != cudaSuccess) {
            mexPrintf("\nmemcpyCudaOpenGL: ERROR in cudaMemcpyFromArray(): %s\n", cudaGetErrorString(cudastatus));
            goto err_unmap;
        }
    }

    if (direction == 1) {
        /* CUDA -> OpenGL copy: */
        cudastatus = cudaMemcpyToArrayAsync((struct cudaArray*) mappedArray, 0, 0, (const void*) gpuptr, (size_t) nrbytes, cudaMemcpyDeviceToDevice, 0);
        if (cudastatus != cudaSuccess) {
            mexPrintf("\nmemcpyCudaOpenGL: ERROR in cudaMemcpyToArray(): %s\n", cudaGetErrorString(cudastatus));
            goto err_unmap;
        }
    }
  }
  
  /* Return of pointers to mapped resource requested? */
  if (cmd == 4) {
    /* Yes: This implies we must not unmap the resource now, as otherwise the
     * returned pointers would be dead on arrival.
     */
    keepmapped = 1;
    
    /* Cast pointer to void* then store it in a 64-Bit unsigned integer return value: */
    plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    *((unsigned long long*) mxGetData(plhs[0])) = (unsigned long long) (void*) ((ispbo) ? mappedPtr : mappedArray);
  }
  
  /* Keep resource mapped? */
  if (slot && !keepmapped) doCudaUnmap(slot);
  
  /* Successfully completed: */
  return;
  
  /* Error handling -- Unwind in reverse order: */
  
err_unmap:
  
  /* Unmap the 'resource': */
  unmapResource(glhandle, gltarget);

err_unregister:
  
  /* Unregister the 'resource': */
  unregisterResource(glhandle, gltarget);

err_final:      
      
  if (cudastatus != cudaSuccess) mexErrMsgTxt("Error in memcpyCudaOpenGL(), reason see above.");
}

