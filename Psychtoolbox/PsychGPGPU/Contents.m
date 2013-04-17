% Psychtoolbox:PsychGPGPU
%
% General Purpose Graphics Processing Unit (GPGPU) computing suppport.
%
% This folder contains functions which employ massively parallel computing
% hardware, typically modern GPU's, to accelerate computing tasks in the
% context of stimulus presentation. The compute hardware is controlled via
% computing api's like NVidia CUDA or the cross-platform, cross-vendor api
% OpenCL.
%
% The current initial implementation as of April 2013 only supports CUDA
% compute capable NVidia GPUs in combination with the free and open-source
% GPUmat toolkit for 64-Bit Matlab. Support for GNU/Octave and other GPU
% and compute toolkits, e.g., OpenCL based systems, will follow in later
% releases. Currently only 64-Bit Matlab is supported. While it would be
% possible to also provide 32-Bit support, we don't expect much demand for
% it, so will not deliver 32-Bit support for the time being.
%
%
% Requirements:
%
% 1. A NVidia GPU with CUDA compute capabilities. This is any OpenGL3/4
%    capable GPU, starting with the GeForce-8000 series or later.
%
% 2. The freely downloadable CUDA SDK Version 5.0 or later from NVidia:
%    https://developer.nvidia.com/cuda-downloads
%
% 3. The free and open-source GPUmat toolbox for Matlab from SourceForge:
%    http://sourceforge.net/projects/gpumat
%
% 4. A 64-Bit version of Matlab for Linux, Windows or MacOSX.
%
% Demos can be found in the Psychtoolbox/PsychDemos/GPGPUDemos folder.
%
% This folder contains support functions for GPU computing:
%
% GPUTypeFromToGL        - Helper function to transfer data efficiently
%                          between CUDA/GPUmat data types and
%                          OpenGL/Psychtoolbox.
%
% memcpyCudaOpenGL.cpp   - Source code for the memcpyCudaOpenGL mex files.
% makememcpyCudaOpenGL.m - Build script to build mex files against Matlab
%                          and CUDA 5.0 SDK.
% memcpyCudaOpenGL.mex*  - Mex files for 64-Bit Matlab.
%
