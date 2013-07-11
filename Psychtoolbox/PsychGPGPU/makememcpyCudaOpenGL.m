function makememcpyCudaOpenGL
% makememcpyCudaOpenGL - Build CUDA<->OpenGL interop helper mex file.
%
% This needs the NVidia CUDA-5.0 SDK or later installed.
%

if exist(['memcpyCudaOpenGL.' mexext], 'file');
    delete(['memcpyCudaOpenGL.' mexext]);
end

if IsOSX
	% -lcuda is optional at the moment. -framework CUDA is optional as well.
	mex -v -g -largeArrayDims memcpyCudaOpenGL.cpp -I/usr/include/ -I"/Developer/NVIDIA/CUDA-5.0/include" CXXFLAGS="\$CXXFLAGS -F/System/Library/Frameworks" -L"/Developer/NVIDIA/CUDA-5.0/lib" -L/usr/local/cuda/lib LDFLAGS="\$LDFLAGS -Xlinker -rpath -Xlinker /usr/local/cuda/lib " -lcudart
end

if IsLinux(1)
	% 64-Bit Linux build needs lib64 folder instead of lib folder:
	mex -v -g -largeArrayDims memcpyCudaOpenGL.cpp -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcudart
end

if IsLinux(0)
	% 32-Bit Linux build needs lib folder:
	mex -v -g -largeArrayDims memcpyCudaOpenGL.cpp -I/usr/local/cuda/include -L/usr/local/cuda/lib -lcudart
end

if IsWin(1)
	% 64-Bit Windows build needs x64 folder instead of Win32 folder:
	mex -v -largeArrayDims memcpyCudaOpenGL.cpp -I"C:\Program Files\Microsoft SDKs\Windows\v7.1\Include" -I"T:\projects\CUDA\v5.0\include" -LT:\projects\CUDA\v5.0\lib\x64 -lcudart
end

if IsWin(0)
	% 32-Bit Windows build needs Win32 folder:
	mex -v -largeArrayDims memcpyCudaOpenGL.cpp -I"T:\projects\CUDA\v5.0\include" -LT:\projects\CUDA\v5.0\lib\Win32 -lcudart
end

delete *.o
