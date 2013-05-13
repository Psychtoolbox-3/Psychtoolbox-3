function retpixels = glReadPixels( x, y, width, height, format, type, bufferoffset )

% glReadPixels  Interface to OpenGL function glReadPixels
%
% usage:  
% For standard readback into a Matlab or Octave image matrix:
% retpixels = glReadPixels( x, y, width, height, format, type )
%
% For readback into the currently bound OpenGL Pixelbuffer object (PBO):
% glReadPixels( x, y, width, height, format, type, bufferoffset )
% where bufferoffset is the positive integer byte-offset into the PBO.
%
% C function:  void glReadPixels(GLint x, GLint y, GLsizei width, GLsizei height, GLenum format, GLenum type, GLvoid* retpixels)

% 24-Jan-2006 -- created (generated automatically from header files)
% 25-Jan-2006 -- return argument allocated for case format=GL_RGB,
%                type=GL_UNSIGNED_BYTE (RFM)
% 13-Jun-2009 -- Remove special code for Octave -> No longer needed. (MK)

% ---allocate---
% ---protected---

if nargout > 0
    % Standard path: Readback pixels to host memory, return image in return
    % argument 'retpixels':

    if nargin~=6,
        error('Invalid number of arguments: Must be (x, y, width, height, format, type) for readback to host-memory.');
    end
    
    bufferoffset = -1;
else
    % PBO backed path: Readback into a Pixelbuffer object. Nothing is
    % directly returned to Matlab, but instead the read pixels are written
    % to the currently bound pixelbuffer object at the offset provided by
    % 'bufferoffset'
    if nargin~=7,
        error('Invalid number of arguments: Must be (x, y, width, height, format, type, bufferoffset) for readback into PBO.');
    end
    
    if isempty(bufferoffset)
        error('You must provide a valid "bufferoffset" into the PBO!');
    end
    
    if bufferoffset < 0
        error('You must provide a valid, zero or positive integral "bufferoffset" into the PBO!');
    end
end    
 
% check format and type
global GL

if x < 0 || y<0
    error('Invalid (negative) (x,y) offset passed to glReadPixels.');
end;

if width <= 0 || height<=0
    error('Invalid (negative or zero) (width, height) passed to glReadPixels.');
end;

switch(format)
    case {GL.RED, GL.GREEN, GL.BLUE, GL.ALPHA, GL.LUMINANCE, GL.INTENSITY, GL.DEPTH_COMPONENT, GL.STENCIL_INDEX }
        numperpixel = 1;
    case {GL.LUMINANCE_ALPHA }
        numperpixel = 2;
    case {GL.RGB , GL.BGR }
        numperpixel = 3;
    case {GL.RGBA , GL.BGRA }
        numperpixel = 4;
    otherwise
        error('Invalid format passed to glReadPixels.');
end;    

if bufferoffset == -1
    % Readback to host memory, aka Matlab- or Octave- matrix:
    % Allocate memory:
    pixels=zeros(numperpixel, width, height);

    % Tell OpenGL that we accept byte-aligned aka unaligned data.
    glPixelStorei(GL.PACK_ALIGNMENT, 1);

    % Perform proper type-cast:
    switch(type)
        case GL.UNSIGNED_BYTE
            pixels = uint8(pixels);
            pclass = 'uint8';
        case GL.BYTE
            pixels = int8(pixels);
            pclass = 'int8';
        case GL.UNSIGNED_SHORT
            pixels = uint16(pixels);
            pclass = 'uint16';
        case GL.SHORT
            pixels = int16(pixels);
            pclass = 'int16';
        case GL.UNSIGNED_INT
            pixels = uint32(pixels);
            pclass = 'uint32';
        case GL.INT
            pixels = int32(pixels);
            pclass = 'int32';
        case GL.FLOAT
            pixels = moglsingle(pixels);
            pclass = 'double';
        otherwise
            error('Invalid type argument passed to glReadPixels()!');
    end;

    % Execute actual call:
    moglcore( 'glReadPixels', x, y, width, height, format, type, pixels );

    % Rearrange data in Matlab friendly format:
    retpixels = zeros(size(pixels,2), size(pixels,3), size(pixels,1), pclass);
    for i=1:numperpixel
        retpixels(:,:,i) = pixels(i,:,:);
    end;
else
    % Readback into bound Pixelbuffer object PBO:
    moglcore( 'glReadPixels', x, y, width, height, format, type, bufferoffset );
    retpixels = [];
end

return


% ---autocode---
%
% function pixels = glReadPixels( x, y, width, height, format, type )
% 
% % glReadPixels  Interface to OpenGL function glReadPixels
% %
% % usage:  pixels = glReadPixels( x, y, width, height, format, type )
% %
% % C function:  void glReadPixels(GLint x, GLint y, GLsizei width, GLsizei height, GLenum format, GLenum type, GLvoid* pixels)
% 
% % 05-Mar-2006 -- created (generated automatically from header files)
% 
% % ---allocate---
% 
% if nargin~=6,
%     error('invalid number of arguments');
% end
% 
% pixels = (0);
% 
% moglcore( 'glReadPixels', x, y, width, height, format, type, pixels );
% 
% return
%
