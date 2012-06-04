function framebufferTrue = MapIndexColorThroughClut(framebufferIndex,clut)
% framebufferTrue = MapIndexColorThroughClut(framebufferIndex,clut)
%
% Take an index color frame buffer and a clut, produce a true color
% framebuffer.  This is based on integers, not the [0-1] clut model of
% OpenGL and PTB-3.
%
% 3/22/05		dhb		Wrote it.

[nRows,nCols] = size(framebufferIndex);
nColors = size(clut,2);
framebufferTrue = clut(double(framebufferIndex)+1,:);
framebufferTrue = reshape(framebufferTrue,nRows,nCols,nColors);

