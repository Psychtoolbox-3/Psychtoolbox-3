function [raster] = MakeGrid(nrow,ncol,xres,yres,plxres,plyres)
% [raster] = MakeGrid(nrow,ncol,xres,yres,plxres,plyres)
%
% Makes raster of elements centered on screen / in image (leftover space is
% divided equally over the edges)
%
% input:
%   ncol:   number of columns
%   nrow:   number of rows
%   xres:   horizontal resolution of screen or image
%   yres:   vertical resolution of screen or element
%   plxres: hor resolution of raster element (can be smaller than xres, not larger)
%   plyres: vert resolution of raster element (can be smaller than yres, not
%           larger)
%
% output
%   raster: struct with boundaries of raster elements
%           is a matrix of structs of size number of raster elements
%           vertically X number of raster element horizontally
%
% IH    2008
% DN    2008 v1.1 - coordinates sometimes faulty
% DN    2008 v1.2 - coordinates sometimes faulty (actually solved this
%                   time) + input checking

psychassert(plxres<=xres,'Raster element (plxres) bigger than image/screen (xres)');
psychassert(plyres<=yres,'Raster element (plyres) bigger than image/screen (yres)');

xresr = plxres*ncol;
yresr = plyres*nrow;

psychassert(xresr<=xres,'Too many raster elements horizontally (%d)\n%d pixels required, but only space for %d (xres)',ncol,xresr,xres);
psychassert(yresr<=yres,'Too many raster elements vertically (%d)\n%d pixels required, but only space for %d (yres)',nrow,yresr,yres);

for p=1:nrow,
    for q=1:ncol,
        raster(p,q).xmin = round((xres - xresr)/2 + (q-1)*plxres);
        raster(p,q).xmax = round((xres - xresr)/2 + (q)  *plxres);
        raster(p,q).ymin = round((yres - yresr)/2 + (p-1)*plyres);
        raster(p,q).ymax = round((yres - yresr)/2 + (p)  *plyres);
    end
end
