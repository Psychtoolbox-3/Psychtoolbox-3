function SaveAsEps(filename,m,pageRect,resolution)
% SaveAsEps(filename,m,[pageRect],[resolution]);
% Based on WindowToEps.c of the VideoToolbox.
% Copyright © 1996-2003 Denis G. Pelli
% 
% SaveAsEps converts a grayscale image (in a Matlab matrix) to a
% PostScript file that a LaserWriter or Linotype will accurately render on
% paper. The output file is a standard Encapsulated PostScript File, with
% file type 'EPSF', commonly referred to as an "eps" file. Most word
% processors, e.g. Word and PageMaker, know how to import eps files.
% Macintoshes only understand QuickDraw, not PostScript, so they don't
% know what image the PostScript code would produce. The Macintosh file
% "creator" of the eps file is set to Microsoft Word, so double-clicking
% it will open it as a Word document containg the image.
% 
% I suggest you use the $25 shareware Macintosh application EPS¥Factory 
% to add a PICT preview to your file, based on the postscript already 
% there. You can download EPS¥Factory from 
% web http://www.artage.com/pages/products/products.html
% (The old freeware ps2eps is incompatible with Mac OS 8.6 and 9, alas.)
% The PICT preview will make the image will look approximately right on 
% your monitor, e.g. in Word. 
% 
% If you want to convert the EPS file to some other format, we recommend
% using Adobe Illustrator, especially the "Export to Web" feature.
%
% To import into Word, use the Insert:Picture or
% Insert:File command. Once imported, you can get Word to rescale the
% image by dragging the image's lower right corner while holding down the
% shift key. Try the VideoToolbox demo Grating. Bear in mind that
% on-screen you're looking at the PICT, whereas, when you print on a
% LaserWriter, the image is produced by the PostScript code. When in your
% word processor, the on-screen images will look best if you use the
% Monitors control panel to set your monitor to 256 Grays, but this won't
% affect printing.
% 
% The filename, by convention, should end in '.eps' to indicate that it's
% an encapsulated postscript file, but this is not enforced. The file's
% type is set to 'EPSF' with creator 'R*ch'. (This creator corresponds to
% BBEdit, so double-clicking will open it as a text file in BBEdit, which
% includes a Special:SendPostScript command for downloading poscript
% images to the laserwriter. You can change the creator to be anything you
% want.)
% 
% 	The image matrix must have 8 bits per pixel (i.e. values 0 to 255). No
% color tables are used. The raw pixel value (Apple calls it an "index")
% is sent directly to the printer with no transformation. PostScript
% assumes that the number, from 0 to 255 is proportional to desired
% reflectance, from zero to 1.
%  	The pageRect is subtle. It describes, in typographers points (1/72"),
% the rectangle that your image will be mapped onto on the printed page.
% It is essential that you keep in mind that Apple and Adobe use different
% coordinate systems. Both Apple and Adobe increase x from left to right.
% However, Apple has y increasing from top to bottom, whereas Adobe
% increases y from bottom to top. Adobe's origin is the lower left corner
% of the page, even though that point is usually not printable, since most
% printers can only print to within about a half inch of the edge. The
% pageRect, though supplied in Apple's Rect data structure, must be in
% Adobe's coordinates, respecting the names of the Rect structure's
% fields: left, top, right, bottom. So, for an image to fill most of an
% 8.5x11 page, with 0.5" margins, you might use the following:
% 
% 		pageRect=SetRect(0.5,10.5,8,0.5)*72;
% 
% In printing PostScript halftones the halftone cell size determines both
% the spatial and graylevel resolutions of the resulting image. If the
% "resolution" parameter is positive, then it specifies the cellsPerInch;
% if it's negative then it specifies the number of grayLevels. If
% "resolution" is omitted then the printer will be left at its default
% cell size, which is usually a good choice. Note that there need not be
% any particular correspondence between pixels in your image and cells in
% the halftone; the printer automatically resamples your image to produce
% the halftone.
% 	If you specify cellsPerInch (resolution>0) then the printer will be
% asked to print its halftone with that many halftone cells per inch. E.g.
% to produce a halftone original for subsequent one-to-one reproduction in
% a journal, you'll want the cells to be coarse enough for them to
% reproduce without re-screening, e.g. 100 cells per inch.
% 	Alternatively, if you specify grayLevels (resolution<0) then the printer
% will be asked to print its halftone with cells containing grayLevels-1
% printer pixels, yielding the specified number of gray levels. E.g. you
% might want to force your 300 dpi LaserWriter to use big cells yielding
% 256 gray levels.
%  	If resolution==0 or resolution==[], it will be ignored.
% 	Here's a minimal example that prints an image to disk, preserving the
% size and scale of the image,
% 
% 		SaveAsEps('test.eps',m);
% 
% The default pageRect is equivalent to this:
% 
% 		pageRect=SetRect(0,-0,size(m,2),-size(m,1));% Adobe coordinates
% 		pageRect=CenterRect(pageRect,SetRect(0,11,8.5,0)*72);
% 		SaveAsEps('test.eps',m,pageRect);
% 
% Tiling is something we often want to do, creating a huge image by taping
% many pages together. You accomplish this by repeatedly printing the huge
% image--only one pageful appears each time--shifting the image so that
% eventually every bit has been printed. (It's slow, since the whole image
% is transmitted to the printer each time.) Here's an example that creates
% a big "width" by "height" image. We tile onto multiple 8.5"x11" pages,
% using 7.5"x10" of each page, allowing for 0.5" nonprinting margins:
% 
% 		delete 'test.eps';
% 		mosaicRect=SetRect(0,height,width,0);	% In units of "points". 72 points/inch.
% 		mosaicRect=OffsetRect(mosaicRect,0.5*72,0.5*72);	% allow for nonprinting margin
% 		for i=0:7.5*72:width
% 			for j=0:10*72:height
% 				pageRect=mosaicRect;
% 				OffsetRect(&pageRect,-i,-j);
% 				ImageToPs('test.ps',m,pageRect);
% 			end
% 		end
% 
% In case your PostScript manual isn't handy, the easy way to obtain
% multiple copies of each page is to use the #copies variable that's built
% into Postscript. Anywhere in your file before "showpage", set #copies to
% the desired value: '/#copies 3 def\n'
% 
% There are many free PostScript downloaders, but I find them clumsy to
% use, especially the ones that don't inform you of errors, should they
% occur. I like these:
% 
% BBEdit (used to come with an extension for downloading postscript)
% 
% (Don't know if this is still available.)
% LaserStatus, a desk accessory included in the 
% MockPackage Plus Utilities
% from:
% CE Software
% 1854 Fuller Road
% PO Box 65580
% West Des Moines, Iowa 50265
% (515)-224-1995
% 
% You may also want to read:
% 
% Adobe Systems (1985) PostScript Language Reference Manual, Second
% Edition. Reading, MA: Addison-Wesley.
% 
% Pelli, D. G. (1987) Programming in PostScript: Imaging on paper from a
% mathematical description. BYTE, 12 (5), 185-202.
% 
% BUGS:
% The "resolution" parameter is ignored by my Apple LaserWriter 16/600.
% I've looked at the PostScript code in the .eps file and it looks
% perfect. For reasons that I don't understand, the printer is ignoring
% the setscreen command, leaving you at the default printer resolution.

% HISTORY:
% 4/21/91 dgp	wrote original PixmapToPostScript.c
% 1995 dgp	wrote WindowToEPS.c based on above.
% 7/14/96 dgp	wrote SaveAsEps.m, based on above
% 7/17/96 dgp	work around Matlab "feature" of logical Boolean indexing
% 7/18/96 dgp	set Mac file CREATOR to Microsoft Word.
% 5/8/97 dgp	make sure pixel values are real (not complex) and in bounds.
% 9/7/98 dgp	convert image (which might be uint8) to double.
% 3/17/99 dgp	mention ps2eps.
% 1/23/00 dgp	ps2eps is incompatible with Mac OS 8.6 and 9; recommend
% 				EPS¥Factory instead.

if nargin<2
	error('Usage: SaveAsEps(filename,m,[pageRect,resolution])');
end
if ~ischar(filename)
	error('Usage: SaveAsEps(filename,m,[pageRect,resolution])');
end
m=round(double(m));
if ~all(all(m<=255 & m>=0))
	fprintf('Image pixels must be in range 0 to 255, not %d to %d\n',min(min(m)),max(max(m)));
	error('Usage: SaveAsEps(filename,m,[pageRect,resolution])');
end	
rect=SetRect(0,0,size(m,2),size(m,1));
if nargin<3
	pageRect=[];
end
if isempty(pageRect)
	pageRect=ScaleRect(rect,1,-1);% convert Apple to Adobe coordinates
	pageRect=CenterRect(pageRect,SetRect(0,11*72,8.5*72,0));
elseif length(pageRect)~=4
	error('pageRect argument should have 1 row and 4 columns, or none []');
end
if nargin<4
	resolution=[];
end
[file,message]=fopen(filename,'wt');
if file==-1
	error(sprintf('SaveAsEps: Error in opening file ''%s''. %s\n',filename,message))
end
t=clock; % [year month day hour minute seconds]
string=sprintf('%s %02d:%02d',date,t(4),t(5));

fprintf(file, ...
	'%%!PS-Adobe-2.0 EPSF-1.2\n');
fprintf(file, ...
	'%%%%Creator:%s\n','SaveAsEps.m, Copyright 1996 Denis G. Pelli');
fprintf(file, ...
	'%%%%For:%s\n','Matlab');
fprintf(file, ...
	'%%%%Title:%s\n',filename);
fprintf(file, ...
	'%%%%CreationDate:%s\n',string);
fprintf(file, ...
	['%%%%Pages: 0\n' ...
	'%%%%BoundingBox:%g %g %g %g\n'] ...
	,pageRect(RectLeft),pageRect(RectBottom) ...
	,pageRect(RectRight),pageRect(RectTop));
fprintf(file, ...
	['%%%%EndComments\n' ...
	'%%%%BeginProlog\n' ...
	'%%%%EndProlog\n' ...
	'%%%%BeginSetup\n' ...
	'/VideoToolbox dup 25 dict def load begin\n' ...
	'end\n' ...
	'VideoToolbox begin\n' ...
	'gsave\n' ...
	'%%%%EndSetup\n']);
fprintf(file, ...
	'/nx %g def				%% pixels per raster line\n',RectWidth(rect));
fprintf(file, ...
	'/ny %g def				%% lines in image\n',RectHeight(rect));
fprintf(file, ...
	['/hypotenuse {dup mul exch dup mul add sqrt} bind def\n' ...
	'/pixelsPerInch gsave initmatrix 72 0 dtransform hypotenuse grestore def\n' ...
	'%g %g translate		%% locate lower left of image\n'] ...
	,pageRect(RectLeft),pageRect(RectBottom));
fprintf(file, ...
	'%g %g scale			%% print image with these dimensions on page\n' ...
	,RectWidth(pageRect),-RectHeight(pageRect));
if resolution>0
	fprintf(file, ...
		'/cellsPerInch %.4f def	%% halftone dot frequency\n',resolution);
	fprintf(file, ...
		'cellsPerInch currentscreen 4 -2 roll pop 3 1 roll setscreen\n');
end
if resolution<0	% -graylevels
	fprintf(file, ...
		'/cellsPerInch pixelsPerInch %.4f div def	%% halftone dot frequency\n' ...
		,sqrt(-resolution-1));
	fprintf(file, ...
		'cellsPerInch currentscreen 4 -2 roll pop 3 1 roll setscreen\n');
end
fprintf(file, ...
	['/s nx string def		%% string to hold one raster line\n' ...
	'nx ny 8					%% dimensions and bits/pixel of source image\n' ...
	'[nx 0 0 ny neg 0 ny]	%% map unit square to PixMap data\n' ...
	'{currentfile s readhexstring pop}	%% read data\n' ...
	'bind					%% speed up reading\n' ...
	'image\n']);
hex=zeros(2,256);
hex(:)=sprintf('%02x',0:255);
if length(hex)==RectWidth(rect)
	% work around Matlab "feature" of logical Boolean indexing when 
	% subscript array is of same size and all 0's and 1's.
	hex=[hex [0;0]];
end
% one row of image per line in postscript file
clear pix;
for i=1+rect(RectTop):rect(RectBottom)
	pix=m(i,1+rect(RectLeft):rect(RectRight));
	fprintf(file,'%s\n',hex(:,1+pix));
end
clear hex buffer
fprintf(file, ...
	['%%%%Trailer\n' ...
	'grestore\n' ...
	'end\n' ...
	'showpage\n' ...
	'%%%%EOF\n']);
fclose(file);
