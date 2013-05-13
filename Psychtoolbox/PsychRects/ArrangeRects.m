function cellRects=ArrangeRects(n,objectRect,windowRect,rightToLeft);
% cellRects=ArrangeRects(n,objectRect,windowRect,[rightToLeft]);
%
% Returns an array of n contiguous rects that achieve a visually 
% pleasing arrangement of n objects of size objectRect within a 
% window of size windowRect. Each row of the returned matrix is a rect.
% If the rightToLeft argument is supplied and nonzero then
% the cells are ordered right to left (like Hebrew letters). 
% Otherwise they are ordered left to right (like English letters).
% The result depends on the proportions of the objectRect, but is 
% independent of its size and position.
% Also see PsychRects.

% 5/12/96 dgp Wrote it, based on my Alphabet.c in TextInNoise.
% 5/16/96 dgp Updated it to work with 1x4, instead of 4x1, rect.
% 5/27/96 dgp more explanation.
% 7/10/96 dgp PsychRects

% Explanation:
% Distribute n objects of size objectRect in a window of size windowRect.
% We're making a table. The row and column spacings are uniform but need not be equal.
% We want a regular grid with uniform horizontal and vertical spaces between objects.
% All rows are fully occupied except the last, which is centered, but still grid-aligned.
% Here's the rationale for the calculation. 
% ry=ratio of window to object heights.
% rx=ratio of window to object width.
% a=ry/rx
% We need to determine how many rows (ny) and columns (nx) we want.
% To show all the n objects we need nx*ny>=n, but we want as little excess as possible
% because excess cells will be blank, wasting space.
% To achieve similar gaps vertically and horizontally between objects
% we want the ratio ny/nx to be approximately equal to the ratio ry/rx=a.
% Multiplying together the two approximate relations gives us ny*ny=a*n, i.e. ny=sqrt(n*a).
% Dividing them gives us nx=sqrt(n/a).
% We then alternately find the smallest integers nx and ny with a product of at least than n.

if nargin<3 || nargin>4 
    error('Usage: cellRects=ArrangeRects(n,objectRect,windowRect,[rightToLeft]);'); 
end
if size(objectRect,2)~=4 || size(windowRect,2)~=4
	error('Wrong size rect argument. Usage: cellRects=ArrangeRects(n,objectRect,windowRect,[rightToLeft]);'); 
end
if nargin<4
	rightToLeft=0;
end
r=windowRect;
a=(RectHeight(r)/RectWidth(r))/(RectHeight(objectRect)/RectWidth(objectRect));
if(a>1)
	nx=ceil(sqrt(n/a));
	ny=ceil(n/nx);
	nx=ceil(n/ny);
else
	ny=ceil(sqrt(n*a));
	nx=ceil(n/ny);
end
ny=ceil(n/nx);
empties=nx*ny-n;
% create a cell for each object
cellHeight=RectHeight(windowRect)/ny;
cellWidth=RectWidth(windowRect)/nx;
xx=1:nx;
if(rightToLeft)
	xx=fliplr(xx);
end
i=1;
cellRects=zeros(n,4);
for y=1:ny;
	for x=xx;
		cellRect=SetRect(0,0,cellWidth,cellHeight);
		cellRect=OffsetRect(cellRect,windowRect(RectLeft),windowRect(RectTop));
		cellRect=OffsetRect(cellRect,(x-1)*cellWidth,(y-1)*cellHeight);
		if(y<ny || ( x-1>=floor(empties/2) && nx-x>=ceil(empties/2) ) )
			cellRects(i,:)=cellRect;
			i=i+1;
			if(i>n)break;end;
		end
	end
end

