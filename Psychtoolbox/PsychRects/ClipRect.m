function newRect = ClipRect(a,b)
% newRect = ClipRect(a,b)
% 
% Returns the rect that is the intersection of the two rects a and b.
% Returns an empty rect [0 0 0 0] if the two rects are apart.
% Also see PsychRects.

% 7/5/96  dgp  Wrote it.
% 7/6/96  dgp  Return empty matrix [] only if apart.
% 7/10/96 dgp  Return empty rect [0 0 0 0] if apart.

if nargin~=2
	error('Usage:  rect=ClipRect(a,b)');
end
if size(a,2)~=4 || size(b,2)~=4
	error('Wrong size rect argument. Usage:  rect=ClipRect(a,b)');
end
newRect=a;
newRect(RectTop)=max(a(RectTop),b(RectTop));
newRect(RectBottom)=min(a(RectBottom),b(RectBottom));
newRect(RectLeft)=max(a(RectLeft),b(RectLeft));
newRect(RectRight)=min(a(RectRight),b(RectRight));
if RectWidth(newRect)<0 || RectHeight(newRect)<0
	newRect=[0 0 0 0];
end
