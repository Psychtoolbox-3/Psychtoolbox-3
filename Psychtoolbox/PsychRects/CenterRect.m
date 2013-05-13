function [rect,dh,dv] = CenterRect(rect,fixedRect)
% [rect,dh,dv] = CenterRect(rect,fixedRect)
% 
% Center the first rect in the second by adding an integer offset.
% Also see PsychRects.

% 5/16/96 dgp  Updated for new OffsetRect, and use symbolic rect indices.
% 7/10/96 dgp  PsychRects
% 8/5/96  dgp  Renamed the arguments, check rect size.
% 7/23/97 dgp  Round the offset.
% 5/4/00  dhb  Return dh and dv as well as centered rect.

if nargin~=2
	error('Usage:  rect=CenterRect(rect,fixedRect)');
end
if size(rect,2)~=4 || size(fixedRect,2)~=4
	error('Wrong size rect argument. Usage:  [rect,dh,dv] = CenterRect(rect,fixedRect)');
end
dv=(fixedRect(RectTop)+fixedRect(RectBottom)-rect(RectTop)-rect(RectBottom))/2;
dh=(fixedRect(RectLeft)+fixedRect(RectRight)-rect(RectLeft)-rect(RectRight))/2;
dv=round(dv);
dh=round(dh);
rect=OffsetRect(rect,dh,dv);


