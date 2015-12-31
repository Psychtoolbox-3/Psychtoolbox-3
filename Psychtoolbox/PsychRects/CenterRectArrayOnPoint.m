function rects = CenterRectArrayOnPoint(rects, x, y)
% rects = CenterRectArrayOnPoint(rects, x, y)
%
% Centers a group of rectangles 'rects' on a given
% coordinate (x,y). The center of the bounding box
% of all given 'rects' is placed onto position (x,y)
% by shifting the whole set of 'rects' accordingly.
%
% See also PsychRects/Contents.

% History:
% 19-Apr-2015   Written by Natalia Zaretskaya.
%

if nargin ~= 3
  error('Usage: CenterRectArrayOnPoint(rects, x, y)');
end

boundingRect = [min(rects(:,1)) min(rects(:,2)) max(rects(:,3)) max(rects(:,4))];
[bRectx, bRecty] = RectCenter(boundingRect);

offsetx = x - bRectx;
offsety = y - bRecty;

rects = OffsetRect(rects, offsetx, offsety);

return;
