function width = RectWidth(rect)
% width = RectWidth(rect)
%
% Returns the rect's width, or the width of each of the rects in the
% rect-array.
%
% See also PsychRects/Contents.

% 5/12/96 dgp wrote it.
% 7/10/96 dgp PsychRects
% 7/27/15 dcn Vectorized

if nargin~=1
	error('Usage:  width = RectWidth(rect)');
end
if size(rect,2)~=4
	error('Wrong size rect argument. Usage:  width = RectWidth(rect)');
end
width = rect(:,3) - rect(:,1);
