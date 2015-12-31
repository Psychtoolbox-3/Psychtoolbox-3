function [rect,dh,dv] = CenterRectd(rect,fixedRect)
% [rect,dh,dv] = CenterRectd(rect,fixedRect)
% 
% Center the first rect in the second. It does the same as CenterRect but
% without rounding to integer pixel position for the rect corners.
%
% rect and fixedRect can both be Mx4 or 4xM rect arrays, but either for one
% of them M must be 1, or rect and fixedRect must have the same shape. If
% rect contains multiple rects and fixedRect only one, each rect is
% centered in fixedRect. If fixedRect contains multiple rects, but rect
% only one, rect is centered in each of the fixedRects. If both rect and
% fixedRect are arrays with the same number of rects, each rect is centered
% in the corresponding fixedrect.
%
% The output will have the same orientation as rect.
%
% See also PsychRects/Contents.

% 5/16/96 dgp  Updated for new OffsetRect, and use symbolic rect indices.
% 7/10/96 dgp  PsychRects
% 8/5/96  dgp  Renamed the arguments, check rect size.
% 7/23/97 dgp  Round the offset.
% 5/4/00  dhb  Return dh and dv as well as centered rect.
% 7/26/15 dcn  Vectorized
% 7/27/15 dcn  Trivially derived from CenterRect
% 7/29/15 dcn  Now handles 4xM and Mx4 inputs

if nargin~=2
	error('Usage:  rect=CenterRectd(rect,fixedRect)');
end
if (size(rect,2)==4 && size(fixedRect,2)~=4) || (size(rect,2)~=4 && size(fixedRect,2)==4)
	% orientation of first rect(-array) is leading for orientation of
	% output. Make sure second input has same orientation
    fixedRect = fixedRect.';
end
if (size(rect, 1)==4) && (size(rect,2)>=1)
    % 4xM
    dv=(fixedRect(2,:)+fixedRect(4,:)-rect(2,:)-rect(4,:))/2;
    dh=(fixedRect(1,:)+fixedRect(3,:)-rect(1,:)-rect(3,:))/2;
elseif (size(rect, 2)==4) && (size(rect,1)>=1)
    % Mx4
    dv=(fixedRect(:,2)+fixedRect(:,4)-rect(:,2)-rect(:,4))/2;
    dh=(fixedRect(:,1)+fixedRect(:,3)-rect(:,1)-rect(:,3))/2;
else
    % Something weird and unknown:
    error('Given matrix of rects not of required 4-by-n or n-by-4 format.');
end
rect=OffsetRect(rect,dh,dv);


