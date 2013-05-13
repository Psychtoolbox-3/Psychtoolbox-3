function r=AdjoinRect(r,fixedRect,fixedSide)
% rect=AdjoinRect(rect,fixedRect,fixedSide)
%
% Moves rect to be outside and alongside the specified edge of fixedRect.
% The legal values for fixedSide are: RectLeft, RectRight, 
% RectTop, and RectBottom.
% Also see PsychRects.

% Denis Pelli 5/26/96, 7/10/96, 8/6/96

if nargin~=3
	error('Usage: rect=AdjoinRect(rect,fixedRect,fixedSide)');
end
if size(r,2)~=4 || size(fixedRect,2)~=4
	error('Wrong size rect argument. Usage:  rect=AdjoinRect(rect,fixedRect,fixedSide)');
end

if fixedSide==RectLeft
	r=OffsetRect(r,fixedRect(fixedSide)-r(RectRight),0);
	return
end
if fixedSide==RectRight
	r=OffsetRect(r,fixedRect(fixedSide)-r(RectLeft),0);
	return
end
if fixedSide==RectTop
	r=OffsetRect(r,0,fixedRect(fixedSide)-r(RectBottom));
	return
end
if fixedSide==RectBottom
	r=OffsetRect(r,0,fixedRect(fixedSide)-r(RectTop));
	return
end
error(sprintf('Illegal value %f of ''fixedSide''',fixedSide));
