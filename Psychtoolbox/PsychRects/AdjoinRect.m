function r=AdjoinRect(r,fixedRect,fixedSide)
% rect=AdjoinRect(rect,fixedRect,fixedSide)
%
% Moves rect to be outside and alongside the specified edge of fixedRect.
% The legal values for fixedSide are: RectLeft, RectRight, 
% RectTop, and RectBottom.
%
% rect and fixedRect can both be Mx4 or 4xM rect arrays, but either for one
% of them M must be 1, or rect and fixedRect must have the same shape. If
% rect contains multiple rects and fixedRect only one, each rect is moved
% w.r.t. fixedRect. If fixedRect contains multiple rects, but rect only
% one, rect is moved w.r.t each of the fixedRects. If both rect and
% fixedRect are arrays with the same number of rects, each rect is moved
% w.r.t. the corresponding fixedrect. fixedSide should always be scalar:
% the same operation will be executed on all rects in the array. 
%
% The output will have the same orientation as rect.
%
% See also PsychRects/Contents.

% Denis Pelli 5/26/96, 7/10/96, 8/6/96
% dcn 7/26/2015: vectorized, cleaned up
% dcn 7/29/2015: now handles 4xM and Mx4 inputs

if nargin~=3
	error('Usage: rect=AdjoinRect(rect,fixedRect,fixedSide)');
end
if (size(r,2)==4 && size(fixedRect,2)~=4) || (size(r,2)~=4 && size(fixedRect,2)==4)
	% orientation of first rect(-array) is leading for orientation of
	% output. Make sure second input has same orientation
    fixedRect = fixedRect.';
end
if ~isscalar(fixedSide)
	error('Wrong size fixedSide argument, must be scalar. Usage:  rect=AdjoinRect(rect,fixedRect,fixedSide)');
end

if (size(r, 1)==4) && (size(r,2)>=1)
    % 4xM
    switch fixedSide
        case 1
            r=OffsetRect(r,fixedRect(fixedSide,:)-r(3,:),0);
        case 3
            r=OffsetRect(r,fixedRect(fixedSide,:)-r(1,:),0);
        case 2
            r=OffsetRect(r,0,fixedRect(fixedSide,:)-r(4,:));
        case 4
            r=OffsetRect(r,0,fixedRect(fixedSide,:)-r(2,:));
        otherwise
            error('Illegal value %f of ''fixedSide''',fixedSide);
    end
elseif (size(r, 2)==4) && (size(r,1)>=1)
    % Mx4
    switch fixedSide
        case 1
            r=OffsetRect(r,fixedRect(:,fixedSide)-r(:,3),0);
        case 3
            r=OffsetRect(r,fixedRect(:,fixedSide)-r(:,1),0);
        case 2
            r=OffsetRect(r,0,fixedRect(:,fixedSide)-r(:,4));
        case 4
            r=OffsetRect(r,0,fixedRect(:,fixedSide)-r(:,2));
        otherwise
            error('Illegal value %f of ''fixedSide''',fixedSide);
    end
else
    % Something weird and unknown:
    error('Given matrix of rects not of required 4-by-n or n-by-4 format.');
end
