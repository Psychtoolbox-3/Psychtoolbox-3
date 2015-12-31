function height = RectHeight(rect)
% height = RectHeight(rect)
%
% Returns the rect's height, or the height of each of the rects in the
% rect-array.
%
% See also PsychRects/Contents.

% 5/12/96 dgp wrote it.
% 7/10/96 dgp PsychRects
% 8/5/96 dgp check rect size.
% 7/27/15 dcn Vectorized
% 7/29/15 dcn Now handles 4xM and Mx4 inputs

if nargin~=1
	error('Usage:  height = RectHeight(rect)');
end

if (size(rect, 1)==4) && (size(rect,2)>=1)
    % 4xM
    height = rect(4,:)-rect(2,:);
elseif (size(rect, 2)==4) && (size(rect,1)>=1)
    % Mx4
    height = rect(:,4)-rect(:,2);
else
    % Something weird and unknown:
    error('Given matrix of rects not of required 4-by-n or n-by-4 format.');
end
