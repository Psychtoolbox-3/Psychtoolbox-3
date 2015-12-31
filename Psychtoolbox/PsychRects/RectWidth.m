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
% 7/29/15 dcn Now handles 4xM and Mx4 inputs

if nargin~=1
	error('Usage:  width = RectWidth(rect)');
end

if (size(rect, 1)==4) && (size(rect,2)>=1)
    % 4xM
    width = rect(3,:) - rect(1,:);
elseif (size(rect, 2)==4) && (size(rect,1)>=1)
    % Mx4
    width = rect(:,3) - rect(:,1);
else
    % Something weird and unknown:
    error('Given matrix of rects not of required 4-by-n or n-by-4 format.');
end
