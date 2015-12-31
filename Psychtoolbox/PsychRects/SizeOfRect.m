function s=SizeOfRect(rect)
% s=SizeOfRect(rect);
% Accepts a Psychtoolbox rect [left, top, right, bottom] or Mx4 rect-array
% and returns the size [rows columns] of a MATLAB array (i.e. image) just
% big enough to hold all the pixels.
%
% See also PsychRects/Contents, RectOfMatrix, RectHeight, RectWidth, RectSize.

% 7/16/06 dgp Wrote it.
% 7/29/15 dcn Now handles 4xM and Mx4 inputs

if (size(rect, 1)==4) && (size(rect,2)>=1)
    % 4xM
    s=[RectHeight(rect); RectWidth(rect)];
elseif (size(rect, 2)==4) && (size(rect,1)>=1)
    % Mx4
    s=[RectHeight(rect) RectWidth(rect)];
else
    % Something weird and unknown:
    error('Given matrix of rects not of required 4-by-n or n-by-4 format.');
end
