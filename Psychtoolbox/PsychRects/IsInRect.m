function inside = IsInRect(x,y,rect)
% inside = IsInRect(x,y,rect)
%
% Is location x,y inside the rect?
% x and/or y can also be column vectors, and rect can be an Mx4 rect-array.
% This means you can test one point against mutliple rects, multiple points
% against one rect or multiple points each against their own rect
%
% See also PsychRects/Contents.

% 3/5/97  dhb  Wrote it.
% 7/26/15 dcn  Vectorized
% 7/29/15 dcn  Now handles 4xM and Mx4 inputs

if (size(rect, 1)==4) && (size(rect,2)>=1)
    % 4xM
    x = x(:).';
    y = y(:).';
    inside = x >= rect(1,:) & x <= rect(3,:) & ...
             y >= rect(2,:) & y <= rect(4,:);
elseif (size(rect, 2)==4) && (size(rect,1)>=1)
    % Mx4
    x = x(:);
    y = y(:);
    inside = x >= rect(:,1) & x <= rect(:,3) & ...
             y >= rect(:,2) & y <= rect(:,4);
else
    % Something weird and unknown:
    error('Given matrix of rects not of required 4-by-n or n-by-4 format.');
end