function newRect = OffsetRect(oldRect,x,y)
% newRect = OffsetRect(oldRect,x,y)
%
% Offset the passed rect-array by the horizontal (x)
% and vertical (y) shift given. You can also pass in column vectors x
% and/or y with e.g., n different "shifts" and the function will return n
% copies of oldRect, each shifted by one of the shift values in x and y.
%
% Alternatively you can give a single scalar x,y shift value, but apply it
% to a whole matrix of rects 'oldRect' --> Apply a common offset to a large
% number of rects simultaneously.
%
% See also PsychRects/Contents.

% 5/16/96  dhb  Relented to Pelli's request to change calling order
%               from v,h to x,y.
%          dhb  Index using RectTop etc.
% 7/10/96  dgp  Wrote it.
% 8/5/96   dgp  Check rect size.
% 5/18/08  mk   Vectorized.
% 5/31/09  mk   Improved error handling. Add support for offsetting a
%               single rect multiple times, ie., pass x,y as vectors, then
%               create many offset versions of single input rect.
% 7/27/15  dcn  Made offsetting a single rect multiple times more flexible,
%               only x or y being a vector is allowed now. Simplified
% 7/29/15  dcn  Now correctly handles a 4x1 input

if nargin~=3
	error('Usage:  newRect = OffsetRect(oldRect,x,y)');
end

% single rect case: expand if any of the offsets are vectors
if PsychNumel(oldRect) == 4
    szs = [length(x) length(y)];
    if ~all(szs==1)
        if ~(szs(1)==1||szs(1)==szs(2)||szs(2)==1)
            error('x and y should be both one, one can be a vector or both should be vectors of the same length')
        end
        % x, y or both are one column vectors with size(x,1) rows/elements.
        % Replicate oldRect into nrpts identical copies, then add point
        % offsets:
        if size(oldRect, 1)==4
            newRect = repmat(oldRect, 1, max(szs));
        else
            newRect = repmat(oldRect, max(szs), 1);
        end
    end
end

% One or multiple rects:
if size(oldRect, 1)==4 && size(oldRect,2)>=1
    % ensure row vectors
    x = x(:).';
    y = y(:).';
    newRect(2, :) = oldRect(2, :) + y;
    newRect(4, :) = oldRect(4, :) + y;
    newRect(1, :) = oldRect(1, :) + x;
    newRect(3, :) = oldRect(3, :) + x;
elseif (size(oldRect, 2)==4) && size(oldRect,1)>=1
    % ensure column vectors
    x = x(:);
    y = y(:);
    newRect(:, 2) = oldRect(:, 2) + y;
    newRect(:, 4) = oldRect(:, 4) + y;
    newRect(:, 1) = oldRect(:, 1) + x;
    newRect(:, 3) = oldRect(:, 3) + x;
else
    % Something weird and unknown:
    error('Given matrix of rects not of required 4-by-n or n-by-4 format.');
end
