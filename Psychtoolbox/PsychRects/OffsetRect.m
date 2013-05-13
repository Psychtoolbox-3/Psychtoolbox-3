function newRect = OffsetRect(oldRect,x,y)
% newRect = OffsetRect(oldRect,x,y)
%
% Offset the passed rect matrix by the horizontal (x)
% and vertical (y) shift given. You can also pass in column vectors x and y
% with e.g., n different "shifts" and the function will return n copies of
% oldRect, each shifted by one of the shift values in x and y.
%
% Alternatively you can give a single scalar x,y shift value, but apply it
% to a whole matrix of rects 'oldRect' --> Apply a common offset to a large
% number of rects simultaneously.
%
% Also see PsychRects.

% 5/16/96  dhb  Relented to Pelli's request to change calling order
%               from v,h to x,y.
%          dhb  Index using RectTop etc.
% 7/10/96  dgp  Wrote it.
% 8/5/96   dgp  Check rect size.
% 5/18/08  mk   Vectorized.
% 5/31/09  mk   Improved error handling. Add support for offsetting a
%               single rect multiple times, ie., pass x,y as vectors, then
%               create many offset versions of single input rect.

if nargin~=3
	error('Usage:  newRect = OffsetRect(oldRect,x,y)');
end

if PsychNumel(oldRect) == 4
    % Single rect case:
    if size(oldRect,2)~=4
        error('Wrong size rect argument. A single rect must be a 1-row 4 element vector!');
    end
    
    if length(x)==1 && length(y)==1
        % Single rect, offset by a single (x,y) point:
        newRect(RectTop) = oldRect(RectTop) + y;
        newRect(RectBottom) = oldRect(RectBottom) + y;
        newRect(RectLeft) = oldRect(RectLeft) + x;
        newRect(RectRight) = oldRect(RectRight) + x;
    else
        % Single rect, but multiple points:
        if ~all(size(x) == size(y)) || size(x,2)~=1
            error('Wrong format of x or y in multipoint case: x and y must be 1-column vectors of matching size!');
        else
            % x and y are one column vectors with size(x,1) rows/elements.
            % Replicate oldRect into nrpts identical copies, then add point
            % offsets:
            newRect = repmat(oldRect, size(x, 1), 1) + [x, y, x, y];
        end
    end
else
    % Multirect case: Offsets (x,y) can be vectors, or a single point:    
    if size(oldRect, 1)==4
        newRect(RectTop, :) = oldRect(RectTop, :) + y;
        newRect(RectBottom, :) = oldRect(RectBottom, :) + y;
        newRect(RectLeft, :) = oldRect(RectLeft, :) + x;
        newRect(RectRight, :) = oldRect(RectRight, :) + x;
    else
        newRect(:, RectTop) = oldRect(:, RectTop) + y;
        newRect(:, RectBottom) = oldRect(:, RectBottom) + y;
        newRect(:, RectLeft) = oldRect(:, RectLeft) + x;
        newRect(:, RectRight) = oldRect(:, RectRight) + x;
    end
end
