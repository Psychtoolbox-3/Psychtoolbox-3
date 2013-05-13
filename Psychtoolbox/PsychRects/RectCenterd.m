function [x, y] = RectCenterd(r)
%   [x,y] = RectCenterd(rect);
%
%	RectCenterd returns the x,y point at the center of a rect. It does the
%	same as RectCenter, just without rounding the result to integer
%	coordinates.
%
%	See also PsychRects, CenterRectOnPoint.

%	9/13/99	Allen Ingling wrote it.
%	10/6/99	dgp Fixed bug.
%   5/18/08 mk Vectorized.
%   5/18/08 mk Non-rounding version, trivially derived from RectCenter.

if nargout~=2
	error('Usage: [x, y] = RectCenterd(rect);');
end

if PsychNumel(r) == 4
    % Single rect:
    x = (0.5*(r(1)+r(3)));
    y = (0.5*(r(2)+r(4)));
else
    % Multi-rect array:
    if (size(r, 1)==4) && (size(r,2)>=1)
        % Multi-column array with one 4-comp. rect per column:
        x = (0.5*(r(1,:)+r(3,:)));
        y = (0.5*(r(2,:)+r(4,:)));
    else
        if (size(r, 2)==4) && (size(r,1)>=1)
            % Multi-row array with one 4-comp. rect per row:
            x = (0.5*(r(:,1)+r(:,3)));
            y = (0.5*(r(:,2)+r(:,4)));
        else
            % Something weird and unknown:
            error('Given matrix of rects not of required 4-by-n or n-by-4 format.');
        end
    end
end
