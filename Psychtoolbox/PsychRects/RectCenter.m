function [x, y] = RectCenter(r)
%   [x,y] = RectCenter(rect);
%
%	RectCenter returns the integer x,y point closest to the center of a rect.  
%
%	See also PsychRects, CenterRectOnPoint.

%	9/13/99	Allen Ingling wrote it.
%	10/6/99	dgp Fixed bug.
%   5/18/08 mk Vectorized.
%   5/30/09 mk Bugfix in vectorization.

if nargout~=2
	error('Usage: [x, y] = RectCenter(rect);');
end

if PsychNumel(r) == 4
    % Single rect:
    x = round(0.5*(r(1)+r(3)));
    y = round(0.5*(r(2)+r(4)));
else
    % Multi-rect array:
    if (size(r, 1)==4) && (size(r,2)>=1)
        % Multi-column array with one 4-comp. rect per column:
        x = round(0.5*(r(1,:)+r(3,:)));
        y = round(0.5*(r(2,:)+r(4,:)));
    else
        if (size(r, 2)==4) && (size(r,1)>=1)
            % Multi-row array with one 4-comp. rect per row:
            x = round(0.5*(r(:,1)+r(:,3)));
            y = round(0.5*(r(:,2)+r(:,4)));
        else
            % Something weird and unknown:
            error('Given matrix of rects not of required 4-by-n or n-by-4 format.');
        end
    end
end
