function r = InsetRect(r,dh,dv)
% r = InsetRect(r,dh,dv)
%
% Shrinks the rectangle by pushing in the left and right
% sides, each by the amount dh, and the top and bottom by the amount dv.
% Negative values of dh and dv produce expansion.
% Any or all of the arguments can be rect array to, for instance perform
% the same operation of a number of rects, or shrink/expand by a different
% number of pixels for the one input rect, or multipe input rects
%
% See also PsychRects/Contents.

% 5/21/96 Denis Pelli
% 7/10/96 dgp PsychRects
% 8/5/96  dgp check rect size.
% 7/26/15 dcn vectorized
% 7/29/15 dcn Can now handle 4xM and Mx4 rect-arrays

if nargin~=3
	error('Usage:  r = InsetRect(r,dh,dv)');
end
if numel(r) == 4 && isscalar(dh) && isscalar(dv)
    % Single rect, single offsets. Multiply per element so that column- and
    % row-vector rects work
    r(1) = r(1)+dh;
    r(2) = r(2)+dv;
    r(3) = r(3)-dh;
    r(4) = r(4)-dv;
else
    % Multi-rect array:
    if (size(r, 1)==4) && (size(r,2)>=1)
        dh=dh(:).';
        dv=dv(:).';
        szs = [size(r,2) size(dh,2) size(dv,2)];
        maxsz = max(szs);
        if ~all(szs==1 | szs==maxsz)
            error('For all inputs, size of second dimension must either be one or equal to the longest input')
        end
        dh = repmat(dh,1,maxsz/szs(2));
        dv = repmat(dv,1,maxsz/szs(3));
        r  = repmat(r,1,maxsz/szs(1))+SetRect(dh,dv,-dh,-dv);
    elseif (size(r, 2)==4) && (size(r,1)>=1)
        dh=dh(:);
        dv=dv(:);
        szs = [size(r,1) size(dh,1) size(dv,1)];
        maxsz = max(szs);
        if ~all(szs==1 | szs==maxsz)
            error('For all inputs, size of first dimension must either be one or equal to the longest input')
        end
        dh = repmat(dh,maxsz/szs(2),1);
        dv = repmat(dv,maxsz/szs(3),1);
        r  = repmat(r,maxsz/szs(1),1)+SetRect(dh,dv,-dh,-dv);
    else
        % Something weird and unknown:
        error('Given matrix of rects not of required 4-by-n or n-by-4 format.');
    end
end
