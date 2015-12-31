function r = ScaleRect(r,horizontalFactor,verticalFactor)
% r = ScaleRect(r,horizontalFactor,verticalFactor)
%
% Scales a rect, by multiplying the left and right coordinates by the 
% horizontal scale factor and multiplying the top and bottom coordinates
% by the vertical scale factor.
% Input can be multiple rects, concatenated into a Mx4 matrix
%
% See also PsychRects/Contents, Expand, InsetRect.

% 5/27/96 dgp Denis Pelli.
% 7/10/96 dgp PsychRects.
% 8/5/96  dgp check rect size.
% 1/13/10 dcn Input can now be multiple rects, concatenated in row
%             direction.
% 7/26/15 dcn Further vectorized, horizontal- and verticalFactor can now
%             also be vectors.
% 7/29/15 dcn Can now handle 4xM and Mx4 rect-arrays

if nargin~=3
	error('Usage: r=ScaleRect(r,horizontalFactor,verticalFactor)');
end

if numel(r) == 4 && isscalar(horizontalFactor) && isscalar(verticalFactor)
    % Single rect, single offsets. Multiply per element so that column- and
    % row-vector rects work
    r(1) = r(1)*horizontalFactor;
    r(2) = r(2)*verticalFactor;
    r(3) = r(3)*horizontalFactor;
    r(4) = r(4)*verticalFactor;
else
    % Multi-rect array:
    if (size(r, 1)==4) && (size(r,2)>=1)
        % 4xM
        horizontalFactor = horizontalFactor(:).';
        verticalFactor = verticalFactor(:).';
        szs = [size(r,2) size(horizontalFactor,2) size(verticalFactor,2)];
        maxsz = max(szs);
        if ~all(szs==1 | szs==maxsz)
            error('For all inputs, size of second dimension must either be one or equal to the longest input')
        end
        horizontalFactor = repmat(horizontalFactor,1,maxsz/szs(2));
        verticalFactor   = repmat(verticalFactor,1,maxsz/szs(3));
        r  = repmat(r,1,maxsz/szs(1)).*SetRect(horizontalFactor,verticalFactor,horizontalFactor,verticalFactor);
    elseif (size(r, 2)==4) && (size(r,1)>=1)
        % Mx4
        horizontalFactor = horizontalFactor(:);
        verticalFactor = verticalFactor(:);
        szs = [size(r,1) size(horizontalFactor,1) size(verticalFactor,1)];
        maxsz = max(szs);
        if ~all(szs==1 | szs==maxsz)
            error('For all inputs, size of first dimension must either be one or equal to the longest input')
        end
        horizontalFactor = repmat(horizontalFactor,maxsz/szs(2),1);
        verticalFactor   = repmat(verticalFactor,maxsz/szs(3),1);
        r  = repmat(r,maxsz/szs(1),1).*SetRect(horizontalFactor,verticalFactor,horizontalFactor,verticalFactor);
    else
        % Something weird and unknown:
        error('Given matrix of rects not of required 4-by-n or n-by-4 format.');
    end
end
