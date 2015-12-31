function newRect = SetRect(left,top,right,bottom)
% newRect = SetRect(left,top,right,bottom);
%
% Create a rect with the specified coordinates.
% This is equivalent to:
%	newRect=[left,top,right,bottom];
% Each of the input can be either a scalar or
% an M-element vector. For those inputs that are
% vectors, they must all have the same length.
% If any of the inputs is a row vector, an 4xM
% rect-array will be output, otherwise if any of
% the inputs is a column vector, a Mx4 rect-array
% will be output.
%
% See also PsychRects/Contents.

% dgp 5/12/96  Make sure the created rect is a column, not a row, vector
% dgp 5/12/96  Use symbolic indices, added usage check.
% dhb 5/13/96  Put created rect back as row vector.
% 7/10/96 dgp  "help PsychRects"
% 2/25/97 dgp  updated
% 7/26/15 dcn  vectorized
% 7/29/15 dcn  Can now produce 4xM or Mx4 depending on orientation of input

if nargin~=4
    error('Usage: newRect=SetRect(left,top,right,bottom)')
end
szs = [length(left) length(top) length(right) length(bottom)];
if all(szs==1)
    newRect=[left top right bottom];
else
    maxsz = max(szs);
    if ~all(szs==1 | szs==maxsz)
        error('All inputs must either have one element or the same number of elements as the longest input')
    end
    nCol = [size(left,2) size(top,2) size(right,2) size(bottom,2)];
    if all(nCol==1)
        % Mx4
        newRect=zeros(maxsz,4);
        newRect(:,2)=top;
        newRect(:,1)=left;
        newRect(:,4)=bottom;
        newRect(:,3)=right;
    else
        % 4xM
        newRect=zeros(4,maxsz);
        newRect(2,:)=top;
        newRect(1,:)=left;
        newRect(4,:)=bottom;
        newRect(3,:)=right;
    end
end