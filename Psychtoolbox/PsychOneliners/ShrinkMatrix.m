function out = ShrinkMatrix(in,fac)
% out = ShrinkMatrix(IN, FAC)
% Shrinks a 2-D or 3-D matrix IN (an image) by a factor FAC.
% matrix will be truncated horizontally and vertically so that the
% resultant shrunk matrix would have integer sizes in the x an y dimension.
% size of 3rd dimension will not be scaled.
% shrinking is performed by mean computation.

% 05/09/08 DN  Wrote it.
% 13/06/12 DN  Urgh, this only went over the diagonal of each NxN submatrix
%              (where N is the scaling factor)


% input checking
if ndims(in)>3
    error('input is not an image');
end

ys = size(in,1);
xs = size(in,2);

if fac ~= round(fac)
    error('scaling factor must be an integer');
end

hcut = mod(xs,fac);
vcut = mod(ys,fac);
if hcut~=0
    disp(sprintf('Warning: right edge of input will be truncated by %d pixels',hcut));
end
if vcut~=0
    disp(sprintf('Warning: lower edge of input will be truncated by %d pixels',vcut));
end

if hcut~=0 || vcut~=0
    in = in(1:end-vcut,1:end-hcut,:);
    xs = xs-hcut;
    ys = ys-vcut;
end

out = zeros(ys/fac,xs/fac,size(in,3));
for p=1:fac     % image rows
    for q=1:fac     % image columns
        out = out + in(p:fac:ys,q:fac:xs,:);
    end
end
out = out./(fac^2);
