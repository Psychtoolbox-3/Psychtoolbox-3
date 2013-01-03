function outImage = BitsPlusPackColorImage(inImage,SPACESCALE,COLORSCALE)
% [outImage,inImage] = BitsPlusPackColorImage(inImage,[SPACESCALE],[COLORSCALE]))
%
% Convert an image for Bits++, assuming it will
% displayed on an even pixel boundary and that
% it has an even number of columns.
%
% Image is resized if SPACESCALE == 1 (default).
% Otherwise scaling is assumed to have been done in advance.
%
% Assume input is in range 0-1 if COLORSCALE == 1 (default).
% Assume input is in range 0-65535 and be integer if COLORSCALE == 0.
%
% If scaling in space, must scale in color, otherwise interpolation
% can violate assumptions of prescaled color (e.g. integer vals).
%
% 8/9/04	dhb		Wrote it.
% 18/4/05   ejw     Converted it to run with OSX version of Psychtoolbox.
% 26.2.07   mk      Bugfix for LSB conversion: Added modulo operation.
% 3/01/07   mk      Bugfix for MSB conversion: Added floor operation.

if (nargin < 2 || isempty(SPACESCALE))
	SPACESCALE = 1;
end
if (nargin < 3 || isempty(COLORSCALE))
	COLORSCALE = 1;
end

if (SPACESCALE && ~COLORSCALE)
	error('COLORSCALE must be set if SPACESCALE is.');
end

if (SPACESCALE)
	[m,n,p] = size(inImage);
	if (rem(n,2) ~= 0)
		error('Image must have even number of columns');
    end
    % imresize is part of the image processing toolbox
    % I have recoded to avoid the need for this
	%inImage = imresize(inImage,[m n/2],'bilinear');
   
    % halve the horizontal width of the image
    % this compensates for the 2:1 aspect ratio for colour++ pixels
    inImage = inImage(:,[1:2:n],:);
   
	index = find(inImage < 0);
	if (~isempty(index))
		inImage(index) = 0; 
	end
	index = find(inImage > 1);
	if (~isempty(index))
		inImage(index) = 1; 
	end
else
	[m,n,p] = size(inImage);
	n = n*2;
end

if (COLORSCALE)
	inImage = round(65535*inImage);
end

% Make output image and stuff with high and low bytes
outImage = uint8(zeros(m,n,p));

% First byte is high 8 bits
outImage(:,1:2:end-1,:) = uint8(floor(inImage/256));

% Second byte is low 8 bits
outImage(:,2:2:end,:) = uint8(mod(inImage, 256));



