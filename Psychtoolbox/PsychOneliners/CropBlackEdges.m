function coords = CropBlackEdges(plaat)
% coords = croblackedges(image)
% returns coordinates that can be used to cut the black edges from an image
% the trick is to sum the image in the x-direction and find where the sum
% is zero, then you know where there are black edges above and below the
% image; same trick for the edges left and right by summing in the
% y-direction.
%
% COORDS = [xmin xmax ymin ymax]
%
% to cut the edges use:
% cutimage = image(coords(3):coords(4),coords(1):coords(2));
%
% DN 2007       Wrote it
% DN 2008-07-30 Now returns correct indices if one or more or all sides of
%               the image do not have a black edge
% DN 2009-02-02 Now returns error msg if all is black

plaat       = sum(plaat,3);
psychassert(any(plaat(:)),'No image in input matrix');

plaatx      = sum(plaat,2);
plaaty      = sum(plaat,1);

% y coordinates
qdiff       = diff([plaatx])~=0;
if plaatx(1)~=0 && plaatx(end)~=0
    coords(3) = 1;
    coords(4) = size(plaat,1);
else
    if plaatx(1)~=0
        coords(3) = 1;
    else
        coords(3) = find(qdiff,1,'first')+1;
    end
    if plaatx(end)~=0
        coords(4) = size(plaat,2);
    else
        coords(4) = find(qdiff,1,'last');
    end
end

% x coordinates
qdiff       = diff([plaaty])~=0;
if plaaty(1)~=0 && plaaty(end)~=0
    coords(1) = 1;
    coords(2) = size(plaat,2);
else
    if plaaty(1)~=0
        coords(1) = 1;
    else
        coords(1) = find(qdiff,1,'first')+1;
    end
    if plaaty(end)~=0
        coords(2) = size(plaat,2);
    else
        coords(2) = find(qdiff,1,'last');
    end
end
