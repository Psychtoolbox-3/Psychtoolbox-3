% TextureTest
%
% Exercises the new Screen('DrawTexture') command in the OS X Psychtolbox 
% which replaces Screen('CopyWindow').   
%

% HISTORY
%  6/28/04    awi     wrote it.

% Tests which we make of textures:
%
%  That that after drawing a texture onto the screen that GetImage returns
%  the identical matrix


screenNumber=max(Screen('Screens'));
% Open and onscreen window on the second display 
w=Screen('OpenWindow',screenNumber, 0,[],32,2);

%make some textures out of matrices
gradientVector=[0:255];
gradientMat=repmat(gradientVector, 256, 1);
gradientMatColor=gradientMat;
gradientMatColor(:,:,2)=fliplr(gradientMat(:,:,1));
gradientMatColor(:,:,3)=gradientMat(:,:,1)';
gradientMatColorAlpha=gradientMatColor;
gradientMatColorAlpha(:,:,4)=ones(256) * 255;
angleMat=zeros(256);
angleMat(1:128,1)=255;
for im=1:256
    for in=1:256
        if im==in
            angleMat(im,in)=255;
        end
    end
end
tColor=Screen('MakeTexture', w, gradientMatColor);
tColorAlpha=Screen('MakeTexture', w, gradientMatColor);
a=Screen('MakeTexture', w, angleMat);
t=Screen('MakeTexture', w, gradientMat);

% See if drawing works at all
Screen('FillRect', w, 0);
Screen('FillRect', w, [255 0 0 255], [600 600 800 800]);

targetRect=RectOfMatrix(gradientMat);
Screen('DrawTexture', w, a, targetRect, OffsetRect(targetRect, 100, 100));
Screen('Flip', w);

WaitSecs(5);
Screen('CloseAll');
