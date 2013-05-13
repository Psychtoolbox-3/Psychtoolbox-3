% TextBugTest
%
% Tests for the text bug where using textures before drawing texts causes
% 'DrawText' not to work.  Reported by lots of people.
%
% TextBugTest tests both text before textures and textures before text.  
%
% Screen should now pass this test.


% HISTORY
%
% mm/dd/yy
% 1/21/05   awi     Wrote it.
% 1/25/05   awi     Loop over both order conditions clearing the mex file between.  


textFirst=0;
textureFirst=1;

drawTextureFlag=1;    %flag specifying whether to draw the texture. 
drawTextFlag=1;       %flag specifying whether to draw the text.

try
    for drawOrder=textFirst:textureFirst
        
        clear mex;
        
        if drawOrder==textureFirst
            fprintf('Drawing texture before drawing text\n');
        elseif drawOrder==textFirst
            fprintf('Drawing text before drawing texture\n');
        else
            error('unrecognized drawing order flag');
        end

        % Chosing the display with the highest dislay number is 
        % a best guess about where you want the stimulus displayed.  
        screens=Screen('Screens');
        screenNumber=max(screens);

        % Make th 
        w=Screen('OpenWindow', screenNumber,0,[],32,2);
        Screen('FillRect', w, 0);
        Screen('Flip',w);

        for i=1:2
            %create the texture and draw it.
            if drawTextureFlag && ( (i==1 && drawOrder==textureFirst) || (i==2 && ~(drawOrder==textureFirst))) 
                textureRectSize=100;
                textureRectOrigin=[0 0 textureRectSize textureRectSize];
                textureRect=OffsetRect(textureRectOrigin, 0, 0);
                magicSquare=magic(textureRectSize);
                normImage=round(magicSquare/max(max(magicSquare)) * 255);
                tx=Screen('MakeTexture', w, normImage);
                Screen('DrawTexture', w, tx, textureRectOrigin, textureRect);
                fprintf('  Drawing Texture\n');
            end %if

            %create the text and draw it
            if drawTextFlag ( (i==1 && ~(drawOrder==textureFirst)) || (i==2 && drawOrder==textureFirst))
                Screen('TextFont',w, 'Courier');
                Screen('TextSize',w, 50);
                Screen('TextStyle', w, 0);
                Screen('DrawText', w, 'Hello World!', 100, 200, [0, 0, 255, 255]);
                Screen('TextFont',w, 'Times');
                Screen('TextSize',w, 30);
                Screen('DrawText', w, 'Hit the space bar to continue.', 100, 300, [255, 0, 0, 255]);
                fprintf('  Drawing Text\n');
            end %if

        end %for

        Screen('Flip',w);
        %there seems to be a bug in KbWait where it won't proceed if we hit
        %the return key.  So requesst an arbitrary non-return key keypress.
        %
        fprintf('Press the space bar to continue\n');
        KbWait;
        Screen('CloseAll');

    end % for drawOrder
    
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    Screen('CloseAll');
    rethrow(lasterror);
end % try..catch
