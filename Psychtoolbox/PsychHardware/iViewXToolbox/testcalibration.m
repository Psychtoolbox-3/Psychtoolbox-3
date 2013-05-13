clear all;
try
    fprintf('Very simple OSX  iViewX calibration demo\n\n\t');
    %     fprintf('At the end of the demo, press any key to quit\n\n\t');

    %     input('Hit the return key to continue.','s');
    %     fprintf('Thanks.\n');


    screenNumber=max(Screen('Screens'));

    white=WhiteIndex(screenNumber);
    black=BlackIndex(screenNumber);
    gray=round((white+black)/2);
    % Open a double buffered fullscreen window and draw a gray background
    % and front and back buffers.

    [w, screenRect]=Screen('OpenWindow',screenNumber, 0,[],32,2);
    Screen('FillRect',w, gray)       ;
    Screen('Flip',w);
    %     if 0
    %         rect=CenterRect(screenRect/4, screenRect);
    %         Screen('FillOval',w, [255 0 0], rect); % red circle
    %         Screen('Flip',w);
    %     else
    %
    % initialize
    if 0 % if we're testing using the iviewxsimulator
        host='localhost';
        port=6666;
    else
        host=[];
        port=[];
    end
    ivx=iViewXInitDefaults(w, [], host, port );

    % change default positions of calibration points
    offset=[0 0]; % default no offset
    scale=[-40 -40];  % default sta ndard  scaling

    ivx=iViewXSetCalPoints(ivx, ivx.nCalPoints, offset, scale);
  


    %         if 0 ivx.localport=4444; end

    [success, ivx]=iViewX('calibrate', ivx);

    if success~=1
        fprintf([mfilename ': calibration failed\n'])
    end
    %         [result, ivx]=iViewX('driftcorrection', ivx);
    %
    %         [result, ivx]=iViewX('driftcorrection', ivx, [200 200]);
    %
    %     end
    Screen('Flip',w);
    %   WaitSecs(2);
    while KbCheck; end
    tEnd=GetSecs+2;
    while ~KbCheck && GetSecs<tEnd, end

    Screen('CloseAll');
    fprintf('\nEnd of demo.\n');
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    %     pnet('closeall');
    Screen('CloseAll');    rethrow(lasterror);
end %try..catch..


