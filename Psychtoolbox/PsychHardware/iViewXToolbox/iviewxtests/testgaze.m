              clear all;
try
    fprintf('Very simple OSX iViewX gaze plot demo\n\n\t');
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
    Screen('FillRect',w, gray);
    Screen('Flip',w);
    if 0
        rect=CenterRect(screenRect/4, screenRect);
        Screen('FillOval',w, [255  0 0], rect); % red circle
        Screen('Flip',w);
    else
        %
        % initialize
        if 0
            host='localhost';
            port=6666;
        else
            host=[];
            port=[];
        end
        ivx=initIViewXDefaults(w, [], host, port );
        if 0 ivx.localport=4444; end

        [success, ivx]=iViewX('calibrate', ivx);
        stop=0;
        if success~=1
            fprintf([mfilename ': calibration failed\n']);
            stop=1;
        end

    end

    while KbCheck; end

    if stop==0
        % open connection for receiving data
        [success, ivx]=iViewX('openconnection', ivx);
        [success, ivx]=iViewX('datastreamingon', ivx);

        gazeRect=[0 0 20 20];
        oldRect=gazeRect;
        newRect=gazeRect;

        tEnd=GetSecs+120;
        i=0;
        k=0;
        data=[];
        while ~KbCheck && GetSecs<tEnd;

%             if mod(i,500)==0 fprintf('%d\n',i); end
%             i=i+1;
%             [data, ivx]=iViewX('receivedata', ivx);
            [data, ivx]=iViewXComm('receive', ivx);
            % data
            if ~isempty(data) && data~=-1
                if 1==strfind(data, 'ET_SPL') % spooled data
                    mygaze=str2num(data(8:end));
                    %                    fprintf('%d\t%d\t%d\t%d\n', k, mydata(1), mydata(2), mydata(3));
                    newRect = CenterRectOnPoint(gazeRect, mygaze(2), mygaze(3));
                    if length(find(oldRect~=newRect))>0
                        %                 fprintf('Tekenen bij %d %d\n', x, y );
                        Screen('FillOval',w, gray, oldRect); % gray circle
                        Screen('FillOval',w, [255 0 0], newRect); % red circle
                        Screen('Flip',w);
                        oldRect=newRect;
                    end
                end
            end
        end
  
        [success, ivx]=iViewX('datastreamingoff', ivx);

        [success, ivx]=iViewX('closeconnection', ivx);

    end   


    Screen('CloseAll');
    fprintf('\nEnd of demo.\n');
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    %     pnet('closeall');
    Screen('CloseAll');    rethrow(lasterror);
end %try..catch..


