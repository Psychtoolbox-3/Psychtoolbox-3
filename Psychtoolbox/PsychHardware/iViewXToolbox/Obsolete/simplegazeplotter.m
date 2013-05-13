clear all;
try
    fprintf('Very simple OSX  screen plotter demo\n\n\t');
    fprintf('At the end of the demo, press any key to quit\n\n\t');

    input('Hit the return key to continue.','s');
    fprintf('Thanks.\n');


    lport=4444;

    
    ivx=initiviewxdefaults;
    ivx

%     params=['ET_CSZ ' num2str(1024) ' ' num2str(768)];
%     [result, ivx]=iViewX('command', ivx, params);

    
    % Open  udpsocket and bind udp port adress to it.
    udp=pnet('udpsocket',lport);
    pnet(udp,'setreadtimeout',.1);
    stat=pnet(udp,'status');
    screenNumber=max(Screen('Screens'));

    white=WhiteIndex(screenNumber);
    black=BlackIndex(screenNumber);
    gray=round((white+black)/2);
    % Open a double buffered fullscreen window and draw a gray background
    % and front and back buffers.
    if 1
        [w, screenRect]=Screen('OpenWindow',screenNumber, 0,[],32,2);
        rect=CenterRect(screenRect/4, screenRect);
        Screen('FillRect',w, gray);
        %   Screen('FillOval',w, [255 0 0], rect); % red circle

        gazeRect=[0 0 10 10];
        oldRect=gazeRect;
        newRect=gazeRect;

        %   WaitSecs(2);
        Screen('Flip',w);
    end
    while KbCheck; end
    tEnd=GetSecs+120;
    
    % open connection for receiving data
    [success, ivx]=iViewX('openconnection', ivx);

    if success==-1
        fprintf([mfilename ': cannot open a connection to iViewX host.\n']);
        return
    end

    
    
    i=0;
    k=0;
    data=[];
    while ~KbCheck && GetSecs<tEnd;
        %         stat=pnet(udp,'status');
        %if stat>=0 fprintf('Stat: %d\n', stat); end
        % Wait/Read udp packet to read buffer
        len=pnet(udp,'readpacket');
        %len=pnet(udp,'readpacket',[],'noblock');
        %         if len>0 fprintf('Len: %d\n', len); end

        if mod(i,500)==0 fprintf('%d\n',i); end
        i=i+1;
        if len>0,
            [ip,port]=pnet(udp,'gethost');
            % if packet larger then 1 byte then read maximum of 1000 doubles in network byte order
            %  data=pnet(udp,'read',1000,'double');
            k=k+1;
            data=pnet(udp,'read',1000);
            % data
            if 1
                if 1==strfind(data, 'ET_SPL') % spooled data
                    mydata=str2num(data(8:end));
%                    fprintf('%d\t%d\t%d\t%d\n', k, mydata(1), mydata(2), mydata(3));
                end
            end
            %data
            if 1
%                 data=str2num(data);
%                 data
                x=mydata(2);
                y=mydata(3);
%                 x
%                 y
                newRect = CenterRectOnPoint(gazeRect, x, y);
                tekenen=oldRect==newRect;
                if length(find(oldRect~=newRect))>0
                    %                 fprintf('Tekenen bij %d %d\n', x, y );
                    Screen('FillRect',w, gray, oldRect); % gray circle
                    Screen('FillRect',w, [255 0 0], newRect); % red circle
                    Screen('Flip',w);
                    oldRect=newRect;
                end
            end
        end
    end

    [success, ivx]=iViewX('closeconnection', ivx);

    Screen('CloseAll');
    fprintf('\nEnd of demo.\n');
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    pnet('closeall');

    Screen('CloseAll');
    rethrow(lasterror);
end %try..catch..


