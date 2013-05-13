% iviewx

    %   14-12-04 volgende zaken uitgeprobeerd.
    %   ET_SIM: start eye image
    %   ET_EIM: end eye image
    %   ET_IMG : image size
    
    %   er komt wel informatie over, maar het is onduidelijk wat deze
    %   precies betekent. (stukjes van 1336 char, en af en toe ook 15 char
    %   (die laatste lijken alleen spaces te zijn). Eye image is 384x288
    %   groot.
    
    
    %   ET_MOV ET_POS
    %   worden niet herkend door iviewx 1.6 build 10




fprintf([mfilename ' start\n']);
clear all;
try

    red.right=KbName('RightArrow');
    red.left=KbName('LeftArrow');
    red.up=KbName('UpArrow');
    red.down=KbName('DownArrow');
    red.pos=KbName('p');
    
    
    pnet('closeall');

    host='192.168.1.2';
    port=4444;
    ivx=iViewXInitDefaults; %([],[], host, port);
    ivx;
    [result, ivx]=iViewXComm('open', ivx);
    result
    [result, ivx]=iViewXComm('check', ivx);
    result
    [success, ivx]=iViewX('setscreensize', ivx);
    
    result=iViewXComm('send', ivx, ['ET_SIM']); % start eye_image
%     [result, ivx]=iViewX('stoprecording', ivx);

%     params='ET_REC';
%     [result, ivx]=iViewX('command', ivx, params);
%     params='ET_CAL 9';
%     [result, ivx]=iViewX('command', ivx, params);

    tEnd=getsecs+300;
    while getsecs<tEnd
        [keyIsDown,secs,keyCode] = KbCheck;
        if 1==keyCode(ivx.modifierKey) && 1==keyCode(ivx.breakKey)
            result=-1;
            return;
            break;
        end
% duration=10;
%         if 1==keyCode(red.left)
%             result=iViewXComm('send', ivx, ['ET_MOV' ' LEFT ' num2str(duration)]);
%         elseif 1==keyCode(red.right)
%             result=iViewXComm('send', ivx, ['ET_MOV' ' RIGHT ' num2str(duration)]);
%         elseif 1==keyCode(red.up)
%             result=iViewXComm('send', ivx, ['ET_MOV' ' UP ' num2str(duration)]);
%         elseif 1==keyCode(red.down)
%             result=iViewXComm('send', ivx, ['ET_MOV' ' DOWN ' num2str(duration)]);
%         end
% 
        if 1==keyCode(red.left)
            result=iViewXComm('send', ivx, ['ET_MOV' ' LEFT ']);
        elseif 1==keyCode(red.right)
            result=iViewXComm('send', ivx, ['ET_MOV' ' RIGHT ']);
        elseif 1==keyCode(red.up)
            result=iViewXComm('send', ivx, ['ET_MOV' ' UP ' ]);
        elseif 1==keyCode(red.down)
            result=iViewXComm('send', ivx, ['ET_MOV' ' DOWN ' ]);
        end
        
        
        if 1==keyCode(red.down)
            result=iViewXComm('send', ivx, ['ET_MOV' ' POS ' num2str(100) ' ' num2str(100) ]);
        end

        [data, ivx]=iViewXComm('receive', ivx);

        if ~isempty(data)
           data;
            if strfind(data, 'ET_SPL')
            elseif strfind(data, 'ET_CHG')
                pt=str2num(data(8:end));
                fprintf('Change to point: %d\n', pt);
            elseif strfind(data, 'ET_PNT')
                coords=str2num(data(8:end));
                fprintf('Coords for point %d:\t%d\t%d\n', coords(1), coords(2), coords(3));
            elseif strfind(data, 'ET_CSZ')
                scr_size=str2num(data(8:end));
                fprintf('Screen: %d\t%d\n', scr_size(1), scr_size(2));
            elseif strfind(data, 'ET_IMG')
                eye_size=str2num(data(8:end));
                fprintf('Eye image: %d\t%d (total: %d)\n', eye_size(1), eye_size(2), eye_size(1)* eye_size(2));
                ivx.udpmaxread=eye_size(1)* eye_size(2);
                ivx.udpreadtimeout=1.0;
                [data, ivx]=iViewXComm('receive', ivx);
                ivx.udpmaxread=1000;
                ivx.udpreadtimeout=0.1;
                size(data)
                decdata=base64decode(data)
%                 if length(data)==15
%                     data
%                 end
                
            end
        end
       % result
        i=i+1;
        drawnow;
    end
    [result, ivx]=iViewXComm('close', ivx);
    
    result=iViewXComm('send', ivx, ['ET_EIM']);
    

    result
    
    [result, ivx]=iViewX('stoprecording', ivx);

    result
    pnet('closeall');
    fprintf([mfilename ' end\n']);

catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    pnet('closeall');
    Screen('CloseAll');    rethrow(lasterror);
end %try..catch..

