% iviewx simulator test
fprintf([mfilename ' start\n']);
clear all;
try

    
    
    
    
    
    pnet('closeall');

    host='localhost';
    port=6666;
    ivx=initIViewXDefaults([],[], host, port);
    ivx.localport=4444;
    ivx;

    [result, ivx]=iViewX('stoprecording', ivx);

    %     fprintf([mfilename ' hier\n']);

    [result, ivx]=iViewX('calibrate', ivx);
    if result==-1
        fprintf([mfilename ': calibration failed.\n']);
    end


    if 0

        [result, ivx]=iViewX('openconnection', ivx);
        result;
        while KbCheck; end

        tEnd=GetSecs+60;

        i=0;
        while GetSecs<tEnd

            [keyIsDown,secs,keyCode] = KbCheck;
            if 1==keyCode(ivx.modifierKey) && 1==keyCode(ivx.breakKey)
                break;
            end
            [data, ivx]=iViewXComm('receive', ivx);

            if ~isempty(data)
                % data
                if strfind(data, 'ET_SPL')
                    mydata=str2num(data(8:end));
                    fprintf('Gaze:\t%d\t%d\t%d\n', mydata(1), mydata(2), mydata(3));

                elseif strfind(data, 'ET_CHG')
                    pt=str2num(data(8:end));
                    fprintf('Change to point: %d\n', pt);
                elseif strfind(data, 'ET_PNT')
                    coords=str2num(data(8:end));
                    fprintf('Coords for point %d:\t%d\t%d\n', coords(1), coords(2), coords(3));
                elseif strfind(data, 'ET_CSZ')
                    scr_size=str2num(data(8:end));
                    fprintf('Screen: %d\t%d\n', scr_size(1), scr_size(2));
                end
            end
            % result
            i=i+1;
            drawnow;
        end
    end
    [success, ivx]=iViewX('closeconnection', ivx);
    [success, ivx]=iViewX('stoprecording', ivx);
    pnet('closeall');
    fprintf([mfilename ' end\n']);

catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    pnet('closeall');
    Screen('CloseAll');    rethrow(lasterror);
end %try..catch..

