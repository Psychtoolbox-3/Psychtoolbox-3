% iviewxsimulator
% simulates iviewx host behaviour


fprintf([mfilename ' start\n']);
clear all;
try
    streaming=0;
    %     [result, ivx]=iViewX('clear', ivx);

    pnet('closeall');

    host='localhost';
    port=4444;
    ivx=iViewXInitDefaults([],[], host, port);
    ivx.localport=6666;
    %     ivx.stopSimulatorKey=KbName('ESC');
    ivx;
    [result, ivx]=iViewX('openconnection', ivx);

    while KbCheck; end

    % start a loop, waiting for commands from client
    tEnd=GetSecs+300;
    i=0;
    while GetSecs<tEnd

        [keyIsDown,secs,keyCode] = KbCheck;
        if 0 && 1==keyCode(ivx.modifierKey) && 1==keyCode(ivx.breakKey)
            break;
        end

        if streaming==1

            message=['ET_SPL ' num2str(round(GetSecs*100000)) ' ' num2str(round(rand*ivx.screenHSize)) ' ' num2str(round(rand*ivx.screenVSize)) ];
            [result, ivx]=iViewX('command', ivx, message);
        end

        [data, ivx]=iViewX('receivedata', ivx);

        if ~isempty(data) && data ~=-1
            % data
            if strfind(data, 'ET_STR')
                fprintf('Start Streaming\n');
                streaming=1;
            elseif strfind(data, 'ET_ACC')
                fprintf('Accept point\n');
            elseif strfind(data, 'ET_EST')
                fprintf('End Streaming\n');
                streaming=0;
            elseif strfind(data, 'ET_REC')
                fprintf('Start Recording\n');
            elseif strfind(data, 'ET_STP')
                fprintf('Stop Recording\n');
            elseif strfind(data, 'ET_PNT')
                coords=str2num(data(8:end));
                fprintf('Coords for point %d:\t%d\t%d\n', coords(1), coords(2), coords(3));
                
                ivx.absCalPos(coords(1),1)=coords(2);
                ivx.absCalPos(coords(1),2)=coords(3);
                
            elseif strfind(data, 'ET_CAL')
                ct=str2num(data(8:end));
                fprintf('Calibration type: %d\n', ct);
                ivx.nCalPoints=ct;

                % send screen info
                WaitSecs(0.01);

                message=['ET_CSZ ' num2str(ivx.screenHSize) ' ' num2str(ivx.screenVSize) ];
                [result, ivx]=iViewX('command', ivx, message);
                fprintf('Sending screen info \n');
                % now start sending calibration point info

                for i=1:ivx.nCalPoints
                    message=['ET_PNT ' num2str(i) ' ' num2str(ivx.absCalPos(i,1)) ' ' num2str(ivx.absCalPos(i,2))];
                    message;
                    [result, ivx]=iViewX('command', ivx, message);
                    fprintf('Sending point info %d\n', i);

                    WaitSecs(0.01);
                end

                for i=[1:ivx.nCalPoints 1]
                    message=['ET_CHG ' num2str(i)];

                    [result, ivx]=iViewX('command', ivx, message);
                    fprintf('Sending point change info %d\n', i);

                    WaitSecs(1+rand*1);
                end
                message=['ET_FIN ' num2str(i)];
                [result, ivx]=iViewX('command', ivx, message);

            elseif strfind(data, 'ET_CHG')
                pt=str2num(data(8:end));
                fprintf('Change to point: %d\n', pt);
            elseif strfind(data, 'ET_PT')
                coords=str2num(data(8:end));
                fprintf('Coords for point %d:\t%d\t%d\n', coords(1), coords(2), coords(3));
            elseif strfind(data, 'ET_CSZ')
                scr_size=str2num(data(8:end));
                ivx.screenHSize=scr_size(1);
                ivx.screenVSize=scr_size(2);


                fprintf('Screen: %d\t%d\n', scr_size(1), scr_size(2));
            else
                fprintf('Unknown command.\n');

            end
            data=[];
        end
    end
    %     params='ET_REC';
    %     [result, ivx]=iViewX('command', ivx, params);
    %     params='ET_CAL 9';
    %     [result, ivx]=iViewX('command', ivx, params);
    %
    [result, ivx]=iViewX('closeconnection', ivx);

    %     result

    %     [result, ivx]=iViewX('stoprecording', ivx);

    %     result
    pnet('closeall');
    fprintf([mfilename ' end\n']);

catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    pnet('closeall');
    Screen('CloseAll');    rethrow(lasterror);
end %try..catch..

