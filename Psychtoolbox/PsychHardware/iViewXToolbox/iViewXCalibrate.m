function [result, ivx]=iViewXCalibrate(ivx)

% calibration routine of iViewX toolbox
% could also double as driftcorrection with a single point?
result=1;

try
    if ~exist('ivx', 'var') || isempty(ivx)
        txt=[ mfilename ' requires a structure with iViewX default values as input.'];
        error(txt);
    end

    % open connection for receiving data
    [success, ivx]=iViewX('openconnection', ivx);

    if success==-1
        fprintf([mfilename ': cannot open a connection to iViewX host.\n']);
        return
    end
    % send screen info
    [success, ivx]=iViewX('setscreensize', ivx);
    [success, ivx]=iViewX('startrecording', ivx);
    
    % send (possibly adapted) calibration point positions
%     [success, ivx]=iViewX('setcalibration', ivx, ivx.nCalPoints); %??
%     [success, ivx]=iViewX('sendcalibrationpoints', ivx);

    [success, ivx]=iViewX('setcalibration', ivx, ivx.nCalPoints);
   % now we should wait for iView to send
    % information on calibration screen size
    % and calibration point positions
    tBreak=GetSecs+ivx.commTimeOut;
    npts=0;
    klaar=0;
    while klaar==0
        [keyIsDown,secs,keyCode] = KbCheck;
        if 1==keyCode(ivx.modifierKey) && 1==keyCode(ivx.breakKey)
            result=-1;
            return;
            break;
        end

        if GetSecs>tBreak
            fprintf([mfilename ': timed out waiting for data.\n']);
            result=-1;
            return
        end
        [data, ivx]=iViewX('receivedata', ivx);

        if ~isempty(data) && data ~=-1
            
%             data
            
            tBreak=GetSecs+ivx.commTimeOut; % reset time out

            if strfind(data, 'ET_CHG')
                pt=str2num(data(8:end));
                fprintf('Change to point: %d\n', pt);
                % this signals we have all the points
                klaar=1;
                break;
            elseif strfind(data, 'ET_PNT')
                coords=str2num(data(8:end));
                %                 fprintf('Coords for point %d:\t%d\t%d\n', coords(1), coords(2), coords(3));
                npts=npts+1;
                if coords(1) ~=npts
                    fprintf('Point order mismatch\n');
                end
                ivx.absCalPos(coords(1),:)=coords(2:3); % 1 is point nr.
            elseif strfind(data, 'ET_CSZ')
                % we should receive this message first
                scr_size=str2num(data(8:end));
                fprintf('Screen: %d\t%d\n', scr_size(1), scr_size(2));
                % we could check if the screen size reported matches the one
                % expected
            end
        end
    end


    fprintf('Calibrating using %d points (%d received).\n', size(ivx.absCalPos,1), npts);

    iViewX('erasescreen', ivx);

    nPointsShown=0;
    allowManAccept=1;
    % we now draw first calibration point and then
    % start a loop in which we wait for iView to tell us
    % when to draw the next calibration point

    tBreak=GetSecs+ivx.commTimeOut;

    klaar=0;
    while klaar==0
        % draw calibration point

        iViewX('drawcalibrationpoint', ivx, [ivx.absCalPos(pt,1), ivx.absCalPos(pt,2)]);

        nPointsShown=nPointsShown+1;
        nextPt=0;
        %     After displaying the first point, the software must monitor
        %     the serial port for further ÒET_CHGÓ commands or  ÒET_BRKÓ,
        %     which cancels the calibration procedure.
        % we also detect manual acceptance of calibration point
        while nextPt==0
            if GetSecs>tBreak
                fprintf([mfilename ': timed out waiting for data.\n']);
                result=-1;
                return
            end
            [keyIsDown,secs,keyCode] = KbCheck;
            if 1==keyCode(ivx.modifierKey) && 1==keyCode(ivx.breakKey)
                % stop program altogether
                result=-1;
                iViewX('break', ivx);
                return;
            end
            if keyCode(ivx.breakKey), % stop calibration
                result=-1;
                iViewX('break', ivx);
                return;
            end
            if keyCode(ivx.nextCalPointKey), % force next calibration point
                iViewX('accept', ivx);
                nextPt=1;
            end

            % If you want the subject to manually accept the point with a
            % key press once they fixate, the ÒET_ACCÓ  command must then
            % be sent to the serial port.  This instructs iView to change
            % to the next calibration point  rather than having it detect
            % fixations itself.
            if keyIsDown==0
                allowManAccept=1; % we only accept manually if key has been released in between
            end

            if keyCode(ivx.calAcceptKey) && allowManAccept==1, % manually accept calibration point
                iViewX('accept', ivx);
                allowManAccept=0;
            end

            [data, ivx]=iViewX('receivedata', ivx);

            data;

            if ~isempty(data) && data ~=-1
                % data
                if strfind(data, 'ET_CHG')
                    pt=str2num(data(8:end));
                    fprintf('Change to point: %d\n', pt);
                    % this signals the point was fixated correctly so we can
                    % go the next point
                    nextPt=1;
                    break;
                elseif strfind(data, 'ET_BRK')
                    result=-1;
                    return
                elseif strfind(data, 'ET_FIN')
                    %                     fprintf('Point: %d (pos %d) finished.\n', nPointsShown, pt);

                    result=1;
                    klaar=1;
                    fprintf('Calibration Finished\n');
                    break
                end
            end
        end
    end

    [success, ivx]=iViewX('stoprecording', ivx); % ?? necessary

    iViewX('erasescreen', ivx);
    WaitSecs(ivx.calExtraTimeSecs);

    [success, ivx]=iViewX('closeconnection', ivx);
    if success~=1
        fprintf([mfilename ': could not close connection.\n']);
    end
    fprintf('\nEnd of calibration.\n');
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.

    Screen('CloseAll');
    rethrow(lasterror);
end %try..catch..


