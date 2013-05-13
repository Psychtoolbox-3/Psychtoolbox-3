clear all;
try
    fprintf('Very simple OSX  iViewX calibration demo\n\n\t');

    w=[];
    ivx=initIViewXDefaults(w, 5 );

    
    [result, ivx]=iViewX('initialize', ivx); % opens connection
    
    result
    
    ivx
     fprintf([mfilename ' 1\n']);
   
    % send calibration command
    [result, ivx]=iViewX('startrecording', ivx);
    fprintf([mfilename ' 2\n']);
    [result, ivx]=iViewX('setcalibration', ivx, ivx.nCalPoints);
    fprintf([mfilename ' 3\n']);

    % now we should wait for iView sending
    % information on calibration screen size
    % and calibration point positions
    [result, ivx]=iViewX('getcalibrationpoints', ivx);

    if result==-1
        return
    end

    
    while 1

        % monitor port for ET_CHG or ET_BRK commands
        ontvangen=[];
        [ontvangen, ivx]=iViewX('receive', ivx);
        
        ontvangen
        % result can now be struct with commands and params
        if 1==isstruct(ontvangen)
            switch ontvangen(1).command
                case 'ET_CHG',
                    % go to the next point
                    i=ontvangen(1).param;
                    break; % does this get us out of the loop?
                case 'ET_BRK', % cancel calibration
                    result=-1;
                    % iViewX('break', ivx);
                    return;
                case 'ET_FIN', % finished calibration
                    result=-1;
                    iViewX('stop', ivx); % ?? necessary
                    return;
                otherwise,
                    % do nothing
            end
        end
        % check if subject decides to break off calibration or experiment
        [keyIsDown,secs,keyCode] = KbCheck;
        if keyCode(ivx.modifierKey) && keyCode(ivx.breakKey),
            % break off experiment
            result=-1;
            ivx.stop=1;

            iViewX('break', ivx);

            return;
        end
        if keyCode(ivx.breakKey), % break off calibration
            result=-1;
            iViewX('break', ivx);
            return;
        end
        if ivx.dummymode==1 && keyCode(ivx.modifierKey) && keyCode(ivx.nextCalPointKey), % break off calibration
            i=i+1;
            if i>ivx.nCalPoints;
                return
            else
                break; % does this get us out of the loop?
            end
        end
        % If you want the subject to manually accept the point with a
        % key press once they fixate, the ÒET_ACCÓ  command must then
        % be sent to the serial port.  This instructs iView to change
        % to the next calibration point  rather than having it detect
        % fixations itself.

        if keyCode(ivx.calAcceptKey),
            result=-1;
            iViewX('accept', ivx);
            % if this is not followed by a change command from iView
            % then we should increment the position ourselves
            % and break out of the loop
        end

    end


    [result, ivx]=iViewX('shutdown', ivx); % opens connection
    
    result


    fprintf('\nEnd of demo.\n');
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    pnet('closeall');
    Screen('CloseAll');    rethrow(lasterror);
end %try..catch..
