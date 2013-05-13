function [result, ivx]=iViewXCalibrate(ivx)

% calibration routine of iViewX toolbox
% could also double as driftcorrection with a single point?
result=1;

try
    if ~exist('ivx', 'var') || isempty(ivx)
        txt=[ mfilename ' requires a structure with iViewX default values as input.'];
        error(txt);
    end

    if ~isfield(ivx, 'absCalPos') || isempty(ivx.absCalPos) || size(ivx.absCalPos,1)<1 || size(ivx.absCalPos,2)~=2
        txt=[ mfilename ': something is wrong with calibration point pixel values.'];
        error(txt);
    end

%     Screen('FillRect',ivx.window, ivx.backgroundColour);
%     Screen('Flip',ivx.window);

    iViewXEraseBackground(ivx);
    
    fprintf('Calibrating using %d points.\n', size(ivx.absCalPos,1));

    pointShown=zeros(1,size(ivx.absCalPos,1));
    nPointsShown=0;

    i=0;
    while 1
        [keyIsDown,secs,keyCode] = KbCheck;
        if keyCode(ivx.modifierKey) && keyCode(ivx.breakKey),
            result=-1;
            ivx.stop=1;
            break;
        end
        if keyCode(ivx.breakKey),
            result=-1;
            break;
        end
        % determine next point
        if ivx.orderedCalibration==1
            i=i+1;
        else
            while 1
                k=randperm(size(ivx.absCalPos,1));
                if pointShown(k(1))==0
                    i=k(1);
                    break;
                end
                drawnow; % enables ctrl-c break?
            end
        end
        iViewXDrawCalibrationPoint(ivx, ivx.absCalPos(i,1), ivx.absCalPos(i,2));
        pointShown(i)=1;
        nPointsShown=nPointsShown+1;

        % here we should probably send some information to the tracker about the
        % present point
        tEnd=GetSecs+ivx.calPaceSecs;
        while GetSecs<tEnd
            [keyIsDown,secs,keyCode] = KbCheck;
            if keyCode(ivx.modifierKey) && keyCode(ivx.breakKey),
                ivx.stop=1;
                result=-1;
                break;
            end
            if keyCode(ivx.breakKey),
                result=-1;
                break;
            end
        end

        if nPointsShown==size(ivx.absCalPos,1)
            break;
        end
    end
    iViewXEraseBackground(ivx);
%     Screen('FillRect',ivx.window, ivx.backgroundColour);
%     Screen('Flip',ivx.window);
    waitSecs(ivx.calExtraTimeSecs);

    fprintf('\nEnd of calibration.\n');
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.

    Screen('CloseAll');
    rethrow(lasterror);
end %try..catch..


