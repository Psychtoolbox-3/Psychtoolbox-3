function [result, ivx]=iViewXDriftCorrection(ivx, pos)

% driftcorrection routine of iViewX toolbox
% could also double as driftcorrection with a single point?
result=1;
absdcpos=[];
try

    if ~exist('ivx', 'var') || isempty(ivx)
        txt=[ mfilename ' requires a structure with iViewX default values as input.'];
        error(txt);
    end

    if 1==exist('pos', 'var') && ~isempty(pos)
        % if pos is supplied, we take the default value, but we first do
        % some checks
        if size(pos,1)==1 && size(pos,2)==2
            if 1==IsInRect(pos(1,1), pos(1,2), Screen('Rect', ivx.window))
                absdcpos=pos;
            else
                fprintf([mfilename ':supplied driftcorrection position is outside window, using default instead.\n']);
            end
        else
            fprintf([mfilename ':something is wrong with the supplied driftcorrection position, using default instead.\n']);
        end
    end

    if isempty(absdcpos)
        % we take default value if no value is found yet
        if ~isfield(ivx, 'absDCPos') || isempty(ivx.absDCPos) || size(ivx.absDCPos,1)~=1 || size(ivx.absCalPos,2)~=2
            txt=[ mfilename ': something is wrong with default driftcorrection point pixel values.'];
            error(txt);
        else
            absdcpos=ivx.absDCPos(1,:);
        end
    end

%     Screen('FillRect',ivx.window, ivx.backgroundColour);
%     Screen('Flip',ivx.window);
    iViewX('erasescreen', ivx);
    fprintf('Driftcorrection at position (%d, %d).\n', absdcpos(1,1), absdcpos(1,2));

    pointShown=zeros(1,size(absdcpos,1));
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
                k=randperm(size(absdcpos,1));
                if pointShown(k(1))==0
                    i=k(1);
                    break;
                end
                drawnow; % enables ctrl-c break?
            end
        end
        iViewX('drawcalibrationpoint', ivx, [absdcpos(i,1), absdcpos(i,2)]);
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

        if nPointsShown==size(absdcpos,1)
            break;
        end
    end
%     Screen('FillRect',ivx.window, ivx.backgroundColour);
%     Screen('Flip',ivx.window);
    iViewX('erasescreen', ivx);
    WaitSecs(ivx.calExtraTimeSecs);

    fprintf('\nEnd of driftcorrection.\n');
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.

    Screen('CloseAll');
    rethrow(lasterror);
end %try..catch..


