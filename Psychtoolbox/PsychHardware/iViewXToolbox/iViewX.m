function [result, ivx]=iViewX(cstr, ivx, params)

% USAGE: [result, ivx]=iViewX(cstr, ivx [, params])
%
% iViewX requires as input:
% 1. a command string
% 2. a structure with iViewX default values.
% cstr:  command to execute
% ivx: structure holding default information
% params: optional parameters to pass on with the command
%
% returns:
% result: any result produced by the command
% ivx: structure with default information (that may have been modified)

% History
% 011204    fwc added additional commands based on manual
% 020105    fwc added sendcalibrationpoints

result=-1;

if ~exist('cstr', 'var') || isempty(cstr)
    %     txt=[ mfilename ' requires a command string as input.'];
    %     error(txt);
    help iViewX
    return
end

switch lower(cstr)
    case {'help', '?'},
        help iViewX
        return
end

if ~exist('ivx', 'var') || isempty(ivx)
    txt=[ mfilename ' requires a structure with iViewX default values as input.'];
    error(txt);
end

if ~exist('params', 'var') || isempty(params)
    params=[];
end

switch lower(cstr)
    case 'command',
        % --> check that input is a text string
        % send string as command
        result=iViewXComm('send', ivx, params);
        %     case 'initialize',
        %         [result, ivx]=iViewXComm('open', ivx);
    % connecting
    case 'openconnection',
        [result, ivx]=iViewXComm('open', ivx);
    case 'closeconnection',
        [result, ivx]=iViewXComm('close', ivx);
    case 'checkconnection',
        [result, ivx]=iViewXComm('check', ivx);
    case 'receivedata',
        % check port for commands or gaze data from iView
        [result, ivx]=iViewXComm('receive', ivx);   
    case 'datafile',
        % saves current data buffer content to a datafile. Data is saved in
        % binary format. Use the IDF Converter to convert to ASCII format. 
        % If the given filename already exists, the existing
        % file will not be overwritten and the command is not executed.
        % params should contain a datafile
        result=iViewXComm('send', ivx, ['ET_SAV ' params]);
    % calibration
    case 'calibrate',
        [result, ivx]=iViewXCalibrate(ivx);
    case {'cancelcalibration', 'break'},
        result=iViewXComm('send', ivx, 'ET_BRK');
    case 'defaultcalibrationpoints',
        result=iViewXComm('send', ivx, 'ET_DEF'); % restores default calibration positions
    case 'sendcalibrationpoints',
        % we could check for the presence of the ivx.absCalPos array
        for i=1:ivx.nCalPoints
            message=['ET_PNT ' num2str(i) ' ' num2str(ivx.absCalPos(i,1)) ' ' num2str(ivx.absCalPos(i,2))];
%             message
            result=iViewXComm('send', ivx, message);
        end
    case 'setpoint',
        % params should contain three digits, 
        % point number, x and y locations
        result=iViewXComm('send', ivx, ['ET_PNT ' num2str(params(1)) ' ' num2str(params(2)) ' ' num2str(params(3))]);

    case 'setcalibration',
        % send calibration type
        % check params contains a single valid digit
        result=iViewXComm('send', ivx, ['ET_CAL ' num2str(params)]);
        %     case 'getcalibrationpoints',
        %         [result, ivx]=iViewXGetCalibrationPoints(ivx); 
    case {'setscreensize', 'setcalibrationarea'},
        % send screen size
        if ~isempty(ivx.window)
            [ivx.screenHSize, ivx.screenVSize]=WindowSize(ivx.window);
        end
        message=['ET_CSZ ' num2str(ivx.screenHSize) ' ' num2str(ivx.screenVSize) ];
        % check params contains a single valid digit
        result=iViewXComm('send', ivx, message);
    case 'accept',
        result=iViewXComm('send', ivx, 'ET_ACC');
    case 'setchecklevel',
        % params should contain a single digit, the new level
        % runs from 0 (weak) to 3 (strong)
        result=iViewXComm('send', ivx, ['ET_LEV ' num2str(params)]);
    case 'driftcorrection',
        result=iViewXComm('send', ivx, 'ET_RCL');
        % this must require more, probably a routine similar to calibration
%         % params can hold position for drift correction
%         [result, ivx]=iViewXDriftCorrection(ivx, params);
    case 'validation',
        fprintf('%s: Sorry, command ''%s'' not implemented yet.\n', mfilename, command );
% recording
    case 'startrecording',
        result=iViewXComm('send', ivx, 'ET_REC');
    case 'pauserecording',
        result=iViewXComm('send', ivx, 'ET_PSE');
    case 'continuerecording',
        result=iViewXComm('send', ivx, 'ET_CNT');
    case 'stoprecording',
        result=iViewXComm('send', ivx, 'ET_STP');
    case 'clearbuffer',
        result=iViewXComm('send', ivx, 'ET_CLR');
    case 'incrementsetnumber',
        % params should contain a single digit, the new set number
        result=iViewXComm('send', ivx, ['ET_INC ' num2str(params)]); 
% streaming
    case 'datastreamingon',
        result=iViewXComm('send', ivx, 'ET_STR');
        % we could/should include ability to specify subsampling
    case 'datastreamingoff',
        result=iViewXComm('send', ivx, 'ET_EST');
    case 'setdataformat',
        % set data format for streaming output
        result=iViewXComm('send', ivx, ['ET_FRM ' params ]);
% drawing
    case 'erasescreen',
        iViewXEraseScreen(ivx);
    case 'drawcalibrationpoint',
        % check params should contain two valid ints
        % that are within the screen bounds
        iViewXDrawCalibrationPoint(ivx, params);
 % miscellaneous
    case 'sendconfiguration',
        result=iViewXComm('send', ivx, 'ET_CFG');
        % not very useful in this manner, as we need to read it out of
        % course
    % message
    case 'message',
        % --> check that input is single text string
        result=iViewXComm('send', ivx, ['ET_REM ' params ]);
    case 'auxdata',
        % --> check that input is single text string
        result=iViewXComm('send', ivx, ['ET_AUX ' params ]);
    case 'loadbitmap',
        % --> check that input is a file location on iviewx
        result=iViewXComm('send', ivx, ['ET_BMP ' params ]);
otherwise,

        fprintf('%s: unknown command ''%s''\n', mfilename, command );
end

%---------------------------

% function result=iViewXSend(ivx, sendstr)
%
% switch ivx.connection
%     case 1, % ethernet
%
%         % open udp connection, and send command string
%         % if this is too time consuming, we should open a port and keep it open
%         % for the time of the experiment
%         % code stolen from udp_send_demo.m
%         result=-1;
%         udp=pnet('udpsocket',ivx.socket);
%         if udp~=-1,
%             try, % Failsafe
%                 pnet(udp,'udpconnect',ivx.host,ivx.port);
%                 %         [ip,port]=pnet(udp,'gethost')
%                 %         stat=pnet(udp,'status')
%                 pnet(udp,'write',[sendstr char(10)]);        % Write to write buffer
%                 pnet(udp,'writepacket',ivx.host,ivx.port);   % Send buffer as UDP packet
%
%             catch,
%                 pnet('closeall');
%                 disp(lasterr)
%             end
%             %     [ip,port]=pnet(udp,'gethost')
%             %     stat=pnet(udp,'status')
%             pnet(udp,'close');
%             result=1;
%         end
%
%     case 2,
%         % assume we use serial communication
%         error([mfilename ': serial communication not implemented yet']);
%
%
%     otherwise,
%         error([mfilename ': unknown connection type']);
% end
%

% ET_ACC  accepts calibration point
% ET_BRK  cancels calibration
% ET_CAL  starts calibration
% ET_CHG  indicates calibration point change
% ET_CSZ  sets size of calibration area
% ET_DEF  resets calibration points to default positions
% ET_FIN  indicates end of calibration
% ET_LEV  sets check level for calibration
% ET_PNT  sets position of calibration point
% ET_RCL  starts drift correction





