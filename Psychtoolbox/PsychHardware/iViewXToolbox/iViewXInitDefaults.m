function ivx=iviewxinitdefaults(window, points, host, port)

% set a number of iviewx default parameters
% we should be able to reset a number of parameters in a convenient way
% such as absCalPos, or calPacSecs
%

% 020105    fwc changed name to 'iviewxinitdefaults' for consistency
%               setting calibrationpoints is now done in seperate function

if exist('window', 'var') && ~isempty(window)
    ivx.window=window;
else
    ivx.window=[];
end

if exist('points', 'var') && ~isempty(points)
    ivx.nCalPoints=points;
else
    ivx.nCalPoints=9;
end

if exist('host', 'var') && ~isempty(host)
    ivx.host=host;
else
    % ivx.host='localhost';
    ivx.host='192.168.1.2';
end

if exist('port', 'var') && ~isempty(port)
    ivx.port=port;
else
    ivx.port=4444;
end

ivx.isconnected=1;
ivx.notconnected=0;
ivx.stop=0;
ivx.dummymode=0;
ivx.udp=[];
ivx.connection=1; % 0=no, 1=ethernet, 2=serial

% defaults for udp communication
ivx.socket=1111; % used when sending messages
ivx.localport=4444; % for receiving messages
ivx.udpreadtimeout=0.1; % time out for reading
ivx.udpmaxread=1000; % maximum number of elements to read in

% set defaults for calibration
ivx.orderedCalibration=0;
ivx.calAcceptManually=0;
ivx.calPace=1000;
ivx.calPaceSecs=ivx.calPace/1000;
ivx.calExtraTimeSecs=0; % time to wait before returning after calibration
ivx.calPointSize=3; % of screen
ivx.calPointSizeCenter=.5;
ivx.modifierKey=KbName('LeftGUI');
ivx.breakKey=KbName('ESCAPE');
ivx.calibrateKey=KbName('c');
ivx.validateKey=KbName('v');
ivx.driftcorrectKey=KbName('d');
ivx.outputKey=KbName('o');
ivx.calAcceptKey=KbName('space');
ivx.nextCalPointKey=KbName('p');
ivx.commTimeOut=120; % how long do we wait for data from iViewX host?

if isempty(ivx.window) % default to something
    h=1024;
    v=768;
else
    [h v]=WindowSize(ivx.window);
end

ivx.screenHSize=h;
ivx.screenVSize=v;


offset=[0 0]; % default no offset
scale=[0 0];    % default standard scaling

ivx=iViewXSetCalPoints(ivx, ivx.nCalPoints, offset, scale);

% switch(ivx.nCalPoints)
%     case 1,
%         ivx.relCalPos=[ 50 50]; % 1 point in % of screen
%     case 3,
%         ivx.relCalPos=[ 5 50; 50 50; 95 50]; % 3 points in % of screen
%     case 5,
%         ivx.relCalPos=[ 5 5; 95 5; 50 50; 5 95; 95 95]; % 5 points in % of screen
%     case 9,
%         ivx.relCalPos=[ 5 5; 50 5; 95 5; 5 50; 50 50; 95 50; 5 95; 50 95; 95 95]; % 9 points in % of screen
%     case 13,
%         ivx.relCalPos=[ 5 5; 50 5; 95 5; 5 50; 50 50; 95 50; 5 95; 50 95; 95 95]; % 9 points in % of screen
%         ivx.relCalPos=[ ivx.relCalPos; 33.3 25; 66.7 25; 33.3 75; 66.7 75]; % 13 points in % of screen
%     otherwise,
%         ivx.relCalPos=[ 5 5; 50 5; 95 5; 5 50; 50 50; 95 50; 5 95; 50 95; 95 95]; % 9 points in % of screen
% end
% 


ivx.relCalPos=ivx.relCalPos/100;

% set default point for driftcorrection
ivx.relDCPos=[ 50 50]; % 1 point in % of screen
ivx.relDCPos=ivx.relDCPos/100;



if 1
%     ivx.absCalPos(:,1)=ivx.relCalPos(:,1)*h;
%     ivx.absCalPos(:,2)=ivx.relCalPos(:,2)*v;
%     ivx.absCalPos=round(ivx.absCalPos);
% 
    ivx.absDCPos(:,1)=ivx.relDCPos(:,1)*h;
    ivx.absDCPos(:,2)=ivx.relDCPos(:,2)*v;
    ivx.absDCPos=round(ivx.absDCPos);

    calPointSize=round(ivx.calPointSize/100*v);
    calPointSize=calPointSize+(1-mod(calPointSize,2));
    if calPointSize<3 calPointSize=3; end

    ivx.calRect=[0 0 calPointSize calPointSize];
    calPointSize=round(ivx.calPointSizeCenter/100*v);
    calPointSize=calPointSize+(1-mod(calPointSize,2));
    if calPointSize<1 calPointSize=1; end
    ivx.calRectCenter=[0 0 calPointSize calPointSize];
end

if ~isempty(ivx.window)
    ivx.calPointColour=WhiteIndex(ivx.window);
    ivx.calPointColourCenter=BlackIndex(ivx.window);
    ivx.backgroundColour=round((WhiteIndex(ivx.window)+BlackIndex(ivx.window))/2);
else
    ivx.calPointColour=[255 255 255];
    ivx.calPointColourCenter=[0 0 0];
    ivx.backgroundColour=[128 128 128];
end



