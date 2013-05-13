function ivx=iviewxsetcalpoints(ivx, points, offset, scale)

% set a number of iviewx default parameters

if ~exist('ivx', 'var') || isempty(ivx)
    error('iviewxsetcalpoints: please specify ''ivx'' structure holding defaults');
    return;
end


if exist('points', 'var') && ~isempty(points)
    ivx.nCalPoints=points;
elseif ~isfield(ivx, 'nCalPoints') || isempty(ivx.nCalPoints)
    ivx.nCalPoints=9;
end

if exist('offset', 'var') && ~isempty(offset)
    ivx.calPointOffset=offset;
else
    ivx.calPointOffset=[0 0];
end


if exist('scale', 'var') && ~isempty(scale)
    ivx.calPointScale=scale;
else
    ivx.calPointScale=[100 100];
end


ivx.calPointInset=5;

dt=ivx.calPointInset;
switch(ivx.nCalPoints)
    case 1,
        ivx.relCalPos=[ 50 50]; % 1 point in % of screen
    case 3,
        ivx.relCalPos=[ 50 50; 5 50; 95 50]; % 3 points in % of screen
    case 5,
        ivx.relCalPos=[ 50 50; 5 5; 95 5; 5 95; 95 95]; % 5 points in % of screen
    case 9,
        ivx.relCalPos=[ 50 50; 5 5; 95 5; 5 95; 95 95; 50 5; 5 50; 95 50; 50 95]; % 9 points in % of screen
    case 13,
        ivx.relCalPos=[ 50 50; 5 5; 50 5; 95 5; 5 50; 95 50; 5 95; 50 95; 95 95]; % 9 points in % of screen
        ivx.relCalPos=[ ivx.relCalPos; 33.3 25; 66.7 25; 33.3 75; 66.7 75]; % 13 points in % of screen
    otherwise,
        ivx.relCalPos=[ 50 50; 5 5; 50 5; 95 5; 5 50; 95 50; 5 95; 50 95; 95 95]; % 9 points in % of screen
        ivx.nCalPoints=length(ivx.relCalPos);
end


if isempty(ivx.window) % default to something
    h=1024;
    v=768;
else
    [h v]=WindowSize(ivx.window);
end
ivx.screenHSize=h;
ivx.screenVSize=v;

% center around 0,0
ivx.relCalPos=ivx.relCalPos-50;
% scale relative to center
ivx.relCalPos(:,1)=ivx.relCalPos(:,1)*(ivx.calPointScale(:,1)/100);
ivx.relCalPos(:,2)=ivx.relCalPos(:,2)*(ivx.calPointScale(:,2)/100);

% apply any offset
ivx.relCalPos(:,1)=ivx.relCalPos(:,1)+50+ivx.calPointOffset(:,1);
ivx.relCalPos(:,2)=ivx.relCalPos(:,2)+50+ivx.calPointOffset(:,2);

ivx.absCalPos=[]; % clear any existing points
ivx.relCalPos=ivx.relCalPos/100;
% if 1
ivx.absCalPos(:,1)=ivx.relCalPos(:,1)*h;
ivx.absCalPos(:,2)=ivx.relCalPos(:,2)*v;
ivx.absCalPos=round(ivx.absCalPos);

%     ivx.absDCPos(:,1)=ivx.relDCPos(:,1)*h;
%     ivx.absDCPos(:,2)=ivx.relDCPos(:,2)*v;
%     ivx.absDCPos=round(ivx.absDCPos);

%     calPointSize=round(ivx.calPointSize/100*v);
%     calPointSize=calPointSize+(1-mod(calPointSize,2));
%     if calPointSize<3 calPointSize=3; end
%
%     ivx.calRect=[0 0 calPointSize calPointSize];
%     calPointSize=round(ivx.calPointSizeCenter/100*v);
%     calPointSize=calPointSize+(1-mod(calPointSize,2));
%     if calPointSize<1 calPointSize=1; end
%     ivx.calRectCenter=[0 0 calPointSize calPointSize];
% end
