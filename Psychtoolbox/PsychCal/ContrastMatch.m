function weight=ContrastMatch(device,dimWeight,foreColor,backColor)
% weight=ContrastMatch(device,dimWeight,[foreColor],[backColor])
%
% Displays two gratings. The bright grating alternates white and weight*white,
% producing luminances LMax and LBright.
% The dim grating alternates dimWeight*white and black, producing luminances
% LDim and LMin.
% When the grating contrasts match, LMax/LBright=LDim/LMin.
% In a more compact notation, L1/Lb=Ld/L0. We set the weight wd that produces
% Ld, and the observer adjusts the weight wb that produces Lb. They are related
% by the gamma function, y=g(w), where y is the normalized luminance,
% y=(L-L0)/(L1-L0). We know wb and wd, and we have access to the gamma function,
% so we can compute yb and yd. Solving for L, we have
%	L=y*(L1-L0)+L0
% Our relation at match can be rewritten as
%	L1*L0=Ld*Lb
%	0=Ld*Lb-L1*L0
%	=(yd*(L1-L0)+L0)*(yb*(L1-L0)+L0)-L1*L0
%	Let's divide by L1^2, and define r=L0/L1.
%	=(yd*(1-r)+r)*(yb*(1-r)+r)-r
%	This is a quadratic equation in r
%	=(yd+(1-yd)*r)*(yb+(1-yb)*r)-r
%	=yd*yb+(yb*(1-yd)+yd*(1-yb)-1)*r+(1-yd)*(1-yb)*r^2
%	c=[(1-yd)*(1-yb) yb*(1-yd)+yd*(1-yb)-1 yd*yb]; % coefficients of quadratic polynomial
%	r=roots(c)
% The answer, r, is the desired ratio of L0/L1.
%
% 5/28/96 dgp  Wrote it.
% 5/28/96 dgp  Updated to use new GetMouse, that flushes mouse events.
% 8/4/96  dhb  Changed name to ContrastMatch.
% 8/16/97 dgp  Changed "text" to "theText" to avoid conflict with TEXT function.
% 7/19/98 dgp  Removed obsolete TIMER.
% 6/30/03 dgp Updated Screen OpenScreen to Screen OpenWindow.

PsychDefaultSetup(1);

if nargin<2 || nargout>1
	error('Usage: weight=GratingMatch(device,dimWeight,[foreColor],[backColor])');
end
if nargin<4
	backColor=[0 0 0];
end
if nargin<3
	foreColor=[255 255 255];
end
dpi=67;
distanceM=0.57;
pixelDeg=57/(dpi*distanceM/0.0254);
[win, screenRect] = Screen('OpenWindow', device);
screenRect=reshape(screenRect,1,4);
white=[255 255 255];
black=[0 0 0];
clut=(0:255)';
clut=clut(:,[1 1 1]);
clut(1,:)=white;
clut(256,:)=black;
clut(2,:)=black;
clut(3,:)=white;
clut(4,:)=black;
clut(5,:)=dimWeight*white+(1-dimWeight)*black;
Screen('SetClut', win, clut);

if 0
	blend=1:RectHeight(barRect);
	blend=blend-blendPeriodPix*floor(blend/blendPeriodPix); % modulo the period
	blend=2+(blend >= round(forePart*blendPeriodPix));
	blend=Expand(blend',RectWidth(barRect),1);
	pure=ones(RectHeight(barRect),RectWidth(barRect));
	Screen('BlitImage255',blend,barRect);
	barRect=OffsetRect(barRect,RectWidth(barRect),0);
	Screen('BlitImage255',pure,barRect);
	barRect=OffsetRect(barRect,RectWidth(barRect),0);
end

% Create two 1 c/deg squarewave gratings at different luminances.
% CLUT entries: 1=adjustable, 2=foreColor, 3=backColor
barWidth=max(1,round(0.5/(1*pixelDeg))); % 1 c/deg grating
testRect=ScaleRect(screenRect,0.5,0.5);
testRect=round(AlignRect(testRect,screenRect,RectLeft,RectTop));
barRect=SetRect(0,0,barWidth,RectHeight(testRect));
for color=0:2:2
	barRect=AlignRect(barRect,testRect,RectLeft,RectBottom);
	for i=0:2:ceil(RectWidth(testRect)/RectWidth(barRect))
		Screen('FillRect', win, color+1,barRect);
		barRect=OffsetRect(barRect,RectWidth(barRect),0);
		Screen('FillRect', win, color+2,barRect);
		barRect=OffsetRect(barRect,RectWidth(barRect),0);
	end
	Screen('FillRect', win, 0,AdjoinRect(testRect,testRect,RectRight))
	testRect=AdjoinRect(testRect,testRect,RectBottom);
end

% Print instructions
theText=sprintf(     'Move the mouse up and');
theText=char(theText,'down to match the grating');
theText=char(theText,'contrasts. Click when ');
theText=char(theText,'you see one long grating,');
theText=char(theText,'partly occluded by a dark');
theText=char(theText,'filter.');
s=24;
%Screen('TextFont', win, 'Chicago');
Screen('TextSize', win, s);
s=s+8;
% textRect=SetRect(0,0,Screen('TextWidth',theText(2,:)),size(theText,1)*s);
% textRect=CenterRect(textRect,screenRect);
% textRect=OffsetRect(textRect,RectWidth(screenRect)/4+20/4,0);
% for i=1:size(theText,1)
% 	Screen('DrawText',textRect(RectLeft),textRect(RectTop)+s*i,255,theText(i,:));
% end
DrawFormattedText(win, theText, 0, 'center', 'center', white);

% animate
% track vertical mouse position with vertical slider knob.
sliderRect=SetRect(0,0,20,RectHeight(screenRect));
sliderRect=CenterRect(sliderRect,screenRect);
knobRect=SetRect(0,0,RectWidth(sliderRect),RectWidth(sliderRect));
knobRect=InsetRect(CenterRect(knobRect,sliderRect),1,0);
top=RectTop;
bottom=RectBottom;
Screen('FillRect', win, 0, sliderRect);
Screen('FrameRect', win, 255, sliderRect);
while 1
	[x,y,button]=GetMouse;
	weight=(sliderRect(bottom)-y)/RectHeight(sliderRect);
	Screen('SetClut', win, weight*foreColor+(1-weight)*backColor,1);
	dy=y-(knobRect(top)+knobRect(bottom))/2;
	residue=knobRect;
	if dy>0
		residue(bottom)=residue(top)+dy;
	else
		residue(top)=residue(bottom)+dy;
	end
	knobRect=OffsetRect(knobRect,0,dy);
	Screen('FillRect', win, 0,residue);
	Screen('FillRect', win, 255,knobRect);
	if(button) break; end;
	WaitSecs(.01); % make sure we miss some frames, so mouse gets updated
end
Screen('CloseAll');
