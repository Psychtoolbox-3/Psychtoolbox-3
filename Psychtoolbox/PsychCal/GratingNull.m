function foreWeight=GratingNull(device,message,forePart,blendPeriodPix,foreColor,backColor,distanceM,dpi)
% foreWeight=GratingNull(device,message,forePart,blendPeriodPix,[foreColor],[backColor],[viewingDistanceM],[dpi])
% GratingNull obtains visual matches that allow calibration of a video monitor's gamma.
% Finds the "matchColor" (dac triplet) that produces a visual match to a two-color optical mixture
% of foreColor and backColor of which the proportion forePart is the 
% "foreground" color. The matchColor is a weighted average, with the weight adjusted by the viewer,
% of foreColor and backColor. 
% matchColor=foreWeight*foreColor+(1-foreWeight)*backColor;
% The program returns foreWeight.
% The blend grating is very fine (normally 60 c/deg).
% The spatial periodicity of the grating for the visual match is
% optimal (i.e. 3 c/deg) at the specified viewing distance.
% If dpi is omitted it is assumed to be 66.
% If viewing distance is omitted it is assumed to be the distance at which
% the mixture spatial frequency is 60 c/deg.
%
% Denis Pelli 6/1/97

% 5/21/96 dgp wrote it
% 5/26/96 dgp enhanced to allow arbitrary foreColor and backColor
% 5/28/96 dgp updated to use new GetMouse, that flushes mouse events
% 6/1/97  dgp updated to use FlushEvents and CopyWindows, which is MUCH faster.
% 6/1/97  dgp Added 'message' argument, to allow indication of progress, e.g. '1 of 6'.
% 8/14/97 dhb Make call to ischar conditional on Version 5, since it fails in version 4. 
% 8/16/97	dgp	Changed "text" to "theText" to avoid conflict with TEXT function.
% 4/10/98 dhb Change TIMER to WaitSecs.
% 7/14/99 dgp At least some parts of VisualGammaDemo assume 8-bit pixelsize, so i'm forcing it here.
% 4/06/02 awi  -Check all elements of the new multi-element button vector returned   
%               by GetMouse on Windows.
%              -Replace Chicago font with Arial because it's available on both Mac and Windows

PsychDefaultSetup(1);

if sscanf(version,'%f',1)<5
	if nargin<4 || nargin>8 || nargout>1 || ~ischar(message)
		error('Usage:	foreWeight=GratingNull(device,message,forePart,blendPeriodPix,[foreColor],[backColor],[viewingDistanceM],[dpi])');
	end
else
	if nargin<4 || nargin>8 || nargout>1 || ~ischar(message)
		error('Usage:	foreWeight=GratingNull(device,message,forePart,blendPeriodPix,[foreColor],[backColor],[viewingDistanceM],[dpi])');
	end
end
if nargin<8
	dpi=66;
end
if nargin<7
	distanceM=0.0254*blendPeriodPix*57*dpi/60;
end
if nargin<6
	backColor=[0 0 0];
end
if nargin<5
	foreColor=[255 255 255];
end
pixelDeg=57/(dpi*distanceM/0.0254);
% 7/14/99 dgp At least some parts of VisualGammaDemo assume 8-bit pixelsize, so i'm forcing it here.
[window,screenRect]=Screen(device,'OpenWindow');
white=[255 255 255];
black=[0 0 0];
clut=(0:255)';
clut=clut(:,[1 1 1]);
clut(1,:)=white;
clut(256,:)=black;
clut(2,:)=(foreColor+backColor)/2;
clut(3,:)=foreColor;
clut(4,:)=backColor;
Screen(window,'SetClut',clut);

% Create a 3 c/deg squarewave grating in which even bars are a blend of
% foreColor and backColor, and
% odd bars are uniform and viewer-adjustable, to match the even bars.
% The spatial frequency of the grating (3 c/deg) is chosen to 
% optimize visual contrast sensitivity.
% The blend is produced at the highest possible spatial 
% frequency, alternating lines vertically, rather than pixels
% horizontally, to minimize retinal contrast, while staying well
% below the display's video bandwidth.
% CLUT entries: 1=adjustable, 2=foreColor, 3=backColor
barWidth=max(1,round(0.5/(3*pixelDeg))); % 3 c/deg grating
testRect=ScaleRect(screenRect,0.5,1);
testRect=round(AlignRect(testRect,screenRect,RectLeft));
barRect=SetRect(0,0,barWidth,RectHeight(testRect));
blend=1:RectHeight(barRect);
blend=blend-blendPeriodPix*floor(blend/blendPeriodPix); % modulo the period
blend=2+(blend >= round(forePart*blendPeriodPix));
blend=Expand(blend',RectWidth(barRect),1);
pure=ones(RectHeight(barRect),RectWidth(barRect));
barRect=AlignRect(barRect,testRect,RectLeft,RectBottom);
Screen(window,'PutImage',blend,barRect);
barRect2=OffsetRect(barRect,RectWidth(barRect),0);
Screen(window,'PutImage',pure,barRect2);
barRect=UnionRect(barRect,barRect2);
barRect2=OffsetRect(barRect,RectWidth(barRect),0);
while RectWidth(barRect)<RectWidth(screenRect)
	Screen('CopyWindow',window,window,barRect,barRect2);
	barRect=UnionRect(barRect,barRect2);
	barRect2=OffsetRect(barRect,RectWidth(barRect),0);
end
Screen(window,'FillRect',0,OffsetRect(testRect,RectWidth(testRect),0))

% Print instructions
theText=message;
theText=char(theText,sprintf('View from %.1f meters.',distanceM));
theText=char(theText,'Move the mouse up and');
theText=char(theText,'down to null out the');
theText=char(theText,'grating. Click when the');
theText=char(theText,'screen appears uniform.');
s=24;
Screen(window,'TextFont','Arial');
Screen(window,'TextSize',s);
s=s+8;
textRect=SetRect(0,0,Screen(window,'TextWidth',theText(2,:)),size(theText,1)*s);
textRect=CenterRect(textRect,screenRect);
textRect=OffsetRect(textRect,RectWidth(screenRect)/4+20/4,0);
for i=1:size(theText,1)
	Screen(window,'DrawText',theText(i,:),textRect(RectLeft),textRect(RectTop)+s*i,255);
end

% animate
% track vertical mouse position with vertical slider knob.
sliderRect=SetRect(0,0,20,RectHeight(screenRect));
sliderRect=CenterRect(sliderRect,screenRect);
knobRect=SetRect(0,0,RectWidth(sliderRect),RectWidth(sliderRect));
knobRect=InsetRect(CenterRect(knobRect,sliderRect),1,0);
top=RectTop;
bottom=RectBottom;
Screen(window,'FillRect',0,sliderRect);
Screen(window,'FrameRect',255,sliderRect);
Screen('Flip', window, [], 1);
while 1
	[x,y,button]=GetMouse;
	foreWeight=(sliderRect(bottom)-y)/RectHeight(sliderRect);
	Screen(window,'SetClut',foreWeight*foreColor+(1-foreWeight)*backColor,1);
	dy=y-(knobRect(top)+knobRect(bottom))/2;
	residue=knobRect;
	if dy>0
		residue(bottom)=residue(top)+dy;
	else
		residue(top)=residue(bottom)+dy;
	end
	knobRect=OffsetRect(knobRect,0,dy);
	Screen(window,'FillRect',0,residue);
	Screen(window,'FillRect',255,knobRect);
    Screen('Flip', window, [], 1);

	if any(button)break;end;
	WaitSecs(.01); % make sure we miss some frames, so mouse gets updated
end
Screen('CloseAll');
