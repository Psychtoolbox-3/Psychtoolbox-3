function UserResponse = TwoStateQuery(TheQuestion,TheChoices)
% Syntax: UserResponse = TwoStateQuery(TheQuestion,[TheChoices])
% 
% Purpose: Create simple dialog box asking user an either/or question
%
% History: 5/5/04		mpr		decided whether or not to celebrate cinco de mayo
%					 10/13/04	mpr		set Yes Button automatically to be enlarged for large
%                           strings 
%						1/21/05	mpr		allowed fontsize to shrink for longer questions
%					  1/27/05	mpr		allowed question to break at file separator instead of
%                           just spaces 
%						6/30/05	mpr		forced underscores to remain as underscores instead of
%                           indicating subscript
%						3/27/06	mpr		changed returned value to -1 when user deletes window;
%                           note that this is dangerous because how a closed
%                           window affects subsequent program behavior is not
%                           implemented the same way in all code calling this
%                           function... 
%           5/3/06  mpr   set menubar to none for Version 2006a
%           8/25/06 mpr   resized No Button as needed; could be done better, but  
%                           the plan is to keep the 'Yes' and 'No' strings short
%                           so that the code doesn't need to be pushed hard
%           9/29/06 mpr   repainted buttons in version 7
%           1/14/07 mpr   undid previoius fix because button color bug is almost
%                             fixed
%          5/20/13    mk Add text only fallback for Octave and non-GUI.

if nargin < 2 || isempty(TheChoices)
	TheChoices{1} = 'Yes';
	TheChoices{2} = 'No';
elseif length(TheChoices) ~= 2 || ~iscell(TheChoices)
	warning('I don''t like the choices I received.  Am substituting unadulterated originals.');
	disp('But for kicks, this is what I was passed:');
	TheChoices
	TheChoices{1} = 'Yes';
	TheChoices{2} = 'No';
end	

UnderScoreSpots = find(TheQuestion == '_');
for k=length(UnderScoreSpots):-1:1
	TheQuestion = [TheQuestion(1:(UnderScoreSpots(k)-1)) '\' TheQuestion(UnderScoreSpots(k):end)];
end

% Provide text fallback for non-GUI mode or Octave:
if ~IsGUI || IsOctave
  fprintf('%s\n', TheQuestion);
  fprintf('0 = %s\n', TheChoices{2});
  fprintf('1 = %s\n', TheChoices{1});
  answer = input('Choose your answer: ','s');
  if isempty(answer) || (answer~='0' && answer~='1')
    UserResponse = -1;
  else
      UserResponse = sscanf(answer, '%i');
  end
  
  return;
end

FigPos = [300 600];
FigWidth = 500;
FigHeight = 200;
YesBPos = [125 30];
YesBWidth = 80;
YesBHeight = 50;
NoBPos = [295 30];

PrimeFontSize = 36;

QueryFigh = figure('Position',[FigPos FigWidth FigHeight],'Resize','off','NumberTitle','off', ...
									'Name','Cross the Rubicon','MenuBar','none','Color','k');
set(QueryFigh,'Tag',['Questionable Figure ' int2str(QueryFigh)]);
axes('Position',[0 0 1 1],'XLim',[0 FigWidth],'YLim',[0 FigHeight],'XColor','k','YColor','k');
th=text(FigWidth/2,0.65*FigHeight,TheQuestion,'FontSize',36,'Color',get_color('hot pink'), ...
				'HorizontalAlignment','Center');
TextExt = get(th,'Extent');
if TextExt(1) < 0
	SpaceLocs = find(TheQuestion == 32 | TheQuestion == abs(filesep));
	[MinDist, MinInd] = min(abs(floor(length(TheQuestion)/2)-SpaceLocs));
	if TheQuestion(SpaceLocs(MinInd)) == abs(filesep)
		set(th,'String',sprintf('%s\n%s',TheQuestion(1:SpaceLocs(MinInd)),TheQuestion((SpaceLocs(MinInd)+1):end)));
	else
		set(th,'String',sprintf('%s\n%s',TheQuestion(1:(SpaceLocs(MinInd)-1)),TheQuestion((SpaceLocs(MinInd)+1):end)));
	end
	while TextExt(1) < 0 && get(th,'FontSize') > 8
		set(th,'FontSize',get(th,'FontSize')-1);
		TextExt = get(th,'Extent');
		% I should create an error if font size does get down to 8, but I don't expect it ever to happen, so I'm not
		% going to deal with it.
	end
end
YesBh = uicontrol('Style','PushButton','Position',[YesBPos YesBWidth YesBHeight],'String',TheChoices{1}, ...
									'ForeGroundColor',get_color('midnight blue'),'BackgroundColor', ...
									get_color('wheat'),'FontSize',36,'Tag','YesButton','UserData',-1);
set(YesBh,'CallBack',['set(findobj(' int2str(QueryFigh) ',''Tag'',''YesButton''),''UserData'',1);']);
NoBWidth = YesBWidth;
NoBHeight = YesBHeight; 
NoBh = uicontrol('Style','PushButton','Position',[NoBPos YesBWidth YesBHeight],'String',TheChoices{2}, ...
                'ForeGroundColor',get_color('midnight blue'),'BackgroundColor', ...
                get_color('wheat'),'FontSize',36,'CallBack', ...
                ['set(findobj(' int2str(QueryFigh) ',''Tag'',''YesButton''),''UserData'',0);']);
YesBSSize = get(YesBh,'Extent');
while YesBSSize(3) > YesBWidth
	if YesBPos(1)+YesBWidth > NoBPos(1) - YesBWidth/2
		set(YesBh,'FontSize',get(YesBh,'FontSize')-12);
		if ~get(uih,'FontSize')
			error('Yes string is too big!');
		end
	else
		YesBPos(1) = YesBPos(1) - (YesBSSize(3)-YesBWidth)/2;
		YesBWidth = 1.1*YesBSSize(3);
		set(YesBh,'Position',[YesBPos YesBWidth YesBHeight]);
	end
	YesBSSize = get(YesBh,'Extent');
end

NoBSSize = get(NoBh,'Extent');
while NoBSSize(3) > NoBWidth
  if NoBWidth < YesBWidth
    NoBWidth = YesBWidth;
    set(NoBh,'Position',[NoBPos NoBWidth NoBHeight]);
  elseif NoBPos(1)+NoBWidth < 0.9*FigWidth
    NoBWidth = NoBWidth*1.1;
    set(NoBh,'Position',[NoBPos NoBWidth NoBHeight]);
  else
		set(NoBh,'FontSize',get(NoBh,'FontSize')-12);
		if ~get(NoBh,'FontSize')
			error('No string is too big!');
		end
	end
	NoBSSize = get(NoBh,'Extent');
end

drawnow;
waitfor(YesBh,'UserData');

% just in case user deletes figure rather than clicking button
if isempty(findobj('Tag',['Questionable Figure ' int2str(QueryFigh) ]))
	UserResponse = -1;
else
	UserResponse = get(YesBh,'UserData');
	delete(QueryFigh);
end

return;
