function ReturnFigh = ConfirmInfo(TheQuestion,ButtonString,HowLongToWait)
% Syntax: [FigureHandle] = ConfirmInfo(TheQuestion,[ButtonString],[TimeoutPeriod])
% 
% Purpose: Create simple dialog box asking user to verify receipt of information
%
% History: 1/11/05		mpr		modified this from TwoStateQuery
%					 1/21/05		mpr		allowed font size to be reduced for longer strings
%           5/3/06    mpr   added third argument to limit time for response;
%                           defaults to Inf; works only in version 7 (requires
%                           timer object class)
%					 5/15/06		mpr		allowed figure to be non-modal; under that
%                           condition, clicking button deletes figure; to make
%                           nonmodal, in version 5, pass non-zero non-Inf third
%                           argument; in version 7, pass NaN as third argument;
%                           also added return argument to allow calling function
%                           to delete figure if user fails to delete figure and
%                           program proceeds 
%          8/26/06    mpr   repainted buttons in version 7
%         12/14/07    mpr   commented call to RepaintAllPushButtons because
%                           the Matlab bug is now almost fixed 
%         12/29/07    mpr   made three-line questions possible
%					 1/3/08			mpr		introduced NLChar for cross-platform compatibility
%          3/5/08     mpr   small cosmetic fix for new case where a third line
%                           was needed even though Extent(3) < 2...
%          3/21/08    mpr   allowed return value to signal if user closed window
%                           instead of clicking uicontrol button (works only if
%                           figure is modal)
%          3/31/08    mpr   added check to make sure search for appropriate
%                           space does not overrun number of spaces
%          5/20/13    mk Add text only fallback for Octave and non-GUI.

DoModal = 1;
if nargin < 3 || isempty(HowLongToWait)
  HowLongToWait = Inf;
end

% Provide text fallback for non-GUI mode or Octave:
if ~IsGUI || IsOctave
  fprintf('%s\nPress any key to continue.\n', TheQuestion);
  [secs, keyCode] = KbStrokeWait (-1, GetSecs + HowLongToWait);
  if any(keyCode)
    ReturnFigh = 1;
  else
    ReturnFigh = 0;
  end
  return;
end

ThisVersion = version;
if str2num(ThisVersion(1)) < 7
	if ~isinf(HowLongToWait)
		DoModal = 0;
	end
  HowLongToWait = Inf;
elseif isnan(HowLongToWait)
	DoModal = 0;
	HowLongToWait = Inf;
end

if nargin < 2 || isempty(ButtonString)
	ButtonString = 'Okay';
end	

NLChar = abs(sprintf('\n'));

FigPos = [300 600];
FigWidth = 500;
FigHeight = 200;
OkayBPos = [210 25];
OkayBWidth = 80;
OkayBHeight = 50;

PrimeFontSize = 36;

QueryFigh = figure('Position',[FigPos FigWidth FigHeight],'Resize','off','NumberTitle','off', ...
									'Name','Acknowledge','Menubar','none','Color','k');
set(QueryFigh,'Tag',['Confirmation Figure ' int2str(QueryFigh)]);
axes('Position',[0 0 1 1],'XLim',[0 FigWidth],'YLim',[0 FigHeight],'XColor','k','YColor','k');
th=text(FigWidth/2,0.7*FigHeight,TheQuestion,'FontSize',PrimeFontSize,'Color',get_color('hot pink'), ...
				'HorizontalAlignment','Center');
set(th,'Units','normalized');
TextExt = get(th,'Extent');

% There must be better ways to do what I do in this next condition...  At
% first this function was used only with short questions... questions that
% fit on one line.  When I started using it with longer questions, I
% started reducing the font size.  Then I decided I wanted to wrap lines so
% that I could put longer questions with larger font sizes.  Then I found
% there were times I wanted to wrap again because two lines required a font
% size too small...  The basic idea is to add carriage returns and shrink
% fonts to get the text to fit in a box of size pre-defined by the figure
% size.  It would probably be better at times to manipulate the size of the
% box as well as font sizes and number of lines.  But I'm not trying to
% accommodate all strings... I don't want to create a word processor!  I
% just want to get the aesthetics acceptable for questions of lengths within
% a reasonably restricted range.  What follows is ugly and inelegant, but
% it suits my needs.  Good enough to pass a cursory inspection anyways.
% Although in some versions of Matlab, the visibility box doesn't extend all the
% way out to the edges and bits of text at the margins are clipped.  Doesn't
% seem to be anything I can do about that except add blank spaces, but I don't
% know where the margins will be a priori and the problem doesn't bother me
% enough to make me want to write a fix...
% -- mpr

if TextExt(1) < 0
	SpaceLocs = find(TheQuestion == 32);
  
  if TextExt(3) > 3
    while TextExt(3) > 3 && get(th,'FontSize') > 8
      set(th,'FontSize',get(th,'FontSize')-1);
      TextExt = get(th,'Extent');
    end
    set(th,'String',TheQuestion(1:(SpaceLocs(1)-1)));
    SpaceNum = 1;
    PTextExt = get(th,'Extent');
    while PTextExt(3) < 0.95
      SpaceNum = SpaceNum+1;
      set(th,'String',TheQuestion(1:(SpaceLocs(SpaceNum)-1)));
      PTextExt = get(th,'Extent');
    end
    SpaceNum = SpaceNum-1;
    TheQuestion(SpaceLocs(SpaceNum)) = NLChar;
    SpaceEnd = SpaceNum;
    set(th,'String',TheQuestion((SpaceLocs(SpaceNum)+1):(SpaceLocs(SpaceEnd)-1)));
    PTextExt = get(th,'Extent');
    while PTextExt(3) < 0.95
      SpaceEnd = SpaceEnd+1;
      set(th,'String',TheQuestion((SpaceLocs(SpaceNum)+1):(SpaceLocs(SpaceEnd)-1)));
      PTextExt = get(th,'Extent');
    end
    TheQuestion(SpaceLocs(SpaceEnd-1)) = NLChar;
    set(th,'String',TheQuestion);
    TextExt = get(th,'Extent');
    while TextExt(3) > 0.95 && get(th,'FontSize') > 8
      set(th,'FontSize',get(th,'FontSize')-1);
      TextExt = get(th,'Extent');
    end
    if TextExt(3) > 1
      error('Question too large!');
    end
  elseif TextExt(3) > 2
    while TextExt(3) > 2 && get(th,'FontSize') > 8
      set(th,'FontSize',get(th,'FontSize')-1)
      TextExt = get(th,'Extent');
    end
    SpaceEnd = 1;
    set(th,'String',TheQuestion(1:(SpaceLocs(1)-1)));
    TextExt = get(th,'Extent');
    while TextExt(3) < 0.95
      SpaceEnd = SpaceEnd+1;
      if SpaceEnd > length(SpaceLocs)
        SpaceEnd = SpaceEnd-1;
        break;
      end
      set(th,'String',TheQuestion(1:(SpaceLocs(SpaceEnd)-1)));
      TextExt = get(th,'Extent');
    end
    TheQuestion(SpaceLocs(SpaceEnd-1)) = NLChar;
    set(th,'String',TheQuestion);
    TextExt = get(th,'Extent');
    while TextExt(3) > 0.95 && get(th,'FontSize') > 8
      set(th,'FontSize',get(th,'FontSize')-1);
      TextExt = get(th,'Extent');
    end
  elseif TextExt(3) > 0.95
    SpaceEnd = 1;
    set(th,'String',TheQuestion(1:(SpaceLocs(1)-1)));
    TextExt = get(th,'Extent');
    while TextExt(3) < 0.95
      SpaceEnd = SpaceEnd+1;
      if SpaceEnd > length(SpaceLocs)
        SpaceEnd = SpaceEnd-1;
        break;
      end
      set(th,'String',TheQuestion(1:(SpaceLocs(SpaceEnd)-1)));
      TextExt = get(th,'Extent');
    end
    TheQuestion(SpaceLocs(SpaceEnd-1)) = NLChar;
    set(th,'String',TheQuestion);
		TextExt = get(th,'Extent');
		if TextExt(3) > 0.95
			TheQuestion(SpaceLocs(end)) = NLChar;
			set(th,'String',TheQuestion);
			TextExt = get(th,'Extent');
			EndSpace = length(SpaceLocs);
			while TextExt(3) > 0.95 && EndSpace > SpaceEnd
				TheQuestion(SpaceLocs(EndSpace)) = ' ';
				EndSpace = EndSpace-1;
				TheQuestion(SpaceLocs(EndSpace)) = NLChar;
				set(th,'String',TheQuestion);
				TextExt = get(th,'Extent');
			end
		end
  end
end
OkayBh = uicontrol('Style','PushButton','Position',[OkayBPos OkayBWidth OkayBHeight],'String',ButtonString, ...
									'ForeGroundColor',get_color('midnight blue'),'BackgroundColor', ...
									get_color('wheat'),'FontSize',36,'Tag','OkayButton','UserData',-1);
set(OkayBh,'CallBack',['set(findobj(' int2str(QueryFigh) ',''Tag'',''OkayButton''),''UserData'',1);']);
OkayBSSize = get(OkayBh,'Extent');
while OkayBSSize(3) > OkayBWidth
	if OkayBPos(1)+OkayBWidth > FigWidth - OkayBWidth/2
		set(uih,'FontSize',get(uih,'FontSize')-12);
		if ~get(uih,'FontSize')
			error('Okay string is too big!');
		end
	else
		OkayBPos(1) = OkayBPos(1) - (OkayBSSize(3)-OkayBWidth)/2;
		OkayBWidth = 1.1*OkayBSSize(3);
		set(OkayBh,'Position',[OkayBPos OkayBWidth OkayBHeight]);
	end
	OkayBSSize = get(OkayBh,'Extent');
end

drawnow;
if ~isinf(HowLongToWait)
  QFh = QueryFigh; % Just to shorten callback lines...
  TimerHandle = timer('TimerFcn',['if ishandle(' int2str(QFh) '), ' ...
                                  '  if strcmp(get(' int2str(QFh) ',''Tag''),''Confirmation Figure ' int2str(QFh) '''), ' ...
                                  '    set(findobj(' int2str(QFh) ',''Tag'',''OkayButton''),''UserData'',1); ' ...
                                  '  end; ' ...
                                  'end; ' ...
                                  'stop(timerfind(''Tag'',''Timer for ' int2str(QFh) ''')); ' ...
                                  'delete(timerfind(''Tag'',''Timer for ' int2str(QFh) '''));'], ...
                      'StartDelay',HowLongToWait,'Tag',['Timer for ' int2str(QFh)],'Period',HowLongToWait);
  start(TimerHandle);
end
if DoModal
	waitfor(OkayBh,'UserData');
	
	if ~isinf(HowLongToWait)
	  TimerHandle = timerfind('Tag',['Timer for ' int2str(QFh)]);
	  if ~isempty(TimerHandle)
	    stop(TimerHandle);
	    delete(TimerHandle);
	  end
  end
	
  UserHitOkay = ishandle(OkayBh);
  
	delete(findobj('Tag',['Confirmation Figure ' int2str(QueryFigh)]));
else
  UserHitOkay = 1;
	set(OkayBh,'CallBack',['delete(' int2str(QueryFigh) ');']);
end

if nargout
  % Variable name only makes sense for modal figure; which is ironic because if
  % user hits okay, then QueryFigh is no longer a valid handle.  This logic will
  % work for me, though.  Even with a bit of illogic...  as long as the
  % functions calling this know what the return value really signifies which
  % should always be clear because a zero return makes no sense for a non-modal
  % window, and a non-zero value makes no sense for a modal window.  Aside from
  % signifying how the window was closed, that is...
  if UserHitOkay
    ReturnFigh = QueryFigh;
  else
    ReturnFigh = 0;
  end
end

return;
