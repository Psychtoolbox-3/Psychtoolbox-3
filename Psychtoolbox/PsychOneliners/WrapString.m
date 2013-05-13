function wrappedString=WrapString(string,maxLineLength)
% wrappedString=WrapString(string,[maxLineLength])
% 
% Wraps text by changing spaces into linebreaks '\n', making each line as
% long as possible without exceeding maxLineLength (default 74
% characters). WrapString does not break words, even if you have a word
% that exceeds maxLineLength. The returned "wrappedString" is identical to
% the supplied "string" except for the conversion of some spaces into
% linebreaks. Besides making the text look pretty, wrapping the text will
% make the printout narrow enough that it can be sent by email and
% received as sent, not made hard to read by mindless breaking of every
% line.
% 
% Note that this schemes is based on counting characters, not pixels, so
% it will give a fairly even right margin only for monospaced fonts, not
% proportionally spaced fonts. The more general solution would be based on
% counting pixels, not characters, using either Screen 'TextWidth' or
% TextBounds.

% 6/30/02 dgp Wrote it.
% 10/2/02 dgp Make it clear that maxLineLength is in characters, not pixels.
% 09/20/09 mk Improve argument handling as per suggestion of Peter April.

if nargin>2 || nargout>1 
	error('Usage: wrappedString=WrapString(string,[maxLineLength])\n');
end

if nargin<2
    maxLineLength=[];
end

if isempty(maxLineLength) || isnan(maxLineLength)
	maxLineLength=74;
end

eol=sprintf('\n');
wrapped='';
while length(string)>maxLineLength
	l=min([strfind(char(string),eol) length(string)+1]);
	if l<maxLineLength
		% line is already short enough
		[wrapped,string]=onewrap(wrapped,string,l);
	else
		s=strfind(char(string),' ');
		n=find(s<maxLineLength);
		if ~isempty(n)
			% ignore spaces before the furthest one before maxLineLength
			s=s(max(n):end);
		end
		% break at nearest space, linebreak, or end.
		s=sort([s l]);
		[wrapped,string]=onewrap(wrapped,string,s(1));
	end
end
wrappedString=[wrapped string];
return

function [wrapped,string]=onewrap(wrapped,string,n)
if n>length(string)
	wrapped=[wrapped string];
	string=[];
	return
end
wrapped=[wrapped string(1:n-1) sprintf('\n')];
string=string(n+1:end);
return
