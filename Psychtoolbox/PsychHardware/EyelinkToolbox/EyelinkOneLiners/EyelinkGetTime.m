function eltime=EyelinkGetTime(el, maxwait)
% USAGE: eltime=EyelinkGetTime(el [, maxwait])
% Get tracker time
%
% el: is eyelinkdefault structure
% maxwait: specifies maximum time to wait for a tracker time event
% if omitted, will wait indefinetely

% history
% 28-10-02	fwc created it, based on 'EXGetEyeLinkTime.m' of jp
% 10-12-02	fwc	no longer quits on errors, but returns -1;
%				return getsecs time in dummymode
% 23-12-02	fwc	fixed small error in setting usemaxwait
% 17-02-03	fwc fixed serious bug due to what very often eltime would be set to 0
% 10-04-09  mk  Deuglified.

eltime=-1;
time=0;

if ~exist('maxwait', 'var')
    maxwait = [];
end

if isempty(maxwait) || maxwait <= 0
	usemaxwait=0; % do not timeout
else
	usemaxwait=1;
end

% get Eyelink time
start=GetSecs;
if Eyelink('isconnected') == el.connected   % readtime doesn't work in dummymode
	if Eyelink('requesttime') ~= 0			% request time over network
		fprintf('EyelinkGetTime: requesttime error');
		return;
	end;
	while time == 0
		if Eyelink('isconnected') == el.notconnected
			fprintf('EyelinkGetTime: lost connection during readtime request');
			return;
		end
	
		if (usemaxwait==1) && (GetSecs-start>maxwait) % repeat until available or timeout
			fprintf('EyelinkGetTime: timeout of time request');
			return;
		end
		time = Eyelink('readtime');	
	end
	eltime=time;
elseif Eyelink('isconnected') == el.dummyconnected
	eltime=GetSecs;
end
