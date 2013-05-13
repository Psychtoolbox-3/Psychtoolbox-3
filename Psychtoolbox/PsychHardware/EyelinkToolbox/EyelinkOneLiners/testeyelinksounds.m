function err=testeyelinksounds(el)

% test for the presence of sounds used by eyelink toolbox
%
% USAGE: err=testeyelinksounds(el)
%
%		el: if not supplied, we use 'initeyelinkdefaults'
%			to open a structure with eyelink default values
% returns 0 when no errors, -1 if one or more are missing
%
err=0;

if eyelink( 'isconnected' ) == 0
	Eyelink('initializedummy');
	initializedhere=1;
else
	initializedhere=0;
end
% if no default values are supplied open
% a structure here
if nargin==0 || (nargin==1 && isempty(el) )
	el=initeyelinkdefaults;
end

% test if the default sounds are available
% by playing them
error1=caltargetbeep(el);
SND('Wait');
if error1~=0
	fprintf( '%s appears to be missing.\n', el.targetdisplaysound );
end
error2=caldonebeep(el, 0);
SND('Wait');
if error2~=0
	fprintf( '%s appears to be missing.\n', el.calibrationfailedsound );
end
error3=caldonebeep(el, 1);
SND('Wait');
if error3~=0
	fprintf( '%s appears to be missing.\n', el.calibrationsuccesssound );
end

if error1~=0 || error2~=0 || error3~=0
	fprintf( '\nOne or more of the sounds\n');
	fprintf( '(specified in "initeyelinkdefaults.m")\n');
	fprintf( 'appear to be missing.\n');
	fprintf( 'Please put these sounds in the system folder.\n\n');
	err=-1;
else
	% fprintf( '\nAll sounds required by EyelinkToolbox are present.\n\n');
end
if initializedhere==1;
	Eyelink('shutdown');
end
