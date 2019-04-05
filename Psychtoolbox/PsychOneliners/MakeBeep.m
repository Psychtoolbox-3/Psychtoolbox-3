function [beep,samplingRate] = MakeBeep(freq,duration,samplingRate)
% [beep,samplingRate] = MakeBeep(freq,duration,[samplingRate])
%
% Compute array that can be used by Snd to produce a pure tone of specified
% "freq" (Hz) and "duration" (s). The "samplingRate" defaults to
% Snd('DefaultRate').
% 
% 	beep = MakeBeep(freq,duration);
% 	Snd('Open');
% 	.... do some stuff ....
% 	Snd('Play',beep);
%
% See Snd.

% 6/21/95		dhb, ccc	PlayBeep: Wrote it.
% 3/29/97		dgp 			Updated
% 4/2/97		dgp				Expanded comments above.
% 11/25/97  dhb     	Fixed comment to correctly indicate milliseconds.
% 12/10/97  dhb				Add samplingRate and DONTPLAY args, snd return.
% 2/9/98		dgp				Updated to use Snd instead of SndPlay.
% 2/13/98   dhb       Return sampling rate.
% 2/16/98   dgp       MakeBeep: Based on PlayBeep, but "duration" is now
%											in s, not ms, and default sampling rate is now same as Snd.
% 11/1/99   dgp       Cosmetic.
% 4/13/02   dgp       Make the default samplingRate platform dependent, to match Snd.
% 4/13/02   dgp       Get the default samplingRate from Snd.
% 2/07/19   dcn       Off by one fixed.

if nargin<2 || isempty(duration)
	error('Usage: beep=MakeBeep(freq,duration,[samplingRate]);')
end
if nargin<3 || isempty(samplingRate)
	samplingRate = Snd('DefaultRate');
end
beep = sin(2*pi*freq*(0:duration*samplingRate-1)/samplingRate);
