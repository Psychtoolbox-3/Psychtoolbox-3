function S = MakeItS(wls)
% S = MakeItS(wls)
%
% If argument is a [start delta n] description, it is
% left alone.
%
% If passed length is not a [start delta n] description,
% convert it to one.  Formats handled are a list of evenly
% spaced wavelengths or a struct with fields start, step, numberSamples.
%
% Format error checking could be more agressive.
%
% 7/26/02  dhb  Allow struct format too.

% Force passed description to S format.
[m,n] = size(wls);
if (isstruct(wls))
	S = [wls.start wls.step wls.numberSamples];
elseif (m == 1 && n == 3)
  if (wls(1) >= 0 && wls(3) > 0)
    S = wls;
  else
    error('Passed wls is not interpretable');
  end
else
  S = WlsToS(wls);
end


