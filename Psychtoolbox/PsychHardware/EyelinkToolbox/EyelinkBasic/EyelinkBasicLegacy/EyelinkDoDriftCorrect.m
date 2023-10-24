function result=EyelinkDoDriftCorrect(el, x, y, draw, allowsetup)
    warning('EyelinkToolbox:LegacyDoDriftCorrect',['The function EyelinkDoDriftCorrect() is deprecated. Please update ', ...
	'your script to use the current method for handling camera setup mode callbacks with PsychEyelinkDispatchCallback.m.']);
	warning('off', 'EyelinkToolbox:LegacyDoDriftCorrect');