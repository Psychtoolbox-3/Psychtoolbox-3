function result=EyelinkLegacyDoTrackerSetup(el, sendkey)
    warning('EyelinkToolbox:LegacyDoTrackerSetup',['Use of the function EyelinkDoTrackerSetup() without providing a callback handler ', ...
        '(such as the included PsychEyelinkDispatchCallback) is deprecated. Please update your script to use the currently supported conventions.']);
	warning('off', 'EyelinkToolbox:LegacyDoTrackerSetup');