function EyelinkClearCalDisplay(el)
warning('EyelinkToolbox:LegacyClearCalDisplay',['The function EyelinkClearCalDisplay() is deprecated. Please update ', ...
    'your script to use the current method for handling camera setup mode callbacks with PsychEyelinkDispatchCallback.m.']);
warning('off', 'EyelinkToolbox:LegacyClearCalDisplay');

Screen( 'FillRect',  el.window, el.backgroundcolour );    % clear_cal_display()
Screen( 'Flip',  el.window);