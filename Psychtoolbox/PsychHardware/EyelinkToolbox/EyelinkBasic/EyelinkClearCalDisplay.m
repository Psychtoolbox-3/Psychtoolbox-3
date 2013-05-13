function EyelinkClearCalDisplay(el)

Screen( 'FillRect',  el.window, el.backgroundcolour );	% clear_cal_display()
Screen( 'Flip',  el.window);
