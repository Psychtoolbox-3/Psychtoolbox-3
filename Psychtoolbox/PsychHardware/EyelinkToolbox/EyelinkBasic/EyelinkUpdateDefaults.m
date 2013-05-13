function EyelinkUpdateDefaults(el)
% USAGE = EyelinkUpdateDefaults(el) 
%
% This function passes changes to the calibration defaults structure
% to the callback function. To be called after EyelinkInitDefaults if any
% changes are made to the output structure of that function.
%
% el=EyelinkInitDefaults(window);
% el.backgroundColour = BlackIndex(window);
% EyelinkUpdateDefaults(el);
%
% 27-1-2011 NJ created 
% 19-12-2012 IA Fix hardcoded callback 

if ~isempty(el.callback) && exist(el.callback,'file')
	
	%evaluate the callback function with the new el defaults
	feval(el.callback, el); 
		
end

end
