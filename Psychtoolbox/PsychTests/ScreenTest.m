% ScreenTest
%
% Test hardware/software configuration for performance.
%
% 10/14/06  dhb  Wrote stub version.

% Get all screens
screens = Screen('Screens');

% Test each screen
for screen = screens
    fprintf('\n***** ScreenTest: Testing Screen %d *****\n',screen);
    
    % Open window.  Prints some diagnostics all on its own.
    w = Screen('OpenWindow',screen);
    
    % Close window
    Screen('Close',w);
    
    fprintf('\n***** ScreenTest: Done With Screen %d *****\n',screen);
end
