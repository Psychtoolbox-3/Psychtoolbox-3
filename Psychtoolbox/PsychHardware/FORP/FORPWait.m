function [KeyPressed,EventTime] = FORPWait(TimeToWait)
% FORPWait     Checks for the specified amount of time, if a button of a FORP device
%              (tested for HH-5-CYL) is pressed , returns if a button has been
%              pressed or specified amount of time has passed.
% 
% Usage:    
%   
%    [KeyPressed,EventTime] = FORPWait([Seconds])  
% 
% Arguments:
%
%    Seconds             Maximum time to check for buttonpresses in seconds.
%                        Default's to 'wait forever' if not provided.
%
%
%
% Returns the keycode (KeyPressed) of the pressed button and the 
% time (EventTime) of the status check.
%  
%
%    KeyPressed          KeyCode of the Pressed Button or empty value if
%                        waiting timed out without any key press.
%
%
%    EventTime           Time of keypress as returned by GetSecs.
% 
% 
% IMPORTANT NOTE:
%       
%    Current-Designs FORP Device (HH-5-CYL) does not return any values for
%    manufacturer or product, so i used the VendorID returned by 
%    PsychHID('Devices') for the HH-5-CYL.(ATM i do no really know if the 
%    VendorID has unique values).
%    Another issue i had to solve was a ?bug? using 'GetReport'. I had to 
%    close the Callbackhandlers to the current device by calling 
%    'ReceiveReportsStop' before calling 'GetReport' on another device.
%    If you have got any advice for a better way solve those problems, feel 
%    free to let me know:
%
%           Florian Stendel 
%           Visual Processing Lab
%           Universitaets - Augenklinik Magdeburg
%           Leipziger Strasse 44
%           39120 Magdeburg
%           Tel:    0049 (0)391 67 21723
%           Email:  vincentdhs@gmx.de
%
%
%    09/10/06   fs   Wrote it.
%    19/10/06   fs   Added some further improvements suggested by Mario
%                    Kleiner.

    KeyPressed  =   '';
    keydata     =   []; 

    if nargin < 1
        TimeToWait=inf;
    end

    deadline = GetSecs + TimeToWait;
    while (deadline > GetSecs) && isempty(KeyPressed)
        [KeyPressed,EventTime] = FORPCheck;
    end
    return;
end
