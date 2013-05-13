function [KeyPressed,EventTime] = FORPCheck()
% FORPCheck Checks if a button of a FORP device (HH-5-CYL) is pressed.
% 
% Usage:    
%   
%    [KeyPressed,EventTime] = FORPCheck()  
% 
%
% Return the key name (KeyPressed) of the pressed button and the 
% time (EventTime) of the status check.
%  
%
%    KeyPressed          Key name of the Pressed Button, empty string
%                        if none pressed.
%
%
%    EventTime           Time of keypress check, as returned by GetSecs.
% 
% IMPORTANT NOTE:
%       
%    
%    See FORPWait.
%
%    Going through each device can be very time consuming. If would advise
%    to unplug each unnecessary device, so less devices has to be checked.
%    If you have got any advice for a better way to solve those problems, feel 
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
%    10/10/06   fs   Wrote it.
%    10/16/06   mk   Add caching of device index.
%    10/17/06   fs   Done some restructuring and testing. 
%    02/08/07   mk   New vendor id 6171 added to valid device lists.

    persistent psychtoolbox_forp_id;

    % List of vendor IDs for valid FORP devices:
    vendorIDs = [1240 6171];
    
    KeyPressed  =   '';
    keydata     =   [];

    % Try to detect first FORP device at first invocation:
    if isempty(psychtoolbox_forp_id)
        Devices = PsychHID('Devices');
	% Loop through all KEYBOARD devices with the vendorID of FORP's vendor:
        for i=1:size(Devices,2)
            if strcmp(Devices(i).usageName,'Keyboard') && ismember(Devices(i).vendorID, vendorIDs)
                psychtoolbox_forp_id=i;
                break;
            end
        end
    end

    if isempty(psychtoolbox_forp_id)
        error('FORPCheck: No FORP-Device detected on your system');
    end

    % Needed for getting the very LAST Buttonpress:
    FORPQueueClear(psychtoolbox_forp_id);

    [keydata,err] = PsychHID('GetReport',psychtoolbox_forp_id,1,0,8); % Get Report from device and save pressed.
    EventTime=GetSecs;
    PsychHID('ReceiveReportsStop',psychtoolbox_forp_id);	      % Stop receiving reports from that device

    if any(keydata)
        index = find(keydata ~=0,1);                                  % get the (1st) index of the not-null element in keydata. Only one index can be processed by KbName
        if ~isempty(index)	                                      % if there is a index in keydata then
            KeyPressed=KbName(double(keydata(index)));	              % translate element at index into keyname
        end
    end

    return;
end
