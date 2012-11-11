% OptiCAL - Psychtoolbox IOPort interface to the CRS OptiCAL luminance
%           meter device
%
% Usage:
%   >> handle = OptiCAL('Open', port)
%   >> lum = OptiCAL('Read', handle)
%   >> OptiCAL('Close', handle)
%   >> OptiCAL('CloseAll')
%
% Inputs:
%   command   - string 'Open', 'Read', 'Close', or 'CloseAll'
%   port      - string device (e.g. '/dev/ttyS0')
%   handle    - double handle
%
% Outputs:
%   handle    - double handle
%   lum       - luminance in cd/m^2
%
% Note:
%   Requires read and write access to the serial port (Ubuntu: add user to
%   the dialout group). Not yet tested with USB-to-serial adapters.
%
% References:
%   [1] Haenel, V. (2010). Pure python interface to CRS OptiCAL. Retrieved
%       October 18, 2012 from https://github.com/esc/pyoptical.
%   [2] Cambridge Research Systems. (1995). OptiCAL user's guide v4.02.
%       Retrieved October 18, 2012 from
%   	http://support.crsltd.com/ics/support/DLRedirect.asp?fileID=63194
%   [3] Cambridge Research Systems. (2009). OptiCAL user's guide v4.02
%       errata. Retrieved October 18, 2012 from
%   	http://support.crsltd.com/ics/support/DLRedirect.asp?fileID=63194
%
% Author: Andreas Widmann, University of Leipzig, 2012

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2012 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% MIT license:
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.

function [ ret ] = OptiCAL(command, handle)

if nargin < 1
    error('Not enough input arguments.')
end

persistent OptiCALCfg;

switch command

    case 'Open'

        if nargin < 2
            error('Not enough input arguments.')
        end

        % Open device
        port = handle;
        handle = OptiCALOpen(port);

        % Define configuration structure
        OptiCALCfg = [OptiCALCfg OptiCALConfiguration(handle)];

        ret = handle;

    case {'Close', 'CloseAll'}
        
        if nargin < 2 || isempty(handle) || strcmp(command, 'CloseAll')
 
            % Close all devices if used without handle argument or requested CloseAll        
            while ~isempty(OptiCALCfg)
               OptiCALClose(OptiCALCfg(1).handle);
               OptiCALCfg(1) = [];
            end

        else
            
            % Close handle
            idx = [];
            if ~isempty(OptiCALCfg)
                idx = find([OptiCALCfg.handle] == handle);
            end
            if isempty(idx)
                error('Invalid OptiCAL device handle.')
            end
            
            OptiCALClose(OptiCALCfg(idx).handle);
            OptiCALCfg(idx) = [];
 
        end

    case 'Read'

        if nargin < 2 || isempty(handle)
            error('Not enough input arguments.')
        end

        idx = [];
        if ~isempty(OptiCALCfg)
            idx = find([OptiCALCfg.handle] == handle);
        end
        if isempty(idx)
            error('Invalid OptiCAL device handle.')
        end

        % Read raw value
        adc = OptiCALRead(OptiCALCfg(idx).handle);

        % Convert to luminance
        adc_adj = adc - OptiCALCfg(idx).Z_count - 524288; % From OptiCAL manual
        lum = ((adc_adj / 524288) * OptiCALCfg(idx).V_ref * 1e-6) / (OptiCALCfg(idx).R_feed * OptiCALCfg(idx).K_cal * 1e-15); % From pyoptical! See errata!!!
        ret = max([0 lum]);

    otherwise
        error('Unknown command.')

end

end

function handle = OptiCALOpen(port)

    % Open serial port
    OptiCALCfg.port = port;
    handle = IOPort('OpenSerialPort', OptiCALCfg.port);
    WaitSecs('YieldSecs', 0.1); % Why needed? Does not reliably work w/o.

    % Calibrate
    IOPort('Write', handle, 'C');
    data = IOPort('Read', handle, 1, 1);
    if isempty(data) || data ~= 6
        IOPort('Close', handle);
        error('Could not calibrate OptiCAL device.')
    end

    % Set in current Mode
    IOPort('Write', handle, 'I');
    data = IOPort('Read', handle, 1, 1);
    if isempty(data) || data ~= 6
        IOPort('Close', handle);
        error('Could not set OptiCAL device into current mode.')
    end

end

function OptiCALClose(handle)

     IOPort('Close', handle);

end

function OptiCALCfg = OptiCALConfiguration(handle)

    % Defaults
    OptiCALCfg.handle = handle;

    % Read configuration
    OptiCALCfg.prod_type  = OptiCALReadEEPROM(handle,  0,  2, 'int');
    OptiCALCfg.optical_sn = OptiCALReadEEPROM(handle,  2,  4, 'int');
    OptiCALCfg.fw_vers    = OptiCALReadEEPROM(handle,  6,  2, 'int');
    OptiCALCfg.V_ref      = OptiCALReadEEPROM(handle, 16,  4, 'int');
    OptiCALCfg.Z_count    = OptiCALReadEEPROM(handle, 32,  4, 'int');
    OptiCALCfg.R_feed     = OptiCALReadEEPROM(handle, 48,  4, 'int');
    OptiCALCfg.R_gain     = OptiCALReadEEPROM(handle, 64,  4, 'int');
    OptiCALCfg.probe_sn   = OptiCALReadEEPROM(handle, 80, 16, 'char');
    OptiCALCfg.K_cal      = OptiCALReadEEPROM(handle, 96,  4, 'int');

end

function val = OptiCALReadEEPROM(handle, address, nBytes, toType)

    offset = 128;

    % Read field from EEPROM
    data = zeros(1, nBytes);
    for iByte = 1:nBytes;
        IOPort('Write', handle, uint8(offset + address + iByte - 1));
        temp = IOPort('Read', handle, 1, 2);
        data(iByte) = temp(1);
    end

    % Cast to type
    if strcmp(toType, 'char')
        val = char(data);
    else
        val = sum(pow2(0:8:8 * (nBytes - 1)) .* data);
    end

end

function adc = OptiCALRead(handle)

    IOPort('Write', handle, 'L'); % See errata!!!
    temp = IOPort('Read', handle, 1, 4);
    adc = sum(pow2([0 8 16]) .* temp(1:3));

end
