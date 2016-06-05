function oldClut = LoadIdentityClut(windowPtr, loadOnNextFlip, lutType, disableDithering)
% oldClut = LoadIdentityClut(windowPtr [, loadOnNextFlip=0][, lutType=auto][, disableDithering=1])
%
% Loads an identity clut on the window specified by windowPtr. If
% loadOnNextFlip is set to 1, then the clut will be loaded on the next call
% to Screen('Flip'). By default, the clut will be loaded immediately or on
% the next vertical retrace, depending on OS and graphics card.
% By default, this also tries to disable digital output dithering on supported
% hardware, setting the optional flag 'disableDithering' to zero will
% leave dithering alone.
%
% If you use Linux and have low level hardware access enabled via a call
% to PsychLinuxConfiguration during installation or later, or if you
% use OS/X and have the PsychtoolboxKernelDriver loaded
% ("help PsychtoolboxKernelDriver") and the graphics card is a GPU of the
% AMD Radeon X1000 series or later card, or a equivalent model of the FireGL
% or FirePro series, then this routine will try to use special low-level setup
% code for optimal identity mapping.
% It will also disable digital display dithering if disableDithering == 1.
% On Windows with AMD cards, digital display dithering will also get disabled
% automatically. On other graphics cards, digital display dithering will not
% be affected and needs to be manually controlled by you in some vendor specific way.
%
% On other than AMD cards under Linux or OSX, the function will make
% a best effort to upload a suitable clut, as follows:
%
% This mechanism relies on heuristics to detect the exact type of LUT to
% upload into your graphics card. These heuristics may go wrong, thanks to
% the ever increasing amount of graphics hardware driver bugs in both
% Windows and MacOS/X. For that reason you can also use the function
% SaveIdentityClut() to manually specify either the type of LUT to use
% (overriding the automatic choice), or to specify the complete LUT, or to
% capture the current identity LUT of a display that works. That function
% will store the override Clut in a per-user, per-GPU, per-Screen configuration
% file. If LoadIdentityClut finds such a matching configuration file, it
% will use the LUT specified there, instead of performing an automatic
% selection.
%
% The routine optionally returns the old clut in 'oldClut'.
%
% You can restore the old "original" LUT at any time by calling
% RestoreCluts, or sca, but only until you call clear all! The original
% LUT's are backed up in a global variable for this to work.
%
% If you use a Cambridge Research systems Bits+ or Bits# box or a VPixx Inc.
% DataPixx, ViewPixx or ProPixx device, use the script BitsPlusIdentityClutTest
% for advanced diagnostic and troubleshooting wrt. identity clut's, display
% dithering and other evils which could spoil your day for high bit depth
% visual stimulus display.

% History:
% ??/??/??   mk  Written.
% 05/31/08   mk  Add code to save backup copy of LUT's for later restore.
% 09/09/09   mk  Add more LUT workarounds for Apple's latest f$%#%$#.
% 09/16/09   mk  Even more workarounds for ATI's latest %#$%#$ on Windows.
% 09/20/09   mk  Add option to load LUT type or LUT from config file to
%                override the auto-detection here.
% 05/30/11   mk  Add option to use Screen's built in low-level GPU setup
%                methods to achieve an identity mapping. Fallback to old
%                heuristics if that is unsupported, disabled or failed.
% 03/05/14   mk  Update help texts and some diagnostic output.
% 03/06/14   mk  Allow control if dithering should be touched or not.
% 09/20/14   mk  Silence "unknown gpu" warning for Intel on Linux.
% 04/04/15   mk  Only use PsychGPUControl() dither disable on Windows by
%                default. Only use on Linux + AMD as fallback if low-level
%                dither disable fails.
% 05/18/16   mk  Add detection for Broadcom VC4 SoC gpu in RaspberryPi.

global ptb_original_gfx_cluts;


if nargin < 1
    error('Invalid number of arguments to LoadIdentityClut.');
end

% If not specified, we'll set the clut to load right away.
if nargin < 2
    loadOnNextFlip = [];
end

if isempty(loadOnNextFlip)
    loadOnNextFlip = 0;
end

if nargin < 3
    lutType = [];
end

if nargin < 4 || isempty(disableDithering)
    disableDithering = 1;
end

% Get screen id for this window:
screenid = Screen('WindowScreenNumber', windowPtr);

% Query what kind of graphics hardware is
% installed to drive the display corresponding to 'windowPtr':

% Query vendor of associated graphics hardware:
winfo = Screen('GetWindowInfo', windowPtr);

% Get current clut for use as backup copy later on:
oldClut = Screen('ReadNormalizedGammaTable', windowPtr);

if disableDithering && IsWin
    fprintf('LoadIdentityClut: Info: Trying to disable digital display dithering.\n');
    % Try to use PsychGPUControl() method to disable display dithering
    % globally on all connected GPUs and displays. We only use this
    % on MS-Windows and it only works with AMD/ATI GPU's,
    % It will no-op silently on other system configurations.
    %
    % We don't use the success status return code of the function, because
    % we don't know how trustworthy it is. Also this only affects dithering,
    % not gamma table identity setup, so the code-pathes below must run anyway
    % for proper setup, even if their dithering disable effect may be redundant.
    PsychGPUControl('SetDitheringEnabled', 0);
end

% Ask Screen to use low-level setup code to configure the GPU for
% untampered framebuffer pixel paththrough. This is only supported on a
% subset of GPU's under certain conditions if the PsychtoolboxKernelDriver
% is loaded, but should be the most reliable solution if it works:
if ~IsWin
    % Low level control possible for some GPU's (e.g., AMD).
    % [] will enable full passthrough and force dithering off.
    % -1 will enable passthrough except for dithering, which is left at the OS default.
    if disableDithering
        % Also disable dithering:
        passthroughrc = Screen('LoadNormalizedGammatable', windowPtr, []);
    else
        % Only identity LUTs, no color conversion, degamma etc., but leave the
        % dither settings untouched, so the OS + graphics driver stays in control:
        passthroughrc = Screen('LoadNormalizedGammatable', windowPtr, -1);        
    end
else
    passthroughrc = intmax;
end

if ismember(passthroughrc, [1, 2])
    % Success. How well did we do?
    if passthroughrc == 2
        fprintf('LoadIdentityClut: Info: Used GPU low-level setup code to configure (hopefully) perfect identity pixel passthrough.\n');
    else
        fprintf('LoadIdentityClut: Warning: Used GPU low-level setup code to configure a perfect identity gamma table, but failed at\n');
        fprintf('LoadIdentityClut: Warning: configuring rest of color conversion hardware. This may or may not work.\n');
    end
else
    % No success. Either the low-level setup is unsupported, or it failed:
    if passthroughrc == 0
        fprintf('LoadIdentityClut: Warning: GPU low-level setup code for pixel passthrough failed for some reason! Using fallback...\n');
    elseif IsOSX || IsLinux
        fprintf('LoadIdentityClut: Info: Could not use GPU low-level setup for setup of pixel passthrough. Will use fallback method.\n');
        % AMD GPU, aka GPU core of format 'R100', 'R500', ... starting with a 'R'?
        if winfo.GPUCoreId(1) == 'R'
            % Some advice for AMD users on Linux and OSX:
            fprintf('LoadIdentityClut: Info: On your AMD/ATI GPU, you may get this working by loading the PsychtoolboxKernelDriver\n');
            fprintf('LoadIdentityClut: Info: on OS/X or using a Linux system properly setup with PsychLinuxConfiguration.\n');

            % On AMD try to use PsychGPUControl to force dithering off. This
            % will only work on Catalyst, and has the side effect of disabling
            % dithering on all displays:
            if disableDithering
                PsychGPUControl('SetDitheringEnabled', 0);
            end
        end
    end

    % Carry on with our bag of clut heuristics and other tricks...

    % Raw renderer string, with leading or trailing whitespace trimmed:
    gpuname = strtrim(winfo.GLRenderer);

    % Replace blanks and '/' with underscores:
    gpuname(isspace(gpuname)) = '_';
    gpuname = regexprep( gpuname, '/', '_' );

    % Same game for version string:
    glversion = strtrim(winfo.GLVersion);

    % Replace blanks with underscores:
    glversion(isspace(glversion)) = '_';
    glversion(glversion == '.') = '_';

    % Is there a predefined configuration file for the proper type of LUT to
    % load? Build kind'a unique path and filename for our system config:
    lpath = sprintf('%sIdentityLUT_%s_%s_%s_Screen%i.mat', PsychtoolboxConfigDir, OSName, gpuname, glversion, screenid);

    if exist(lpath, 'file')
        % Configfile! Load it: This will define a variable lutconfig in the
        % workspace:
        load(lpath);

        if isscalar(lutconfig)
            % Single scalar id code given to select among our hard-coded LUT's
            % below:
            gfxhwtype = lutconfig;
            fprintf('Applying the gamma lookup table of type %i, as specified in the following configuration file:\n%s\n', gfxhwtype, lpath);
        else
            % Full blown LUT given:
            lut = lutconfig;
            if ~isnumeric(lut) || size(lut,1) < 1 || size(lut,2) ~= 3
                sca;
                error('LoadIdentityClut: Loaded data from config file is not a valid LUT! Not a numeric matrix or less than 1 row, or not 3 columns!');
            end
            fprintf('Applying the gamma lookup table for identity mapping from the following configuration file:\n%s\n', lpath);
            gfxhwtype = -1;
        end

    else
        % No specific config file. Try our best at auto-detection of flawed
        % systems:

        % Query OS type and version:
        if IsOSX
            % Which OS/X version?
            compinfo = Screen('Computer');
            osxversion = sscanf(compinfo.system, '%*s %*s %i.%i.%i');
        end

        % We derive type of hardware and thereby our strategy from the vendor name:
        gfxhwtype = winfo.GLVendor;

        if ~isempty(strfind(gfxhwtype, 'NVIDIA')) || ~isempty(strfind(gfxhwtype, 'nouveau'))
            % NVidia card:

            % We start with assumption that it is a "normal" one:
            gfxhwtype = 0;

            % GeForce 250 (e.g., GeForce GTS 250) under MS-Windows?
            if IsWin && ~isempty(strfind(winfo.GPUCoreId, 'G80')) && ~isempty(strfind(winfo.GLRenderer, '250'))
                % GeForce 250 on MS-Windows. Needs LUT type 1, according to Jon Peirce:
                gfxhwtype = 1;
            end
            
            % Is it a Geforce-8000 or later (G80 core or later) and is this OS/X?
            if ~isempty(strfind(winfo.GPUCoreId, 'G80')) && IsOSX
                % 10.5.x Leopard?
                if (osxversion(1) == 10) && (osxversion(2) == 5) && (osxversion(3) >=0)
                    % Yes. One of the releases with an embarassing amount of bugs,
                    % brought to you by Apple. Need to apply an especially ugly
                    % clut to force these cards into an identity mapping:
                    fprintf('LoadIdentityClut: NVidia Geforce 8000 or later on OS/X 10.5.x detected. Enabling special type-I LUT hacks for totally broken operating systems.\n');
                    gfxhwtype = 2;
                end

                % 10.6.x Snow Leopard or later?
                if (osxversion(1) == 10) && (osxversion(2) >= 6) && (osxversion(3) >=0)
                    % Yes. One of the releases with an embarassing amount of bugs,
                    % brought to you by Apple. Need to apply an especially ugly
                    % clut to force these cards into an identity mapping:
                    
                    % Peter April reports his GeForce GT 330M on 10.6.8
                    % needs a type 0 LUT, so we handle this case:
                    if ~isempty(strfind(winfo.GLRenderer, '330'))
                        fprintf('LoadIdentityClut: NVidia Geforce GT 330 or later on OS/X 10.6.x or later detected. Enabling type-0 LUT.\n');
                        gfxhwtype = 0;
                    else
                        % Something else: LUT-III is so far correct.
                        fprintf('LoadIdentityClut: NVidia Geforce 8000 or later on OS/X 10.6.x or later detected. Enabling special type-II LUT hacks for totally broken operating systems.\n');
                        gfxhwtype = 3;
                    end
                end
            end
            
            % Is this a latest generation NVidia QuadroFX card under Linux with NVidia binary blob?
            [dummy, dummy2, nrslots] = Screen('ReadNormalizedGammatable', windowPtr); %#ok<ASGLU>
            if IsLinux && ~isempty(strfind(winfo.GLRenderer, 'Quadro')) && (nrslots > 256)
                % LUT type 3:
                gfxhwtype = 3;
                fprintf('LoadIdentityClut: NVidia Quadro or later on Linux with binary Blob detected. Enabling special type-III LUT.\n');
            end
        else
            if ~isempty(strfind(gfxhwtype, 'ATI')) || ~isempty(strfind(gfxhwtype, 'AMD')) || ~isempty(strfind(gfxhwtype, 'Advanced Micro Devices')) || ...
                    ~isempty(strfind(winfo.GLRenderer, 'DRI R')) || ~isempty(strfind(winfo.GLRenderer, 'on ATI R')) || ~isempty(strfind(winfo.GLRenderer, 'on AMD'))
                % AMD/ATI card:

                % A good default at least on OS/X is type 1:
                gfxhwtype = 1;

                if IsWin && ~isempty(strfind(winfo.GPUCoreId, 'R600'))
                    % At least the Radeon HD 3470 under Windows Vista and Linux needs type 0
                    % LUT's. Let's assume for the moment this is true for all R600
                    % cores, ie., all Radeon HD series cards.
                    fprintf('LoadIdentityClut: ATI Radeon HD-2000 or later detected. Using type-0 LUT.\n');
                    gfxhwtype = 0;
                elseif (IsLinux) && (~isempty(strfind(winfo.GLRenderer, 'DRI R')) || ~isempty(strfind(winfo.GLRenderer, 'on ATI R')) || ~isempty(strfind(winfo.GLRenderer, 'on AMD')))
                    if ~isempty(strfind(winfo.GPUCoreId, 'R600'))
                        % At least the Radeon R9 380 Tonga Pro with DCE10 display engine under Linux with DRI2 Mesa needs type 3
                        % LUT's. Let's assume for the moment this is true for all such R600+ cores, ie., all Radeon HD series cards.
                        fprintf('LoadIdentityClut: ATI Radeon HD-2000 or later on Linux DRI2/DRI3 detected. Using type-3 LUT.\n');
                        gfxhwtype = 3;
                    else
                        % At least the Radeon R3xx/4xx/5xx under Linux with DRI2 Mesa needs type 0
                        % LUT's. Let's assume for the moment this is true for all R600
                        % cores, ie., all Radeon HD series cards.
                        fprintf('LoadIdentityClut: ATI Radeon on Linux DRI2 detected. Using type-0 LUT.\n');
                        gfxhwtype = 0;
                    end
                elseif IsOSX && (~isempty(strfind(winfo.GLRenderer, 'Radeon HD 5')))
                    % At least on OS/X 10.6 with ATI Radeon HD-5000 series,
                    % type 2 seems to be the correct choice, according to Cesar
                    % Ramirez (NASA):
                    fprintf('LoadIdentityClut: ATI Radeon HD-5000 Evergreen on OS/X detected. Using type-2 LUT.\n');
                    gfxhwtype = 2;
                end
            elseif ~isempty(strfind(gfxhwtype, 'Intel'))
                % Intel card: Type 0 LUT is correct at least on Linux:
                gfxhwtype = 0;
                fprintf('LoadIdentityClut: Intel integrated graphics chip detected. Using type-0 LUT.\n');
            elseif ~isempty(strfind(gfxhwtype, 'Broadcom')) && ~isempty(strfind(winfo.GLRenderer, 'VC4'))
                % VC4 in RaspberryPi: Type 0 LUT is correct at least on Linux:
                gfxhwtype = 0;
                fprintf('LoadIdentityClut: VideoCore-4 SoC graphics chip detected. Using type-0 LUT.\n');
            else
                % Unknown card: Default to NVidia behaviour:
                gfxhwtype = 0;
                warning('LoadIdentityClut: Warning! Unknown graphics hardware detected. Set up identity CLUT may be wrong!'); %#ok<WNTAG>
            end
        end
    end

    % Use LUT type override from command line, if any is given:
    if ~isempty(lutType)
        gfxhwtype = lutType;
        fprintf('LoadIdentityClut: OVERRIDE: Will use identity LUT of type %i, as specified in function call argument.\n', gfxhwtype);
    end

    % We have different CLUT setup code for the different gfxhw-vendors:
    if gfxhwtype == -1;
        % Upload the lut defined in 'lut':
        oldClut = Screen('LoadNormalizedGammaTable', windowPtr, lut, loadOnNextFlip);
    end

    if gfxhwtype == 0
        % This works on WindowsXP with NVidia GeForce 7800 and OS/X 10.4.9 PPC
        % with GeForceFX-5200. Both are NVidia cards, so we assume this is a
        % correct setup for NVidia hardware:
        oldClut = Screen('LoadNormalizedGammaTable', windowPtr, (0:1/255:1)' * ones(1, 3), loadOnNextFlip);
    end

    if gfxhwtype == 1
        % This works on OS/X 10.4.9 on Intel MacBookPro with ATI Mobility
        % Radeon X1600, also on 10.5.8 with ATI Radeon HD-2600: We assume this
        % is the correct setup for ATI hardware:
        oldClut = Screen('LoadNormalizedGammaTable', windowPtr, ((1/256:1/256:1)' * ones(1, 3)), loadOnNextFlip);
    end

    if gfxhwtype == 2 && IsOSX
        % This works on OS/X 10.5 with NVidia GeForce 8800. It is an ugly hack
        % to compensate for the absolutely horrible, embarassing bugs in Apple's NVidia
        % graphics drivers and their clut handling. Does this company still
        % have a functioning Quality control department, or are they all
        % fully occupied fiddling with their iPhones?!?

        % This works on 10.5.8 with Geforce-9200M in the MacMini according to
        % CRS, and with 10.5.8 with Geforce-8800  in the MacPro according to
        % me. We assume this works on all G80 et al. GPU's:
        loadlut = (linspace(0,(1023/1024),1024)' * ones(1, 3));
        oldClut = Screen('LoadNormalizedGammaTable', windowPtr, loadlut, loadOnNextFlip);
    end

    if gfxhwtype == 3 && IsOSX
        % This works on OS/X 10.6.0 with NVidia Geforce-9200M according to CRS
        % and with 10.6.0 with Geforce-8800  in the MacPro according to
        % me. We assume this works on all G80 et al. GPU's:
        for i=0:1023
            if ((i / 1024.0) < 0.5)
                loadlut(i+1) = (i / 1024.0);
            else
                loadlut(i+1) = (i / 1024.0) - (1.0/256.0);
            end
        end
        loadlut = loadlut' * ones(1, 3);
        oldClut = Screen('LoadNormalizedGammaTable', windowPtr, loadlut, loadOnNextFlip);
    end

    if gfxhwtype == 3 && IsWin
        % This is an experimental variant of the OS/X type 3 lut, but with 256
        % slots. It is supposed for WindowsXP, assuming some NVidia GPU's,
        % e.g., some QuadroFX 3700 GPU's have similar problems:
        for i=0:255
            if ((i / 255.0) < 0.5)
                loadlut(i+1) = (i / 255.0);
            else
                loadlut(i+1) = (i / 255.0) - (1.0/256.0);
            end
        end
        loadlut = loadlut' * ones(1, 3);
        oldClut = Screen('LoadNormalizedGammaTable', windowPtr, loadlut, loadOnNextFlip);
    end

    if gfxhwtype == 3 && IsLinux
        % NVidia QuadroFX cards with binary blob and lut's with more than
        % 256 slots. Upload a standard linear lut with matching number of
        % slots:
        [dummy, dummy2, nrslots] = Screen('ReadNormalizedGammatable', windowPtr); %#ok<ASGLU>
        loadlut = (linspace(0,((nrslots-1) / nrslots), nrslots)' * ones(1, 3));
        oldClut = Screen('LoadNormalizedGammaTable', windowPtr, loadlut, loadOnNextFlip);
    end

    if ~ismember(gfxhwtype, [-1,0,1,2,3])
        sca;
        error('Could not upload identity CLUT to GPU! Invalid LUT or invalid LUT id or other invalid arguments passed?!?');
    end

    % End of high level lut setup path.
end

% Store backup copies of clut's for later restoration by RestoreCluts():

% Create global clut backup cell array if it does not exist yet:
if isempty(ptb_original_gfx_cluts)
    % Create 10 slots for out up to 10 screens:
    ptb_original_gfx_cluts = cell(10,1);
end

% Do we have already a backed up original clut for 'screenid'?
% If so, we don't store this clut, as an earlier invocation of a clut
% manipulation command will have stored the really real original lut:
if isempty(ptb_original_gfx_cluts{screenid + 1})
    % Nope. Store backup:
    ptb_original_gfx_cluts{screenid + 1} = oldClut;
end

return;
