function CMClose(meterType)
% CMClose([meterType])
%
% Close color meter interface.
%
% meterType 1 is the PR650 (default)
% meterType 2 is the CVI (need CVIToolbox)
% meterType 3 is the CRS Colorimeter
% meterType 4 is the PR655
% meterType 5 is the PR670
% meterType 6 is the PR705
%
% 2/15/02  dhb  Wrote it.
% 4/13/02  dgp	Cosmetic.
% 2/26/03  dhb  Added more meter types. 
% 3/27/03  dhb, jmh  Fix up default argument.
% 2/07/09  mk, tbc  Add PR-655 support.
% 12/04/12 zlb  Add PR-705 support.

% Set default meterType.
if nargin < 1 || isempty(meterType)
  meterType = 1;
end

switch meterType
	case 1
		% PR-650
        PR650close;
	case 2
		if ~exist('CVIClose') %#ok<EXIST>
			error('Need CVIToolbox to call CVIClose (meterType 2)');
		end
		CVIClose;
	case 3
		CRSColorClose;
	case 4
		% PR-655
        PR655close;
    case 5
		% PR-670
        PR670close;
    case 6
        % PR-705
        PR705close;
	otherwise,
		error('Unknown meter type');
end

