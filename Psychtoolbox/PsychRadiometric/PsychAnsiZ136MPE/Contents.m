% Psychtoolbox:PsychAnsiZ136MPE.
%
% Partial and in progress implementation of Ansi Z136.1-2007 standard for computing
% maximum permissable exposure.
%
% *****************************************************************
% IMPORTANT:
%   a) Individuals using these routines must accept full responsibility 
%   for light exposures they implement. We recommend that values computed
%   with these routines be carefully checked against independent calculations.
%   We have done our best to follow the standard, but it is very complex and
%   there may be errors.
%
%   b) As of March, 2013, these routines are still very much a work
%   in progress and should thus be treated with special caution.  See
%   AnsiZ136MPEBasicTest and AnsiZ136MPEDeloriTest, both for test 
%   code, and for comments about points of uncertainty.
%
%   c) Particularly obscure to me at present is the limiting cone aperture
%   section described with an asterisk in Table 2 of the 2007 standard 
%   document.
%
%   d) The field's knowledge of safe light levels is evolving, and it is important
%   to keep abreast of current research and not to rely soley on the
%   Ansi standard.  In particular, note that there
%   are reports of disrubption at light levels below the limits specified
%   in the Ansi 2007 stanard. See for example:
%     Morgan et al., (2008), IOVS, 49, 3715-3729.
%     Morgan et al., (2009), IOVS, 50, 6015-6022.
%     Hunter et al., (2012), Prog. Ret. Eye. Res., 31, 28-42.
%
%   e) Apparently a new version of the Ansi Z136 standard is forthcoming,
%   but I do not have information as to what changed.
%
%   f) The Z136 standard is for laser light, which is coherent and monochromatic.
%   Different standards apply for broadband lights.  Two relevant references are:
%     ISO 15004-2:2007, Ophthalmic instruments -- Fundamental requirements and
%     test methods -- Part 2: Light hazard protection.
%
%     ICNIRP (International Commission on Non-Ionizing Radiation Protection) Guidelines,
%     (1997), Guidelines on exposure to broad-band inchoherent optical radiation (0.38 to 3 uM),
%     Health Physics, 73, 539-554.
%   The second seems a little old, but is the most recent ICNIRP guideline for broadband
%   that I could easily locate.
%
%  - David Brainard, March 6, 2013.
% *****************************************************************
%
% REFERENCES.
%   Ansi Z136.1-2007.  The standard document. The more recent Ansi Z136.8-2012
%   refers back to the 2007 document for the MPE calculations.  But see
%   point d) in the notes above.
%
%   Delori et al., 2007, JOSA A, 24, 1250-1265.  Provides explanation
%   about many of the calculations in the 2000 version of the standard.
%
%   AnsiZ136MPEBasicTest - Test the suite of routines.  Generates many figures that should match those in the standard.
%   AnsiZ136MPEComputeCa - Compute constant Ca, Table 6
%   AnsiZ136MPEComputeCb - Compute constant Cb, Table 6
%   AnsiZ136MPEComputeCc - Compute constant Cc, Table 6
%   AnsiZ136MPEComputeCe - Compute constant Ce, Table 6
%   AnsiZ136MPEComputeExtendedSourceLimit - Compute overall MPE limit for extended sources, Table 5b
%   AnsiZ136MPEComputeExtendedSourcePhotochemicalLimit - Compute photochemical MPE for extended sources, Table 5b
%   AnsiZ136MPEComputeLimitingConeAngle - Compute limiting cone angle, Table 6
%   AnsiZ136MPEComputeT2 - Compute constant T2, Table 6
%   AnsiZ136MPEDeloriTest - Tests our computations against those from Delori's spreadsheet.

% Copyright (c) 2013 by David Brainard