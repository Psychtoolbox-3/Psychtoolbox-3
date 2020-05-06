function [seed,whichGen] = ClockRandSeed(seed,whichGen)
% [seed,whichGen] = ClockRandSeed([seed],[whichGen])
% ClockRandSeed seeds the random number from the time-of-day clock or another
% "random enough" random source, e.g., high resolution time, or a true
% random number source. If 'seed' is passed this is used instead of a
% random source or clock time.
%
% 'whichGen' optionally allows selection of the random number generator.
% If omitted, the recommended Mersenne-Twister will be used on both
% Matlab and Octave.
%
% If supported, the rng() function will be used, with rng('shuffle') if
% 'seed' is omitted. Otherwise the old rand()/randn() setup method is
% used, with 'reset' on Octave for truly random init from /dev/urandom
% or high resolution clock time. On Matlab the old seed = fix(1e6*sum(clock))
% will be used.
%
% NOTE: There is a page in the 2017 Matlab docs describing the issues around
% random number generator selection and seeding. Type "help rand" and then
% follow to the note on "Replace Discouraged Syntaxes of rand and randn" that
% is at the end of the Description section.
%
% For reseeding on old Matlab versions, note the following:
%
% The multiplier 1e6 in fix(1e6*sum(clock)) is bigger than the 100
% suggested by Mathworks (see help RAND). They suggest 100*sum(clock), but
% this would only change the seed once every 1/100 sec, which might
% correspond to many iterations of a loop. With the bigger multiplier it's
% pretty sure (depending only on clock grain) that two successive calls
% will generate different seeds.
% 
% Also unlike the Mathworks suggestion (see help RAND), we call FIX. This
% has no effect on RAND and RANDN, but makes it easier to correctly print
% the seed to a file. When a random seed generates an interesting result
% one would like to be able to recreate the conditions that generated the
% run. If the seed 1e6*sum(clock) is printed to a file, '%.0f' will often
% round it differently than RAND would.
%
% See also: rng, rand, randn.

% 10/15/93  dhb Made this a function.
% 4/11/94   dhb Added randn seeding.
% 8/9/97    dgp Update for more reliable seed.
%           dgp Gutted it to make version-dependence explicit.
% 8/9/97    dhb Ungutted it to hide version-dependence.
% 8/13/97   dhb Optional return of which generator used and passing of seed/generator.
% 7/24/04   awi Cosmetic.
% 1/29/05   dgp Cosmetic.
% 06/15/17  dhb Added warning and help text that this is obsolete.
% 9-Aug-17  mk  Undid dhb's warning, as it promotes writing non-portable code
%               that only works on specific versions of Matlab, and thereby
%               creates vendor lock-in. Instead rewrote the function to call
%               equivalent methods to rng('shuffle') in current Octave and older
%               versions of Matlab, and rng('shuffle') in current Matlab and
%               future Octave. That was the whole point of having ClockRandSeed()
%               in the first place: Something that doesn't break each time that
%               Mathworks incompetent or careless software engineering breaks
%               backward compatibility. At least i don't want this to degrade to
%               Mathworks or Apples standards if possible.

if nargin < 2 || isempty(whichGen)
  % Recommended generator on Matlab and Octave is 'twister' for Mersenne-Twister.
  % This is also the default startup generator of Octave and modern Matlab:
  whichGen = 'twister';
end

if nargin < 1 || isempty(seed)
  % No seed specified. Let Matlab/Octave do the job:
  if exist('rng')
    % Reinit from system time or other appropriate source:
    rng('shuffle', whichGen);
    % Retrieve new seed after reinit:
    seed = rng;
  else
    if IsOctave
      % Octave 'reset' will reinit from system time, or from the
      % truly random random number generator /dev/urandom if it
      % exists, e.g. on Linux:
      rand(whichGen, 'reset');
      % Retrieve new seed after init:
      seed = rand(whichGen);
      % Apply it to the other randn() generator, so both operate
      % off the same 'seed':
      randn(whichGen, seed);
    else
      % Matlab: Manually reinit with clock time:
      seed = fix(1e6*sum(clock));
      rand(whichGen, seed);
      randn(whichGen, seed);
    end
  end

  % Got our reinit and whichGen and seed contain the new selected
  % generator type and seed, so we are done:
  return;
end

% If we get to this point, then user wants total control:
if exist('rng')
  if isstruct(seed) % later Matlab return struct by rng()
    rng(seed);
  else
    rng(seed, whichGen);
  end
else
  rand(whichGen,seed);
  randn(whichGen,seed);
end

return;
