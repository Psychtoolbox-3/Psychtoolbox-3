function [key, el]=EyelinkGetKey(el)
% USAGE: [key, el]=EyelinkGetKey(el)
%
%  el: eyelink default values, also stores
%      getkeytime (last time this function was used)

% matlab version of eyelink supplied getkey() function
% we can't call it getkey() as there is already a similar,
% function in PsychToolbox

% 02-06-01  fwc changed to accept el structure and tickcount value
%           which should hold the time at which getkeyforeyelink was used
%           last time.
%           key definitions are now set in 'initeyelinkdefaults.m'
% 12-01-04  fwc adapted for use with OS X version of psychToolbox
%                 replaced tickcount with el.getkeytime, that is
%                 stored in el structure. Also max key repetition time is
%                 now a parameter (el.getkeyrepeat) set in initeyelinkdefaults.m
%                 there's now also a variant that tests
%                 for key presses that only uses KbCheck
%                 CharAvail fails on OS X so we only test for keypresses using keyIsDown
%                 parameter from KbCheck function
% 280606    fwc changed name to EyelinkGetKey
%               now actually returning char codes rather than true
%               characters
% 270309    edf added function and modifier keys
% 201213    lj  added parameters for function modifyKey, for octave to get key responses from display side



% excerpt from "exptsppt.h"
% ******** KEY SCANNING ********/
% some useful keys returned by getkey()                    */
% These keys allow remote control of tracker during setup. */
% on non-DOS platforms, you should produce these codes and */
% all printable (0x20..0x7F) keys codes as well.           */
% Return JUNK_KEY (not 0) if untranslatable key pressed.   */
% TERMINATE_KEY can be to break out of EXPTSPPT loops.     */

% Returns 0 if no key pressed          */
% returns a single UINT16 integer code */
% for both standard and extended keys  */
% Standard keys == ascii value.        */
% MSBy is set for extended codes       */

  key=0;

  if nargin < 1  || ~isstruct(el)
    error( 'USAGE: [key, el]=EyelinkGetKey(el)' );
  end

  if GetSecs-el.getkeytime<el.getkeyrepeat  % to prevent too rapid repeats of this function?
    return;
  end

  % in original getkey() there was a test for cmd . (call to UserAbort()).
  % Obsolete within matlab environment as this stops matlab execution completely (on OS 9).
  % Here we provide an alternative (default is apple-esc).
  % you can change this by setting el.modifierkey en el.quitkey
  % specific quitkey defined in 'EyelinkInitDefaults.m' file.

  [keyIsDown,secs,keyCodes] = KbCheck(el.devicenumber);

  if 1==isequal(keyCodes, el.lastKeyCodes)
    return;
  end
  el.lastKeyCodes=keyCodes;

  if keyIsDown==0 % no keypress available
    return;
  end

  if keyCodes(el.modifierkey) && keyCodes( el.quitkey )
    key=el.TERMINATE_KEY;
    return;
  end

  el.getkeytime=secs;
  if ~isfield( el, 'keysCached') || isempty(el.keysCached) || el.keysCached==0
    el;
    error( 'Please call el=EyelinkInitDefaults before calling EyelinkGetKey' );
  end

  % first see if any of the special eyelink keys have been pressed
  % most of these need to be translated into a special value
  % These keys allow remote control of tracker during setup. */
  % on non-DOS platforms, you should produce these codes */

  %     key=find(keyCodes);

  if any(keyCodes(el.uparrow))==1
    key=el.CURS_UP;
  elseif any(keyCodes(el.downarrow))==1
    key=el.CURS_DOWN;
  elseif any(keyCodes(el.leftarrow))==1
    key=el.CURS_LEFT;
  elseif any(keyCodes(el.rightarrow))==1
    key=el.CURS_RIGHT;
  elseif any(keyCodes(el.return))==1
    key=el.ENTER_KEY;
  elseif any(keyCodes(el.enter))==1
    key=el.ENTER_KEY;
  elseif keyCodes(158)==1 % Return
    key=el.ENTER_KEY;
  elseif keyCodes(40)==1 % Return
    key=el.ENTER_KEY;
  elseif keyCodes(187)==1 % KeypadBackspace
    key=el.ENTER_KEY;
  elseif any(keyCodes(el.backspace))==1
    key=el.BACKSPACE;
  elseif any(keyCodes(el.escape))==1
    key=el.ESC_KEY;
  elseif any(keyCodes(el.pageup))==1
    key=el.PAGE_UP;
  elseif any(keyCodes(el.pagedown))==1
    key=el.PAGE_DOWN;
  elseif any(keyCodes(el.space))==1
    key=el.SPACE_BAR; % returns 32
  elseif any(keyCodes(el.f1))
    key=el.F1_KEY;
  elseif any(keyCodes(el.f2))
    key=el.F2_KEY;
  elseif any(keyCodes(el.f3))
    key=el.F3_KEY;
  elseif any(keyCodes(el.f4))
    key=el.F4_KEY;
  elseif any(keyCodes(el.f5))
    key=el.F5_KEY;
  elseif any(keyCodes(el.f6))
    key=el.F6_KEY;
  elseif any(keyCodes(el.f7))
    key=el.F7_KEY;
  elseif any(keyCodes(el.f8))
    key=el.F8_KEY;
  elseif any(keyCodes(el.f9))
    key=el.F9_KEY;
  elseif any(keyCodes(el.f10))
    key=el.F10_KEY;
  end

  function modifyKey(el, keyCodes, key)
    if any(keyCodes(el.left_shift))
      key=bitand(key,el.ELKMOD_LSHIFT); 
    end
    if any(keyCodes(el.right_shift)) 
      key=bitand(key,el.ELKMOD_RSHIFT); 
    end
    if any(keyCodes(el.left_control)) 
      key=bitand(key,el.ELKMOD_LCTRL); 
    end
    if any(keyCodes(el.right_control)) 
      key=bitand(key,el.ELKMOD_RCTRL); 
    end
    if any(keyCodes(el.lalt))
      key=bitand(key,el.ELKMOD_LALT); 
    end
    if any(keyCodes(el.ralt))
      key=bitand(key,el.ELKMOD_RALT); 
    end
    if any(keyCodes(el.num))
      key=bitand(key,el.ELKMOD_NUM); 
    end
    if any(keyCodes(el.caps)) 
      key=bitand(key,el.ELKMOD_CAPS); 
    end
  end

  % already found a key response,
  if key~=0
    modifyKey(el, keyCodes, key);
    return;
  end

  % we should produce all printable (0x20..0x7F) keys codes as well.           */
  % this is probably enough

  key=find(keyCodes); % find keyCodes>0
  %     key=key(1); % take only first one in case more keys were pressed
  key=KbName(key(1)); % convert to a character, potential problem: KbName does not translate all keys to actual ascii characters.
  key=double(key(1)); % some keyNames have two characters (e.g. 8*);
  if key >=hex2dec('20') && key < hex2dec('7F') % only return  printable chars (0x20..0x7F)
    modifyKey(el, keyCodes, key);
    return;
  end

  key=el.JUNK_KEY; % return code for untranslatable key
  return;
end
