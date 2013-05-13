function strs=getmodestrs(el)
strs={};
mode = Eyelink('CurrentMode');
switch mode
    case -1
        strs{end+1}='disconnected';
    case 0
        strs{end+1}='unknown';
    otherwise
        if bitand(mode,el.IN_DISCONNECT_MODE)
            strs{end+1}='IN_DISCONNECT_MODE';
        end

        if bitand(mode,el.IN_UNKNOWN_MODE)
            strs{end+1}='IN_UNKNOWN_MODE';
        end

        if bitand(mode,el.IN_IDLE_MODE)
            strs{end+1}='IN_OFFLINE_MODE';
        end

        if bitand(mode,el.IN_SETUP_MODE)
            strs{end+1}='IN_SETUP_OR_CAL_VAL_DCORR_MODE';
        end

        if bitand(mode,el.IN_RECORD_MODE)
            strs{end+1}='IN_RECORD_MODE';
        end

        if bitand(mode,el.IN_TARGET_MODE)
            strs{end+1}='IN_FIXATION_TARGETS_NEEDED_MODE';
        end

        if bitand(mode,el.IN_DRIFTCORR_MODE)
            strs{end+1}='IN_DRIFTCORR_MODE';
        end

        if bitand(mode,el.IN_USER_MENU)
            strs{end+1}='IN_USER_MENU_MODE';
        end

        if bitand(mode,el.IN_PLAYBACK_MODE)
            strs{end+1}='IN_PLAYBACK_MODE';
        end

        if bitand(mode,el.IN_IMAGE_MODE)
            strs{end+1}='IN_IMAGE_MODE';
        end
end
end
