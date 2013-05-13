function WaitSetMouse(newX, newY, windowPtrOrScreenNumber)

% set and wait for new cursor position to take effect

while 1 % wait for new cursor position to be set
    SetMouse(newX,newY,windowPtrOrScreenNumber);
    [mx, my, buttons]=GetMouse(windowPtrOrScreenNumber);
    if mx==newX && my==newY
        break;
    end
end
