function err=caldonebeep(el, error)

if 0
if error<=0
	err=SND('Play', el.calibrationfailedsound);
else
	err=SND('Play', el.calibrationsuccesssound);
end
SND('Wait');

end
