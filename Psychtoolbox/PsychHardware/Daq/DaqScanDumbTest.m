% DaqScanDumbTest - Minimal test of DaqAInScan() functionality.
% Tries 1 channel acquisition of samples from channel 0 with
% default sampling parameters (1000 samples at 1 KhZ typically).

 daq = DaqDeviceIndex
 options.FirstChannel=0
 options.LastChannel=0
 params=DaqAInScanBegin(daq, options )
 tic;
 while toc < 2
 DaqAInScanContinue (daq, options);
 end
 v = DaqAInScanEnd (daq, options);
