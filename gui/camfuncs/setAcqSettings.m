function out = setAcqSettings(acq)

disp(' ');
disp('Updating acquisition settings ...');

if ~isValidAcq(acq)
    warning('Bad acquisition settings provided. Aborting.');
    return;
end
      

% Set Acquisition Mode 
% 1: Single Scan, 2: Accumulate, 3: Kinetics, 4: Fast Kinetics, 5: Run till abort
fprintf([' SetAcquisitionMode       : ' num2str(acq.AcquisitionMode) ' ... ']);
[ret]=SetAcquisitionMode(acq.AcquisitionMode);
disp(error_code(ret))

% Set Read Mode
% 0: Full Vertical Binning 1: Multi-Track, 2: Random-Track, 3: Single-Track, 4: Image
fprintf([' SetReadMode              : ' num2str(acq.ReadMode) ' ... ']);
[ret]=SetReadMode(acq.ReadMode);
disp(error_code(ret))

% Set Shutter Mode
% Type : (0: LOW to open, 1: HIGH to open), Mode: (0: Auto, 1: Open, 2:Close), Times in ms
fprintf([' SetShutter               : ' num2str(acq.ShutterType) ',' ...
    num2str(acq.ShutterMode) ',' num2str(acq.ClosingTime) ',' ...
    num2str(acq.OpeningTime) ' ... ']);
[ret]=SetShutter(acq.ShutterType, acq.ShutterMode,acq.ClosingTime, acq.OpeningTime);
disp(error_code(ret))

% Set Fan Mode
% 0: on Full, 1: on low, 2: off
fprintf([' SetFanMode               : ' num2str(acq.FanMode) ' ... ']);
[ret]=SetFanMode(acq.FanMode);
disp(error_code(ret))

% Set AD Channel
fprintf([' SetADChannel             : ' num2str(acq.ADChannelIndex) ' ... ']);
[ret]=SetADChannel(acq.ADChannelIndex);
disp(error_code(ret))

% Set PreAmp Gain 
fprintf([' SetPreAmpGain            : ' num2str(acq.PreAmpGainIndex) ' ... ']);
[ret]=SetPreAmpGain(acq.PreAmpGainIndex);
disp(error_code(ret))

% DOES THIS NEED TO BE SET?
% % Set Output Amplfier Gain Mode 
% % 0: DAC 0-255, 1: DAC 0-4095, 2: Linear, 3: Real EM
% fprintf([' SetOutputAmplifier       : ' num2str(0) ' ... ']);
% [ret]=SetOutputAmplifier(0);
% disp(error_code(ret))

% Set EM Gain Mode 
% 0: DAC 0-255, 1: DAC 0-4095, 2: Linear, 3: Real EM
fprintf([' SetEMGainMode            : ' num2str(acq.EMGainMode) ' ... ']);
[ret]=SetEMGainMode(acq.EMGainMode);
disp(error_code(ret))

% Enable/Disable High EM Gains
fprintf([' SetEMAdvanced            : ' num2str(acq.EMAdvanced) ' ... ']);
[ret]=SetEMAdvanced(acq.EMAdvanced);
disp(error_code(ret))

% Set EM Gain
fprintf([' SetEMCCDGain             : ' num2str(acq.EMCCDGain) ' ... ']);
[ret]=SetEMCCDGain(acq.EMCCDGain);
disp(error_code(ret))

% Set Vertical Shift Speed
fprintf([' SetVSSpeed               : ' num2str(acq.VSSpeedIndex) ' ... ']);
[ret]=SetVSSpeed(acq.VSSpeedIndex);
disp(error_code(ret))

% Set Horizontal Shift Speed
fprintf([' SetHSSpeed               : ' ...
    num2str(acq.AmpTypeIndex) ',' num2str(acq.HSSpeedIndex) ' ... ']);
[ret]=SetHSSpeed(acq.AmpTypeIndex,acq.HSSpeedIndex);
disp(error_code(ret))

% Set Exposure Time 
fprintf([' SetExposureTime          : ' num2str(acq.ExposureTime) ' ... ']);
[ret]=SetExposureTime(acq.ExposureTime);
disp(error_code(ret))

% Set Trigger Mode
fprintf([' SetTriggerMode           : ' num2str(acq.TriggerMode) ' ... ']);
[ret]=SetTriggerMode(acq.TriggerMode);
disp(error_code(ret))

% Set Accumulation Cycle Time
fprintf([' SetAccumulationCycleTime : ' num2str(acq.AccCycleTime) ' ... ']);
[ret]=SetAccumulationCycleTime(acq.AccCycleTime);
disp(error_code(ret))

% Set Number of Accumulations
fprintf([' SetNumberAccumulations   : ' num2str(acq.NumAcc) ' ... ']);
[ret]=SetNumberAccumulations(acq.NumAcc);
disp(error_code(ret))

% Set Number of Kinetics
fprintf([' SetNumberKinetics        : ' num2str(acq.NumKin) ' ... ']);
[ret]=SetNumberKinetics(acq.NumKin);
disp(error_code(ret))

% Set Kinetic Cycle Time
fprintf([' SetKineticCycleTime      : ' num2str(acq.KinCycleTime) ' ... ']);
[ret]=SetKineticCycleTime(acq.KinCycleTime);
disp(error_code(ret))

% Set binning and image acquisition area
fprintf([' SetImage                : ' ...
    num2str(acq.xbin) ',' num2str(acq.ybin) ',' ...
    num2str(acq.xstart) ',' num2str(acq.xend) ',' ...
    num2str(acq.ystart) ',' num2str(acq.yend) ' ... ']);
[ret]=SetImage(acq.xbin,acq.ybin,...
    acq.xstart,acq.xend,acq.ystart,acq.yend);
disp(error_code(ret))

% acq.validExposureTime = 0;
% acq.validKineticsCycleTime = 0;
% acq.validAccumulationCycleTime = 0;
% [ret,acq.validExposureTime,acq.validAccumulationCycleTime,...
%     acq.validKineticCycleTime]=GetAcquisitionTimings;
end

