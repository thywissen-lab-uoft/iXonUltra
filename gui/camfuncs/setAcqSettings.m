function out = setAcqSettings(acq)

if ~isValidAcq(acq)
    warning('Bad acquisition settings provided. Aborting.');
    return;
end
      

% Set Acquisition Mode 
% 1: Single Scan, 2: Accumulate, 3: Kinetics, 4: Fast Kinetics, 5: Run till abort
[ret]=SetAcquisitionMode(acq.AcquisitionMode);

% Set Read Mode
% 0: Full Vertical Binning 1: Multi-Track, 2: Random-Track, 3: Single-Track, 4: Image
[ret]=SetReadMode(acq.ReadMode);

% Set Shutter Mode
% Type : (0: LOW to open, 1: HIGH to open), Mode: (0: Auto, 1: Open, 2:Close), Times in ms
[ret]=SetShutter(acq.ShutterType, acq.ShutterMode,acq.ClosingTime, acq.OpeningTime);

% Set Fan Mode
% 0: on Full, 1: on low, 2: off
[ret]=SetFanMode(acq.FanMode);

% Set AD Channel
[ret]=SetADChannel(acq.ADChannelIndex);

% Set PreAmp Gain 
[ret]=SetPreAmpGain(acq.PreAmpGainIndex);

% Set EM Gain Mode 
% 0: DAC 0-255, 1: DAC 0-4095, 2: Linear, 3: Real EM
[ret]=SetEMGainMode(acq.EMGainMode);

% Enable/Disable High EM Gains
[ret]=SetEMAdvanced(acq.EMAdvanced);

% Set EM Gain
[ret]=SetEMCCDGain(acq.EMCCDGain);

% Set Vertical Shift Speed
[ret]=SetVSSpeed(acq.VSSpeedIndex);

% Set Horizontal Shift Speed
[ret]=SetHSSpeed(acq.AmpTypeIndex,acq.HSSpeedIndex);

% Set Exposure Time 
[ret]=SetExposureTime(acq.ExposureTime);

% Set Trigger Mode
[ret]=SetTriggerMode(acq.TriggerMode);

% Set Accumulation Cycle Time
[ret]=SetAccumulationCycleTime(acq.AccCycleTime);

% Set Number of Accumulations
[ret]=SetNumberAccumulations(acq.NumAcc);

% Set Number of Kinetics
[ret]=SetNumberKinetics(acq.NumKin);

% Set Kinetic Cycle Time
[ret]=SetKineticCycleTime(acq.KinCycleTime);

% Set binning and image acquisition area
[ret]=SetImage(acq.xbin,acq.ybin,...
    acq.xstart,acq.xend,acq.ystart,acq.yend);

% acq.validExposureTime = 0;
% acq.validKineticsCycleTime = 0;
% acq.validAccumulationCycleTime = 0;
% [ret,acq.validExposureTime,acq.validAccumulationCycleTime,...
%     acq.validKineticCycleTime]=GetAcquisitionTimings;
end

