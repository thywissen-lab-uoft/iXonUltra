function acq = defaultLiveAcqSettings
% THIS WILL NOT WORK AS I AM STILL CODING IT UP
acq=struct;

% Acquisition Modes
acq.AcquisitionMode=0;  % Kinetics
acq.ReadMode=0;         % Full Image

% Dont allow changing of the shutter from here
% Shutter
% acq.ShutterType=1;      % HIGH to open
% acq.ShutterMode=0;      % Auto
% acq.OpeningTime=5;      % 5 ms
% acq.ClosingTime=5;      % 5 ms

% Fan
acq.FanMode=0;          % Always on

% Gains
acq.ADChannelIndex=0;   % AD 0 (only one)
acq.PreAmpGainIndex=0;  % x1
acq.EMGainMode=0;       % Real EM
acq.EMAdvanced=0;       % Dont allow over 300
acq.EMCCDGain=0;      % 300 Gain

% Read Out
acq.AmpTypeIndex=0;     % Amp Type (EM)
acq.VSSpeedIndex=0;     % Speed Index 3: 1.7 us
acq.HSSpeedIndex=0;     % Speed Index 3: 1 MHz

% Timings
acq.TriggerMode=0;      % External
acq.ExposureTime=0;     % 2 seconds
acq.AccCycleTime=0;     % none
acq.KinCycleTime=0;     % none
acq.NumAcc=1;           % 1 accumulation
acq.NumKin=2;           % two images (PWA,bkgd)
acq.xbin=1;             % no x binning
acq.ybin=1;             % no y binning
acq.xstart=1;           % hardware ROI is full 1:512
acq.xend=512;           % hardware ROI is full 1:512
acq.ystart=1;           % hardware ROI is full 1:512
acq.yend=512;           % hardware ROI is full 1:512
end

