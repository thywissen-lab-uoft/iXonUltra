function desc= acqDescription(acq)

% Copy acuisition settings
desc=acq;

% Clear the value (to be populated with a descriptor)
fnames=fieldnames(desc);
for k=1:numel(fnames)
  desc.(fnames{k})='';
end

%% Go through each

% AcquisitionMode
strs={'Single Scan','Accumulate','Kinetics','Fast Kinetics','Run til abort'};
desc.AcquisitionMode=strs{acq.AcquisitionMode};

% ReadMode
strs={'Full Vertical Binning','Multi-Track','Random-Track','Single-Track','Image'};
desc.ReadMode=strs{acq.ReadMode+1};


% Disable shutter manipulation as setting that can be changed (always open
% or closed)
%{
% ShutterType
strs={'open low','open high'};
desc.ShutterType=strs{acq.ShutterType+1};

% ShutterMode
strs={'Auto','Open','Close'};
desc.ShutterMode=strs{acq.ShutterMode+1};

% ClosingTime
desc.ClosingTime='ms';

% OpeningTime
desc.OpeningTime='ms';
%}

% FanMode
strs={'on full','on low','off'};
desc.FanMode=strs{acq.FanMode+1};

% ADChannelIndex
desc.ADChannelIndex='amp 0';

% PreAmpGainIndex
strs={'x1','x2','x3'};
desc.PreAmpGainIndex=strs{acq.PreAmpGainIndex+1};

% EMGainMode
strs={'DAC 255','DAC 4095','Linear','Real EM'};
desc.EMGainMode=strs{acq.EMGainMode+1};

% EMAdvanced
strs={'yes >x300','no >x300'};
desc.EMAdvanced=strs{acq.EMAdvanced+1};

% EMCCDGain
desc.EMCCDGain='xN gain';

% AmpTypeIndex
strs={'EM/convential','convential/NIR'};
desc.AmpTypeIndex=strs{acq.AmpTypeIndex+1};

% VSSpeedIndex

% HSSpeedIndex

% TriggerMode
switch acq.TriggerMode
    case 0
        desc.TriggerMode='Internal';
    case 1
        desc.TriggerMode='External';
    case 6
        desc.TriggerMode='ExternalStart';
    case 7
        desc.TriggerMode='External Exposure';
    case 9
        desc.TriggerMOde='External FVB EM';
    case 10
        desc.TriggerMode='Software Trigger';
    case 12
        desc.TriggerMode='External Charge Shifting';
    otherwise
        desc.TriggerMode='??';
end

% ExposureTime
desc.ExposureTime='seconds';

% AccCycleTime
desc.AccCycleTime='seconds';

% KinCycleTime
desc.KinCycleTime='seconds';

% NumAcc
desc.NumAcc='img/scan';

% NumKin
desc.NumKin='scan/trig';

% xbin
desc.xbin='pixel(s)';

% ybin
desc.ybin='pixel(s)';

% xstart
desc.xstart='x px start';

% xend
desc.xend='x px end';

% ystart
desc.ystart='y px start';

% yend
desc.yend='y px end';
end

