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
strs={'no >x300','yes >x300'};
desc.EMAdvanced=strs{acq.EMAdvanced+1};

% EMCCDGain
desc.EMCCDGain='xN gain';

% AmpTypeIndex
strs={'EM/conventonial','conventional/NIR'};
desc.AmpTypeIndex=strs{acq.AmpTypeIndex+1};

% VSSpeedIndex
[~,sp] = GetVSSpeed(acq.VSSpeedIndex);
desc.VSSpeedIndex = [num2str(sp) ' us'];

% HSSpeedIndex

[~,sp]=GetHSSpeed(acq.ADChannelIndex,acq.AmpTypeIndex,acq.HSSpeedIndex);
desc.HSSpeedIndex = [num2str(sp) ' MHz'];
% keyboard
% [ret, nfo.NumAmp] = GetNumberAmp;
% [ret, nfo.NumADChannels] = GetNumberADChannels;
% for amp = 0:nfo.NumAmp-1
%     for adc = 0:nfo.NumADChannels-1
%         [ret, nfo.NumHSSpeeds{amp+1,adc+1}] = GetNumberHSSpeeds(adc,amp);
%         nfo.AvailableHSSpeeds{amp+1,adc+1} = zeros(1,nfo.NumHSSpeeds{amp+1,adc+1});
%         for j = 0:nfo.NumHSSpeeds{amp+1,adc+1}-1
%             [ret,nfo.AvailableHSSpeeds{amp+1,adc+1}(j+1)] = GetHSSpeed(adc,amp,j);
%         end
%     end
% end


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

