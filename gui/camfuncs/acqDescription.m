function desc= acqDescription(acq)


fnames={...
    'AcquisitionMode', [1 2 3 4 5],{'Single Scan','Accumulate','Kinetics','Fast Kinetics','Run til abort'};
    'ReadMode',[0 1 2 3 4], {'Full Vertical Binning','Multi-Track','Random-Track','Single-Track'};
    'ShutterType', [0 1], {'low is open','high is open'};
    'ShutterMode', [0 1 2], {'Auto','Open','Close'};
    'ClosingTime', 'int','int';
    'OpeningTime', 'int','int';
    'FanMode', [0 1 2],{'on full','on low','off'};
    'ADChannelIndex', 0, {'0'};
    'PreAmpGainIndex', [0 1 2];
    'EMGainMode';
    'EMAdvanced';
    'EMCCDGain';
    'VSSpeedIndex';
    'AmpTypeIndex';
    'HSSpeedIndex';
    'ExposureTime';
    'TriggerMode';
    'AccCycleTime';
    'NumAcc';
    'NumKin';
    'KinCycleTime';
    'xbin';
    'ybin';
    'xstart';
    'xend';
    'ystart';
    'yend'};


f = fieldnames(A)';
f{2,1} = {};
B = struct(f{:});
end

