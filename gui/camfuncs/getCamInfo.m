function nfo = getCamInfo
% Acquire information of the camera.
%
% Note that this is independent of the current operating conditions and
% only needs to be acquired once upon connecting to the camera


nfo=struct;

%%%%%% Get basica camera info %%%%%%
[ret,nfo.HeadModel]=GetHeadModel;
[ret,nfo.SerialNumber]=GetCameraSerialNumber;
[ret,nfo.xpixels,nfo.ypixels]=GetDetector;
[ret,nfo.PCBVersion,nfo.FlexVersion,~,~,nfo.FirmwareVersion,nfo.FirmwareBuild]=GetHardwareVersion;
[ret,nfo.EPROMVersion,nfo.COFVersion,nfo.DriverRev,nfo.DriverVersion,nfo.DLLRev,nfo.DLLVersion]=GetSoftwareVersion;
[ret,nfo.PixelSize(1),nfo.PixelSize(2)]=GetPixelSize;
%%%%%% Get Settable Settings %%%%%%

% Shutter minimum times
[ret,nfo.MinShutterTime(1),nfo.MinShutterTime(2)]=GetShutterMinTimes;

% Sensor temperature range
[ret,nfo.TempRange(1),nfo.TempRange(2)]=GetTemperatureRange;

% EM Gain range
[ret,nfo.EMGainRange(1),nfo.EMGainRange(2)]=GetEMGainRange;

% Availabe vertical shifts speeds (us) (0.3, 0.5, 0.9, 1.7, 3.3 µs)
[ret, nfo.NumVSSpeeds] = GetNumberVSSpeeds;
nfo.AvailableVSSpeeds = zeros(1,nfo.NumVSSpeeds);
[ret, nfo.NumVSSpeeds] = GetNumberVSSpeeds;
for j = 0:nfo.NumVSSpeeds-1
    [ret,nfo.AvailableVSSpeeds(j+1)] = GetVSSpeed(j);
end

% Get available horizontal shift speeds (for each available combination af
% Amplfier (0: EMCCD, 1: conventional) and ADChannel (only one channel for
% our iXon Ultra). HSSpeeds are given in MHz (see Performance sheet)
% For our iXon Ultra camera (single ADC available)
% EMCCD: HSS = 17, 10, 5, 1MHz
% conventional: 3, 1, 0.08 MHz
[ret, nfo.NumAmp] = GetNumberAmp;
[ret, nfo.NumADChannels] = GetNumberADChannels;
for amp = 0:nfo.NumAmp-1
    for adc = 0:nfo.NumADChannels-1
        [ret, nfo.NumHSSpeeds{amp+1,adc+1}] = GetNumberHSSpeeds(adc,amp);
        nfo.AvailableHSSpeeds{amp+1,adc+1} = zeros(1,nfo.NumHSSpeeds{amp+1,adc+1});
        for j = 0:nfo.NumHSSpeeds{amp+1,adc+1}-1
            [ret,nfo.AvailableHSSpeeds{amp+1,adc+1}(j+1)] = GetHSSpeed(adc,amp,j);
        end
    end
end


% PreAmp Gains (all available for both amplifier types). (1,2,3)
[ret, nfo.NumPreAmpGains] = GetNumberPreAmpGains; 
nfo.AvailablePreAmpGains = zeros(1,nfo.NumPreAmpGains);
for j = 0:nfo.NumPreAmpGains-1
    [ret,nfo.AvailablePreAmpGains(j+1)] = GetPreAmpGain(j);
end



% [ret,cam.currentTemp]=GetTemperature;
% % [ret,cam.PreAmpGain]=GetPreAmpGain(cam.PreAmpGainIndex);
% [ret,cam.VSSpeed]=GetVSSpeed(cam.VSSpeedIndex);
% [ret,cam.HSSpeed]=GetHSSpeed(cam.ADChannelIndex,cam.AmpTypeIndex,cam.HSSpeedIndex);
%  [ret,cam.EMCCDGain]=GetEMCCDGain;


end