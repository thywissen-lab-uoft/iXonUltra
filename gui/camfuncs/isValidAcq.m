function isValid = isValidAcq(acq)
isValid=1;

% Field Names
fnames={'AcquisitionMode';
    'ReadMode';
%     'ShutterType'; % remove changing of shutter mode
%     'ShutterMode'
%     'ClosingTime';
%     'OpeningTime';
    'FanMode';
    'ADChannelIndex';
    'PreAmpGainIndex';
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

allFields=1;

% Check if every field is present
for kk=1:length(fnames)
    if ~ismember(fnames{kk},fieldnames(acq))
       allFields=0;
       warning([fnames{kk} ' not provided']);
    end    
end

if ~allFields
    isValid=0; 
end

%% Check if valid value for each field

try

    % Acquisition mode
    if isempty(find(acq.AcquisitionMode==[1 2 3 4 5],1)); warning('Bad AcquisitionMode'); isValid=0; end

    % Read mode
    if isempty(find(acq.ReadMode==[1 2 3 4 5],1)); warning('Bad ReadMode');isValid=0; end

    % Remove shutter checking as we disable this as a user option
    %{
    
    % Shutter Type
    if isempty(find(acq.ShutterType==[0 1],1)); warning('Bad ShutterType'); isValid=0; end

    % Shutter Mode
    if isempty(find(acq.ShutterMode==[0 1 2],1)); warning('Bad ShutterMode'); isValid=0; end

    % Closing Time
    if ~(isnumeric(acq.ClosingTime) && isreal(acq.ClosingTime) && ...
            rem(acq.ClosingTime,1)==0 && acq.ClosingTime>=0)
        warning('Bad ClosingTime'); 
        isValid=0;
    end

    % Opening Time
    if ~(isnumeric(acq.OpeningTime) && isreal(acq.OpeningTime) && ...
            rem(acq.OpeningTime,1)==0 && acq.OpeningTime>=0)
        warning('Bad ClosingTime'); 
        isValid=0;
    end
    
    %}

    % FanMode
    if isempty(find(acq.FanMode==[0 1 2],1)); warning('Bad FanMode'); isValid=0; end

    % ADChannelIndex (hardcode it to 0 but this could be read)
    if ~isequal(acq.ADChannelIndex,0); warning('Bad ADChannelIndex'); isValid=0; end

    % PreAmpGainIndex (hardcode it to 0,1,2 but this could be read)
    if isempty(find(acq.PreAmpGainIndex==[0 1 2],1)); warning('Bad PreAmpGainIndex'); isValid=0; end

    % EMGainMode
    if isempty(find(acq.EMGainMode==[0 1 2 3],1)); warning('Bad EMGainMode'); isValid=0; end

    % EMAdvanced
    if isempty(find(acq.EMAdvanced==[0 1],1)); warning('Bad EMAdvanced'); isValid=0; end

    % EMCCDGain
    if ~(isnumeric(acq.EMCCDGain) && isreal(acq.EMCCDGain) && ...
            rem(acq.EMCCDGain,1)==0 && acq.EMCCDGain>=1 && acq.EMCCDGain<=1000)
        warning('Bad EMCCDGain');
        return;
    end

    % Check if gain is sound
    if ~acq.EMAdvanced && acq.EMCCDGain>300
        warning('Gain greater than x300 need EMAdvanced engaged.');
        isValid=0;
    end

    % VSSpeedIndex (hardcode to [0 1 2], but this could be read)
    if isempty(find(acq.VSSpeedIndex==[0 1 2 3],1)); warning('Bad VSSpeedIndex'); isValid=0; end

    % AmpTypeIndex
    if isempty(find(acq.AmpTypeIndex==[0 1],1)); warning('Bad AmpTypeIndex'); isValid=0; end

    % HSSpeedIndex
    if isempty(find(acq.HSSpeedIndex==[0 1 2 3],1)); warning('Bad HSSpeedIndex'); isValid=0; end

    % ExposureTime
    if ~(isnumeric(acq.ExposureTime) && isreal(acq.ExposureTime) && ...
            acq.ExposureTime>0)
        warning('Bad ExposureTime');
        isValid=0;
    end

    % TriggerMode
    if isempty(find(acq.TriggerMode==[0 1 6 7 9 10 12],1)); warning('Bad TriggerMode'); isValid=0; end

    % AccCycleTime
    if ~(isnumeric(acq.AccCycleTime) && isreal(acq.AccCycleTime) && ...
            acq.AccCycleTime>=0)
        warning('Bad AccCycleTime');
        isValid=0;
    end

    % NumAcc
    if ~(isnumeric(acq.NumAcc) && isreal(acq.NumAcc) && ...
            rem(acq.NumAcc,1)==0 && acq.NumAcc>=0)
        warning('Bad NumAcc'); 
        isValid=0;
    end

    % NumKin
    if ~(isnumeric(acq.NumKin) && isreal(acq.NumKin) && ...
            rem(acq.NumKin,1)==0 && acq.NumKin>=0)
        warning('Bad NumKin'); 
        isValid=0;
    end

    % KinCycleTime
    if ~(isnumeric(acq.KinCycleTime) && isreal(acq.KinCycleTime) && ...
            acq.KinCycleTime>=0)
        warning('Bad KinCycleTime');
        isValid=0;
    end

    % xbin arbitrarily restricting to [1 2 3 4]
    if isempty(find(acq.xbin==1:4,1)); warning('Bad xbin'); isValid=0; end

    % ybin arbitrarily restricting to [1 2 3 4]
    if isempty(find(acq.ybin==1:4,1)); warning('Bad xbin'); isValid=0; end

    % x limits (harcode to 512x512 sensor)
    if ~(acq.xend>acq.xstart && acq.xend<=512 && acq.xstart>=1); warning('Bad xlimits'); isValid=0; end

    % x limits (harcode to 512x512 sensor)
    if ~(acq.yend>acq.ystart && acq.yend<=512 && acq.ystart>=1); warning('Bad ylimits'); isValid=0; end

catch ME
    disp(ME)
    warning('Error on checking acquisition settings');
    isValid=0;
end
end
