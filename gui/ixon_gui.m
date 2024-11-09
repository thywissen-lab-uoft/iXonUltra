function ixon_gui(doDebug)
% ixon_gui.m
%
% Author      : C. Fujiwara
%
% This code operates the iXon Ultra camera that the lattice experiment
% uses to take fluorescence images of the quantum gas microscope.

if nargin == 0;doDebug = 0;end

if doDebug
   warning(['You are operating in DEBUG MODE. This removes ' ...
       'certain safety precautions. If not intended set the doDebug ' ...
       'flag to 0.']); 
end

% Manification (as measured in Graham's thesis)
mag=[82.6 83.2];

% Load the mask file
analysis_path = fullfile(fileparts(fileparts(mfilename('fullpath'))),'analysis');
maskname=fullfile(analysis_path,'ixon_mask.mat');

ixon_mask=load(maskname);
ixon_mask=ixon_mask.BW;

% Add analysis paths
addpath(analysis_path);addpath(genpath(analysis_path))

camera_control_file='Y:\_communication\camera_control.mat';

%% Other Settings

cmap=purplemap;                                 % default colormap
guiname='iXon GUI';                             % Figure Name
defaultDir=['C:' filesep 'IxonImageHistory'];   % Temporary save directory

doSaveGUIAnalysis = 1;
GUIAnalysisSaveDir = 'X:\IxonGUIAnalysisHistory';
GUIAnalysisLocalDir = 'C:\IxonGUIAnalysisHistory';

currDir=defaultDir;                             % Current directory of navigator is the default one
if ~exist(defaultDir,'dir'); mkdir(defaultDir);end % Make temp directory if necessary

% Dummy file to load on startup
fname='example_data_EIT_RAMAN.mat';
data=load(fname);
data=data.data;
data.Z=data.RawImages(:,:,2)-data.RawImages(:,:,1);
Z=data.Z;




%% Initialize Drivers and GUI

% Add the Andor MATLAB drivers to the MATLAB path. You need to have
% installed the MATLAB drivers in order for this to work.
dir_andor = fullfile(matlabroot,'toolbox','Andor');
if ~exist(dir_andor,'dir')
    warning(['You have not installed the Andor toolbox ' ...
        ' for running the camera. Certain capabilities will not ' ...
        'work.']);
else
    addpath(genpath(dir_andor));
end

% Add all subdirectories for this m file
mpath = fileparts(mfilename('fullpath'));
addpath(mpath);addpath(genpath(mpath))
% Find any instances of the GUI and bring it to focus, this is tof avoid
% restarting the GUI which may leave the shutter open.
h = findall(0,'tag','GUI');
for kk=1:length(h)    
    if isequal(h(kk).Name,guiname)        
        warning(['iXon GUI instance detected.  Bringing into focus. ' ...
            ' If you want to start a new instance, close the original iXon GUI.']); 
       figure(h(kk));
       return;
    end    
end
disp('Initializing iXon GUI...');

%% Load In Images
img_Snap    = imresize(imread(fullfile(mpath,'icons','snapLim.png')),[15 15]);
img_Select  = imresize(imread(fullfile(mpath,'icons','target.jpg')),[15 15]);
img_Full    = imresize(imread(fullfile(mpath,'icons','fullLim.png')),[15 15]);


%% Camera Settings
% Initialization camera status structure

% Declare cam_info struct;
cam_info=struct;

% Initialize Camera Status
cam_status=struct;
cam_status.isConnected=0;       % Are you connected to camera?
cam_status.Temperature=NaN;     % sensor temperature
cam_status.TemperatureSP=-80;   % temperature set point
cam_status.isTempStable=0;      % is sensor tempreature stable?
cam_status.isCooling=0;         % is the TEC active?4
cam_status.isAcquiring=0;       % is the camera acquiring frames?

acq=defaultNormalAcqSettings;   % load default settings
desc=acqDescription(acq);       % Get description of settings
%% Initialize Figure

% Initialize the primary figure
hF=figure;clf
set(hF,'Color','w','units','pixels','Name',guiname,...
    'Tag','GUI','CloseRequestFcn',@closeGUI,'NumberTitle','off',...
    'Position',[50 50 1200 630],'SizeChangedFcn',@SizeChangedFcn,...
    'toolbar','none');

topBarHeight = 40;

% Callback for when the GUI is requested to be closed.
    function closeGUI(fig,~)
        answer = questdlg('Close the iXon GUI?','Close iXon?',...
            'Yes','No','No') ;        
        switch answer
            case 'Yes'
                disp('Closing the iXon GUI ...')
                doClose = 1;
            case 'No'
                disp('Not closing the iXon GUI.')
                doClose = 0;
            otherwise
                doClose = 0;
        end        
        if doClose
            disp('Closing iXon GUI...');      
            stop(statusTimer);
            if cam_status.isConnected;ixon_disconnectCamera;end
            delete(statusTimer);
            delete(fig);      % Delete the figure          
        end
    end

% Change the figure picture (will be depreciated)
warning off
javaFrame = get(hF,'JavaFrame');
javaFrame.setFigureIcon(javax.swing.ImageIcon(fullfile(mpath,'icons','ixon_pic.PNG')));
warning on

function SizeChangedFcn(~,~)
        % This resize fucntion ensures that the X and Y cut/sum plot has
        % commenserate positioning with respect the actual image shown
        x0=hpControl.Position(1)+hpControl.Position(3);
        Wfig=hF.Position(3);Hfig=hF.Position(4);                          % Figure size
        
        if (Wfig>360 && Hfig>55)
            hp.Position=[x0 1 max([hF.Position(3)-x0 5]) Hfig]; 
        end             
        resizePlots;                                                % Resize plots                          
        drawnow;

        % Control Panel
        hpControl.Position(4) = hF.Position(4);

        hpCam.Position(2) = hpControl.Position(4)-hpCam.Position(4);
        hpAcq.Position(2) = hpCam.Position(2)-hpAcq.Position(4);
        % hpSave.Position(2) = hpAcq.Position(2)-hpSave.Position(4);
        hpNav.Position(2)       = hpAcq.Position(2)-hpNav.Position(4);

 
        hpProcess.Position(2)           = hpNav.Position(2)-hpProcess.Position(4);

        hpFeedback.Position(2)      = hpProcess.Position(2)-hpFeedback.Position(4);
        hpPosition.Position(2)       = hpFeedback.Position(2)-hpPosition.Position(4);   

        % hpKspace.Position(2)    = hpPosition.Position(2)-hpKspace.Position(4);

        hpKspace.Position(2)   = hpNav.Position(2)-hpKspace.Position(4);
      hpBin.Position(2)       = hpKspace.Position(2) - hpBin.Position(4);   

        % hpBin.Position(2)       = hpNav.Position(2) - hpBin.Position(4);   
        hpDig.Position(2)       = hpBin.Position(2) - hpDig.Position(4);   

        hpDisp_Select.Position(2)=hpNav.Position(2)-hpDisp_Select.Position(4);
        hpDispOpt.Position(2)=hpDisp_Select.Position(2)-hpDispOpt.Position(4);
        if hpDispOpt.Position(2)>0
            hpFit.Position(4)=hpDispOpt.Position(2);
        else
            hpFit.Position(4) = 1;
        end


        hpFit.Position(4) = max([hpDispOpt.Position(2) 5]);

        % strstatus.Position(1)   = hpCam.Position(3)-strstatus.Position(3)-2;        
        drawnow;       
end
%% Control Panel Container
hpControl = uipanel(hF,'units','pixels');
hpControl.Position = [0 0 600 hF.Position(4)];

%% Camera Panel

hpCam=uipanel(hpControl,'units','pixels','backgroundcolor','w','title','camera',...
    'fontsize',6);
hpCam.Position(3) = hpControl.Position(3);
hpCam.Position(4) = 30;
hpCam.Position(1) = 0;
hpCam.Position(2) = hpControl.Position(4) - hpCam.Position(4);

% Connect camera
ttstr='Connect to the iXon camera.';
hbConnect=uicontrol(hpCam,'style','pushbutton','string','connect','units','pixels',...
    'fontsize',7,'Position',[2 1 45 18],'backgroundcolor',[80 200 120]/255,...
    'Callback',@connectCB,'ToolTipString',ttstr);

% Callback for the connect button
    function connectCB(~,~)
        strstatus.String='CONNECTING';
        drawnow;
        out=ixon_connectCam;                 % Connect to the camera       
        if ~out && ~doDebug             % Give warning if connection fails
           warning('Unable to connect to camera');
           return;
        end     
        strstatus.String='CONNECTED';
        drawnow;       
        cam_status.isConnected=1;       
        ixon_setCameraShutter(0);            % Close the shutter
        loadAcquisitionSettings;        % Load default acquisition settings      
        ixon_setCameraShutter(0);            % Close the shutter (again to be safe)
        hbOpenShutter.Enable='on';      % allow shutter to be opened
        cam_info=getCamInfo;            % Get the camera information
        disp(' ');
        disp(cam_info);        
        
        % Set the temperature set point
        ixon_setTemperature(tblTemp.Data);
        hbCamInfo.Enable='on';       
        hbCamAbilities.Enable='on';
        
        % Enable/Disable connect
        hbDisconnect.Enable='on';
        hbConnect.Enable='off';   

        % Enable/Disable temperature
        hbCool.Enable='on';
        tblTemp.Enable='on';
        strtemp.Enable='on';
        start(statusTimer);

        % Enable acquisition
        hbstart.Enable='on';
        rbSingle.Enable='on';
        rbLive.Enable='on';
        tbl_acq.Enable='on';        
    end

% Disconnect camera
ttstr='Disconnect to the iXon camera.';
hbDisconnect=uicontrol(hpCam,'style','pushbutton','string','disconnect','units','pixels',...
    'fontsize',7,'Position',[hbConnect.Position(1)+hbConnect.Position(3) 1 55 18],'backgroundcolor',[255 102 120]/255,...
    'Callback',@disconnectCB,'ToolTipString',ttstr,'enable','off');

% Callback for the disconnect button
    function disconnectCB(~,~)
        % Stop temperature monitor
        stop(statusTimer);
        
        strstatus.String='DISCONNECTING';
        drawnow;
        % Disconnect from camera
        out=ixon_disconnectCamera; 
       
        if ~out && ~doDebug
           return;
        end
        
        strstatus.String='DRV_NOT_INITIALIZED';
        drawnow;
          
        cam_status.isConnected=0;
        cam_status.isCooling=0;
        cam_status.isTempStable=0;
        cam_status.Temperature=NaN;
        
        hbDisconnect.Enable='off';
        hbConnect.Enable='on'; 
        hbCool.Enable='off';
        hbCoolOff.Enable='off';
        tblTemp.Enable='off';
        hbCamInfo.Enable='off';
        hbCamAbilities.Enable='off';
        
        hbOpenShutter.Enable='off';
        hbCloseShutter.Enable='off';
        rbSingle.Enable='off';
        rbLive.Enable='off';
        
        hbstart.Enable='off';
        hbstop.Enable='off';
        tbl_acq.Enable='off';

        set(strtemp,'ForegroundColor','r','String','NaN','Enable','off');
    end

% Info button
ttstr='Display camera info';
cdata=imresize(imread(fullfile(mpath,'icons','info.jpg')),[16 16]);
hbCamInfo=uicontrol(hpCam,'style','pushbutton','CData',cdata,'callback',@infoCB,...
    'enable','on','backgroundcolor','w','position',[130 1 18 18],...
    'ToolTipString',ttstr,'enable','off');
hbCamInfo.Position(1)= hbDisconnect.Position(1)+hbDisconnect.Position(3);

    function infoCB(~,~)        
        cam_info=getCamInfo;        
        disp(cam_info); 
        disp(cam_info.NumHSSpeeds)
        disp(cam_info.AvailableHSSpeeds{1});
        disp(cam_info.AvailableHSSpeeds{2});
    end

% Capabilities button
ttstr='Display camera capabilities';
cdata=imresize(imread(fullfile(mpath,'icons','infob.jpg')),[16 16]);
hbCamAbilities=uicontrol(hpCam,'style','pushbutton','CData',cdata,'callback',@abilitiesCB,...
    'enable','on','backgroundcolor','w','position',[152 1 18 18],...
    'ToolTipString',ttstr,'enable','off');
hbCamAbilities.Position(1)= hbCamInfo.Position(1)+hbCamInfo.Position(3);

    function abilitiesCB(~,~)  
        % Read Camera Capabilities
        cam_skills=getCameraCapabilities;
        % Dont display camera capabilities
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        disp('Displaying Camera Capabilities');        
        fnames=fieldnames(cam_skills);
        for nn=1:length(fnames)
           disp(fnames{nn});
           disp(cam_skills.(fnames{nn})); 
           disp(' ');
        end
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    end

% Start Cooling
ttstr='Begin cooling the sensor to set point.';
hbCool=uicontrol(hpCam,'style','pushbutton','string','cooler on',...
    'units','pixels','fontsize',7,'Position',...
    [175 1 50 18],'enable','off',...
    'backgroundcolor',[173 216 230]/255,'callback',{@coolCB ,1},...
    'ToolTipString',ttstr);
hbCool.Position(1)= hbCamAbilities.Position(1)+hbCamAbilities.Position(3);


% Stop Cooling
ttstr='Stop cooling the sensor to set point.';
hbCoolOff=uicontrol(hpCam,'style','pushbutton','string','cooler off',...
    'units','pixels','fontsize',7,'Position',...
    [238 1 50 18],'enable','off',...
    'backgroundcolor',[.8 .8 .8],'callback',{@coolCB, 0},...
    'ToolTipString',ttstr);
hbCoolOff.Position(1)= hbCool.Position(1)+hbCool.Position(3);

    function coolCB(~,~,state)  
        out=ixon_setTEC(state);        
        if ~out && ~doDebug
           return; 
        end        
        % Enable/Disable buttons
        if state
            hbCool.Enable='off';
            hbCoolOff.Enable='on';
            cam_status.isCooling=1;
        else
            hbCool.Enable='on';
            hbCoolOff.Enable='off';
            cam_status.isCooling=0;
        end        
    end

% Table set point
ttstr='Change the sensor temperature set point.';
tblTemp=uitable(hpCam,'units','pixels','ColumnWidth',{30},...
    'ColumnEditable',[true],'ColumnName',{},...
    'Data',cam_status.TemperatureSP,'FontSize',8,...
    'CellEditCallback',@chTempCB,'RowName',{},'ToolTipString',ttstr,...
    'enable','off');
% tblTemp.Position(3:4)=tblTemp.Extent(3:4);
tblTemp.Position(3:4)=[35 20];

tblTemp.Position(1:2)=[300 1];
tblTemp.Position(1)= hbCoolOff.Position(1)+hbCoolOff.Position(3);

    function chTempCB(src,evt)
        disp('Changing temperature set point');
        Told=evt.PreviousData;
        Tnew=evt.NewData;        
        % Check that the data is numeric
        if ~isnumeric(Tnew) || isinf(Tnew) || isnan(Tnew)
            warning('Incorrect data type provided for temperature.');
            src.Data=Told;
            return;
        end            
        Tnew=round(Tnew);        
        if Tnew<-120 || Tnew>20
            warning('Temperature set point out of bounds. Resetting.');
            src.Data=Told; 
           return; 
        end
        out=ixon_setTemperature(Tnew);        
        if ~out  && ~doDebug
           src.Data=Told; 
        end            
        src.Data=Tnew;        
    end

% Text
ttstr='Camera sensor temperature. (green: stable; yellow: unstable; red: set point not reached)';
strtemp=uicontrol(hpCam,'style','text','string','NaN','units','pixels',...
    'backgroundcolor','w','fontsize',12,'horizontalalignment','center',...
    'foregroundcolor','r','enable','off','fontweight','bold',...
    'ToolTipString',ttstr);
strtemp.Position(3:4)=[40 20];
strtemp.Position(1:2)=[340 1];
strtemp.Position(1)= tblTemp.Position(1)+tblTemp.Position(3);


% Open camera shutter
ttstr='Open camera shutter.';
hbOpenShutter=uicontrol(hpCam,'style','pushbutton','string','open shutter',...
    'units','pixels','fontsize',7,'Position',[385 1 65 18],'enable','off',...
    'backgroundcolor',[255 204 0]/255,'callback',{@shutterCB,1},...
    'ToolTipString',ttstr);
hbOpenShutter.Position(1)= strtemp.Position(1)+strtemp.Position(3);

ttstr='Close camera shutter.';
hbCloseShutter=uicontrol(hpCam,'style','pushbutton','string','close shutter',...
    'units','pixels','fontsize',7,'Position',[465 1 65 18],'enable','off',...
    'backgroundcolor',[255 102 120]/255,'callback',{@shutterCB,0},...
    'ToolTipString',ttstr);
hbCloseShutter.Position(1)= hbOpenShutter.Position(1)+hbOpenShutter.Position(3);


    function shutterCB(~,~,state)        
        if state && cam_status.Temperature>-60
            warning('Denying your request to open the shutter above -60C.');
            return;
        end        
        out=ixon_setCameraShutter(state);        
        % Exit if bad return
        if ~out && ~doDebug
           return; 
        end        
        % Enable/Disable buttons
        if state
            hbOpenShutter.Enable='off';
            hbCloseShutter.Enable='on';
        else
            hbOpenShutter.Enable='on';
            hbCloseShutter.Enable='off';
        end
    end


% Text
ttstr='Camera status.';
strstatus=uicontrol(hpCam,'style','text','string','DRV_NOT_INITIALIZED','units','pixels',...
    'backgroundcolor','w','fontsize',8,'horizontalalignment','left',...
    'foregroundcolor','k','enable','on','fontweight','bold',...
    'ToolTipString',ttstr);
strstatus.Position(3:4)=[157 16];
strstatus.Position(1) = hbCloseShutter.Position(1) +hbCloseShutter.Position(3);
strstatus.Position(2)=1;

% Timer to update temperature
statusTimer=timer('Name','iXonTemperatureTimer','Period',1,...
    'TimerFcn',@statusTimerFcn,'ExecutionMode','FixedSpacing');

    function statusTimerFcn(~,~)
        try
        % Get the temperature
        [out,temp,outstr]=ixon_getTemperature;
        strtemp.String=[num2str(temp) ' C'];        
        cam_status.Temperature=temp;   
        
        if cam_status.Temperature>-60 && isequal(hbCloseShutter.Enable,'on')
            warning('Shutter is open and temperature above -60. Closing shutter');
            shutterCB([],[],0);
        end
        
        switch outstr
            case 'DRV_TEMPERATURE_STABILIZED'
                cam_status.isTempStable=1;
                strtemp.ForegroundColor=[80 200 120]/255;
            case 'DRV_TEMP_NOT_STABILIZED'
                cam_status.isTempStable=0;
                strtemp.ForegroundColor=[255 204 0]/255;
            otherwise
                cam_status.isTempStable=0;
                strtemp.ForegroundColor='r';
        end     
        
        % Camera Status
        [out,outstr]=ixon_getCameraStatus;
        strstatus.String=outstr;
        drawnow;
        catch ME
            warning('status timer failed');
            warning(ME.message);
        end
    end


%% Image Acqusition Panel
% Panel for image acquisition controls and settings.
hpAcq=uipanel(hpControl,'units','pixels','backgroundcolor','w',...
    'title','acquisition','fontsize', 6);
hpAcq.Position(3) = hpControl.Position(3);
hpAcq.Position(4) = 50;
hpAcq.Position(1) = 0;
hpAcq.Position(2) = hpCam.Position(2) - hpAcq.Position(4);

% Start acquisition button
ttstr='Start acquisition.';
hbstart=uicontrol(hpAcq,'style','pushbutton','string','start',...
    'units','pixels','fontsize',8,'backgroundcolor',[80 200 120]/255,...
    'Position',[1 1 40 18],...    
    'Callback',@startCamCB,'ToolTipString',ttstr,'enable','off');
hbstart.Position(2) = hpAcq.Position(4)-hbstart.Position(4)-10;
    function startCamCB(src,evt)
        disp('Starting acquisition');        
        % Send acquistion start command
        out=ixon_startCamera;
        tic
        start(acqTimer);        
        % Enable/Disable Button/Tables
        rbSingle.Enable='off';
        rbLive.Enable='off';
        hbstart.Enable='off';
        hbstop.Enable='on';        
        tbl_acq.ColumnEditable(2)=false;
        tbl_acq.Enable='off';
    end

% Stop acquisition button
ttstr='Abort acquisition.';
hbstop=uicontrol(hpAcq,'style','pushbutton','string','stop',...
    'units','pixels','fontsize',7,'backgroundcolor',[250 102 120]/255,...
    'Position',[45 hbstart.Position(2) 40 18], ...
    'Callback',@stopCamCB,'ToolTipString',ttstr,'enable','off');
hbstop.Position(1) = hbstart.Position(1)+hbstart.Position(3);

    function stopCamCB(src,evt)
        disp('Stopping acquisition');     
        stop(acqTimer);
        ixon_stopCamera;        
        % Enable/Disable Button/Tables
        hbstart.Enable='on';
        hbstop.Enable='off';        
        rbSingle.Enable='on';
        rbLive.Enable='on';        
        tbl_acq.ColumnEditable(2)=true;
        tbl_acq.Enable='on';
    end

% Acquisition help button
ttstr='Help on acquisition.';
cdata=imresize(imread(fullfile(mpath,'icons','help.jpg')),[15 15]);
hbAcqInfo=uicontrol(hpAcq,'style','pushbutton','CData',cdata,'callback',@helpCB,...
    'enable','on','backgroundcolor','w','position',[88 hbstart.Position(2) 18 18],...
    'ToolTipString',ttstr,'enable','on');
hbAcqInfo.Position(1) =hbstop.Position(1) + hbstop.Position(3);

    function helpCB(~,~)
       disp('Im helping you'); 
    end
% Continuous Acquisition checkbox
ttstr='Reinitialize camera acquisition after image acquisition.';
hcAcqRpt=uicontrol(hpAcq,'style','checkbox','string','repeat acquisition?','fontsize',6,...
    'backgroundcolor','w','Position',[hbAcqInfo.Position(1)+hbAcqInfo.Position(3) hbstart.Position(2) 85 18],...
    'ToolTipString',ttstr,'enable','on','value',1);

% Button group for acquisition mode
bgAcq = uibuttongroup(hpAcq,'units','pixels','backgroundcolor','w','BorderType','None',...
    'SelectionChangeFcn',@chAcqCB);  
bgAcq.Position(3:4)=[120 18];
bgAcq.Position(1:2)=[hcAcqRpt.Position(1)+hcAcqRpt.Position(3) hbstart.Position(2)];    

% Radio buttons for triggered versus software acquisition
rbSingle=uicontrol(bgAcq,'Style','radiobutton','String','triggered',...
    'Position',[0 0 60 18],'units','pixels','backgroundcolor','w','Value',1,...
    'UserData','Normal','Enable','off','fontsize',7);
rbLive=uicontrol(bgAcq,'Style','radiobutton','String','software',...
    'Position',[60 0 60 18],'units','pixels','backgroundcolor','w',...
    'UserData','Live','Enable','off','fontsize',7);

% Change acqusition mode callback
    function chAcqCB(~,evt)        
        oldStr=evt.OldValue.UserData;
        newStr=evt.NewValue.UserData;     
        
        disp(['Changing camera acquisition mode to "' newStr '"']);

        if isequal(oldStr,newStr)
           disp('The old and new setting are the same.  How did you get here?'); 
           return
        end  
        
        % Close the shutter
        disp('Closing shutter in case you forgot.');                
        ixon_setCameraShutter(0);        
        hbCloseShutter.Enable='off';
        hbOpenShutter.Enable='on';
        
        switch newStr
            case 'Live'
                msg=['Entering live mode. Verify your settings before ' ...
                    'starting acquisition. You can break the camera.'];
                msgbox(msg,'Live Mode','warn','modal');  
                acq=defaultLiveAcqSettings;    
                loadAcquisitionSettings;
            case 'Normal'
                acq=defaultNormalAcqSettings;
                loadAcquisitionSettings;
            otherwise
                warning('Unexpected acqusition mode. What happened?');
        end        
        
    end

% Auto Save check box
ttstr=['Use the camera control file to set the camera settings?'];
hcAdwinCamera=uicontrol(hpAcq,'style','checkbox','string','auto camera config?','fontsize',7,...
    'backgroundcolor','w','Position',[bgAcq.Position(1)+bgAcq.Position(3) bgAcq.Position(2) 110 18],...
    'ToolTipString',ttstr);




% Timer to check on acquisition
% boop = struct;
% boop.ReadTime = [];
% boop.
acqTimer.UserData=[0 0];
acqTimer=timer('Name','iXonAcquisitionWatchTimer','Period',1,...
    'TimerFcn',@acqTimerFcn,'ExecutionMode','FixedSpacing','StartFcn',@acqTimerStartFcn);

    function acqTimerStartFcn(src,evt)
        src.UserData = [0 0];
    end

    function autoCameraConfig
        if exist(camera_control_file,'file') && cam_status.isConnected 
            try
                CameraControl = load(camera_control_file);
                if ~isfield(CameraControl,'IxonMultiExposures')
                    return;
                end
                if acq.NumKin~=length(CameraControl.IxonMultiExposures)    
                    disp('incompatible number of exposures detected. Automatically changing')
                    stopCamCB;
                    pause(0.1);
                    acq.NumKin = length(CameraControl.IxonMultiExposures);
                    loadAcquisitionSettings;
                    pause(0.1);
                    startCamCB;
                end           
            end
        end  
    end

    function autoSaveDir
        if exist(camera_control_file,'file')
            try 
                CameraControl = load(camera_control_file);
                if ~isfield(CameraControl,'SaveDir')
                    return; 
                end

                if isempty(CameraControl.SaveDir)
                    return;
                end
                dirToday=ixon_getDayDir;
                
                SaveDir = fullfile(dirToday,['ixon_' CameraControl.SaveDir]);
                tSaveDir.UserData = SaveDir;

                str=strsplit(SaveDir,filesep);
                str=[str{end-1} filesep str{end}];
                tSaveDir.String=str; % display string
            end
        end

    end

% Timer function checks if acquisition is over and restarts it
    function acqTimerFcn(src,evt)      
        try
          
        % Camera Status
        [out,outstr]=ixon_getCameraStatus;  
        % Check number of available images and post and upate
        [ret,first,last] = GetNumberAvailableImages;                    
        if last~=src.UserData(1)
            src.UserData(1)=last;
            disp([datestr(now) ' ixon image ' num2str(last)])
        end
                    
        switch outstr
            case 'DRV_IDLE'      
                % Grab the images from the camera
                imgs=ixon_grabRawImages; 
                toc                
                % Assign images metadata
                mydata=makeImgDataStruct(imgs);
                % Restart Acquisition if desired (auto-stopts)
                if hcAcqRpt.Value
                    ixon_startCamera;  
                else
                    stop(src);
                    % Enable/Disable Button/Tables
                    hbstart.Enable='on';
                    hbstop.Enable='off';        
                    rbSingle.Enable='on';
                    rbLive.Enable='on';        
                    tbl_acq.ColumnEditable(2)=true;
                    tbl_acq.Enable='on';                    
                end                  
                % Save data to image history
                % only save to history if not in live mode
                if ~rbLive.Value
                    saveData(mydata);
                end
                % Save images to save directory
                if hcSave.Value
                   saveData(mydata,tSaveDir.UserData); 
                end              
                % Update live preview if new                
                if ~cAutoUpdate.Value
                    currDir=defaultDir;
                    data=mydata;                       
                    newDataCallback;                     

                else                    
                    % Just update index
                    updateHistoryInd(data);   
                end                  
            case 'DRV_ACQUIRING'
                % Acquisition is still going. 
                if hcAdwinCamera.Value
                    autoCameraConfig;
                end

                if hcAdwinSaveDir.Value
                    autoSaveDir;
                end
                
                 % if h
            otherwise
                warning('Acuisition timer has unexpected result');
                stopCamCB;
        end        
        catch ME
            warning('Acqtimer failed.');
            for tt=1:length(ME.stack)
                disp([ME.stack(tt).file ' (' num2str(ME.stack(tt).line) ']']);        
            end
        end
    end

% Auto Save check box
ttstr=['Enable/Disable automatic saving to external directory. Does ' ...
    'not override saving to image history.'];
hcSave=uicontrol(hpAcq,'style','checkbox','string','save?','fontsize',7,...
    'backgroundcolor','w','Position',[0 hbstart.Position(2)-18 50 18],...
    'ToolTipString',ttstr);

ttstr=['Use the camera control file to set the save directory?'];
hcAdwinSaveDir=uicontrol(hpAcq,'style','checkbox','string','CurrentJob SaveDir','fontsize',7,...
    'backgroundcolor','w','Position',[hcSave.Position(1)+hcSave.Position(3) hcSave.Position(2) 110 18],...
    'ToolTipString',ttstr);

% Browse button
ttstr='Select directory to save images.';
cdata=imresize(imread(fullfile(mpath,'icons','browse.jpg')),[15 15]);
bBrowse=uicontrol(hpAcq,'style','pushbutton','CData',cdata,'callback',@browseCB,...
    'enable','on','backgroundcolor','w','position',[95 hcSave.Position(2) 18 18],...
    'tooltipstring',ttstr);
bBrowse.Position(1) = hcAdwinSaveDir.Position(1) + hcAdwinSaveDir.Position(3);
% String for current save directory
ttstr='The current save directory.';
tSaveDir=uicontrol(hpAcq,'style','text','string','save directory','fontsize',8,...
    'backgroundcolor','w','units','pixels','horizontalalignment','left',...
    'enable','on','UserData','','Position',[bBrowse.Position(1)+bBrowse.Position(3) hcSave.Position(2)-5 hF.Position(3)-135 20],...
    'tooltipstring',ttstr);

% Browse button callback
    function browseCB(~,~)
        str=ixon_getDayDir;
        str=uigetdir(str);        
        if str
            tSaveDir.UserData=str; % Full directory to save
            str=strsplit(str,filesep);
            str=[str{end-1} filesep str{end}];
            tSaveDir.String=str; % display string
        else
            disp('no directory chosen!');
        end
    end



%% Navigator Panel 

hpNav=uipanel(hpControl,'units','pixels','backgroundcolor','w',...
    'title','GUI Image Source','fontsize',6);
hpNav.Position(3) = hpControl.Position(3);
hpNav.Position(4) = 50;
hpNav.Position(1) = 0;
hpNav.Position(2) = hpAcq.Position(2) - hpNav.Position(4);

% Button to change navigator directory to default
ttstr='Revert previewer source directory to default location.';
cdata=imresize(imread('icons/home.jpg'),[17 17]);
hbhome=uicontrol(hpNav,'style','pushbutton','CData',cdata,...
    'callback',@defaultDirCB,'enable','on','backgroundcolor','w',...
    'position',[1 18 20 20],'ToolTipString',ttstr);

% Change directory to default
    function defaultDirCB(~,~)
        disp(['Changing previwer directory to ' defaultDir]);
        currDir=defaultDir;
        chData([],[],0);        
    end

% Button to change preview source directory
ttstr='Change previwer source directory.';
cdata=imresize(imread('icons/browse.jpg'),[20 20]);
hbchdir=uicontrol(hpNav,'style','pushbutton','CData',cdata,'callback',@chDirCB,...
    'enable','on','backgroundcolor','w','position',[hbhome.Position(1)+hbhome.Position(3) 18 20 20],...
    'ToolTipString',ttstr);

% Get directory from user and load first image in the folder
    function chDirCB(~,~)
        str=ixon_getDayDir;
        str=uigetdir(str);        
        if ~isequal(str,0) && ~isequal(str,currDir)       
            disp(['Changing directory to ' str]);
            currDir=str;
            chData([],[],0);   
        end
    end

% Button to load an image into the acquisition
ttstr='Load an image into the previer and change the source directory.';
cdata=imresize(imread('icons/file.jpg'),[17 17]);
hbload=uicontrol(hpNav,'style','pushbutton','CData',cdata,...
    'callback',@browseImageCB,'enable','on','backgroundcolor','w',...
    'position',[hbchdir.Position(1)+hbchdir.Position(3) 18 20 20],'ToolTipString',ttstr);

    function browseImageCB(~,~)
       loadImage; 
    end

% Button to delete image
ttstr='Delete this image from the source directory';
cdata=imresize(imread('icons/garbage.jpg'),[17 17]);
hbdelete=uicontrol(hpNav,'style','pushbutton','CData',cdata,...
    'callback',@deleteImageCB,'enable','on','backgroundcolor','w',...
    'position',[hbload.Position(1)+hbload.Position(3) 18 20 20],'ToolTipString',ttstr);

    function deleteImageCB(src,evt)
        filename = tNavName.String;
        str = ['delete ' filename '?'];
         answer = questdlg('Send image to recyle bin?',str,...
            'Yes','No','No') ;        
        switch answer
            case 'Yes'
                disp(['deleting' filename  ' ...'])
                
                if exist(filename,'file')                    
                    previousState = recycle('on');
                    delete(filename);
                    recycle(previousState);
                    disp('deleted');
                    chData([],[],tNavInd.Data);
                else
                    disp('no file to delete?')
                end
            case 'No'
                disp('cancel delete')
            otherwise
                disp('cancel delete')
        end        
     

    end
      

ttstr='Jump to most recent image acquired.';
hbNavNow=uicontrol(hpNav,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',[char(10094) char(10094)],'fontsize',10,...
    'callback',{@chData, 1},'ToolTipString',ttstr,'UserData','absolute');
hbNavNow.Position=[hbdelete.Position(1)+hbdelete.Position(3) 18 24 20];

ttstr='Step to next more recent image';
hbNavLeft=uicontrol(hpNav,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',char(10094),'fontsize',10,...
    'callback',{@chData, -1},'ToolTipString',ttstr,'UserData','increment');
hbNavLeft.Position=[hbNavNow.Position(1)+hbNavNow.Position(3) 18 12 20];

tNavInd=uitable('parent',hpNav,'units','pixels','columnformat',{'numeric'},...
    'fontsize',8,'columnwidth',{25},'rowname',{},'columnname',{},...
    'columneditable',[true],'CellEditCallback',@tblch,'UserData','absolute');
tNavInd.Data=[200];
tNavInd.Position(1:2)=[hbNavLeft.Position(1)+hbNavLeft.Position(3)  18];
tNavInd.Position(3:4)=tNavInd.Extent(3:4);

    function tblch(src,evt)
        n=evt.NewData;        
        if isnumeric(n) && n > 0 && floor(n) == n
            chData(src,evt,n)% n is a natural number
        else
            src.Data = evt.PreviousData;
        end
    end

% Text for total number of images in directory
ttstr='Total number of images in the directory';
tNavMax=uicontrol(hpNav,'style','text','string','of 2001','fontsize',7,...
    'backgroundcolor','w','units','pixels','horizontalalignment','center',...
    'Position',[tNavInd.Position(1)+tNavInd.Position(3) 18 35 14],'tooltipstring',ttstr);

ttstr='Step to later image.';
hbNavRight=uicontrol(hpNav,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',char(10095),'fontsize',10,...
    'callback',{@chData, 1},'ToolTipString',ttstr,'UserData','increment');
hbNavRight.Position=[tNavMax.Position(1)+tNavMax.Position(3) 18 12 20];

ttstr='Jump to first image acquired.';
hbNavLast=uicontrol(hpNav,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',[char(10095) char(10095)],'fontsize',10,...
    'callback',{@chData, inf},'ToolTipString',ttstr,'UserData','absolute');
hbNavLast.Position=[hbNavRight.Position(1)+hbNavRight.Position(3) 18 24 20];


    function loadImage(filename)
        if nargin<1
            [filename,pathname]=uigetfile([defaultDir filesep '*.mat']);
            if ~filename
                disp('No mat file chosen!');
                return;
            end
            filename=[pathname filename];
            currDir=pathname;
        end          
        disp(['     Loading ' filename]);        
        olddata=data;
        try
            newdata=load(filename);            
            data=newdata.data;
            newDataCallback;      
        catch ME                   
            warning(getReport(ME,'extended'));
            errordlg('Unable to load image, reverting to old data');
            beep            
            disp(['FileName : ' filename]);            
            data=olddata;
            newDataCallback;      
        end
    end

% Callback function for changing number of ROIs
    function chData(src,evt,state)     
        if isempty(src)
            index_type='absolute';
        else
            index_type = src.UserData;
        end
       % Get mat files in history directory          
       filenames=dir([currDir  filesep '*.mat']);
       filenames={filenames.name};       
       filenames=sort(filenames);
       filenames=flip(filenames);
       
       if isempty(filenames)
          warning('No data in this folder. Aborting loading file.');
          return;
       end
        
        % Get list of .mat files in the directory, buggy (this is faster)
%         f0=java.io.File([currDir filesep '*.mat']);    
%         if isempty(f0.list)
%             filenames={};
%         else
%             filenames=transpose(sort(cell(f0.list)));
%         end        
    
        % Current data mat  
       myname=[data.Name '.mat'];        

      % Find current filename in directory
       i0=find(strcmp(filenames,myname),1);
       if isempty(i0)
          i0=1; 
       end

       switch index_type
           case 'increment'
               if isequal(state,-1)
                    i1=max([i0-1 1]); 
               elseif isequal(state,1)
                    i1=min([i0+1 length(filenames)]);
               end
           case  'absolute'               
               i1 = max([min([state length(filenames)]) 1]); 
       end

 

        % switch state
        %     case -1
        %         i1=max([i0-1 1]);            
        %     case 1
        %         i1=min([i0+1 length(filenames)]);
        %     case 0
        %         i1=1;    
        %     case nan
        %         i1 = length(filename);
        % end   
        
        newfilename=fullfile(currDir,filenames{i1});
        % tNavInd.String=sprintf('%03d',i1);
        tNavInd.Data(1)=i1;


        % [a,b,ext]=fileparts(newfilename);
        % tNavName.String=fullfile(a,[b ext]); 
        tNavName.String = newfilename;
        tNavMax.String=['of ' num2str(length(filenames))];        

        drawnow;   
        loadImage(newfilename);
    end


% Checkbox for auto updating when new images are taken
ttstr='Automatically refresh to most recent image upon new image acquisition.';
cAutoUpdate=uicontrol('parent',hpNav,'units','pixels','string',...
    'hold image?','value',0,'fontsize',8,'backgroundcolor','w',...
    'Style','checkbox','ToolTipString',ttstr);
cAutoUpdate.Position=[hbNavLast.Position(1)+hbNavLast.Position(3) 20 90 14];

% Text for string of full file name
ttstr='Name of current image';
tNavName=uicontrol(hpNav,'style','text','string','FILENAME','fontsize',7,...
    'backgroundcolor','w','units','pixels','horizontalalignment','left',...
    'Position',[1 1 hpNav.Position(3) 14],'tooltipstring',ttstr);
tNavName.String=data.Name;


%% Image Process Panel

hpProcess=uipanel(hpControl,'units','pixels','backgroundcolor','w',...
    'Position',[0 hpControl.Position(4)-270 160 270],'title','processing');
hpProcess.Position(3) = 160;
hpProcess.Position(4) = 270;
hpProcess.Position(1) = 0;
hpProcess.Position(2) = hpNav.Position(2) - hpProcess.Position(4);

tbl_process_1 = uitable('parent',hpProcess,'columnformat',{'logical','char'},...
    'columnname',{},'rowname',{},'columnwidth',{15,135},'columneditable',[true false],...
    'fontsize',7);
tbl_process_1.Data={true 'subtract bias';
    true 'subtract background';
    false 'apply image mask';
    false 'deconvolve psf Rich-Lucy';
    true 'auto-detect PSF noise?'};
tbl_process_1.Position(3:4)=tbl_process_1.Extent(3:4);
tbl_process_1.Position(1:2)=[1 hpProcess.Position(4)-tbl_process_1.Position(4)-15];



% 
% 
% 
% 
% 
% ttstr='Apply gaussian filter to smooth image';
% cKGaussFilter=uicontrol('style','checkbox','string','fft filter (px)',...
%     'units','pixels','parent',hpADV,'backgroundcolor','w',...
%     'value',1,'ToolTipString',ttstr,'fontsize',7);
% cKGaussFilter.Position=[5 1 80 15];
% 
% % Fitlering the FFT
% tblKGaussFilter=uitable('parent',hpADV,'units','pixels',...
%     'rowname',{},'columnname',{},'Data',1,'columneditable',[true],...
%     'columnwidth',{25},'fontsize',7,'ColumnFormat',{'numeric'});
% tblKGaussFilter.Position=[80 cKGaussFilter.Position(2)-3 30 20];
% 
% % Mask IR Checkbox
% hcIRMask=uicontrol(hpADV,'style','checkbox','string','mask IR (1/px)','fontsize',7,...
%     'backgroundcolor','w','Position',[5 18 100 15],...
%     'ToolTipString',ttstr,'enable','on','value',1);
% 
% % IR Mask Value
% tblIRMask=uitable('parent',hpADV,'units','pixels',...
%     'rowname',{},'columnname',{},'Data',.01,'columneditable',[true],...
%     'columnwidth',{45},'fontsize',7,'ColumnFormat',{'numeric'});
% tblIRMask.Position=[hpADV.Position(3)-55 hcIRMask.Position(2)+1 50 20];
% 
% 
% % Checkbox for computing FFT
% ttstr='compute fft';
% hcFFT=uicontrol(hpADV,'style','checkbox','string','compute FFT','fontsize',7,...
%     'backgroundcolor','w','Position',[5 33 80 15],...
%     'ToolTipString',ttstr,'enable','off','Value',1);
% 
% % Checkbox for enabling 2D gauss fitting
% ttstr='Apply gaussian filter to smooth image';
% cGaussFilter=uicontrol('style','checkbox','string','gauss filter (px)',...
%     'units','pixels','parent',hpADV,'backgroundcolor','w',...
%     'value',0,'ToolTipString',ttstr,'fontsize',7);
% cGaussFilter.Position=[5 hcFFT.Position(2)+16 100 15];
% 
% tblGaussFilter=uitable('parent',hpADV,'units','pixels',...
%     'rowname',{},'columnname',{},'Data',.25,'columneditable',[true],...
%     'columnwidth',{40},'fontsize',8,'ColumnFormat',{'numeric'});
% tblGaussFilter.Position=[hpADV.Position(3)-50 cGaussFilter.Position(2)-2 45 20];
% 
% % Rotate data (makes it easier for lattice grid stuff)
% ttstr='Rotate data';
% cRotate=uicontrol('style','checkbox','string','rotate (deg.)',...
%     'units','pixels','parent',hpADV,'backgroundcolor','w',...
%     'value',0,'ToolTipString',ttstr,'fontsize',7);
% cRotate.Position=[5 cGaussFilter.Position(2)+15 80 15];
% 
% tblTheta=uitable('parent',hpADV,'units','pixels',...
%     'rowname',{},'columnname',{},'Data',59.81,'columneditable',[true],...
%     'columnwidth',{60},'fontsize',8,'ColumnFormat',{'numeric'});
% tblTheta.Position=[hpADV.Position(3)-70 cRotate.Position(2)+3 65 20];
% 
% 
% % Rotate data (makes it easier for lattice grid stuff)
% ttstr='scale up image';
% cScale=uicontrol('style','checkbox','string','scale factor',...
%     'units','pixels','parent',hpADV,'backgroundcolor','w',...
%     'value',0,'ToolTipString',ttstr,'fontsize',7);
% cScale.Position=[5 cRotate.Position(2)+18 80 15];
% 
% tblScale=uitable('parent',hpADV,'units','pixels',...
%     'rowname',{},'columnname',{},'Data',2,'columneditable',[true],...
%     'columnwidth',{40},'fontsize',8,'ColumnFormat',{'numeric'});
% tblScale.Position=[hpADV.Position(3)-70 cScale.Position(2) 50 20];
% 
% 
% % Checkbox for applying point spread function 
% ttstr='Deconvolve data with point spread function using the Richardson-Lucy algorithm.';
% hcPSF=uicontrol(hpADV,'style','checkbox','string','denconvolve psf Rich-Lucy','fontsize',7,...
%     'backgroundcolor','w','Position',[5 tblScale.Position(2)+70 155 20],...
%     'ToolTipString',ttstr,'enable','on');
% 
% ttstr='Calculate noise statistics from image?';
% hcPSF_noise=uicontrol(hpADV,'style','checkbox','string','auto detect noise?','fontsize',6,...
%     'backgroundcolor','w','Position',[15 hcPSF.Position(2)-10 155 14],'Value',1,...
%     'ToolTipString',ttstr,'enable','on');

tblPSF=uitable('parent',hpProcess,'units','pixels',...
    'columnname',{char(963),'size','iter','noise'},'rowname',{},'Data',[1.32 11 31 25],'columneditable',[true],...
    'columnwidth',{35 25 25 35},'fontsize',7,'ColumnFormat',{'numeric'},'CellEditCallback',{@(src,evt) updatePSFKGraphic});
tblPSF.Position(3:4) = tblPSF.Extent(3:4);
% tblPSF.Position(1:2)=[20 hcPSF_noise.Position(2)-tblPSF.Extent(4)];  
tblPSF.Position(1:2)=[10 tbl_process_1.Position(2)-tblPSF.Extent(4)-2];  

tbl_process_2 = uitable('parent',hpProcess,'columnformat',{'logical','char'},...
    'columnname',{},'rowname',{},'columnwidth',{15,95, 40},'columneditable',[true false true],...
    'fontsize',7);
tbl_process_2.Data={
    false 'gauss filter (px)' 1;
     false 'scale (factor)' 2;
    false 'rotate (deg)' 60;
    false 'fft filter (px)' 1;
    false 'fft ir mask (1/px)' 1;
   };
tbl_process_2.Position(3:4)=tbl_process_2.Extent(3:4);
tbl_process_2.Position(1:2)=[1 tblPSF.Position(2)-tbl_process_2.Position(4)-5];


% % Subtract background image
% ttstr='Subtract off background image from atoms images.';
% hcSubBG=uicontrol(hpADV,'style','checkbox','string','subtract bgd','fontsize',7,...
%     'backgroundcolor','w','Position',[5 hcPSF.Position(2)+15 80 15],...
%     'ToolTipString',ttstr,'enable','on','Value',1);
% 
% ttstr='Apply mask to image to eliminate aperture clipping';
% hcMask=uicontrol(hpADV,'style','checkbox','string','apply image mask','fontsize',7,...
%     'backgroundcolor','w','Position',[5 hcSubBG.Position(2)+15 120 13],...
%     'ToolTipString',ttstr,'enable','on','Value',0);
% 
% % Subtract bias
% ttstr='Subtract off electronic/software bias of 200 counts from raw images.';
% hcSubBias=uicontrol(hpADV,'style','checkbox','string','subtract bias','fontsize',7,...
%     'backgroundcolor','w','Position',[5 hcMask.Position(2)+13 80 15],...
%     'ToolTipString',ttstr,'enable','on','Value',1);



% process button
hbprocess=uicontrol(hpProcess,'style','pushbutton','string','process',...
    'units','pixels','callback',@processCB,'parent',hpProcess,'backgroundcolor',[80 200 120]/255);
hbprocess.Position=[hpProcess.Position(3)-45 1 45 15];
% hbprocess.Position=[hpADV.Position(3)-45 hpADV.Position(4)-30 45 15];
hbprocess.Position=[3 1 hpProcess.Position(3)-8 18];
hbprocess.String = 'process images';
    function processCB(src,evt)
        newDataCallback;
    end

%% Feedback Panel

hpFeedback=uipanel(hpControl,'units','pixels','backgroundcolor','w',...
    'Position',[0 hpProcess.Position(2)-50 160 50],'title','feedback');


% ttstr='Open the feedback file';
% cdata=imresize(imread('icons/file_matlab.png'),[14 14]);
open_func=uicontrol(hpFeedback,'style','pushbutton','string','open feedback function',...
    'callback',@(src,evt) open('ixon_gui_feedback'),'backgroundcolor','w',...
    'fontsize',7,...
    'position',[3 hpFeedback.Position(4)-34 hpFeedback.Position(3)-8  16],'ToolTipString',ttstr);

% Feedback function
hbfeedback=uicontrol(hpFeedback,'style','pushbutton','string','evaluate feedback function',...
    'units','pixels','callback',{@(~,~) ixon_gui_feedback},'parent',hpFeedback,...
    'backgroundcolor','w','position',[3 1 hpFeedback.Position(3)-8 16],'backgroundcolor',[80 200 120]/255,...
    'fontsize',7);

%% Analysis Panel
% hpAnl=uipanel(hF,'units','pixels','backgroundcolor','w','title','position analysis');
% hpAnl.Position=[0 hpADV.Position(2)-130 160 180];

hpPosition=uipanel(hpControl,'units','pixels','backgroundcolor','w',...
    'Position',[0 hpFeedback.Position(2)-170 160 170],'title','position analysis');


% Checkbox for center of mass and sigma 
ttstr='Automatically perform analysis on new image';
hc_anlX_auto=uicontrol(hpPosition,'style','checkbox','string','auto-analyze on new image?','fontsize',7,...
    'backgroundcolor','w','Position',[1 hpPosition.Position(4)-35 hpPosition.Position(3)-1 15],...
    'ToolTipString',ttstr,'enable','off','Value',1);

% Table of ROIs
tblROI=uitable(hpPosition,'units','pixels','ColumnWidth',{30 30 30 30},...
    'ColumnEditable',true(ones(1,4)),'ColumnName',{'X1','X2','Y1','Y2'},...
    'Data',[1 512 1 512],'FontSize',8,...
    'CellEditCallback',@chROI,'RowName',{});
tblROI.Position(3:4)=tblROI.Extent(3:4)+0*[18 0];
tblROI.Position(1:2)=[5 hc_anlX_auto.Position(2)-tblROI.Position(4)];

% Callback function for changing ROI via table
    function chROI(src,evt)
        m=evt.Indices(1); n=evt.Indices(2);
        
        ROI=src.Data(m,:);
        % Check that the data is numeric
        if sum(~isnumeric(ROI)) || sum(isinf(ROI)) || sum(isnan(ROI))
            warning('Incorrect data type provided for ROI.');
            src.Data(m,n)=evt.PreviousData;
            return;
        end
        
        ROI=round(ROI);      % Make sure this ROI are integers   
        % Check that limits go from low to high
        if ROI(2)<=ROI(1) || ROI(4)<=ROI(3)
           warning('Bad ROI specification given.');
           ROI(evt.Indices(2))=evt.PreviousData;
        end               
        % Check that ROI is within image bounds
        if ROI(1)<1; ROI(1)=1; end       
        if ROI(3)<1; ROI(3)=1; end   
        if ROI(4)>512; ROI(4)=512; end       
        if ROI(2)>512; ROI(2)=512; end         
        % Reassign the ROI
        src.Data(m,:)=ROI;      
        % Try to update ROI graphics
        try
            pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];
            set(pROI(m),'Position',pos);
        catch
           warning('Unable to change display ROI.');
           src.Data(m,n)=evt.PreviousData;
        end
    end

% Button to enable GUI selection of analysis ROI
ttstr='Select the analysis ROI.';
uicontrol(hpPosition,'style','pushbutton','Cdata',img_Select,'Fontsize',10,...
    'Backgroundcolor','w','Position',[130 tblROI.Position(2)+25 18 18],...
    'Callback',@slctROICB,'ToolTipString',ttstr);

% % Button to snap display ROI to the data ROI
ttstr='Snap analysis ROI to display ROI.';
uicontrol(hpPosition,'style','pushbutton','Cdata',img_Snap,'Fontsize',10,...
    'Backgroundcolor','w','Position',[130 tblROI.Position(2)+7 18 18],...
    'Callback',@snapPosAnalysisROI,'ToolTipString',ttstr);

    function snapPosAnalysisROI(src,evt)
        ROI = tbl_dROI_X.Data;
        set(tblROI,'Data',ROI);
        pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];
        set(pROI,'Position',pos);
    end

% Button for maximizing the display limits
ttstr='Max analysis ROI to full image size.';
uicontrol(hpPosition,'style','pushbutton','Cdata',img_Full,'Fontsize',10,...
    'Backgroundcolor','w','Callback',@maxPosAnalysisROI,...
    'ToolTipString',ttstr,'Position',[130 tblROI.Position(2)-11 18 18]);

    function maxPosAnalysisROI(~,evt)
        ROI = [1 512 1 512];
        tblROI.Data = ROI;
        pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];
        set(pROI,'Position',pos);
    end


% Callback for selecting an ROI based upon mouse click input.
    function slctROICB(~,~)
        disp(['Selecting display ROI .' ...
            ' Click two points that form the rectangle ROI.']);
        set(hp,'SelectedTab',tabX);
        [x1,y1]=ginputMe(1);          % Get a mouse click
        x1=round(x1);y1=round(y1);  % Round to interger       
        p1=plot(x1,y1,'+','color','r','linewidth',1,'markersize',10,'hittest','off'); % Plot it        
        [x2,y2]=ginputMe(1);          % Get a mouse click
        enableInteractivity;

        x2=round(x2);y2=round(y2);  % Round it        
        p2=plot(x2,y2,'+','color','r','linewidth',1,'markersize',10,'hittest','off');  % Plot it

        
        % Create the ROI
        ROI=[min([x1 x2]) max([x1 x2]) min([y1 y2]) max([y1 y2])];
        % Constrain ROI to image
        if ROI(1)<1; ROI(1)=1; end       
        if ROI(3)<1; ROI(3)=1; end   
        if ROI(4)>512; ROI(4)=512; end       
        if ROI(2)>512; ROI(2)=512; end           
        % Try to update ROI graphics
        tblROI.Data=ROI;           
        try
            pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];
            set(pROI(1),'Position',pos);
            drawnow;        
        catch
           warning('Unable to change ROI.');
        end
        delete(p1);delete(p2);                   % Delete markers
    end


% Checkbox for center of mass and sigma 
ttstr='Find 1st and 2nd moments in the ROI box.';
hc_anlX_Box=uicontrol(hpPosition,'style','checkbox','string','box','fontsize',7,...
    'backgroundcolor','w','Position',[5 tblROI.Position(2)-15 120 15],...
    'ToolTipString',ttstr,'enable','off','Value',1);

% Checkbox for image histogram
ttstr='Calculate histogram';
hc_anlX_Histogram=uicontrol(hpPosition,'style','checkbox','string','histogram','fontsize',7,...
    'backgroundcolor','w','Position',[5 hc_anlX_Box.Position(2)-15 120 15],...
    'ToolTipString',ttstr,'enable','off','Value',1);

% Checkbox for principal component analysis
ttstr='Principal component analysis to determine cloud axes..';
hc_anlX_PCA=uicontrol(hpPosition,'style','checkbox','string','find principal axes','fontsize',7,...
    'backgroundcolor','w','Position',[5 hc_anlX_Histogram.Position(2)-15 120 15],...
    'ToolTipString',ttstr,'enable','off','callback',@hcpcaCB);

    function hcpcaCB(src,~)
       for n=1:length(pPCA)
          if ~src.Value
             pPCA(n).Visible='off'; 
          end 
       end  
       
       if src.Value
          hc_anlX_GaussRot.Enable='on';
       else
            hc_anlX_GaussRot.Enable='off';
            hc_anlX_GaussRot.Value=0;
       end
    end


hc_anlX_Gauss=uicontrol(hpPosition,'style','checkbox','string','gaussian','fontsize',7,...
    'backgroundcolor','w','Position',[5 hc_anlX_PCA.Position(2)-15 60 15],...
    'ToolTipString',ttstr,'enable','on','callback',@hcgaussCB);

    function hcgaussCB(src,~)
       for n=1:length(pXF)
          pXF(n).Visible=src.Value;
          pYF(n).Visible=src.Value;
          
          if ~src.Value
             pGaussRet(n).Visible='off';
             cGaussRet.Value=0;
             cGaussRet.Enable='off';
          end
       end       
       if src.Value;hc_anlX_GaussRot.Value=0;end
    end

hc_anlX_GaussRot=uicontrol(hpPosition,'style','checkbox','string','rotatable?','fontsize',7,...
    'backgroundcolor','w','Position',[65 hc_anlX_Gauss.Position(2) 100 15],'callback',@(~,~) disp('hi'),...
    'ToolTipString',ttstr,'enable','off');

ttstr='Analyze stripe pattern in image to measure field stability';
hcStripe=uicontrol(hpPosition,'style','checkbox','string','stripe pattern','fontsize',7,...
    'backgroundcolor','w','Position',[5 hc_anlX_Gauss.Position(2)-15 100 15],...
    'ToolTipString',ttstr);

% Do Position space analysis
hbposition=uicontrol(hpPosition,'style','pushbutton','string','position analysis',...
    'units','pixels','callback',{@(~,~) updatePositionAnalysis},'parent',hpPosition,...
    'backgroundcolor','w','position',[3 1 hpPosition.Position(3)-8 18],'backgroundcolor',[80 200 120]/255);



%% Momentum Panel
% hpKspace=uipanel(hF,'units','pixels','backgroundcolor','w','title','momentum analysis');
% hpKspace.Position=[0 hpAnl.Position(2)-90 160 105];

hpKspace=uipanel(hpControl,'units','pixels','backgroundcolor','w',...
    'Position',[0 hpPosition.Position(2)-110 160 130],'title','momentum analysis');
hpKspace.Position(1)=hpProcess.Position(1)+hpProcess.Position(3);
hpKspace.Position(2)=hpNav.Position(2)-hpKspace.Position(4);

% Checkbox for center of mass and sigma 
ttstr='Automatically perform analysis on new image';
hc_anlK_auto=uicontrol(hpKspace,'style','checkbox','string','auto-analyze on new image?','fontsize',7,...
    'backgroundcolor','w','Position',[1 hpKspace.Position(4)-35 hpKspace.Position(3)-1 15],...
    'ToolTipString',ttstr,'enable','on','Value',0);

% Table of ROIs
tblROIK=uitable(hpKspace,'units','pixels','ColumnWidth',{30 30 30 30},...
    'ColumnEditable',true(ones(1,4)),'ColumnName',{'kx1','kx2','ky1','ky2'},...
    'Data',[-.5 .5 -.5 .5],'FontSize',6,...
    'CellEditCallback',@chROIK,'RowName',{});
tblROIK.Position(3:4)=tblROIK.Extent(3:4)+0*[18 0];
tblROIK.Position(1:2)=[5 hc_anlK_auto.Position(2)-tblROIK.Position(4)];

% Callback function for changing ROI via table
    function chROIK(src,evt)
        m=evt.Indices(1); n=evt.Indices(2);
        
        ROI=src.Data(m,:);
        % Check that the data is numeric
        if sum(~isnumeric(ROI)) || sum(isinf(ROI)) || sum(isnan(ROI))
            warning('Incorrect data type provided for ROI.');
            src.Data(m,n)=evt.PreviousData;
            return;
        end
        
        ROI=round(ROI);      % Make sure this ROI are integers   
        % Check that limits go from low to high
        if ROI(2)<=ROI(1) || ROI(4)<=ROI(3)
           warning('Bad ROI specification given.');
           ROI(evt.Indices(2))=evt.PreviousData;
        end               
        % Check that ROI is within image bounds
        if ROI(1)<-.5; ROI(1)=-.5; end       
        if ROI(3)<-.5; ROI(3)=-.5; end   
        if ROI(4)>.5; ROI(4)=.5; end       
        if ROI(2)>.5; ROI(2)=.5; end         
        % Reassign the ROI
        src.Data(m,:)=ROI;      
        % Try to update ROI graphics
        try
            pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];
            set(pROI_K(m),'Position',pos);
        catch
           warning('Unable to change display ROI.');
           src.Data(m,n)=evt.PreviousData;
        end
    end


% Mask IR Checkbox
hcFindLattice=uicontrol(hpKspace,'style','checkbox','string','lattice basis and phase','fontsize',7,...
    'backgroundcolor','w','Position',[5 20 120 15],...
    'ToolTipString',ttstr,'enable','on','value',1);
hcFindLattice.Position(2) = tblROIK.Position(2) - 15;

% Mask IR Checkbox
hcKFocus=uicontrol(hpKspace,'style','checkbox','string','multi-shot focusing','fontsize',7,...
    'backgroundcolor','w','Position',[5 5 120 15],...
    'ToolTipString',ttstr,'enable','on','value',1);
hcKFocus.Position(2) = hcFindLattice.Position(2) - 15;

% Refit button
hb_Kanalyze=uicontrol(hpKspace,'style','pushbutton','string','momentum analysis',...
    'units','pixels','callback',@analyze_k,'parent',hpKspace,'backgroundcolor',[80 200 120]/255);
hb_Kanalyze.Position=[3 1 hpKspace.Position(3)-8 18];



% Callback function for redoing fits button
    function analyze_k(~,~)
        if hcFindLattice.Value && isfield(data,'Zf')
            try
            hb_Kanalyze.BackgroundColor=[255 219 88]/255;
            drawnow;
            opts = data.ProcessOptions;  

            if isfield(data,'LatticeK')
                data = rmfield(data,'LatticeK');
            end
            for kk=1:size(data.Zf,3)
                tic;
                fprintf(['(' num2str(kk) '/' num2str(size(data.Zf,3)) ') Fitting reciprocal lattice ...']);
                data.LatticeK(kk) = findLatticeK(data.f,data.f,data.Zf(:,:,kk),opts);                
                k1 = data.LatticeK(kk).k1;
                k2 = data.LatticeK(kk).k2;          
                t2 = toc;
                fprintf([' done (' num2str(t2,2) ' sec.)' ' phase ...']);
                tic;
                data.LatticePhase(kk) = findLatticePhase(data.X,data.Y,data.Z,k1,k2);              
                t=toc;
                disp([' done (' num2str(t,2) ' sec.)']);                
            end                        
            updateGridGraphics;
            latticeGridCB(cDrawLattice);
            
            updateReciprocalReticleGraphics;
            reciprocalLatticeReticleCB(cLat_KReticle);
            
            updateReciprocalTextGraphics;
            reciprocalLatticeTextCB(cLat_K_text);
            
            hb_Kanalyze.BackgroundColor=[80 200 120]/255;
            drawnow;

            if cAutoColor_K.Value;setClim('K');end  
            end
        end
                
        if hcKFocus.Value &&  ...
            isfield(data,'Zf') ...
                && isfield(data.Flags,'lattice_fluor_multi_mode') ...
                &&  (data.Flags.lattice_fluor_multi_mode==2)
            

            try
                opts=struct;
                opts.doDebug=1;
                opts.Parent = tabFocus;
                opts.ROI = tbl_dROI_focus.Data;
                if isfield(data,'KFocusing')
                    data=rmfield(data,'KFocusing');
                end             

                data=ixon_fft_multi_shot_focusing(data,opts); 
            catch ME
                warning('OH NO');                
            end
        end
    end

%% Binning Panel
hpBin=uipanel(hpControl,'units','pixels','backgroundcolor','w',...
    'title','binned analysis');
hpBin.Position(3)=160;
hpBin.Position(4)=200;
hpBin.Position(1)=hpProcess.Position(1)+hpProcess.Position(3);
hpBin.Position(2)=hpKspace.Position(2)-hpBin.Position(4);

ttstr='auto do bin analysis';
hc_anlB_auto=uicontrol(hpBin,'style','checkbox','string','auto-analyze on new image?','fontsize',7,...
    'backgroundcolor','w','Position',[5 hpBin.Position(4)-35 hpBin.Position(3)-1 15],...
    'ToolTipString',ttstr,'enable','on','Value',0);

hb_BinOptions=uicontrol('style','pushbutton','string',...
    'open(ixon_gui_bin_options.m)','units','pixels','parent',hpBin,...
    'fontsize',7,'Callback',@(src,evt) open('ixon_gui_bin_options.m'));
hb_BinOptions.Position(3:4)=[150 18];
hb_BinOptions.Position(1:2)=hc_anlB_auto.Position(1:2) + ...
    [0 -hb_BinOptions.Position(4)-2];   

% Checkbox for image sharpness
ttstr='Calculate histgoram';
hc_anlB_Histogram=uicontrol(hpBin,'style','checkbox','string','histogram','fontsize',7,...
    'backgroundcolor','w','Position',[1 hb_BinOptions.Position(2)-15 100 15],...
    'ToolTipString',ttstr,'enable','off','Value',1);

ttstr='radial binned histgoram analysis';
hc_anlB_Radial=uicontrol(hpBin,'style','checkbox','string','radial analysis','fontsize',7,...
    'backgroundcolor','w','Position',[1 hc_anlB_Histogram.Position(2)-15 hpBin.Position(3)-1 15],...
    'ToolTipString',ttstr,'enable','on','Value',1);


ttstr='bin stripe focus';
hc_anlB_stripe=uicontrol(hpBin,'style','checkbox','string','stripe','fontsize',7,...
    'backgroundcolor','w','Position',[1 hc_anlB_Radial.Position(2)-15 hpBin.Position(3)-1 15],...
    'ToolTipString',ttstr,'enable','on','Value',0);

ttstr='  focus';
hc_anlB_multishotfocus=uicontrol(hpBin,'style','checkbox','string','multi-shot focusing','fontsize',7,...
    'backgroundcolor','w','Position',[1 hc_anlB_stripe.Position(2)-15 hpBin.Position(3)-1 15],...
    'ToolTipString',ttstr,'enable','on','Value',0);

ttstr='circle stripe phase';
hc_anlB_circlestripe=uicontrol(hpBin,'style','checkbox','string','circle stripe','fontsize',7,...
    'backgroundcolor','w','Position',[1 hc_anlB_multishotfocus.Position(2)-15 hpBin.Position(3)-1 15],...
    'ToolTipString',ttstr,'enable','on','Value',0);


% Refit button
hb_Binanalyze=uicontrol(hpBin,'style','pushbutton','string','binned analysis',...
    'units','pixels','callback',@analyze_bin,'parent',hpBin,'backgroundcolor',[80 200 120]/255);
hb_Binanalyze.Position=[3 1 hpBin.Position(3)-8 18];

    function [a1, a2, p1, p2] = getLattice        
        switch bgBasis.SelectedObject.UserData
            case 'fft'
                if isfield(data,'LatticePhase')
                    a1 = data.LatticePhase(kk).a1;
                    a2 = data.LatticePhase(kk).a2;                        
                    p1 = data.LatticePhase(kk).p1;
                    p2 = data.LatticePhase(kk).p2;    
                    
                    a1 = data.LatticePhase(1).a1;
                    a2 = data.LatticePhase(1).a2;                        
                    p1 = data.LatticePhase(1).p1;
                    p2 = data.LatticePhase(1).p2;   
                else
                    errordlg('Lattice phase not calculated yet.');
                    beep
                    return;                        
                end
            case 'manual'
                a1 = tblBasis.Data(1,(1:2))';
                a2 = tblBasis.Data(2,(1:2))';
                p1 = tblBasis.Data(1,3);
                p2 = tblBasis.data(2,3);
        end
    end

function analyze_bin(src,evt)
    if isfield(data,'LatticeBin')
        data = rmfield(data,'LatticeBin');
    end
        
     try
        hb_Binanalyze.BackgroundColor=[255 219 88]/255;
        drawnow;       
        

        [data] = ixon_ProcessImagesToBin(data,ixon_gui_bin_options());
        
        updateBinnedGraphics;     
        updateBinnedHistogramGraphics;    

         if hc_anlB_Radial.Value
            ropts=struct;
            ropts.FigureNumber = 3001;
            data=bin_radialHistogram(data);
            bin_showRadialHistogram(data,ropts);
         end                
         
        if (hc_anlB_stripe.Value && ~isfield(data.Flags,'plane_selection_dotilt')) || ...
            (hc_anlB_stripe.Value && isfield(data.Flags,'plane_selection_dotilt') &&  data.Flags.plane_selection_dotilt)
            opts_stripe = struct;
            opts_stripe.FigureNumber=3000;
            opts_stripe.Threshold = stripe_threshold_tbl.Data;
            ll=1;
            n1 = data.LatticeBin(ll).n1;
            n2 = data.LatticeBin(ll).n2;
            Zb = data.LatticeBin(ll).Zbin;    
            opts_stripe.LGuess = 26.62;
            opts_stripe.SumIndex = 1; % 1 for trees % 2 for fallen trees                
            out = bin_StripeFit(n1,n2,Zb,opts_stripe);                  
            data.BinStripe = out;     
            bin_showStripeBin(data,[],opts_stripe);
        end  
        
        if (hc_anlB_multishotfocus.Value && ~isfield(data.Flags,'lattice_fluor_multi_mode')) || ...
            (hc_anlB_multishotfocus.Value && isfield(data.Flags,'lattice_fluor_multi_mode') &&  data.Flags.lattice_fluor_multi_mode==2)
            opts_focusing = struct;
            opts_focusing.FigureNumber=4000;
            opts_focusing.Sites = [25 25];
            opts_focusing.XVal = data.Params.qgm_MultiPiezos(2:5);
            data = bin_recenter(data);
            n1 = data.LatticeBin(1).n1;
            n2 = data.LatticeBin(1).n2;             

            Zb = zeros(length(n2),length(n1),length(data.LatticeBin));
            for jj=1:length(data.LatticeBin)
                Zb(:,:,jj) = data.LatticeBin(jj).Zbin;
            end
            Focusing = bin_multiShotFocusing(n1,n2,Zb,opts_focusing);      
            data.MultiShotFocusing = Focusing;
        end 


        if (hc_anlB_circlestripe.Value && ~isfield(data.Flags,'plane_selection_dotilt')) || ...
            (hc_anlB_circlestripe.Value && isfield(data.Flags,'plane_selection_dotilt') &&  data.Flags.plane_selection_dotilt)
            opts_stripe = struct;
            opts_stripe.doDebug=1; 

            ll=1;
            n1 = data.LatticeBin(ll).n1;
            n2 = data.LatticeBin(ll).n2;
            Zb = zeros(length(n2),length(n1));
            for ll=1:length(data.LatticeBin)
                myZ = data.LatticeBin(ll).Zbin;
                myZ(isnan(myZ))=0;
                myZ(isinf(myZ))=0;
                Zb = Zb + myZ;                 
            end
            data.BinStripeCircular = ...
                    bin_StripeFitCircular(n1,n2,Zb,opts_stripe);
        end  

            hb_Binanalyze.BackgroundColor=[80 200 120]/255;
            drawnow;   

     catch ME
            warning(getReport(ME));
         end
    end    
        
    
    
%% Digitization Panel
% hpDig=uipanel(hF,'units','pixels','backgroundcolor','w','title','digitization');
% hpDig.Position=[0 hpBin.Position(2)-80 160 100];


hpDig=uipanel(hpControl,'units','pixels','backgroundcolor','w',...
    'title','digital analysis');
hpDig.Position(3:4) = [160 140];
hpDig.Position(1) = hpBin.Position(1);
hpDig.Position(2) = hpBin.Position(2)-hpDig.Position(4);

ttstr='auto do digital analysis';
hc_anlD_auto=uicontrol(hpDig,'style','checkbox','string','auto-analyze on new image?','fontsize',7,...
    'backgroundcolor','w','Position',[1 hpDig.Position(4)-35 hpDig.Position(3)-1 15],...
    'ToolTipString',ttstr,'enable','on','Value',0);

hb_DigOptions=uicontrol('style','pushbutton','string',...
    'open(ixon_gui_dig_options.m)','units','pixels','parent',hpDig,...
    'fontsize',7,'Callback',@(src,evt) open('ixon_gui_dig_options.m'));
hb_DigOptions.Position(3:4)=[150 18];
hb_DigOptions.Position(1:2)=hc_anlD_auto.Position(1:2) + ...
    [0 -hb_DigOptions.Position(4)-2];   

% Thresholding Fidelity
hcDigStandard=uicontrol(hpDig,'style','checkbox','string','standard','fontsize',7,...
    'backgroundcolor','w','Position',[1 hb_DigOptions.Position(2)-20 100 15],'Value',1,...
    'enable','on');

% Thresholding Fidelity
hcDigRadial=uicontrol(hpDig,'style','checkbox','string','radial','fontsize',7,...
    'backgroundcolor','w','Position',[1 hcDigStandard.Position(2)-16 100 15],'Value',0,...
    'enable','on');

% Thresholding Fidelity
hcDigCorrelator=uicontrol(hpDig,'style','checkbox','string','correlators','fontsize',7,...
    'backgroundcolor','w','Position',[1 hcDigRadial.Position(2)-16 100 15],'Value',1,...
    'enable','on');

% Thresholding Fidelity
hcDigFidelity=uicontrol(hpDig,'style','checkbox','string','two-shot fidelity','fontsize',7,...
    'backgroundcolor','w','Position',[1 hcDigCorrelator.Position(2)-16 100 15],'Value',0);

% Refit button
hb_Diganalyze=uicontrol(hpDig,'style','pushbutton','string','digital analysis',...
    'units','pixels','callback',@analyze_dig,'parent',hpDig,'backgroundcolor',[80 200 120]/255);
hb_Diganalyze.Position=[3 1 hpDig.Position(3)-8 18];

    function chBasis(src,evt)
       defaultBasis = src.Data;
    end

    function analyze_dig(src,evt)
        if ~isfield(data,'LatticeBin')
            warning('No binnined data to digitize');
            return;
        end
        
        hb_Diganalyze.BackgroundColor=[255 219 88]/255;
        drawnow;
        
        dig_threshold = tblDig.Data;
        data = ixon_digitize(data,dig_threshold);

         if hcDigFidelity.Value             
            opts=struct;
            opts.threshold = dig_threshold;            
            opts.Label = data.Name;
            n1 = data.LatticeDig(1).n1;
            n2 = data.LatticeDig(1).n2;
            Zdig = zeros(length(n2),length(n1),length(data.LatticeDig));            
            for jj = 1 :length(data.LatticeDig)
                Zdig(:,:,jj) = data.LatticeDig(jj).Zdig;
            end
            opts.FigureNumber = 4001; 
            bin_Fidelity(data,opts);            
            opts.FigureNumber = 4002; 
            out = dig_Fidelity(Zdig,n1,n2,opts);  
            data.DigFideltiy = out; 
        end          

        updateDigitalGraphics;        
        hb_Diganalyze.BackgroundColor=[80 200 120]/255;
        drawnow;
    end

    function updateDigitalGraphics    
         if ~isfield(data,'LatticeDig')
            return;
         end  
        
        imgnum = menuSelectImg.Value; % get image to analyze
        dig = data.LatticeDig(imgnum);  
        set(hImg_D,'XData',dig.n1,'YData',dig.n2,'CData',dig.Zdig);  
        % Update box count string
        str=[ num2str(dig.Natoms) ' atoms' newline ...
            '$(X_\mathrm{c},Y_\mathrm{c}) = ' '('  num2str(round(dig.Xc_site,1)) ',' ...
            num2str(round(dig.Yc_site,1)) ')$' newline ...
            '$(\sigma_X,\sigma_Y) = ' '('  num2str(round(dig.Xs_site,1)) ',' ...
            num2str(round(dig.Ys_site,1)) ')$' newline ...
            'thresh = ' num2str(dig.Threshold)]; 
        %Update box count string object
        set(tD_SE,'String',str); 
        set(tD_SW,'String',[data.Name ' (' num2str(menuSelectImg.Value) ')']);  
        
        drawAtomsCB;
        set(pAtoms,'XData',dig.Rn(1,:),'YData',dig.Rn(2,:));
        
        pAtoms.DataTipTemplate.DataTipRows(3:end)=[];

        row = dataTipTextRow('Counts',dig.Counts);
        pAtoms.DataTipTemplate.DataTipRows(3) = row;
        
        row = dataTipTextRow('n1',dig.N(1,:));
        pAtoms.DataTipTemplate.DataTipRows(4) = row;
        
         row = dataTipTextRow('n2',dig.N(2,:));
        pAtoms.DataTipTemplate.DataTipRows(5) = row;
        
        drawnow;
        
    end

%% Image Number Selector
% hpDisp_Select = uipanel(hF,'units','pixels','backgroundcolor','w','title','image selector');
% hpDisp_Select.Position=[160 500 160 80];

hpDisp_Select = uipanel(hpControl,'units','pixels','backgroundcolor','w',...
    'title','image selector');
hpDisp_Select.Position=[320 hpNav.Position(2)-50 280 50];

menuSelectCMAP=uicontrol('style','popupmenu','string',...
    {'black-purple','black-purple-white','white-purple'},'units','pixels','parent',hpDisp_Select,...
    'Callback',@updateCMAP,'fontsize',12,'Value',3);
menuSelectCMAP.Position(3:4)=[125 18];
menuSelectCMAP.Position(1:2)=[2 15];   

menuSelectImg=uicontrol('style','popupmenu','string',...
    {'image 1 of 1','image 1 of 2'},'units','pixels','parent',hpDisp_Select,...
    'Callback',{@(a,b) updateGraphics},'fontsize',12);
menuSelectImg.Position(3:4)=[125 18];
% menuSelectImg.Position(1:2)=[2 15];   
menuSelectImg.Position(1)=menuSelectCMAP.Position(1)+menuSelectCMAP.Position(3)+2;
menuSelectImg.Position(2)=15;   

    function updateCMAP(src,evt)
        switch src.Value
            case 1
                ca = [0 0 0];       
                cb = [0.7 .1 .6];
                cc = [linspace(ca(1),cb(1),1000)' ...
                    linspace(ca(2),cb(2),1000)' linspace(ca(3),cb(3),1000)'];
                colormap(hF,cc);
                pGrid.Color = [.5 .5 .5 .5];
                pAtoms.Color=[255,215,0]/255;
                pAtoms.MarkerFaceColor=pAtoms.Color;
                pAtoms.MarkerSize=3;

                % stripe_pBar.Color='w';
                % stripe_pAngleCirc.Color='w';  
                % stripe_pCloudEllipse.Color='w';  

            case 2
                colormap(hF,purplemap);
                pGrid.Color = [.5 .5 .5 .5];
                % stripe_pBar.Color='w';
                % stripe_pAngleCirc.Color='w';  
                % stripe_pCloudEllipse.Color='w';
                pAtoms.Color=[255,215,0]/255;
     pAtoms.MarkerFaceColor=pAtoms.Color;
                pAtoms.MarkerSize=3;
            case 3
                pGrid.Color = [.3 .3 .3 .3];
                colormap(hF,purplemap);
                ca = [1 1 1];
                cb = [0.6 0 .5];
                cc = [linspace(ca(1),cb(1),1000)' ...
                    linspace(ca(2),cb(2),1000)' linspace(ca(3),cb(3),1000)'];
                colormap(hF,cc);
                % stripe_pBar.Color='k';
                % stripe_pAngleCirc.Color='k';  
                % stripe_pCloudEllipse.Color='k';
                pAtoms.Color=[218,165,32]/255;
%                         pAtoms.Color=[0,128,128]/255;

                pAtoms.MarkerFaceColor=pAtoms.Color;
                pAtoms.Color='k';
                pAtoms.MarkerFaceColor='y';

                pAtoms.MarkerSize=4;
        end
    end

    function updateImageDataLists
        menuSelectImg.String = {};
        for nn = 1 :size(data.Z,3)
            menuSelectImg.String{nn} = ['image ' num2str(nn) ' of ' num2str(size(data.Z,3))];
        end        
        menuSelectImg.Value = min([menuSelectImg.Value length(menuSelectImg.String)]);
    end

%% Display Options Panel

% hpDispOpt=uitabgroup(hF,'units','pixels');
% hpDispOpt=uipanel(hF,'units','pixels','title','display');

hpDispOpt=uitabgroup(hpControl,'units','pixels');
hpDispOpt.Position(1) = hpDisp_Select.Position(1);
hpDispOpt.Position(3) = hpDisp_Select.Position(3);
hpDispOpt.Position(4) = 250;
hpDispOpt.Position(2) = hpDisp_Select.Position(2)-hpDispOpt.Position(4);

disp_opt_tabs(1)=uitab(hpDispOpt,'Title','position','units','pixels');
disp_opt_tabs(2)=uitab(hpDispOpt,'Title','stripe','units','pixels');
disp_opt_tabs(3)=uitab(hpDispOpt,'Title','focus','units','pixels');

disp_opt_tabs(4)=uitab(hpDispOpt,'Title','histogram','units','pixels');
disp_opt_tabs(5)=uitab(hpDispOpt,'Title','momentum','units','pixels');
disp_opt_tabs(6)=uitab(hpDispOpt,'Title','binned','units','pixels');
disp_opt_tabs(7)=uitab(hpDispOpt,'Title','binned histogram','units','pixels');
disp_opt_tabs(8)=uitab(hpDispOpt,'Title','digitized','units','pixels');


%% Display Options Panel

hpDisp_X = uipanel(hF,'units','pixels','backgroundcolor','w');
hpDisp_X.Position=[160 500 160 190];

set(hpDisp_X,'parent',hpDispOpt.Children(1))

hpDisp_X.Position=[1 1 hpDispOpt.Position(3) hpDispOpt.Position(4)];


menuSelectImgType=uicontrol('style','popupmenu','string',...
    {'processed','no filter'},'units','pixels','parent',hpDisp_X,...
    'Callback',{@ (src,evt) updateDispPosImg},'fontsize',8,'Value',2);
menuSelectImgType.Position(3:4)=[140 18];
menuSelectImgType.Position(1:2)=[2 hpDisp_X.Position(4)-menuSelectImgType.Position(4)-35];   

    function updateDispPosImg
        imgnum = menuSelectImg.Value;
        switch menuSelectImgType.Value
            case 1
                set(hImg,'XData',data.X,'YData',data.Y,'CData',data.Z(:,:,imgnum));                
                set(hImg_K,'XData',data.f,'YData',data.f,'CData',data.ZfNorm(:,:,imgnum));
            case 2
                set(hImg,'XData',data.X,'YData',data.Y,'CData',data.ZNoFilter(:,:,imgnum));
                set(hImg_K,'XData',data.f,'YData',data.f,'CData',data.ZfNorm(:,:,imgnum));

        end
        if cAutoColor_X.Value;setClim('X');end  
        foo;
    end

% Table for changing display limits
tbl_dROI_X=uitable('parent',hpDisp_X,'units','pixels','RowName',{},...
    'columnname',{'x1','x2','y1','y2'},'UserData','X',...
    'ColumnEditable',[true true true true],'CellEditCallback',@tbl_dispROICB,...
    'ColumnWidth',{30 30 30 30},'FontSize',8,'Data',[1 size(Z,2) 1 size(Z,1)]);
tbl_dROI_X.Position(3:4)=tbl_dROI_X.Extent(3:4);
tbl_dROI_X.Position(1:2)=[2 menuSelectImgType.Position(2)-tbl_dROI_X.Position(4)-7];

% Button for maximizing the display limits
ttstr='Maximize display ROI to full image size.';
hbFullLim_X=uicontrol(hpDisp_X,'style','pushbutton','Cdata',img_Full,'Fontsize',10,...
    'Backgroundcolor','w','Callback',{@(~,~) chDispROI('max','X');},'ToolTipString',ttstr);
hbFullLim_X.Position = [tbl_dROI_X.Position(1)+tbl_dROI_X.Position(3) ...
    tbl_dROI_X.Position(2)-12 18 18];

% Button to snap display ROI to the data ROI
ttstr='Snap display ROI to data ROI.';
hbSnapLim_X=uicontrol(hpDisp_X,'style','pushbutton','Cdata',img_Snap,'Fontsize',10,...
    'Backgroundcolor','w','Position',hbFullLim_X.Position,...
    'Callback',{@(~,~) chDispROI('min','X');},'ToolTipString',ttstr);
hbSnapLim_X.Position(2) = [hbSnapLim_X.Position(2)+18];

% Button to enable GUI selection of display limits
ttstr='Select the display ROI.';
hbSlctLim_X=uicontrol(hpDisp_X,'style','pushbutton','Cdata',img_Select,'Fontsize',10,...
    'Backgroundcolor','w','Position',hbFullLim_X.Position,...
    'Callback',{@(src,evt) slctDispCB(tbl_dROI_X,'X',1)},'ToolTipString',ttstr);
hbSlctLim_X.Position(2) = [hbSnapLim_X.Position(2)+18];

% Table to adjust color limits on image
climtbl_X=uitable('parent',hpDisp_X,'units','pixels','RowName',{},'ColumnName',{'c1','c2'},...
    'Data',[0 12000],'ColumnWidth',{40,65},'ColumnEditable',[true true],...
    'CellEditCallback',{@climCB 'X'});
climtbl_X.Position(3:4)=climtbl_X.Extent(3:4);
climtbl_X.Position(1:2)=[2 tbl_dROI_X.Position(2)-climtbl_X.Position(4)-2];

cAutoColor_X=uicontrol(hpDisp_X,'style','checkbox','string','auto?',...
    'units','pixels','fontsize',8,'backgroundcolor','w','callback',{@(src,evt) cAutoCLIMCB(src,evt,'X')},...
    'enable','on','value',1);
cAutoColor_X.Position=[climtbl_X.Position(1)+climtbl_X.Position(3)+1 climtbl_X.Position(2) 50 20];


% Button group for deciding what the X/Y plots show
bgPlot_X = uibuttongroup(hpDisp_X,'units','pixels','backgroundcolor','w','BorderType','None',...
    'SelectionChangeFcn',@foo);  
bgPlot_X.Position(3:4)=[125 15];
bgPlot_X.Position(1:2)=[2 60];
    
% Radio buttons for cuts vs sum
rbCut_X=uicontrol(bgPlot_X,'Style','radiobutton','String','plot cut',...
    'Position',[0 0 60 15],'units','pixels','backgroundcolor','w','Value',1,...
    'fontsize',7);
rbSum_X=uicontrol(bgPlot_X,'Style','radiobutton','String','plot sum',...
    'Position',[60 0 60 15],'units','pixels','backgroundcolor','w','Value',0,...
    'fontsize',7);

% Checkbox for enabling display of the CoM analysis
cCoMStr_X=uicontrol(hpDisp_X,'style','checkbox','string','center of mass text',...
    'units','pixels','fontsize',7,'backgroundcolor','w','callback',@cCoMCB,...
    'enable','on','value',1);
cCoMStr_X.Position=[2 47 125 15];

% Checkbox for showing/hiding crosshair
cCross_X=uicontrol(hpDisp_X,'style','checkbox','string','cross hair',...
    'units','pixels','fontsize',7,'backgroundcolor','w','callback',@cCrossCB,...
    'enable','on','value',1,'UserData','X');
cCross_X.Position=[2 32 120 15];

% Checkbox for showing/hiding lattice
cDrawLattice=uicontrol(hpDisp_X,'style','checkbox','string','lattice grid',...
    'units','pixels','fontsize',7,'backgroundcolor','w','callback',@latticeGridCB,...
    'enable','on','value',1);
cDrawLattice.Position=[2 17 80 15];

% Checkbox for showing/hiding lattice
cTextLattice=uicontrol(hpDisp_X,'style','checkbox','string','lattice text',...
    'units','pixels','fontsize',7,'backgroundcolor','w','callback',@latticeTextCB,...
    'enable','on','value',1);
cTextLattice.Position=[cDrawLattice.Position(3)+cDrawLattice.Position(1)+5 17 120 15];


% Checkbox for showing/hiding digitization
cDrawAtoms=uicontrol(hpDisp_X,'style','checkbox','string','draw atoms',...
    'units','pixels','fontsize',7,'backgroundcolor','w','callback',@drawAtomsCB,...
    'enable','on','value',1);
cDrawAtoms.Position=[2 2 80 15];

%% Display Options Panel Focus

hpFocus = uipanel(hF,'units','pixels','backgroundcolor','w');
hpFocus.Position=[160 500 160 190];
set(hpFocus,'parent',hpDispOpt.Children(3))
hpFocus.Position=[1 1 hpDispOpt.Position(3) hpDispOpt.Position(4)];


% Table for changing display limits
tbl_dROI_focus=uitable('parent',hpFocus,'units','pixels','RowName',{},...
    'columnname',{'x1','x2','y1','y2'},'UserData','X',...
    'ColumnEditable',[true true true true],...
    'ColumnWidth',{30 30 30 30},'FontSize',8,'Data',[260 300 150 300],...
    'ColumnFormat',{'numeric', 'numeric','numeric','numeric'});
tbl_dROI_focus.Position(3:4)=tbl_dROI_focus.Extent(3:4);
tbl_dROI_focus.Position(1:2)=[2 5];

%% Display Options Panel

hpDisp_K = uipanel(hF,'units','pixels','backgroundcolor','w');
hpDisp_K.Position=[160 hpDisp_X.Position(2)-180 160 180];

set(hpDisp_K,'parent',hpDispOpt.Children(4));
hpDisp_K.Position=[1 1 hpDispOpt.Position(3) hpDispOpt.Position(4)];

% Table for changing display limits
tbl_dROI_K=uitable('parent',hpDisp_K,'units','pixels','RowName',{},...
    'columnname',{'kx1','kx2','ky1','ky2'},'UserData','K',...
    'ColumnEditable',[true true true true],'CellEditCallback',@tbl_dispROICB,...
    'ColumnWidth',{30 30 30 30},'FontSize',6,'Data',[-.5 .5 -.5 .5]);
tbl_dROI_K.Position(3:4)=tbl_dROI_K.Extent(3:4);
tbl_dROI_K.Position(1:2)=[2 hpDisp_K.Position(4)-tbl_dROI_K.Position(4)-15];

% Button for maximizing the display limits
ttstr='Maximize display ROI to full image size.';
hbFullLim_K=uicontrol(hpDisp_K,'style','pushbutton','Cdata',img_Full,'Fontsize',10,...
    'Backgroundcolor','w','Callback',{@(~,~) chDispROI('max','K');},'ToolTipString',ttstr);
hbFullLim_K.Position = [tbl_dROI_K.Position(1)+tbl_dROI_K.Position(3) ...
    tbl_dROI_K.Position(2)-15 18 18];

% Button to snap display ROI to the data ROI
ttstr='Snap display ROI to data ROI.';
hbSnapLim_K=uicontrol(hpDisp_K,'style','pushbutton','Cdata',img_Snap,'Fontsize',10,...
    'Backgroundcolor','w','Position',hbFullLim_K.Position,...
    'Callback',{@(~,~) chDispROI('min','K');},'ToolTipString',ttstr);
hbSnapLim_K.Position(2) = [hbSnapLim_K.Position(2)+18];

% Button to enable GUI selection of display limits
ttstr='Select the display ROI.';
hbSlctLim_K=uicontrol(hpDisp_K,'style','pushbutton','Cdata',img_Select,'Fontsize',10,...
    'Backgroundcolor','w','Position',hbFullLim_K.Position,...
    'Callback',{@(src,evt) slctDispCB(tbl_dROI_K,'K',1);},'ToolTipString',ttstr);
hbSlctLim_K.Position(2) = [hbSnapLim_K.Position(2)+18];

% Table to adjust color limits on image
climtbl_K=uitable('parent',hpDisp_K,'units','pixels','RowName',{},'ColumnName',{'c1','c2'},...
    'Data',[0 3e5],'ColumnWidth',{40,65},'ColumnEditable',[true true],...
    'CellEditCallback',{@climCB 'K'});
climtbl_K.Position(3:4)=climtbl_K.Extent(3:4);
climtbl_K.Position(1:2)=[2 tbl_dROI_K.Position(2)-climtbl_K.Position(4)-2];

cAutoColor_K=uicontrol(hpDisp_K,'style','checkbox','string','auto?',...
    'units','pixels','fontsize',7,'backgroundcolor','w','callback',{@(src,evt) cAutoCLIMCB(src,evt,'K')},...
    'enable','on','value',1);
cAutoColor_K.Position=[climtbl_K.Position(1)+climtbl_K.Position(3)+1 climtbl_K.Position(2) 50 20];


% Button group for deciding what the X/Y plots show
bgPlot_K = uibuttongroup(hpDisp_K,'units','pixels','backgroundcolor','w','BorderType','None',...
    'SelectionChangeFcn',@foo3);  
bgPlot_K.Position(3:4)=[125 15];
bgPlot_K.Position(1:2)=[2 climtbl_K.Position(2)-bgPlot_K.Position(4)-2];
    
% Radio buttons for cuts vs sum
rbCut_K=uicontrol(bgPlot_K,'Style','radiobutton','String','plot cut',...
    'Position',[0 0 60 20],'units','pixels','backgroundcolor','w','Value',1,...
    'fontsize',7);
rbSum_K=uicontrol(bgPlot_K,'Style','radiobutton','String','plot sum',...
    'Position',[60 0 60 20],'units','pixels','backgroundcolor','w','Value',0,...
    'fontsize',7);

% Checkbox for showing/hiding crosshair
cCross_K=uicontrol(hpDisp_K,'style','checkbox','string','cross hair',...
    'units','pixels','fontsize',7,'backgroundcolor','w','callback',@cCrossCB,...
    'enable','on','value',1,'UserData','K');
cCross_K.Position=[2 47 120 15];


% Checkbox for showing/hiding PSF in fourier domain
cPSF_K=uicontrol(hpDisp_K,'style','checkbox','string','1/e^2 point spread function',...
    'units','pixels','fontsize',7,'backgroundcolor','w','callback',@PSFKCB ,...
    'enable','on','value',1);
cPSF_K.Position=[2 32 160 15];

% Checkbox for showing/hiding lattice in fourier domain
cLat_KReticle=uicontrol(hpDisp_K,'style','checkbox','string','reciprocal lattice reticle',...
    'units','pixels','fontsize',7,'backgroundcolor','w','callback',@reciprocalLatticeReticleCB,...
    'enable','on','value',1);
cLat_KReticle.Position=[2 17 120 15];

% Checkbox for reciprocal lattice text
cLat_K_text=uicontrol(hpDisp_K,'style','checkbox','string','reciprocal lattice text',...
    'units','pixels','fontsize',7,'backgroundcolor','w','callback',@reciprocalLatticeTextCB ,...
    'enable','on','value',1);
cLat_K_text.Position=[2 2 120 15];

%% Bin and Digital Plot Display Options

hpDisp_B = uipanel(hF,'units','pixels','backgroundcolor','w');
hpDisp_B.Position=[160 hpDisp_K.Position(2)-120 160 120];

set(hpDisp_B,'parent',hpDispOpt.Children(5));
hpDisp_B.Position=[1 1 hpDispOpt.Position(3) hpDispOpt.Position(4)];


% Table for changing display limits
tbl_dROI_B=uitable('parent',hpDisp_B,'units','pixels','RowName',{},...
    'columnname',{'n1i','n1f','n2i','n2f'},'UserData','B',...
    'ColumnEditable',[true true true true],'CellEditCallback',@tbl_dispROICB,...
    'ColumnWidth',{30 30 30 30},'FontSize',8,'Data',[1 100 1 100]);
tbl_dROI_B.Position(3:4)=tbl_dROI_B.Extent(3:4);
tbl_dROI_B.Position(1:2)=[2 hpDisp_B.Position(4)-tbl_dROI_X.Position(4)-30];


% Button to snap display ROI to the data ROI
ttstr='Snap display ROI to data ROI.';
hbSnapLim_B=uicontrol(hpDisp_B,'style','pushbutton','Cdata',img_Snap,'Fontsize',10,...
    'Backgroundcolor','w',...
    'Callback',{@(~,~) chDispROI('min','B');},'ToolTipString',ttstr);
hbSnapLim_B.Position = [tbl_dROI_B.Position(1)+tbl_dROI_B.Position(3) ...
    tbl_dROI_B.Position(2) 18 18];

% Button to enable GUI selection of display limits
ttstr='Select the display ROI.';
hbSlctLim_B=uicontrol(hpDisp_B,'style','pushbutton','Cdata',img_Select,'Fontsize',10,...
    'Backgroundcolor','w','Position',hbSnapLim_B.Position,...
    'Callback',{@(src,evt) slctDispCB(tbl_dROI_B,'B',1)},'ToolTipString',ttstr);
hbSlctLim_B.Position(2) = [hbSnapLim_B.Position(2)+18];

% Table to adjust color limits on image
climtbl_B=uitable('parent',hpDisp_B,'units','pixels','RowName',{},'ColumnName',{'c1','c2'},...
    'Data',[0 12000],'ColumnWidth',{40,65},'ColumnEditable',[true true],...
    'CellEditCallback',{@climCB 'B'});
climtbl_B.Position(3:4)=climtbl_B.Extent(3:4);
climtbl_B.Position(1:2)=[2 tbl_dROI_B.Position(2)-climtbl_B.Position(4)-2];

cAutoColor_B=uicontrol(hpDisp_B,'style','checkbox','string','auto?',...
    'units','pixels','fontsize',8,'backgroundcolor','w','callback',{@(src,evt) cAutoCLIMCB(src,evt,'B')},...
    'enable','on','value',1);
cAutoColor_B.Position=[climtbl_B.Position(1)+climtbl_B.Position(3)+1 climtbl_B.Position(2) 50 20];


%% Count Histogram

hpDisp_HB = uipanel(hF,'units','pixels','backgroundcolor','w','parent',hpDispOpt.Children(7));
hpDisp_HB.Position=[0 1 hpDispOpt.Position(3) hpDispOpt.Position(4)];


menuSelectBinType=uicontrol('style','popupmenu','string',...
    {'raw','processed','normalized'},'units','pixels','parent',hpDisp_HB,...
    'Callback',{@ (src,evt) updateBinnedHistogramGraphics},'fontsize',8,'Value',2);
menuSelectBinType.Position(3:4)=[140 18];
menuSelectBinType.Position(1:2)=[2 hpDisp_HB.Position(4)-menuSelectBinType.Position(4)-35];   

%% Digitized Display

hpDisp_D = uipanel(hF,'units','pixels','backgroundcolor','w');
hpDisp_D.Position=[160 hpDisp_HB.Position(2)-110 160 110];

set(hpDisp_D,'parent',hpDispOpt.Children(8));
hpDisp_D.Position=[1 1 hpDispOpt.Position(3) hpDispOpt.Position(4)];



% Table for changing display limits
tbl_dROI_D=uitable('parent',hpDisp_D,'units','pixels','RowName',{},...
    'columnname',{'n1i','n1f','n2i','n2f'},'UserData','D',...
    'ColumnEditable',[true true true true],'CellEditCallback',@tbl_dispROICB,...
    'ColumnWidth',{30 30 30 30},'FontSize',8,'Data',[1 size(Z,2) 1 size(Z,1)]);
tbl_dROI_D.Position(3:4)=tbl_dROI_D.Extent(3:4);
tbl_dROI_D.Position(1:2)=[2 hpDisp_D.Position(4)-tbl_dROI_D.Position(4)-20];

% Button to snap display ROI to the data ROI
ttstr='Snap display ROI to data ROI.';
hbSnapLim_B=uicontrol(hpDisp_D,'style','pushbutton','Cdata',img_Snap,'Fontsize',10,...
    'Backgroundcolor','w',...
    'Callback',{@(~,~) chDispROI('min','D');},'ToolTipString',ttstr);
hbSnapLim_B.Position = [tbl_dROI_D.Position(1)+tbl_dROI_D.Position(3) ...
    tbl_dROI_D.Position(2) 18 18];

% Button to enable GUI selection of display limits
ttstr='Select the display ROI.';
hbSlctLim_B=uicontrol(hpDisp_D,'style','pushbutton','Cdata',img_Select,'Fontsize',10,...
    'Backgroundcolor','w','Position',hbSnapLim_B.Position,...
    'Callback',{@(src,evt) slctDispCB(tbl_dROI_D,'D',1)},'ToolTipString',ttstr);
hbSlctLim_B.Position(2) = [hbSnapLim_B.Position(2)+18];

% Checkbox for enabling display of the CoM analysis
cCoMStr_D=uicontrol(hpDisp_D,'style','checkbox','string','center of mass text',...
    'units','pixels','fontsize',7,'backgroundcolor','w','callback',@cCoMCB_D,...
    'enable','on','value',1);
cCoMStr_D.Position=[2 2 125 15];



%% Display Callbacks  

% Callback for changing display table ROI
    function tbl_dispROICB(src,evt)             
        ROI=src.Data;                               % Grab the new ROI             
        [ROI,err] = chDispROI(ROI,src.UserData);    % Update the ROI               
        if err
            src.Data(evt.Indices(2))=evt.PreviousData;
        else            
            src.Data = ROI;                     % Update table in case changed     
        end
    end
   

% Change display ROI on plots
    function [ROI,err] = chDispROI(ROI,img_type)  
        err = 0;        
        % Determine the ROI limits depending on image type
        switch img_type
            case 'X'
                ROI_LIM = [min(data.X) max(data.X) min(data.Y) max(data.Y)];
            case 'K'
                ROI_LIM = [min(data.f) max(data.f) min(data.f) max(data.f)];
            case 'dig'
                ROI_LIM = [min(data.X) max(data.X) min(data.Y) max(data.Y)];
            case 'hop'
                ROI_LIM = [min(data.X) max(data.X) min(data.Y) max(data.Y)];
            case 'B'
                ROI_LIM = [-500 500 -500 500];
            case 'D'
                ROI_LIM = [-500 500 -500 500];
            otherwise
                warning('OH GOD NO');            
        end        
        % Check for string inputs of max and min
        if isa(ROI,'char')
            if isequal(ROI,'max')
               ROI = ROI_LIM;
            else               
                switch img_type
                    case 'X'
                        aROI = tblROI.Data;
                    case 'K'
                        aROI = tblROIK.Data;
                    case 'B'
                        if isfield(data,'LatticeBin')
                           aROI = [min(data.LatticeBin(1).n1) max(data.LatticeBin(1).n1) ...
                               min(data.LatticeBin(1).n2) max(data.LatticeBin(1).n2)];
                        else
                            aROI=[0 100 0 100];
                        end
                    case 'D'
                        if isfield(data,'LatticeDig')
                           aROI = [min(data.LatticeDig(1).n1) max(data.LatticeDig(1).n1) ...
                               min(data.LatticeDig(1).n2) max(data.LatticeDig(1).n2)];
                        else
                            aROI=[0 100 0 100];
                        end                        
                end                      
                ROI=[min(aROI(:,1)) max(aROI(:,2)) min(aROI(:,3)) max(aROI(:,4))];
           end
        end
        
        % Make sure ROI is numeric
         if sum(~isnumeric(ROI)) || sum(isinf(ROI)) || sum(isnan(ROI))
            warning('Incorrect data type provided for ROI.');
            err = 1;
            return
        end        

        % Make sure ROI is increasing order
        if ROI(2)<=ROI(1) || ROI(4)<=ROI(3)
           warning('Bad ROI specification given.');
            err = 1;
            return
        end       
        
        % Keep ROI within the bounds
        if ROI(1)<ROI_LIM(1); ROI(1)=ROI_LIM(1); end       
        if ROI(3)<ROI_LIM(3); ROI(3)=ROI_LIM(3); end   
        if ROI(2)>ROI_LIM(2); ROI(2)=ROI_LIM(2);end       
        if ROI(4)>ROI_LIM(4); ROI(2)=ROI_LIM(4);end       
        
        % Attempt to change the display ROI
        try
            
            switch img_type
                case 'X'            
                    set(axImg,'XLim',ROI(1:2),'YLim',ROI(3:4));
                    % set(hAxX,'XLim',ROI(1:2));
                    % set(hAxY,'YLim',ROI(3:4));
                case 'K'
                    set(axImg_K,'XLim',ROI(1:2),'YLim',ROI(3:4));
                    % set(hAxX_K,'XLim',ROI(1:2));
                    % set(hAxY_K,'YLim',ROI(3:4));  
                case 'B'
                    set(axImg_B,'XLim',ROI(1:2),'YLim',ROI(3:4));
                case 'D'
                    set(axImg_D,'XLim',ROI(1:2),'YLim',ROI(3:4));
            end
            drawnow;
            resizePlots;
        catch ME
            warning('Unable to change display ROI.');            
            err = 1;
        end    
    end

% Callback for selecting an ROI based upon mouse click input.
    function slctDispCB(src,img_type,chDisp)
        disp(['Selecting display ROI.' ...
            ' Click two points that form the rectangle ROI.']);
        switch img_type
            case 'X'
                set(hp,'SelectedTab',tabX);
            case 'K'
                set(hp,'SelectedTab',tabK);
            case 'B'
                set(hp,'SelectedTab',tabB);
        end

        % Get click 1
        [x1,y1]=ginputMe(1);                                    % Get a mouse click
        if isequal(img_type,'X') || isequal(img_type,'B');x1=round(x1);y1=round(y1);end  % Round Pos.      
        p1=plot(x1,y1,'+','color','k','linewidth',1);           % Marker 1
        
        % Get click 2
        [x2,y2]=ginputMe(1);                                    % Get a mouse click
        if isequal(img_type,'X') || isequal(img_type,'B');x2=round(x2);y2=round(y2);end  % Round Pos.      
        p2=plot(x2,y2,'+','color','k','linewidth',1);           % Marker 2
enableInteractivity;

        delete(p1);delete(p2);                                  % Delete markers

        ROI=[min([x1 x2]) max([x1 x2]) ...                      % ROI
            min([y1 y2]) max([y1 y2])]; 

        if nargin == 3 && chDisp        
            [ROI,err] = chDispROI(ROI,img_type);                    % Change ROI

            % Update table
            if ~err                        
                src.Data = ROI;                     
            end  
        end
        
          
    end

% Callback for changing the color limits table
    function climCB(src,evt,img_type)
        try
            switch img_type
                case 'X'
                    ax = axImg;
                case 'K'
                    ax = axImg_K;
                case 'B'
                    ax = axImg_B;
            end                    
            ax.CLim=src.Data;
        catch ME
            warning('Unable to change color limits.');
            src.Data(evt.Indices(2))=evt.PreviousData;
        end
    end

    function cAutoCLIMCB(src,evt,img_type)   
        if src.Value 
            setClim(img_type);
        end      
    end 

    function setClim(img_type,CLIM)
        switch img_type
            case 'X'
                ax = axImg;
                z = hImg.CData;
                tbl = climtbl_X;
                N0 = round(max(max(z))*.95);
            case 'K'
                ax = axImg_K;
                z = hImg_K.CData;
                tbl = climtbl_K;
                N0 = round(max(max(z))*.95);                
                if isfield(data,'LatticeK')
                    imgnum = menuSelectImg.Value;
                    A = max([data.LatticeK(imgnum).A1 data.LatticeK(imgnum).A2]);
                    N0 = round(A);
                end
                
              case 'B'
                ax = axImg_B;
                z = hImg_B.CData;
                tbl = climtbl_B;
                N0 = round(max(max(z))*.95);

        end              
        if nargin == 1
            CLIM = [0 N0];
        end
        
        try
            set(ax,'CLim',CLIM);       
            tbl.Data = CLIM;            
        end
    end

    function chPlotCB(~,~)
        % Update Data Plot
        updateDataPlots(data);

        %   ;
        
        if hc_anlX_Gauss.Value
           % updateGaussPlot(data); 
           updateGaussGraphics
        end
    end

    function cCoMCB(src,~)       
        if ~isfield(data,'BoxCount')
            tCoMAnalysis.Visible='off';
            return;
        end
        set(tCoMAnalysis,'Visible',src.Value);    
    end

    function cCoMCB_D(src,~)       
        if ~isfield(data,'BoxCount')
            tD_SE.Visible='off';
            return;
        end
        set(tD_SE,'Visible',src.Value);    
    end

    function cCrossCB(src,evt)                
        switch src.UserData
            case 'X'
                set(pCrossX,'Visible',src.Value);
                set(pCrossY,'Visible',src.Value);
            case 'K'
                set(pCrossX_K,'Visible',src.Value);
                set(pCrossY_K,'Visible',src.Value);
        end
    end


%% Tabular Data Results Panel
% Panel for parameters and analysis results.

hpFit=uitabgroup(hpControl,'units','pixels');
% hpFit.Position=[320 0 300 ...
    % hF.Position(4)-(hpCam.Position(4)+hpSave.Position(4)+hpNav.Position(4))];
% hpFit.Position=[320 0 300 hpNav.Position(2)];

% hpFit.Position=[320 0 300 hpNav.Position(2)-200];

hpFit.Position(1) = hpDispOpt.Position(1);
hpFit.Position(3) = hpDispOpt.Position(3);
hpFit.Position(2) = 0;
hpFit.Position(4) = hpDispOpt.Position(2);

tabs(1)=uitab(hpFit,'Title','acq','units','pixels');
tabs(2)=uitab(hpFit,'Title','param','units','pixels');
tabs(3)=uitab(hpFit,'Title','flags','units','pixels');

tabs(4)=uitab(hpFit,'Title','pos','units','pixels');
tabs(5)=uitab(hpFit,'Title','fft','units','pixels');
tabs(6)=uitab(hpFit,'Title','hop','units','pixels');

% Table for acquisition
tbl_acq=uitable(tabs(1),'units','normalized','RowName',{},'fontsize',7,...
    'ColumnName',{},'ColumnWidth',{120 30 69},'columneditable',[false true false],...
    'Position',[0 0 1 1],'celleditcallback',@acqTblCB,'ColumnFormat',...
    {'char','numeric','char'},'enable','off');

    function loadAcquisitionSettings
        if ~isValidAcq(acq)            
            warning('Invalid acquisition settings to send to camera.');
            return; 
        end
        
        % Update descriptions
        desc=acqDescription(acq);

        % Update table with settings and descriptions        
        tbl_acq.Data=[fieldnames(acq), ...
            struct2cell(acq), struct2cell(desc)];     
    
        % Send commands to camera
        [out,bNames]=setAcqSettings(acq);
        
        % Process output
        if out
            disp(' All acquisition settings written');
        else
            disp(' ');
            str=strcat(bNames','\n');
            str=[str{:}];
            str(end-1:end)=[];
            warning(sprintf([' Unexpected return result. ' ...
                ' Possible issue with : \n' ...
                str]));
         end   
    end

    function acqTblCB(src,evt)  
        % Temp acquisition setting with new data
        f=src.Data{evt.Indices(1),1};        
        temp_acq=acq;
        temp_acq.(f)=evt.NewData;
        
        if isValidAcq(temp_acq)
            acq=temp_acq;
            desc=acqDescription(acq);

            loadAcquisitionSettings;
        else
            src.Data{evt.Indices(1),evt.Indices(2)}=evt.PreviousData;     
        end   
        
        temp = get(src,'Data');
        set(src,'Data',{ 'dummy' });
        set(src,'Data', temp );
    end


% Table for run parameters
tbl_params=uitable(tabs(2),'units','normalized','RowName',{},'fontsize',7,...
    'ColumnName',{},'ColumnWidth',{215 40},'columneditable',[false false],...
    'Position',[0 0 1 1]);

% Table for run flags
tbl_flags=uitable(tabs(3),'units','normalized','RowName',{},'fontsize',7,...
    'ColumnName',{},'ColumnWidth',{215 40},'columneditable',[false false],...
    'Position',[0 0 1 1]);

% Table for pos analysis outputs
tbl_pos_analysis=uitable(tabs(4),'units','normalized','RowName',{},'ColumnName',{},...
    'fontsize',8,'ColumnWidth',{160 95},'columneditable',false(ones(1,2)),...
    'Position',[0 0 1 1]);

% Table for fft analysis outputs
tbl_fft_analysis=uitable(tabs(5),'units','normalized','RowName',{},'ColumnName',{},...
    'fontsize',8,'ColumnWidth',{160 95},'columneditable',false(ones(1,2)),...
    'Position',[0 0 1 1]);

% Table for hopping analysis outputs
tbl_fidelity_analysis=uitable(tabs(6),'units','normalized','RowName',{},'ColumnName',{},...
    'fontsize',8,'ColumnWidth',{160 95},'columneditable',false(ones(1,2)),...
    'Position',[0 0 1 1]);


%% Initialize the image panel

hp=uitabgroup(hF,'units','pixels','Position',...
    [400 0 hF.Position(3)-200 hF.Position(4)-40],...
    'SelectionChangedFcn',@hp_tab_cb);
hp.Position=[hpControl.Position(3) 0 hF.Position(3)-hpControl.Position(3) hF.Position(4)];

    function hp_tab_cb(src,evt)
        newTitle = evt.NewValue.Title;

        ind = 0;
        for jj=1:length(hpDispOpt.Children)
            if isequal(hpDispOpt.Children(jj).Title,newTitle)
                ind = jj;
            end
        end
        if ind
            hpDispOpt.SelectedTab=hpDispOpt.Children(ind);
        end

       
    end


% Tab Groups for each display
tabX=uitab(hp,'Title','position','units','pixels','backgroundcolor','w');
tabStripe=uitab(hp,'Title','stripe','units','pixels','backgroundcolor','w');
tabFocus=uitab(hp,'Title','focus','units','pixels','backgroundcolor','w');

tabH=uitab(hp,'Title','histogram','units','pixels','backgroundcolor','w');
tabK=uitab(hp,'Title','momentum','units','pixels','backgroundcolor','w');
tabB=uitab(hp,'Title','binned','units','pixels','backgroundcolor','w');

tabHB=uitab(hp,'Title','binned histogram','units','pixels','backgroundcolor','w');
tabD=uitab(hp,'Title','digitized','units','pixels','backgroundcolor','w');
% tabC=uitab(hp,'Title','correlators','units','pixels','backgroundcolor','w');
% tabFidelity=uitab(hp,'Title','fidelity','units','pixels','backgroundcolor','w');

% Define spacing for images, useful for resizing
l=80;   % Left gap for fitting and data analysis summary

% Resize the X/Y plots and images to maximize area in figure and line up
% properly

    function resizePlot(axi,cbar,axx,axy)
        if hp.Position(3)<200 || hp.Position(4)<200
            return;
        end 
        axi.Position = [40 110 hp.Position(3)-200 hp.Position(4)-200];             
        % Get the aspect ratio of plot objects
        Rimg=axi.PlotBoxAspectRatio;Rimg=Rimg(1)/Rimg(2);
        Rax=axi.Position(3:4);Rax=Rax(1)/Rax(2);        
        % Size of plot objects (position is weird in axis equal tight);
        if Rax>Rimg
            h1=axi.Position(4);
            w1=axi.Position(4)*Rimg;   
            axx.Position=[40+(axi.Position(3)-w1)/2 axi.Position(2)-l w1 80-15];
            axy.Position=[40+(axi.Position(3)+w1)/2+15 axi.Position(2) 80 h1];
        else
            w1=axi.Position(3);
            h1=w1/Rimg;            
            axx.Position=[axi.Position(1) 110+(axi.Position(4)-h1)/2-l ...
                w1 80-15];
            axy.Position=[axi.Position(1)+axi.Position(3)+15 ...
                110+(axi.Position(4)-h1)/2 l h1];            
        end        
        % Match cut limits with the images limits (isn't necessary?)
        % set(axx,'XLim',axi.XLim,'XTick',axi.XTick);
        % set(axy,'YLim',axi.YLim,'YTick',axi.YTick);        
        % Move the colorbar
        cbar.Position=[axx.Position(1) axy.Position(2)+axy.Position(4)+23 ...
            axx.Position(3) 15]; 
        drawnow;        
    end

    function resizePlots       
        resizePlot(axImg,cBar,hAxX,hAxY);
        resizePlot(axImg_K,cBar_K,hAxX_K,hAxY_K);        
    end

%% Position Images

% Initialize image axis
axImg=axes('parent',tabX);cla
co=get(gca,'colororder');
hImg=imagesc(data.X,data.Y,data.Z);
set(axImg,'box','on','linewidth',.1,'fontsize',10,'units','pixels',...
    'XAxisLocation','top','colormap',colormap(cmap),...
    'xcolor',co(4,:),'ycolor',co(4,:),'YDir','normal','UserData','X');

hold on
axImg.Position=[50 150 tabX.Position(3)-200 tabX.Position(4)-200];
axis equal tight

% This seems to mess with axis interactivity in certain versions of matlab
% pxinfo = impixelinfo(hF,hImg);
enableDefaultInteractivity(axImg);


% Cross Hair Plots
pCrossX=plot([1 512],[512/2 512/2],'-',...
    'color',[1 0 0 .2],'linewidth',1,'hittest','off','parent',axImg);
pCrossY=plot([512/2 512/2],[1 512],'-',...
    'color',[1 0 0 .2],'linewidth',1,'hittest','off','parent',axImg);

% Initialize Grid Objects
pGrid = line([1 512],[1 512],'linestyle','-',...
    'color',[.5 .5 .5 .5],'linewidth',1,'Visible','off','parent',axImg,...
    'hittest','off');

% Initialize Atom Dots
pAtoms = plot(5,5,'o',...
    'color','w','linewidth',1,'Visible','off','parent',axImg,...
    'hittest','on','markersize',8);


% file name string
tImageFile=text(3,3,'FILENAME','units','pixels','fontsize',8,'fontweight','bold',...
    'horizontalalignment','left','verticalalignment','bottom','margin',1,...
    'interpreter','none','backgroundcolor',[1 1 1 .5]);

% Box Count Analysis String
tCoMAnalysis=text(.99,0.01,'FILENAME','units','normalized','fontsize',9,'fontweight','bold',...
    'horizontalalignment','right','verticalalignment','bottom','margin',1,...
    'interpreter','latex',...
    'color','k','backgroundcolor',[1 1 1 .7]);

% Box Count Analysis String
tTopLeft=text(.01,.99,'FILENAME','units','normalized','fontsize',9,'fontweight','bold',...
    'horizontalalignment','left','verticalalignment','top','margin',1,...
    'interpreter','latex',...
    'color','k','backgroundcolor',[1 1 1 .7],'Visible','off');


% Box for ROI (this will become an array later)
pROI=rectangle('position',[1 1 512 512],'edgecolor',co(1,:),'linewidth',2,'hittest','off');

% Reticle for gaussian fit (this will become an array later)
pGaussRet=plot(0,0,'-','linewidth',1,'Visible','off','color',co(1,:),'hittest','off');
% Color bar
cBar=colorbar('fontsize',8,'units','pixels','location','northoutside');

pPCA(1)=plot(0,0,'-','linewidth',1,'color','r','hittest','off');
pPCA(2)=plot(0,0,'-','linewidth',1,'color','r','hittest','off');

axImg.CLim=climtbl_X.Data;
drawnow;

% X Cut/Sum Axis
hAxX=axes('box','on','linewidth',1,'fontsize',10,...
    'XAxisLocation','Bottom','units','pixels','parent',tabX);
hAxX.Position=[axImg.Position(1) axImg.Position(2)-l axImg.Position(3) l];
hold on
% Add X data data and fit plots
pX=plot(data.X,ones(length(data.X),1),'k.-');
pXF=plot(data.X,ones(length(data.X),1),'-','Visible','off','color',co(1,:),'linewidth',2);


% Y Cut/Sum Axis
hAxY=axes('box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'YAxisLocation','Right','YDir','normal','parent',tabX);
hAxY.Position=[axImg.Position(1)+axImg.Position(3) axImg.Position(2) l axImg.Position(4)];
hold on
% Add Y data data and fit plotsf
pY=plot(ones(length(data.Y),1),data.Y,'k.-'); 
pYF=plot(data.X,ones(length(data.X),1),'-','Visible','off','color',co(1,:),'linewidth',2);

set(axImg,'XLim',tbl_dROI_X.Data(1:2),'YLim',tbl_dROI_X.Data(3:4));
drawnow

linkaxes([axImg hAxY],'y');
linkaxes([axImg hAxX],'x');

%% Stripe
% 
% ax_stripe_img = subplot(2,3,[1 2 4 5],'parent',tabStripe);
% stripe_hImgStripe=imagesc(data.X,data.Y,data.Z,'parent',ax_stripe_img);
% set(ax_stripe_img,'box','on','linewidth',.1,'fontsize',10,...
%     'XAxisLocation','bottom','colormap',colormap(cmap),...
%     'YDir','normal');
% hold on
% axis equal tight
% 
% stripe_pFringe=plot(0,0,'-','color',co(1,:),'linewidth',1);
% stripe_pPerp=plot(0,0,'-','color',co(5,:),'linewidth',1);     
% stripe_pBar=plot(0,0,'-','color','w','linewidth',1);     
% stripe_pAngleCirc=plot(0,0,'-','color','w','linewidth',1);     
% stripe_pCloudEllipse=plot(0,0,'-','color','w','linewidth',1);     
% 
% 
% % clear stripe_lines
% % stripe_lines(1) = plot(0,0,'-','color','r','linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
% % stripe_lines(2) = plot(0,0,'-','color','r','linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
% % stripe_lines(3) = plot(0,0,'-','color','r','linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
% % stripe_lines(4) = plot(0,0,'-','color','r','linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
% % stripe_lines(5) = plot(0,0,'-','color','r','linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
% % stripe_lines(6) = plot(0,0,'-','color','r','linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
% % stripe_lines(7) = plot(0,0,'-','color','r','linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
% % stripe_lines(8) = plot(0,0,'-','color','r','linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
% % stripe_lines(9) = plot(0,0,'-','color','r','linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
% % stripe_lines(10) = plot(0,0,'-','color','r','linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
% % 
% 
% % file name string
% stripe_str_bottom_left=text(3,3,'FILENAME','units','pixels','fontsize',10,'fontweight','bold',...
%     'horizontalalignment','left','verticalalignment','bottom','margin',1,...
%     'interpreter','none','backgroundcolor',[1 1 1 .5]);
% 
% % Box Count Analysis String
% stripe_str_bottom_right=text(.99,0.01,'FILENAME','units','normalized','fontsize',10,'fontweight','normal',...
%     'horizontalalignment','right','verticalalignment','bottom','margin',1,...
%     'interpreter','latex',...
%     'color','k','backgroundcolor',[1 1 1 .7]);
% 
% 
% ax_stripe_fit = subplot(2,3,3,'parent',tabStripe);
% set(ax_stripe_fit,'box','on','fontsize',10,...
%     'XAxisLocation','Bottom');
% 
% % Add X data data and fit plots
% stripe_pSum2_fit=plot(0,0,'k--','linewidth',1,'parent',ax_stripe_fit);
% hold on
% stripe_pSum2_data=plot(0,0,'-','color',co(5,:),'linewidth',1,'parent',ax_stripe_fit);
% 
% stripe_pSum1_fit=plot(0,0,'-','linewidth',2,'color',co(2,:),'parent',ax_stripe_fit);
% hold on
% stripe_pSum1_data=plot(0,0,'-','color',co(1,:),'linewidth',1,'parent',ax_stripe_fit);
% xlabel('rotated position (px)');
% ylabel('sum counts');
% 
% legend([stripe_pSum1_data, stripe_pSum1_fit, stripe_pSum2_data, stripe_pSum2_fit],...
%     {'fringe','fringe fit','perp','perp fit'},'fontsize',6,...
%     'location','northeast')
% text(.01,.98,'projected sum counts','units','normalized',...
%     'verticalalignment','top','fontsize',6,'parent',ax_stripe_fit);
% 
% ax_stripe_focus = subplot(2,3,6,'parent',tabStripe);
% set(ax_stripe_focus,'box','on','fontsize',10,...
%     'XAxisLocation','Bottom');
% 
% stripe_pFocus=plot(0,0,'o-','color','k','linewidth',1);     
% hold on
% stripe_pFocusFit=plot(0,0,'--','color','k','linewidth',1);     
% 
% yyaxis right
% 
% stripe_pFocus2=plot(0,0,'o-','color','r','linewidth',1);     
% yyaxis left
%% Histgoram

ax_h1 = subplot(2,2,1,'parent',tabH);
pHist1 = bar(1:100,1:100,'parent',ax_h1,'linestyle','none');
ylabel('occurences');
xlabel('counts/pixel');
hold on
set(ax_h1,'box','on','linewidth',.1,'fontsize',8,'units','normalized',...
    'XAxisLocation','bottom','YDir','normal','UserData','H1',...
    'YAxisLocation','left');
ylim(ax_h1,'auto');
text(.99,.99,'processed, linear scale','units','normalized','fontsize',10,...
    'verticalalignment','top','horizontalalignment','right');

ax_h2 = subplot(2,2,2,'parent',tabH);
pHist2 = bar(1:100,1:100,'parent',ax_h2,'linestyle','none');
ylabel('occurences');
xlabel('counts/pixel');
hold on
set(ax_h2,'box','on','linewidth',.1,'fontsize',8,'units','normalized',...
    'XAxisLocation','bottom','YDir','normal','UserData','H2',...
    'YAxisLocation','left','yScale','log');
ylim(ax_h2,'auto');
text(.99,.99,'processed, log scale','units','normalized','fontsize',10,...
    'verticalalignment','top','horizontalalignment','right');

ax_h3 = subplot(2,2,3,'parent',tabH);
pHist3 = bar(1:100,1:100,'parent',ax_h3,'linestyle','none');
ylabel('occurences');
xlabel('counts/pixel');
hold on
set(ax_h3,'box','on','linewidth',.1,'fontsize',8,'units','normalized',...
    'XAxisLocation','bottom','YDir','normal','UserData','H1',...
    'YAxisLocation','left');
ylim(ax_h3,'auto');
text(.99,.99,'no filtering, linear scale','units','normalized','fontsize',10,...
    'verticalalignment','top','horizontalalignment','right');

ax_h4 = subplot(2,2,4,'parent',tabH);
pHist4 = bar(1:100,1:100,'parent',ax_h4,'linestyle','none');
ylabel('occurences');
xlabel('counts/pixel');
hold on
set(ax_h4,'box','on','linewidth',.1,'fontsize',8,'units','normalized',...
    'XAxisLocation','bottom','YDir','normal','UserData','H2',...
    'YAxisLocation','left','yScale','log');
ylim(ax_h4,'auto');
text(.99,.99,'no filtering, log scale','units','normalized','fontsize',10,...
    'verticalalignment','top','horizontalalignment','right');

%% Momentum Images

% Initialize image axis
axImg_K=axes('parent',tabK);cla
hImg_K=imagesc(data.X,data.Y,data.Z,'parent',axImg_K);
set(axImg_K,'box','on','linewidth',.1,'fontsize',10,'units','pixels',...
    'XAxisLocation','top','colormap',colormap(cmap),'YDir','normal','UserData','K');
hold on
axImg_K.Position=[50 150 tabX.Position(3)-200 tabX.Position(4)-200];
axis equal tight

clear pKReticles
pKReticles(1) = plot(0,0,'-','color',[1 1 1],'linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
pKReticles(2) = plot(0,0,'-','color',[1 1 1],'linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
pKReticles(3) = plot(0,0,'-','color',[1 1 1],'linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
pKReticles(4) = plot(0,0,'-','color',[1 1 1],'linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
pKReticles(5) = plot(0,0,'-','color',[1 1 1],'linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
pKReticles(6) = plot(0,0,'-','color',[1 1 1],'linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
pKReticles(7) = plot(0,0,'-','color',[1 1 1],'linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');
pKReticles(8) = plot(0,0,'-','color',[1 1 1],'linewidth',1,'parent',axImg_K,'Visible','off','hittest','off');

pKPSF = plot(0,0,'--','color',[1 1 1],'parent',axImg_K,'Visible','off','hittest','off');

% Cross Hair Plots
pCrossX_K=plot([1 512],[512/2 512/2],'-','color',[1 0 0 .2],'linewidth',1,'parent',axImg_K,'hittest','off');
pCrossY_K=plot([512/2 512/2],[1 512],'-','color',[1 0 0 .2],'linewidth',1,'parent',axImg_K,'hittest','off');

% file name string
tImageFile_K=text(3,3,'FILENAME','units','pixels','fontsize',8,'fontweight','bold',...
    'horizontalalignment','left','verticalalignment','bottom','margin',1,...
    'interpreter','none','backgroundcolor',[1 1 1 .5],'parent',axImg_K);

% box count analysis sytring
tTopLeftK=text(.01,.99,'FILENAME','units','normalized','fontsize',9,'fontweight','bold',...
    'horizontalalignment','left','verticalalignment','top','margin',1,...
    'interpreter','latex',...
    'color','k','backgroundcolor',[1 1 1 .3],'Visible','off');

pROI_K=rectangle('position',[-.5 -.5 1 1],'edgecolor',co(1,:),'linewidth',2,'parent',axImg_K,'hittest','off');
cBar_K=colorbar('fontsize',8,'units','pixels','location','northoutside');
axImg_K.CLim=climtbl_K.Data;
drawnow;

% X Cut/Sum Axis
hAxX_K=axes('box','on','linewidth',1,'fontsize',10,...
    'XAxisLocation','Bottom','units','pixels','parent',tabK);
hAxX_K.Position=[axImg_K.Position(1) axImg_K.Position(2)-l axImg_K.Position(3) l];
hold on
% Add X data data and fit plots
pX_K=plot(data.X,ones(length(data.X),1),'k.-');


% Y Cut/Sum Axis
hAxY_K=axes('box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'YAxisLocation','Right','YDir','normal','parent',tabK);
hAxY_K.Position=[axImg_K.Position(1)+axImg_K.Position(3) axImg_K.Position(2) l axImg_K.Position(4)];
hold on
% Add Y data data and fit plots
pY_K=plot(ones(length(data.Y),1),data.Y,'k.-'); 
set(axImg_K,'XLim',tbl_dROI_K.Data(1:2),'YLim',tbl_dROI_K.Data(3:4));
drawnow

% Link Axes
linkaxes([axImg_K hAxY_K],'y');
linkaxes([axImg_K hAxX_K],'x');
set(axImg_K,'XLim',tbl_dROI_K.Data(1:2),'YLim',tbl_dROI_K.Data(3:4));

%% Binned Image

axImg_B=axes('parent',tabB);cla
hImg_B=imagesc(1:500,1:500,zeros(500,500),'parent',axImg_B);
set(axImg_B,'box','on','linewidth',.1,'fontsize',8,'units','normalized',...
    'XAxisLocation','bottom','colormap',colormap(cmap),'YDir','normal','UserData','B',...
    'YAxisLocation','left');
hold on
c=colorbar;
c.Label.String = 'counts/site';
axis equal tight
xlabel('lattice site (a_1)','fontsize',8);
ylabel('lattice site (a_2)','fontsize',8);
caxis([0 1]);

% Box for ROI (this will become an array later)
% pROIB=rectangle('position',[1 1 500 500],'edgecolor',co(1,:),'linewidth',2,'hittest','off');

%% Binned Histgoram

%% Digital Image

axImg_D=axes('parent',tabD);
hImg_D=imagesc(1:500,1:500,zeros(500,500),'parent',axImg_D);
set(axImg_D,'box','on','linewidth',.1,'fontsize',8,'units','normalized',...
    'XAxisLocation','bottom','colormap',colormap(cmap),'YDir','normal','UserData','B',...
    'YAxisLocation','left');
hold on
axis equal tight
xlabel('lattice site (a_1)','fontsize',8);
ylabel('lattice site (a_2)','fontsize',8);
caxis([0 1]);

% Box Count Analysis String
tD_SE=text(.99,0.01,'FILENAME','units','normalized','fontsize',9,'fontweight','bold',...
    'horizontalalignment','right','verticalalignment','bottom','margin',1,...
    'interpreter','latex',...
    'color','k','backgroundcolor',[1 1 1 .7]);


% file name string
tD_SW=text(3,3,'FILENAME','units','pixels','fontsize',8,'fontweight','bold',...
    'horizontalalignment','left','verticalalignment','bottom','margin',1,...
    'interpreter','none','backgroundcolor',[1 1 1 .5]);


%% Graphical Callbacks
    function updateGraphics  
        updateDispPosImg;   
        updatePositionAnalysisGraphics
        updatePositionGraphics;
        updateMomentumGraphics;          
        updateBinnedGraphics;
        updateBinnedHistogramGraphics;               
        updateDigitalGraphics;
    end

%% Lattice Grid Callbacks
  function updateGridGraphics       
        imgnum = menuSelectImg.Value;
        if ~isfield(data,'LatticeBin');return;end              
        
        n1i = min(data.LatticeBin(imgnum).n1);
        n1f = max(data.LatticeBin(imgnum).n1);        
        n2i = min(data.LatticeBin(imgnum).n2);
        n2f = max(data.LatticeBin(imgnum).n2);     
        
        a1 = data.LatticeBin(imgnum).a1;
        a2 = data.LatticeBin(imgnum).a2;
        
        theta=acos(sum(a1.*a2)/(norm(a1)*norm(a2)))*180/pi;
        
        p = data.LatticeBin(imgnum).p;
        tTopLeft.String = ['$\vec{a}_1 = (' num2str(round(a1(1),4)) ',' num2str(round(a1(2),4)) ')$' newline ...
            '$\vec{a}_2 = (' num2str(round(a2(1),4)) ',' num2str(round(a2(2),4)) ')$' newline ...
            '$\vec{a}_1\cdot\vec{a}_2 = a_1a_2\cos(' num2str(theta,4) '^\circ )$' newline ...
            '$\phi = 2\pi\times (' num2str(round(p(1),3)) ',' num2str(round(p(2),3)) ')$'];            
        
        % n1 Grid
        rVec = a1*([n1i:n1f]+p(1)+0.5); % All lattice sites 1        
        R1 = rVec + a2*(n2i-1);
        R2 = rVec + a2*(n2f+1);
        Q = nan(1,size(rVec,2));        
        X1 = [R1(1,:); R2(1,:); Q];X1=X1(:);
        Y1 = [R1(2,:); R2(2,:); Q];Y1=Y1(:);

        % n2 Grid
        rVec = a2*([n2i:n2f]+p(2)+0.5); % All lattice sites 1        
        R1 = rVec + a1*(n1i-1);
        R2 = rVec + a1*(n1f+1);
        Q = nan(1,size(rVec,2));        
        X2 = [R1(1,:); R2(1,:); Q];X2=X2(:);
        Y2 = [R1(2,:); R2(2,:); Q];Y2=Y2(:);
        
        set(pGrid,'XData',[X1; X2],'YData',[Y1; Y2]);
  end           


function drawAtomsCB(~,~)
    if ~isfield(data,'LatticeDig')
       pAtoms.Visible='off';
       return;
    end
    
    
   if cDrawAtoms.Value        
        pAtoms.Visible='on';
   else
        pAtoms.Visible='off';
   end  
        
end

    function latticeGridCB(src,~)      
        if ~isfield(data,'LatticeBin')
            pGrid.Visible='off';
            return;
        end
        
        if src.Value        
            pGrid.Visible='on';
        else
            pGrid.Visible='off';
        end        
    end

   function latticeTextCB(src,evt)          
        if src.Value && isfield(data,'LatticeBin')
            tTopLeft.Visible='on';
        else
            tTopLeft.Visible='off';
        end
   end

%% Position Callbacks

    function updatePositionGraphics  
        cCrossCB(cCross_X);
        updateGridGraphics;
        latticeGridCB(cDrawLattice);
        latticeTextCB(cTextLattice);  
        set(tImageFile,'String',[data.Name ' (' num2str(menuSelectImg.Value) ')']);  

        
    end

    % function updateSharpnessGraphics
    %     if ~isfield(data,'SharpnessScore')
    %        t_pos_TopRight.Visible='off';
    %        return
    %     end
    %     imgnum = menuSelectImg.Value;
    %     sh1 = data.SharpnessScore(imgnum); 
    %     sh2 = data.SharpnessScoreNoFilter(imgnum); 
    %     str = ['sharpness scores $(' num2str(sh1,'%.4e') ',' num2str(sh2,'%.4e') ')$'];
    %     set(t_pos_TopRight,'String',str); 
    % end

    function updateBoxGraphics
        if ~isfield(data,'BoxCount')
           tCoMAnalysis.Visible='off';
           return
        end
        imgnum = menuSelectImg.Value;
        bc = data.BoxCount(imgnum);        
        % Update box count string
        str=[ num2str(bc.Npeak,'%.2e') ' max counts ' newline ...
            num2str(bc.Nraw,'%.2e') ' counts' newline ...
            '$(X_\mathrm{c},Y_\mathrm{c}) = ' '('  num2str(round(bc.Xc,1)) ',' ...
            num2str(round(bc.Yc,1)) ')$' newline ...
            '$(\sigma_X,\sigma_Y) = ' '('  num2str(round(bc.Xs,1)) ',' ...
            num2str(round(bc.Ys,1)) ')$']; 
        %Update box count string object
        set(tCoMAnalysis,'String',str);   

        stranl={'box sum (counts)',bc.Nraw;
                'box peak (counts)',bc.Npeak;
                'box Yc (px)',bc.Yc;
                'box Xc (px)',bc.Xc;
                ['box Y' char(963) ' (px)'],bc.Ys;
                ['box X' char(963) ' (px)'],bc.Xs};       
        tbl_pos_analysis.Data=[tbl_pos_analysis.Data; stranl]; 
    end

%% PCA Stuff
    function updatePCAGraphics
        if ~isfield(data,'PCA')
            pPCA(1).Visible = 'off';
            pPCA(2).Visible = 'off';
            return;
        end

        out=data.PCA;        
        x1=out.Mean(1)+out.Radii(1)*out.PCA(1,1)*[-1 1];
        y1=out.Mean(2)+out.Radii(1)*out.PCA(2,1)*[-1 1];
        x2=out.Mean(1)+out.Radii(2)*out.PCA(1,2)*[-1 1];
        y2=out.Mean(2)+out.Radii(2)*out.PCA(2,2)*[-1 1];        
        % PCA analysis table string
        stranl={'','';
            ['pca ' char(952) '1 (deg)'] ,atan(out.PCA(2,1)/out.PCA(1,1))*180/pi;
            ['pca ' char(952) '2 (deg)'],atan(out.PCA(2,2)/out.PCA(1,2))*180/pi;
            ['pca ' char(963) '1 (px)'],out.Radii(1);
            ['pca ' char(963) '2 (px)'],out.Radii(2);
            ['pca xc (px)'],out.Mean(1);
            ['pca yc (px)'],out.Mean(2);};
        tbl_pos_analysis.Data=[tbl_pos_analysis.Data; stranl];        
        set(pPCA(1),'XData',x1,'YData',y1,'Visible','on');
        set(pPCA(2),'XData',x2,'YData',y2,'Visible','on');
    end
%% Histgoram Callbacks


    function updatePositionHistogramGraphics        
        if ~isfield(data,'Histogram') || ~isfield(data,'HistogramNoFilter')
            return;
        end
       imgnum = menuSelectImg.Value;          
        set(pHist1,'XData',data.Histogram(imgnum).Centers,...
            'YData',data.Histogram(imgnum).N);     
        set(pHist2,'XData',data.Histogram(imgnum).Centers,...
            'YData',data.Histogram(imgnum).N);    
        set(pHist3,'XData',data.HistogramNoFilter(imgnum).Centers,...
            'YData',data.HistogramNoFilter(imgnum).N);     
        set(pHist4,'XData',data.HistogramNoFilter(imgnum).Centers,...
            'YData',data.HistogramNoFilter(imgnum).N);   
    end

%% Momentum Callbacks

    function updateMomentumGraphics                   

        if isfield(data,'ZfNorm')
            set(hImg_K,'XData',data.f,'YData',data.f,'CData',...
                data.ZfNorm(:,:,menuSelectImg.Value));
            if cAutoColor_K.Value;setClim('K');end  
        end    
        set(tImageFile_K,'String',[data.Name ' (' num2str(menuSelectImg.Value) ')']);
        
        updateReciprocalReticleGraphics;
        reciprocalLatticeReticleCB(cLat_KReticle);
        
        updateReciprocalTextGraphics;
        reciprocalLatticeTextCB(cLat_K_text);
        
        updatePSFKGraphic;
        PSFKCB(cPSF_K);
    end

    function reciprocalLatticeReticleCB(src,evt)
        if ~isfield(data,'LatticeK')
            for ii=1:length(pKReticles)          
                pKReticles(ii).Visible='off';        
            end
            return;
        end
        if src.Value        
            for ii=1:length(pKReticles)          
                pKReticles(ii).Visible='on';        
            end            
        else
            for ii=1:length(pKReticles)          
                pKReticles(ii).Visible='off';        
            end
        end        
    end

    function updateReciprocalTextGraphics
        imgnum = menuSelectImg.Value;
        if ~isfield(data,'LatticeK')
            tTopLeftK.Visible = 'off';
            return;
        end     
        k1 = data.LatticeK(imgnum).k1;
        k2 = data.LatticeK(imgnum).k2;      
        A1 = data.LatticeK(imgnum).A1;
        A2 = data.LatticeK(imgnum).A2;
        s1 = data.LatticeK(imgnum).s1;
        s2 = data.LatticeK(imgnum).s2;
        Zfsum = data.LatticeK(imgnum).Zfsum;
        NkScore = data.LatticeK(imgnum).NkScore;
        
        theta=acos(sum(k1.*k2)/(norm(k1)*norm(k2)))*180/pi;        
        tTopLeftK.String = ['$\vec{k}_1 = (' num2str(round(k1(1),4)) ',' num2str(round(k1(2),4)) ')$' newline ...
            '$\vec{k}_2 = (' num2str(round(k2(1),4)) ',' num2str(round(k2(2),4)) ')$' newline ...
            '$\vec{k}_1\cdot\vec{k}_2 = k_1k_2\cos(' num2str(theta,4) '^\circ )$' newline ...
            '$N_\mathrm{k}/N_\mathrm{tot}=' num2str(NkScore,'%.3g') '$'];
    end

    function reciprocalLatticeTextCB(src,evt)
        if ~isfield(data,'LatticeK')
            tTopLeftK.Visible = 'off';
            return;
        end 
        if src.Value
           tTopLeftK.Visible='on';
        else
           tTopLeftK.Visible='off';
        end
    end

    function updateReciprocalReticleGraphics
        imgnum = menuSelectImg.Value;        
        if ~isfield(data,'LatticeK')
            return;
        end              
        tt=linspace(0,2*pi,100);
        k1 = data.LatticeK(imgnum).k1;
        k2 = data.LatticeK(imgnum).k2;
        s1 = data.LatticeK(imgnum).s1;
        s2 = data.LatticeK(imgnum).s2;           
        k12p = k1+k2;
        k12n = k1-k2;
        s12 = sqrt(s1*s2);             
        set(pKReticles(1),'XData',k1(1)+4*s1*cos(tt),...
            'YData',k1(2)+4*s1*sin(tt),'color','r');
        set(pKReticles(2),'XData',-k1(1)+4*s1*cos(tt),...
            'YData',-k1(2)+4*s1*sin(tt),'color','r');
        set(pKReticles(3),'XData',k2(1)+4*s1*cos(tt),...
            'YData',k2(2)+4*s2*sin(tt),'color','b');
        set(pKReticles(4),'XData',-k2(1)+4*s1*cos(tt),...
            'YData',-k2(2)+4*s2*sin(tt),'color','b');   
        set(pKReticles(5),'XData',k12p(1)+4*s12*cos(tt),...
            'YData',k12p(2)+4*s12*sin(tt),'color','g');
        set(pKReticles(6),'XData',-k12p(1)+4*s12*cos(tt),...
            'YData',-k12p(2)+4*s12*sin(tt),'color','g');   
        set(pKReticles(7),'XData',k12n(1)+4*s12*cos(tt),...
            'YData',k12n(2)+4*s12*sin(tt),'color',[255, 215, 0]/255);
        set(pKReticles(8),'XData',-k12n(1)+4*s12*cos(tt),...
            'YData',-k12n(2)+4*s12*sin(tt),'color',[255, 215, 0]/255);                   
    end

    function updatePSFKGraphic
        s=tblPSF.Data(1);
        kR = sqrt(1/(4*pi*s^2));
        tt=linspace(0,2*pi,100);        
        set(pKPSF,'XData',kR*cos(tt),...
            'YData',kR*sin(tt),'color',[0 0.4470 0.7410]);   
    end

    function PSFKCB(src,evt)
        if src.Value
            pKPSF.Visible='on';
        else
            pKPSF.Visible='off';
        end 
    end
%% Binned Callbacks
    function updateBinnedGraphics
        if ~isfield(data,'LatticeBin')
            return;
        end         
        imgnum = menuSelectImg.Value;        
        set(hImg_B,'XData',data.LatticeBin(imgnum).n1,...
            'YData',data.LatticeBin(imgnum).n2,...
            'CData',data.LatticeBin(imgnum).Zbin);  
        



% RL = [data.LatticeBin(imgnum).n1(1) data.LatticeBin(imgnum).n1(end) ...
%            data.LatticeBin(imgnum).n2(1) data.LatticeBin(imgnum).n2(end)];
% 
% 
%         ROI=tblROIB.Data;     
% 
%         ROI=round(ROI);      % Make sure this ROI are integers   
%         % Check that limits go from low to high
%         % if ROI(2)<=ROI(1) || ROI(4)<=ROI(3)
%         %    % warning('Bad ROI specification given.');
%         %    % ROI(evt.Indices(2))=evt.PreviousData;
%         % end               
%         % Check that ROI is within image bounds
%         if ROI(1)<RL(1); ROI(1)=RL(1); end       
%         if ROI(3)<RL(3); ROI(3)=RL(3); end   
%         if ROI(4)>RL(4); ROI(4)=RL(4); end       
%         if ROI(2)>RL(2); ROI(2)=RL(2); end      
% 
%         tblROIB.Data = ROI;
%         try
%             pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];
%             set(pROIB(1),'Position',pos);
%         end
        if cAutoColor_B.Value;setClim('B');end     
        updateGridGraphics;
        latticeGridCB(cDrawLattice);
        latticeTextCB(cTextLattice);   
    end

%% Digital Callbacks
    function updateCoM_D
        if ~isfield(data,'LatticeDig') 
           tD_SE.Visible='off';
           return
        end
        imgnum = menuSelectImg.Value;
        bc = data.LatticeDig(imgnum);  
        
        % Update box count string
        str=[ num2str(bc.Natoms) ' atoms' newline ...
            '$(X_\mathrm{c},Y_\mathrm{c}) = ' '('  num2str(round(bc.Xc_site,1)) ',' ...
            num2str(round(bc.Yc_site,1)) ')$' newline ...
            '$(\sigma_X,\sigma_Y) = ' '('  num2str(round(bc.Xs_site,1)) ',' ...
            num2str(round(bc.Ys_site,1)) ')$']; 
        %Update box count string object
        set(tD_SE,'String',str);          
    end

%% Binned Histgoram Callbacks

    function updateBinnedHistogramGraphics   
         opts=struct;
        switch menuSelectBinType.Value
            case 1
                opts.BinSource = 'ZbinRaw';
            case 2
                opts.BinSource = 'Zbin';
            case 3
                opts.BinSource = 'ZbinNormalized';
        end

        opts.ImageNum= menuSelectImg.Value; % get image to analyze
        opts.Parent = tabHB;
        bin_showHistogram(data,opts);      
    end

%% 
    function updatePositionAnalysis
    hbposition.BackgroundColor=	[255 219 88]/255;
    drawnow;
        
        
        data.ROI=tblROI.Data;     


        tbl_pos_analysis.Data={};
        %% Box Analysis
        if hc_anlX_Box.Value                    
            data=ixon_boxCount(data);
        end
        % %% Sharpness Analysis
        % if hc_anlX_Sharpness.Value    
        %     data=ixon_Sharpness(data);
        % end
        %% Histogram Analysis
        if hc_anlX_Histogram.Value            
            data = ixon_PositionHistogram(data);
        end
        %%  Update Principal Component Analysis
        if hc_anlX_PCA.Value      
            % Finding cloud principal axes
            data = ixon_simple_pca(data);   
        end
    %% Guassian Analysis
        if hc_anlX_Gauss.Value
            opts=struct;
            opts.doRescale=1;
            opts.doMask=hcMask.Value;
            opts.Scale=0.5;
            opts.doRotate=0;
            opts.Mask=ixon_mask;        
            data = ixon_gaussFit(data,opts);            
        end    

    %% Stripe Analysis
    if (hcStripe.Value) && (~isfield(data.Flags,'plane_selection_dotilt') || ...
            (isfield(data.Flags,'plane_selection_dotilt') && data.Flags.plane_selection_dotilt))
        opts = struct;
        opts.Parent=tabStripe;
        opts.Name=data.Name;
        X = data.X;
        Y = data.Y;
        myz=sum(data.Z,3);    
        data.StripeCircular = ...
                StripeCircle(X,Y,myz,opts);
    end          

        updateDispPosImg;   
        updatePositionAnalysisGraphics
        updatePositionGraphics;
        
        
        hbposition.BackgroundColor=[80 200 120]/255;
        drawnow;
    end

    function saveGUIData(saveDir)       
        if doSaveGUIAnalysis
            try
            gui_saveData = struct;     
            gui_saveData.Date = data.Date;
            gui_saveData.FileName = data.Name;
            gui_saveData.Params = data.Params;
            gui_saveData.Flags = data.Flags;            
            
            if isfield(data,'BoxCount')
                BoxCount = data.BoxCount;
                gui_saveData.BoxCount = data.BoxCount; 
            end
            
            if isfield(data,'GaussFit')
               gui_saveData.GaussFit = data.GaussFit; 
            end
            
           if isfield(data,'StripeFit')
               gui_saveData.StripeFit = data.StripeFit; 
           end
           
            if isfield(data,'StripeFocus')
               gui_saveData.StripeFocus = data.StripeFocus; 
            end
            
            if isfield(data,'KFocusing')
               gui_saveData.KFocusing = data.KFocusing; 
            end
            
            if isfield(data,'BinStripe')
               gui_saveData.BinStripe = data.BinStripe; 
            end
            
            if isfield(data,'BinStripeCircular')                
               gui_saveData.BinStripeCircular = data.BinStripeCircular; 
            end

            if isfield(data,'StripeCircular')
                gui_saveData.StripeCircular = data.StripeCircular; 
            end
            
            if isfield(data,'MultiShotFocusing')
               gui_saveData.MultiShotFocusing = data.MultiShotFocusing; 
            end
           
            filenames=dir([saveDir filesep '*.mat']);
            filenames={filenames.name};
            filenames=sort(filenames);

            % Delete old images
            if length(filenames)>200
               f=[saveDir filesep filenames{1}];
               delete(f);
            end               

            filename=[data.Name '.mat']; 
            if ~exist(saveDir,'dir')                
               mkdir(saveDir);
            end         

            filename=fullfile(saveDir,filename);
            fprintf('%s',[filename ' ...']);
            save(filename,"-struct",'gui_saveData');


            disp(' done');      
            catch ME
                warning('Unable to save GUI data.')
            end



        end 
    end

    function updatePositionAnalysisGraphics
        tbl_pos_analysis.Data={};
        updateBoxGraphics;
        % updateSharpnessGraphics;
        updatePositionHistogramGraphics;
        updatePCAGraphics;    
        updateGaussGraphics; 
        % updateStripeGraphics;
    end

%% New Data
% What to do when new data is put into the GUI
    function newDataCallback
        try
        % Grab the ROI
        ROI=tblROI.Data;data.ROI=ROI;       
    
        %% Image Processing
        % Grab the RawImages and process them into usable data
        opt = struct;    
        opt.doSubtractBias     = tbl_process_1.Data{1,1};%hcSubBias.Value;
        opt.doSubtractBG       = tbl_process_1.Data{2,1};%hcSubBG.Value;
        opt.doMask             = tbl_process_1.Data{3,1};%hcMask.Value;
        opt.Mask               = ixon_mask;
        
        opt.doPSF              = tbl_process_1.Data{4,1};%hcPSF.Value ;    
        opt.PSF                = tblPSF.Data;    
        opt.DetectNoise        = tbl_process_1.Data{5,1};%hcPSF_noise.Value;
        opt.Noise              = tblPSF.Data(4);
        
        opt.doGaussFilter      = tbl_process_2.Data{1,1};%cGaussFilter.Value;
        opt.GaussFilterRadius  = tbl_process_2.Data{1,3};%tblGaussFilter.Data;
        
        opt.doScale      = tbl_process_2.Data{2,1};%cScale.Value;
        opt.ScaleFactor  = tbl_process_2.Data{2,3};%tblScale.Data;   
    
        opt.doRotate           = tbl_process_2.Data{3,1};%cRotate.Value;
        opt.Theta              = tbl_process_2.Data{3,3};%tblTheta.Data;
        
        % Momentum Space
        opt.doFFT              = 1;%hcFFT.Value;
        opt.doFFTFilter        = tbl_process_2.Data{4,1};%cKGaussFilter.Value;
        opt.FFTFilterRadius    = tbl_process_2.Data{4,3};%tblKGaussFilter.Data;         
        opt.doMaskIR           = tbl_process_2.Data{5,1};%hcIRMask.Value;
        opt.IRMaskRadius       = tbl_process_2.Data{5,3};%tblIRMask.Data;        
    
        
        % Process the Images
        hbprocess.BackgroundColor=	[255 219 88]/255;
        drawnow;
        data = ixon_ProcessImages(data,opt);  
        
        tblPSF.Data(4) =  round(data.NoiseEstimation(1));
        hbprocess.BackgroundColor=[80 200 120]/255;
        
        pAtoms.Visible='off';
        set(pAtoms,'XData',[],'YData',[]);
    
        % Update history index
        updateHistoryInd(data); 
    
        % Update the record of the process images
        updateImageDataLists;    
    
        updateDispPosImg;
    %% Update Parameter Data
        % Update table parameters (alphebetically)
        [~,inds] = sort(lower(fieldnames(data.Params)));
        params = orderfields(data.Params,inds);    
                    
        tbl_params.Data={};
    
        fnames=fieldnames(params);
        for nn=1:length(fnames)
          tbl_params.Data{nn,1}=fnames{nn};
            val=data.Params.(fnames{nn});
            if isa(val,'double')
                tbl_params.Data{nn,2}=num2str(val);
            end
    
            if isa(val,'struct')
               tbl_params.Data{nn,2}='[struct]'; 
            end  
        end        
    %% Update Flag Data
        tbl_flags.Data={};
        if isfield(data,'Flags')
            flags = data.Flags;
            fnames=fieldnames(flags);
            for nn=1:length(fnames)
              tbl_flags.Data{nn,1}=fnames{nn};
                val=flags.(fnames{nn});
                if isa(val,'double')
                    tbl_flags.Data{nn,2}=num2str(val);
                end
                if isa(val,'struct')
                   tbl_flags.Data{nn,2}='[struct]'; 
                end  
            end     
        end
        %% Position Space Analysis
        updatePositionAnalysis;
        
        %% Momentum Space Analysis
        if hc_anlK_auto.Value        
           analyze_k       
        end
        
        if hc_anlB_auto.Value
            analyze_bin
        end
        
        if hc_anlD_auto.Value
            analyze_dig
        end
      
    %% Update Fit Results (depreciated)
    
        
        drawnow;
        climtbl_X.Data=axImg.CLim;
    
        %% Save Analysis
    
        saveGUIData(GUIAnalysisLocalDir);
        saveGUIData(GUIAnalysisSaveDir);

    catch ME
        warning(getReport(ME,'extended','hyperlinks','on'));
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% updateDataPlots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function graphically updates the x and y data plots

function updateDataPlots(data)
    % Grab the ROI
    ROI=tblROI.Data;
    data.ROI=ROI;       
    x=ROI(1):ROI(2);
    y=ROI(3):ROI(4); 

    % Image over the ROI domain
    Zsub=data.Z(y,x);

    % Update plots if sum
    if rbSum_X.Value
        Zy=sum(Zsub,2);
        Zx=sum(Zsub,1);          
        set(pX,'XData',x,'YData',Zx);
        set(pY,'XData',Zy,'YData',y);
        drawnow;
    end
    
    % Update plots if cut
    if rbCut_X.Value
        indy=find(round(pCrossX.YData(1))==y,1);           % Y center
        indx=find(round(pCrossY.XData(1))==x,1);           % X center 
        Zy=data.Z(:,indx);        
        Zx=data.Z(indy,:);        
        set(pX,'XData',data.X,'YData',Zx);
        set(pY,'XData',Zy,'YData',data.Y);
        drawnow;
    end   
end

  % 
  % 
  % function updateStripeGraphics
  %        if ~isfield(data,'StripeFit') 
  %           return;
  %       end
  %       imgnum = menuSelectImg.Value;   
  % 
  %       x = data.X;
  %       y = data.Y;
  %       z = data.ZNoFilter(:,:,imgnum);
  %       stripe = data.StripeFit(imgnum);
  %       focus = data.StripeFocus(imgnum);
  % 
  %       [xx,yy] = meshgrid(x,y);        
  %       Zfit=feval(stripe.Fit,xx,yy);        
  %       theta=stripe.theta;
  %       L = stripe.L;
  %       phi = stripe.phi;
  %       xC = stripe.xC;
  %       yC = stripe.yC;
  %       s1 = stripe.s1;
  %       s2 = stripe.s2;        
  %       set(stripe_str_bottom_left,'String',[data.Name ' (' num2str(menuSelectImg.Value) ')']);
  % 
  %       str=['$(X_\mathrm{c},Y_\mathrm{c}) = ' '('  num2str(round(xC,1)) ',' ...
  %           num2str(round(yC,1)) ')$' newline ...
  %           '$(\sigma_1,\sigma_2) = ' '('  num2str(round(s1,1)) ',' ...
  %           num2str(round(s2,1)) ')$' newline ...
  %           '$(\theta,\phi,\lambda) = (' num2str(round(theta,2)) '^\circ,' ...
  %           num2str(round(phi/(2*pi),2)) '\cdot 2\pi,' ...
  %           num2str(round(L,1)) ')$']; 
  %       set(stripe_str_bottom_right,'String',str);
  %       % Update Image
  %       set(stripe_hImgStripe,'XData',x,'YData',y,'CData',z);       
  % 
  %       % Update Image markers
  %       set(stripe_pFringe,'XData',xC+[-2 2]*s1*cosd(theta),...
  %           'Ydata',yC+[-2 2]*s1*sind(theta));
  %       set(stripe_pPerp,'XData',xC+[-2 2]*s2*cosd(theta+90),...
  %           'Ydata',yC+[-2 2]*s2*sind(theta+90));
  %        set(stripe_pBar,'XData',xC+[50 0],...
  %           'Ydata',yC*[1 1]);
  % 
  %       tt=linspace(0,2*pi,100);
  %       xell = 2*s1*cosd(theta)*cos(tt)-2*s2*sind(theta)*sin(tt)+xC;
  %       yell = 2*s2*cosd(theta)*sin(tt)+2*s1*sind(theta)*cos(tt)+yC;        
  %       set(stripe_pCloudEllipse,'XData',xell,'YData',yell);
  %       tt=linspace(0,theta,100);
  %       set(stripe_pAngleCirc,'XData',xC+25*cosd(tt),...
  %           'YData',yC+25*sind(tt));
  %       ax_stripe_img.CLim = axImg.CLim;
  %       set(ax_stripe_img,'XLim',[1 512],'YLim',[1 512]);
  % 
  %       % Update fits   
  %       set(stripe_pSum1_fit,'XData',x,'YData',sum(imrotate(Zfit,theta,'crop'),1));
  %       set(stripe_pSum1_data,'XData',x,'YData',sum(imrotate(z,theta,'crop'),1));
  %       % Show the sum counts orthogonal to the stripe axis
  %       set(stripe_pSum2_fit,'XData',y,'YData',sum(imrotate(Zfit,theta,'crop'),2));
  %       set(stripe_pSum2_data,'XData',y,'YData',sum(imrotate(z,theta,'crop'),2));      
  % 
  %         % stripe=ixon_fitStripe(data,opts);   
  %       stranl={'','';
  %           ['stripe A (amp)'] ,stripe.A;
  %           ['stripe xC (px)'],stripe.xC;
  %           ['stripe yC (px)'],stripe.yC;
  %           ['stripe ' char(963) '1 (px)'],stripe.s1;
  %           ['stripe ' char(963) '2 (px)'],stripe.s2;
  %           ['stripe B'],stripe.B;
  %           ['stripe ' char(952) ' (deg)'],stripe.theta;
  %           ['stripe ' char(955) ' (px)'],stripe.L;
  %           ['stripe ' char(966) ' (2pi)'],stripe.phi/(2*pi);};  
  %       tbl_pos_analysis.Data=[tbl_pos_analysis.Data; stranl];     
  % 
  % 
  %       set(stripe_pFocus,'Xdata',focus.y_coms,'YData',focus.scores);        
  %       tt=linspace(min([focus.y_coms]),max([focus.y_coms]),100);        
  %       set(stripe_pFocusFit,'XData',tt,'YData',polyval(focus.poly,tt))        
  %       axes(ax_stripe_focus);
  %       yyaxis left
  %       set(ax_stripe_focus,'YLim',[0 max(focus.scores)*1.1]);
  %       set(stripe_pFocus2,'Xdata',focus.y_coms,'YData',focus.sums);
  % 
  % 
  %   end
  % 


    function updateGaussGraphics
        if ~isfield(data,'GaussFit') 
            return;
        end
        imgnum = menuSelectImg.Value;   
        % Grab fit data
        fout=data.GaussFit{imgnum}; 
        A = fout.A;
        Xc = fout.Xc;
        Yc = fout.Yc;
        nbg = fout.nbg;
        if ismember('theta',coeffnames(fout))
            theta = fout.theta;
            s1 = fout.s1;
            s2 = fout.s2;
        else
            s1 = fout.Xs;
            s2 = fout.Ys;
            theta = 0;
        end
        % Evaluate and plot 1/e^2 gaussian reticle
        t=linspace(0,2*pi,100);          
        if ismember('theta',coeffnames(fout))
            xR=Xc+s1*cos(theta)*cos(t)-s2*sin(theta)*sin(t);
            yR=Yc+s1*sin(theta)*cos(t)+s2*cos(theta)*sin(t); 
        else
            xR=Xc+s1*cos(t);
            yR=Yc+s2*sin(t);    
        end   
        % Might need to differentiate between graphics and plot
        set(pGaussRet,'XData',xR,'YData',yR,'linewidth',2);  
        updateGaussLinePlot;      
        drawnow;
        % Gaussian analysis table string
        stranl={'','';
            ['gauss N (counts)'] ,2*pi*A*s1*s2;
            ['gauss A (counts)'],A;
            ['gauss s1' char(963) ' (px)'],s1;
            ['gauss s2' char(963) ' (px)'],s2;
            ['gauss xc (px)'],Xc;
            ['gauss yc (px)'],Yc;
            ['gauss nbg (counts)'],nbg;
            ['gauss theta'],theta;};    
        tbl_pos_analysis.Data=[tbl_pos_analysis.Data; stranl];  
    end

    function updateGaussLinePlot
        if ~isfield(data,'GaussFit') 
            return;
        end
        imgnum = menuSelectImg.Value;         
        
        % Grab fit data
        fout=data.GaussFit{imgnum};

        [xx,yy]=meshgrid(data.X,data.Y);
        Z=feval(fout,xx,yy); 

        % Find indeces in which correspond to ROI boundary
        ROI=tbl_dROI_X.Data;
        [~,c1] = min(abs(data.X-ROI(1)));
        [~,c2] = min(abs(data.X-ROI(2)));
        [~,r1] = min(abs(data.Y-ROI(3)));
        [~,r2] = min(abs(data.Y-ROI(4)));

        zsub = Z(r1:r2,c1:c2);    
        xsub = data.X(c1:c2);
        ysub = data.Y(r1:r2);

        if rbCut_X.Value
            [~,iy] = min(abs(pCrossX.YData(1)-data.Y));
            [~,ix] = min(abs(pCrossY.XData(1)-data.X));
            ZyF=Z(:,ix);ZxF=Z(iy,:);
            set(pXF,'XData',data.X,'YData',ZxF,'Visible','on');
            set(pYF,'XData',ZyF,'YData',data.Y,'Visible','on');
        else       
            % ZyF=sum(zF(:,c1:c2),2);ZxF=sum(zF(r1:r2,:),1);   
            % set(pXF,'XData',data.X,'YData',ZxF,'Visible','on');
            % set(pYF,'XData',ZyF,'YData',data.Y,'Visible','on');
            set(pXF,'XData',xsub,'YData',sum(zsub,1),'Visible','on');
            set(pYF,'YData',ysub,'XData',sum(zsub,2),'Visible','on');
        end  
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% updateAnalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs analysis and updates graphics as required.
function data=updateAnalysis(data)    
    
    % Update PCA analysis
    if hc_anlX_PCA.Value      
        % Finding cloud principal axes
        data=ixon_simple_pca(data);        
        out=data.PCA;        
        x1=out.Mean(1)+out.Radii(1)*out.PCA(1,1)*[-1 1];
        y1=out.Mean(2)+out.Radii(1)*out.PCA(2,1)*[-1 1];
        x2=out.Mean(1)+out.Radii(2)*out.PCA(1,2)*[-1 1];
        y2=out.Mean(2)+out.Radii(2)*out.PCA(2,2)*[-1 1];        
        % PCA analysis table string
        stranl={'','';
            ['pca ' char(952) '1 (deg)'] ,atan(out.PCA(2,1)/out.PCA(1,1))*180/pi;
            ['pca ' char(952) '2 (deg)'],atan(out.PCA(2,2)/out.PCA(1,2))*180/pi;
            ['pca ' char(963) '1 (px)'],out.Radii(1);
            ['pca ' char(963) '2 (px)'],out.Radii(2);
            ['pca xc (px)'],out.Mean(1);
            ['pca yc (px)'],out.Mean(2);};
        tbl_pos_analysis.Data=[tbl_pos_analysis.Data; stranl];        
        set(pPCA(1),'XData',x1,'YData',y1,'Visible','on');
        set(pPCA(2),'XData',x2,'YData',y2,'Visible','on');
    end
    
    % Update Guassian Analysis
    if hc_anlX_Gauss.Value
        disp('Fitting data to 2D gaussian...')   
        opts=struct;
        opts.doRescale=1;
        opts.doMask=hcMask.Value;
        opts.Scale=0.25;
        opts.doRotate=0;
        opts.Mask=ixon_mask;        
        data=ixon_gaussFit(data,opts);   
        cGaussRet.Enable='on';        
        % Gaussian analysis table string
        stranl={'','';
            ['gauss N (counts)'] ,2*pi*data.GaussFit{1}.A*data.GaussFit{1}.s1*data.GaussFit{1}.s2;
            ['gauss A (counts)'],data.GaussFit{1}.A;
            ['gauss x' char(963) ' (px)'],data.GaussFit{1}.s1;
            ['gauss y' char(963) ' (px)'],data.GaussFit{1}.s2;
            ['gauss xc (px)'],data.GaussFit{1}.Xc;
            ['gauss yc (px)'],data.GaussFit{1}.Yc;
            ['gauss nbg (counts)'],data.GaussFit{1}.nbg;};
        tbl_pos_analysis.Data=[tbl_pos_analysis.Data; stranl];  
        % updateGaussPlot(data);
        updateGaussGraphics;
    end    
    
    % Update Guassian Analysis
    if hcStripe.Value        
        opts.Theta = [10 190];
        stripe=ixon_fitStripe(data,opts);   
        stranl={'','';
            ['stripe A (amp)'] ,stripe.A;
            ['stripe xC (px)'],stripe.xC;
            ['stripe yC (px)'],stripe.yC;
            ['stripe ' char(963) '1 (px)'],stripe.s1;
            ['stripe ' char(963) '2 (px)'],stripe.s2;
            ['stripe B'],stripe.B;
            ['stripe ' char(952) ' (deg)'],stripe.theta;
            ['stripe ' char(955) ' (px)'],stripe.L;
            ['stripe ' char(966) ' (2pi)'],stripe.phi/(2*pi);};  
        tbl_pos_analysis.Data=[tbl_pos_analysis.Data; stranl];  
        % updateStripePlot(data,stripe);
        % focus = ixon_focusStripe(data,stripe);
        % updateFocusPlot(data,focus);
    end  

    % Update Guassian Analysis
    if hc_anlX_GaussRot.Value
        disp('Fitting data to 2D gaussian...')   
        opts=struct;
        opts.doRescale=1;
        opts.doMask=hcMask.Value;
        opts.Scale=0.5;
        opts.doRotate=1;
        opts.Mask=ixon_mask; 
        data=ixon_gaussFit(data,opts);   
        cGaussRet.Enable='on';        
        stranl={'','';
            ['gauss N (counts)'] ,2*pi*data.GaussFit{1}.A*data.GaussFit{1}.s1*data.GaussFit{1}.s2;
            ['gauss A (counts)'],data.GaussFit{1}.A;
            ['gauss ' char(963) '1 (px)'],data.GaussFit{1}.s1;
            ['gauss ' char(963) '2 (px)'],data.GaussFit{1}.s2;
            ['gauss xc (px)'],data.GaussFit{1}.Xc;
            ['gauss yc (px)'],data.GaussFit{1}.Yc;
            ['gauss nbg (counts)'],data.GaussFit{1}.nbg;
            ['gauss ' char(952) ' (deg)'],data.GaussFit{1}.theta*180/pi};
        tbl_pos_analysis.Data=[tbl_pos_analysis.Data; stranl]; 
        % updateGaussPlot(data);
        updateGaussGraphics
    end     
    
    %% Fit Results    
    % fr=tbl_pos_analysis.Data(:,2)';    
    % % Ensure fit results is a number
    % for n=1:length(fr)
    %     % If value is empty, assign a zero
    %     if isempty(fr{n})
    %         fr{n}=0;
    %     end
    % 
    %     % If string conver to number
    %     if isstr(fr{n})
    %         try
    %             fr{n}=str2double(fr{n});
    %         catch ME
    %             fr{n}=NaN;
    %         end
    %     end
    % end
    
    % Get fit results variable
    % frVar=frslct.String{frslct.Value};    
    % val=data.Params.(frVar);
    % 
    % % Convert execution date into a time
    % if isequal(frVar,'ExecutionDate') 
    %    val=datenum(val); 
    %    val=val-floor(val);
    %    val=val*24*60;
    % end
    % 
    % % Create fit results object
    % fr=[data.Name frVar val fr];   
    
    %%%%%% Output to fit results
    % Output some analysis to the main workspace, this is done to be
    % comptaible with old regimens for fitting and analysis
    
    % try
    %     % Read in fitresults
    %     ixon_fitresults=evalin('base','ixon_fitresults');        
    % catch ME
    %     % Error means that it is probably undefined
    %     ixon_fitresults={};
    % end
    % 
    % M=size(ixon_fitresults,1)+1;                         % Find next row
    % ixon_fitresults(M:(M+size(fr,1)-1),1:size(fr,2))=fr; % Append data        
    % assignin('base','ixon_fitresults',ixon_fitresults);  % Rewrite fitresults
    % 
    
    % if doSaveGUIAnalysis
    %     gui_saveData = struct;
    %     gui_saveData.Date = data.Date;
    %     gui_saveData.FileName = data.Name; 
    %     gui_saveData.Params = data.Params;
    %     gui_saveData.BoxCount = data.BoxCount;    
    %     if hcStripe.Value  
    %         gui_saveData.Stripe = stripe;
    %     end
    % 
    %    filenames=dir([GUIAnalysisSaveDir filesep '*.mat']);
    %    filenames={filenames.name};
    %    filenames=sort(filenames);
    % 
    %    % Delete old images
    %    if length(filenames)>200
    %        f=[GUIAnalysisSaveDir filesep filenames{1}];
    %        delete(f);
    %    end               
    % 
    %     filename=[data.Name '.mat']; 
    %     if ~exist(GUIAnalysisSaveDir,'dir')
    %        mkdir(GUIAnalysisSaveDir);
    %     end        
    %     filename=fullfile(GUIAnalysisSaveDir,filename);
    %     fprintf('%s',[filename ' ...']);
    %     save(filename,'gui_saveData');
    %     disp(' done');    
    % end
end

%% OTHER HELPER FUNCTIONS

    function mydata=makeImgDataStruct(imgs)
        mydata=struct;
        % Create the image data structure
        mydata=struct;
        mydata.Date=datevec(now);
        mydata.Name=['iXonUltra_' datestr(mydata.Date,'yyyy-mm-dd_HH-MM-SS')];   

        % Grab the images
        mydata.RawImages=imgs;

        % Add magnification
        mydata.Magnification=mag;        
        
        % Grab the sequence parameters and flags   
        [mydata.Params,mydata.Units,mydata.Flags]=ixon_grabSequenceParams;      
        
        imgs = mydata.RawImages;  
        
        % Remove the first image if there is a buffer
        if isfield(mydata.Flags,'qgm_doClearBufferExposure') 
            if mydata.Flags.qgm_doClearBufferEposure && size(imgs,3)>1
                iD = imgs(:,:,1);
                imgs(:,:,1)=[];    
            end
        else
            if size(imgs,3)>1
                iD = imgs(:,:,1);
                imgs(:,:,1)=[];   
            end
        end       

        % Append acquisition information
        mydata.CameraInformation=cam_info;
        mydata.AcquisitionInformation=acq;
        mydata.AcquisitionDescription=desc;        
    end


    function saveData(data,saveDir)
        if nargin==1
           saveDir=defaultDir;
           filenames=dir([saveDir filesep '*.mat']);
           filenames={filenames.name};
           filenames=sort(filenames);

           % Delete old images
           if length(filenames)>200
               f=[saveDir filesep filenames{1}];
               delete(f);
           end               
        end
        filename=[data.Name '.mat']; 
        if ~exist(saveDir,'dir')
           mkdir(saveDir);
        end        
        filename=fullfile(saveDir,filename);
        fprintf('%s',[filename ' ...']);
        save(filename,'data');
        disp(' done'); 
    end

    function updateHistoryInd(data)       
        % Upate history list
        filenames=dir([currDir  filesep '*.mat']);
        filenames={filenames.name};       
        filenames=sort(filenames);
        filenames=flip(filenames);    

        % Find image in image history
        myname=[data.Name '.mat'];               
        ind=find(ismember(filenames,myname));    % index in filenames        
        if isempty(ind)
          ind=1; 
        end

        % Update string
        % tNavInd.String=sprintf('%03d',ind); 

        tNavMax.String=['of ' num2str(length(filenames))];        

        tNavName.String=fullfile(currDir,[data.Name '.mat']);        

        tNavInd.Data(1)=ind;
    end
    
    %% FINISH
newDataCallback;
% Go to most recent image
chData([],[],1);   
drawnow;
SizeChangedFcn
axes(axImg);
updateCMAP(menuSelectCMAP);
addlistener(axImg,'XLim','PostSet',@foo); 
addlistener(axImg,'YLim','PostSet',@foo); 
addlistener(axImg_K,'XLim','PostSet',@foo3); 
addlistener(axImg_K,'YLim','PostSet',@foo3); 
addlistener(axImg_B,'XLim','PostSet',@foo4); 
addlistener(axImg_B,'YLim','PostSet',@foo4); 
addlistener(axImg_D,'XLim','PostSet',@foo5); 
addlistener(axImg_D,'YLim','PostSet',@foo5); 

% set(hF,'WindowState','maximized');
    function foo5(~,~)
        tbl_dROI_D.Data = round([axImg_D.XLim axImg_D.YLim]);
    end

    function foo4(~,~)
        tbl_dROI_B.Data = round([axImg_B.XLim axImg_B.YLim]);
    end

    function foo3(~,~)
        set(pCrossX_K,'XData',axImg_K.XLim,'YData',[1 1]*mean(axImg_K.YLim));
        set(pCrossY_K,'YData',axImg_K.YLim,'XData',[1 1]*mean(axImg_K.XLim));
        
        tbl_dROI_K.Data = [axImg_K.XLim axImg_K.YLim];

        if rbCut_K.Value
            ind1 = find([data.f>pCrossX_K.YData(1)],1);
            ind2 = find([data.f>pCrossY_K.XData(1)],1);       
            set(pX_K,'XData',data.f,'YData',hImg_K.CData(ind1,:));
            set(pY_K,'YData',data.f,'XData',hImg_K.CData(:,ind2));      
        end     
    end

    function foo(~,~)         
        imgnum = menuSelectImg.Value;
        switch menuSelectImgType.Value
            case 1
                Z = data.Z(:,:,imgnum);
            case 2
                Z = data.ZNoFilter(:,:,imgnum);
        end
    
        if cAutoColor_X.Value;setClim('X');end 
        % Find center of update
        xC = mean(axImg.XLim);yC = mean(axImg.YLim);

        % Update crosshair
        set(pCrossX,'XData',axImg.XLim,'YData',[1 1]*yC);
        set(pCrossY,'YData',axImg.YLim,'XData',[1 1]*xC); 

        % Round the table limits
        tbl_dROI_X.Data = round([axImg.XLim axImg.YLim]); 

        % Get the region of interest
        ROI =  [axImg.XLim axImg.YLim];

        % Find indeces in which correspond to ROI boundary
        [~,c1] = min(abs(data.X-ROI(1)));
        [~,c2] = min(abs(data.X-ROI(2)));
        [~,r1] = min(abs(data.Y-ROI(3)));
        [~,r2] = min(abs(data.Y-ROI(4)));

        % Find indeces corresponding to center of displayed image
        [~,iC] = min(abs(data.X-xC));
        [~,iR] = min(abs(data.Y-yC));

        % Update plots if cut
        if rbCut_X.Value
            set(pX,'XData',data.X,'YData',Z(iR,:));
            set(pY,'YData',data.Y,'XData',Z(:,iC));
        end     

        if rbSum_X.Value    
            zsub = Z(r1:r2,c1:c2);    
            xsub = data.X(c1:c2);
            ysub = data.Y(r1:r2);
            set(pX,'XData',xsub,'YData',sum(zsub,1));
            set(pY,'YData',ysub,'XData',sum(zsub,2));
        end 

        if hc_anlX_Gauss.Value
            updateGaussLinePlot;
        end

    end

axes(axImg);
set(axImg,'XLim',[1 512],'YLim',[ 1 512]); 

enableInteractivity;

function enableInteractivity
    enableDefaultInteractivity(axImg);
    enableDefaultInteractivity(axImg_K);
    enableDefaultInteractivity(axImg_B);
    enableDefaultInteractivity(axImg_D);
    enableDefaultInteractivity(axImg);
end

end



function savePosAnalysis(tbl,execdate,src)
    if nargin == 1
        src = 'X\boop.txt';
    end

%     src = 'C:\Users\coraf\OneDrive\Desktop\bob.csv';
    src = 'Y:\ixon_gui_analysis.csv';

    vals = [execdate; tbl(:,2)]';
    names = ['ExecutionDate'; tbl(:,1)];

    vals = vals(~cellfun('isempty',vals));
    names = names(~cellfun('isempty',names));

    T = cell2table(vals,'VariableNames',names);

    if exist(src,'file')
        try
        writetable(T,src);
        end
    end
end
    
% Connect to the Andor iXon Ultra
function out=ixon_connectCam
    out =0;
    disp('Connecting camera');
    
    % Find available cameras
    [ret,ncams]=GetAvailableCameras();
    disp(['Detected ' num2str(ncams) ' valid cameras.']);
    
    if ncams~=1
        warning('Invalid number of cameras detected.');
        return;
    end 
    
    % Initialize the camera and load DLLs
    currDir = pwd;
    fileDir = fileparts(mfilename('fullpath'));
    cd(fileDir);
    fprintf('Connecting to Andor camera ... ');
    [ret] = AndorInitialize('');
    disp(error_code(ret))
    cd(currDir);
    
    % Return status
    if isequal(error_code(ret),'DRV_SUCCESS')
        out=1;
    else
        warning('Unable to initize the Andor iXon Ultra.');
        out=0;
    end
end
