function ixon_gui
% ixon_gui.m
%
% Author      : C. Fujiwara
% Last Edited : 2021/05/06
%
% This code operates the iXon Ultra camera that the lattice experiment
% uses to take fluorescence images of the quantum gas microscope.
%
% This rewrite is primary meant to improve the user experience running the
% camera and to simplify the code.  Furthermore, this code does not use the
% GUIDE interface, which is a GUI designer from MATLAB.  The GUIDE protocl
% has depreciated as of 2021.  This code does not utilize the uifigure
% interface from MATLAB, but uses the figure interface for desinging GUIs.
% This design choice is because CF is not as familiar with the uifigure
% interface.

% Enable debug mode?
doDebug=1;

% Manification (as measured in Graham's thesis)
mag=[82.6 83.2];

%% Other Settings

% Choose the default colormap
cmap=purplemap;

% Figure Name
guiname='iXon GUI';

% Directory where images are automatically saved.
historyDir=['C:' filesep 'IxonImageHistory'];

if ~exist(historyDir,'dir')
    mkdir(historyDir);
end

% % Dummy file to load on startup
% fname='example_data_plane.mat';
% initdata=load(fname);
% Z=initdata.Z;
% % Initializse image data structure with dummy data
% data=struct;
% data.X=1:size(initdata.Z,2);
% data.Y=1:size(initdata.Z,1);
% data.Z=initdata.Z;
% data.Date=datevec(now);
% data.Name=['iXonUltra_' datestr(data.Date,'yyyy-mm-dd_HH-MM-SS')];


% Dummy file to load on startup
fname='example_data_EIT_RAMAN.mat';
data=load(fname);
data=data.data;
data.Z=data.RawImages(:,:,2)-data.RawImages(:,:,1);
Z=data.Z;
% Initializse image data structure with dummy data


%% Initialize Drivers and GUI

% Add the Andor MATLAB drivers to the MATLAB path. You need to have
% installed the MATLAB drivers in order for this to work.
addpath(genpath(fullfile(matlabroot,'toolbox','Andor')));

% Add all subdirectories for this m file
mpath = fileparts(mfilename('fullpath'));
addpath(mpath);addpath(genpath(mpath))

% Find any instances of the GUI and bring it to focus, this is to avoid
% restarting the GUI which may leave the shutter open.
h = findall(0,'tag','GUI');
for kk=1:length(h)
    if isequal(h(kk).Name,guiname)
        disp(['iXon GUI instance detected.  Bringing into focus. ' ...
            ' If you want to start a new instance, close the original iXon GUI.']); 
       figure(h);
       return;
    end    
end

disp('Initializing iXon GUI...');

%% Camera Settings

% Declare cam_info struct;
cam_info=struct;

% Initialize Camera Status
cam_status=struct;
cam_status.isConnected=0;       % Are you connected to camera?
cam_status.Temperature=NaN;     % sensor temperature
cam_status.TemperatureSP=10;    % Must lie between -120 C and 20 C
cam_status.isTempStable=0;      % is sensor tempreature stable?
cam_status.isCooling=0;         % is the TEC active?
cam_status.isAcquiring=0;       % is the camera acquiring frames?

% Load default acquisition settings
acq=defaultNormalAcqSettings;

% Description of Acquisitino settings
desc=acqDescription(acq);


%% Initialize Figure

% Initialize the primary figure
hF=figure;clf
set(hF,'Color','w','units','pixels','Name',guiname,'toolbar','none',...
    'Tag','GUI','CloseRequestFcn',@closeGUI,'NumberTitle','off',...
    'Position',[50 50 1200 850],'SizeChangedFcn',@SizeChangedFcn);

% Callback for when the GUI is requested to be closed.
    function closeGUI(fig,~)
        answer = questdlg('Close the iXon GUI?','Close iXon?',...
            'Yes','No','No') ;
        
         doClose=0;
        % Handle response
        switch answer
            case 'Yes'
                disp('Closing the iXon GUI ...')
                doClose = 1;
            case 'No'
                disp('Not closing the iXon GUI.')
                doClose=0;
        end
        
        if doClose
            disp('Closing iXon GUI...');

            % Stop acquistion if necessary

            % Stop Cooling if necessary

            % Actually, maybe don't allow closing of GUI. Or have a confirm
            % dialog option


            stop(statusTimer);
            if cam_status.isConnected
                disconnectCam;
            end
            delete(statusTimer);

            delete(fig);                % Delete the figure
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
        
        % Grab figure dimensions
        W=hF.Position(3);H=hF.Position(4);         
        
        % Resize image panel  
        if W>360 && H>55        
            hp.Position=[360 1 W-360 H-55];           
        end
        
        % Resize plots
        resizePlots;                            
        
        % Resize Panels
        hpCam.Position(2:3)=[H-hpCam.Position(4) hF.Position(3)];        
        hpSave.Position(2:3)=[hpCam.Position(2)-hpSave.Position(4) hF.Position(3)-300];
        hpNav.Position(1:2)=[hpSave.Position(1)+hpSave.Position(3) hpSave.Position(2)];
        hpAcq.Position(2)=hpSave.Position(2)-hpAcq.Position(4);
        hpADV.Position(2)=hpAcq.Position(2)-hpADV.Position(4);
        hpAnl.Position(2)=hpADV.Position(2)-hpAnl.Position(4);        
        hpDisp.Position(4)=max([hpAnl.Position(2) 1]);        
        hpFit.Position(4)=H-55;        
        
        % Reposition objects in hpDisp because it has variable height.
        tbl_dispROI.Position(2)=hpDisp.Position(4)-tbl_dispROI.Position(4)-20;
        hbFullLim.Position(2)=tbl_dispROI.Position(2)+2;
        hbSnapLim.Position(2)=tbl_dispROI.Position(2)+2;
        hbSlctLim.Position(2)=tbl_dispROI.Position(2)+2;
        climtbl.Position(2)=tbl_dispROI.Position(2)-climtbl.Position(4)-5;
        climtext.Position(2)=climtbl.Position(2);
        bgPlot.Position(2)=climtbl.Position(2)-bgPlot.Position(4)-2;
        cGaussRet.Position(2)=bgPlot.Position(2)-20;
        cCoMStr.Position(2)=cGaussRet.Position(2)-20;
        cCross.Position(2)=cCoMStr.Position(2)-20;
        cDrag.Position(2)=cCross.Position(2);
        tblcross.Position(2)=cCross.Position(2)-tblcross.Position(4);
        
        
        % Move status string
        strstatus.Position(1)=hpCam.Position(3)-strstatus.Position(3)-2;        
        drawnow;
end

%% Camera Panel
hpCam=uipanel(hF,'units','pixels','backgroundcolor','w',...
    'Position',[0 hF.Position(4)-30 hF.Position(3) 30]);

% Connect camera
ttstr='Connect to the iXon camera.';
hbConnect=uicontrol(hpCam,'style','pushbutton','string','connect','units','pixels',...
    'fontsize',10,'Position',[2 5 55 20],'backgroundcolor',[80 200 120]/255,...
    'Callback',@connectCB,'ToolTipString',ttstr);

% Callback for the connect button
    function connectCB(~,~)
        % Connect to the camera
       out=connectCam; 
       
       % Give warning if connection fails
       if ~out && ~doDebug
           warning('Unable to connect to camera');
           return;
       end
       
       cam_status.isConnected=1;
       
       % Close the shutter
       setCameraShutter(0);

       
       % Load default acquisition settings
       loadAcquisitionSettings;     
       
       % Close the shutter (again to be safe)
        setCameraShutter(0);        
        hbOpenShutter.Enable='on'; % allow shutter to be opened

        % Get the camera information
        cam_info=getCamInfo;
        disp(' ');
        disp(cam_info);
        
        % Set the temperature set point
        setTemperature(tblTemp.Data);
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
    'fontsize',10,'Position',[58 5 70 20],'backgroundcolor',[255 102 120]/255,...
    'Callback',@disconnectCB,'ToolTipString',ttstr,'enable','off');

% Callback for the disconnect button
    function disconnectCB(~,~)
        % Stop temperature monitor
        stop(statusTimer);
        
        % Disconnect from camera
        out=disconnectCam; 
       
        if ~out && ~doDebug
           return;
        end
          
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
cdata=imresize(imread(fullfile(mpath,'icons','info.jpg')),[18 18]);
hbCamInfo=uicontrol(hpCam,'style','pushbutton','CData',cdata,'callback',@infoCB,...
    'enable','on','backgroundcolor','w','position',[130 5 20 20],...
    'ToolTipString',ttstr,'enable','off');

    function infoCB(~,~)        
        cam_info=getCamInfo;        
        disp(cam_info); 
    end

% Capabilities button
ttstr='Display camera capabilities';
cdata=imresize(imread(fullfile(mpath,'icons','infob.jpg')),[18 18]);
hbCamAbilities=uicontrol(hpCam,'style','pushbutton','CData',cdata,'callback',@abilitiesCB,...
    'enable','on','backgroundcolor','w','position',[152 5 20 20],...
    'ToolTipString',ttstr,'enable','off');

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
    'units','pixels','fontsize',10,'Position',...
    [175 5 60 20],'enable','off',...
    'backgroundcolor',[173 216 230]/255,'callback',{@coolCB ,1},...
    'ToolTipString',ttstr);

% Stop Cooling
ttstr='Stop cooling the sensor to set point.';
hbCoolOff=uicontrol(hpCam,'style','pushbutton','string','cooler off',...
    'units','pixels','fontsize',10,'Position',...
    [238 5 60 20],'enable','off',...
    'backgroundcolor',[.8 .8 .8],'callback',{@coolCB, 0},...
    'ToolTipString',ttstr);

    function coolCB(~,~,state)  
        out=coolCamera(state);
        
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
tblTemp.Position(3:4)=tblTemp.Extent(3:4);
tblTemp.Position(1:2)=[300 4];

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
        out=setTemperature(Tnew);
        
        if ~out  && ~doDebug
           src.Data=Told; 
        end
            
        src.Data=Tnew;        
    end

% Text
ttstr='Camera sensor temperature. (green: stable; yellow: unstable; red: set point not reached)';
strtemp=uicontrol(hpCam,'style','text','string','NaN','units','pixels',...
    'backgroundcolor','w','fontsize',12,'horizontalalignment','left',...
    'foregroundcolor','r','enable','off','fontweight','bold',...
    'ToolTipString',ttstr);
strtemp.Position(3:4)=[45 20];
strtemp.Position(1:2)=[340 5];

% Text
ttstr='Camera status.';
strstatus=uicontrol(hpCam,'style','text','string','DRV_NOT_INITIALIZED','units','pixels',...
    'backgroundcolor','w','fontsize',10,'horizontalalignment','right',...
    'foregroundcolor','k','enable','on','fontweight','bold',...
    'ToolTipString',ttstr);
strstatus.Position(3:4)=[157 20];
strstatus.Position(1:2)=[hpCam.Position(3)-155 5];

% Timer to update temperature
statusTimer=timer('Name','iXonTemperatureTimer','Period',1,...
    'TimerFcn',@statusTimerFcn,'ExecutionMode','FixedSpacing');

    function statusTimerFcn(~,~)
        
        % Get the temperature
        [out,temp,outstr]=getTemperature;
        strtemp.String=[num2str(temp) ' C'];        
        cam_status.Temperature=temp;        
        
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
        [out,outstr]=getCameraStatus;
        strstatus.String=outstr;
        
    end

% Open camera shutter
ttstr='Open camera shutter.';
hbOpenShutter=uicontrol(hpCam,'style','pushbutton','string','open shutter',...
    'units','pixels','fontsize',10,'Position',[385 5 80 20],'enable','off',...
    'backgroundcolor',[255 204 0]/255,'callback',{@shutterCB,1},...
    'ToolTipString',ttstr);

ttstr='Close camera shutter.';
hbCloseShutter=uicontrol(hpCam,'style','pushbutton','string','close shutter',...
    'units','pixels','fontsize',10,'Position',[465 5 80 20],'enable','off',...
    'backgroundcolor',[255 102 120]/255,'callback',{@shutterCB,0},...
    'ToolTipString',ttstr);

    function shutterCB(~,~,state)
        out=setCameraShutter(state);
        
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


%% Save Panel

hpSave=uipanel(hF,'units','pixels','backgroundcolor','w',...
    'Position',[0 hpCam.Position(2)-30 hF.Position(3)-150 25]);

% Auto Save check box
ttstr=['Enable/Disable automatic saving to external directory. Does ' ...
    'not override saving to image history.'];
hcauto=uicontrol(hpSave,'style','checkbox','string','save images?','fontsize',10,...
    'backgroundcolor','w','Position',[5 0 100 25],'callback',@saveCheck,...
    'ToolTipString',ttstr);

% Save checkbox callback
    function saveCheck(src,~)
        if src.Value
            tSaveDir.Enable='on';
            bBrowse.Enable='on';
        else
            tSaveDir.Enable='off';
            bBrowse.Enable='off';
        end
    end
% Browse button
ttstr='Select directory to save images.';
cdata=imresize(imread(fullfile(mpath,'icons','browse.jpg')),[20 20]);
bBrowse=uicontrol(hpSave,'style','pushbutton','CData',cdata,'callback',@browseCB,...
    'enable','off','backgroundcolor','w','position',[110 2 size(cdata,[1 2])],...
    'tooltipstring',ttstr);

% String for current save directory
ttstr='The current save directory.';
tSaveDir=uicontrol(hpSave,'style','text','string','directory','fontsize',8,...
    'backgroundcolor','w','units','pixels','horizontalalignment','left',...
    'enable','off','UserData','','Position',[135 0 hF.Position(3)-135 20],...
    'tooltipstring',ttstr);

% Browse button callback
    function browseCB(~,~)
        str=getDayDir;
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

%% History 

hpNav=uipanel(hF,'units','pixels','backgroundcolor','w',...
    'Position',[hpSave.Position(1)+hpSave.Position(3) hpSave.Position(2) 300 25]);

% Button to load an image into the acquisition
ttstr='Load an image into the acquisition GUI.';
cdata=imresize(imread('icons/browse.jpg'),[20 20]);
uicontrol(hpNav,'style','pushbutton','CData',cdata,...
    'callback',@browseImageCB,'enable','on','backgroundcolor','w',...
    'position',[0 2 size(cdata,[1 2])],'ToolTipString',ttstr);

    function browseImageCB(~,~)
       loadImage; 
    end

% Checkbox for auto updating when new images are taken
ttstr='Automatically refresh to most recent image upon new image acquisition.';
cAutoUpdate=uicontrol('parent',hpNav,'units','pixels','string',...
    'hold?','value',0,'fontsize',7,'backgroundcolor','w',...
    'Style','checkbox','ToolTipString',ttstr);
cAutoUpdate.Position=[21 5 45 14];

ttstr='Jump to most recent image acquired.';
hbhistoryNow=uicontrol(hpNav,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',[char(10094) char(10094)],'fontsize',10,...
    'callback',{@chData, '0'},'ToolTipString',ttstr);
hbhistoryNow.Position=[65 2 24 20];

ttstr='Step to next more recent image';
hbhistoryLeft=uicontrol(hpNav,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',char(10094),'fontsize',10,...
    'callback',{@chData, '-'},'ToolTipString',ttstr);
hbhistoryLeft.Position=[89 2 12 20];

thistoryInd=uicontrol(hpNav,'Style','text','units','pixels',...
    'backgroundcolor','w','string','000','fontsize',12);
thistoryInd.Position=[101 2 30 20];


ttstr='Step to later image.';
hbhistoryRight=uicontrol(hpNav,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',char(10095),'fontsize',10,...
    'callback',{@chData, '+'},'ToolTipString',ttstr);
hbhistoryRight.Position=[131 2 12 20];

    function loadImage(filename)
        if nargin<1
            [filename,pathname]=uigetfile([historyDir filesep '*.mat']);
            if ~filename
                disp('No mat file chosen!');
                return;
            end
            filename=[pathname filename];
        end          
        disp(['     Loading ' filename]);        
        olddata=data;
        try
            newdata=load(filename);
            data=newdata.data;
            updateImages(data);      
        catch ME
            warning('Unable to load image, reverting to old data');
            data=olddata;
            updateImages(data);      
        end
    end

% Callback function for changing number of ROIs
    function chData(~,~,state)       
       % Get mat files in history directory
       filenames=dir([historyDir  filesep '*.mat']);
       filenames={filenames.name};       
       filenames=sort(filenames);
       filenames=flip(filenames);
       
       myname=[data.Name '.mat'];           	% Current data mat       
       ind=find(ismember(filenames,myname));    % index in filenames   
       if isempty(ind)
          ind=1; 
       end
       
        switch state
            case '-'
                ind=max([ind-1 1]);            
            case '+'
                ind=min([ind+1 length(filenames)]);
            case '0'
                ind=1;        
        end        
        thistoryInd.String=sprintf('%03d',ind);
        drawnow;        
        filename=filenames{ind};
        
        newfilename=fullfile(historyDir,filename);
        loadImage(newfilename);
    end


ttstr='Name of current image';
tImageFile=uicontrol(hpNav,'style','text','string','FILENAME','fontsize',7,...
    'backgroundcolor','w','units','pixels','horizontalalignment','left',...
    'Position',[145 2 hpNav.Position(3)-133 14],'tooltipstring',ttstr);
tImageFile.String=data.Name;
%% Image Acqusition Panel
% Panel for image acquisition controls and settings.

hpAcq=uipanel(hF,'units','pixels','backgroundcolor','w',...
    'Position',[0 hF.Position(4)-55-100 160 100],'title','acquisition');

% Start acquisition button
ttstr='Start acquisition.';
hbstart=uicontrol(hpAcq,'style','pushbutton','string','start',...
    'units','pixels','fontsize',10,'backgroundcolor',[80 200 120]/255,...
    'Position',[2 hpAcq.Position(4)-35 40 20],...    
    'Callback',@startCamCB,'ToolTipString',ttstr,'enable','off');

    function startCamCB(src,evt)
        disp('Starting acquisition');
        
        % Send acquistion start command
        out=startCamera;
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
    'units','pixels','fontsize',10,'backgroundcolor',[250 102 120]/255,...
    'Position',[45 hpAcq.Position(4)-35 40 20], ...
    'Callback',@stopCamCB,'ToolTipString',ttstr,'enable','off');


    function stopCamCB(src,evt)
        disp('Stopping acquisition');     
        
        %
        stop(acqTimer);
        stopCamera;
        
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
    'enable','on','backgroundcolor','w','position',[88 hpAcq.Position(4)-35 20 20],...
    'ToolTipString',ttstr,'enable','on');

    function helpCB(~,~)
       disp('Im helping you'); 
    end

% Continuous Acquisition checkbox
ttstr='Reinitialize camera acquisition after image acquisition.';
hcAcqRpt=uicontrol(hpAcq,'style','checkbox','string','repeat acquisition?','fontsize',8,...
    'backgroundcolor','w','Position',[5 hpAcq.Position(4)-55 120 20],...
    'ToolTipString',ttstr,'enable','on');

% Button group for acquisition mode
bgAcq = uibuttongroup(hpAcq,'units','pixels','backgroundcolor','w','BorderType','None',...
    'SelectionChangeFcn',@chAcqCB);  
bgAcq.Position(3:4)=[175 40];
bgAcq.Position(1:2)=[5 hpAcq.Position(4)-95];    

% Radio buttons for cuts vs sum
rbSingle=uicontrol(bgAcq,'Style','radiobutton','String','Normal (pwa,pwa,...,bkgd)',...
    'Position',[1 0 200 20],'units','pixels','backgroundcolor','w','Value',1,...
    'UserData','Normal','Enable','off');
rbLive=uicontrol(bgAcq,'Style','radiobutton','String','Live (be careful)',...
    'Position',[1 20 200 20],'units','pixels','backgroundcolor','w',...
    'UserData','Live','Enable','off');

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
        setCameraShutter(0);        
        hbCloseShutter.Enable='off';
        hbOpenShutter.Enable='on';
        
        switch newStr
            case 'Live'
                msg=['Entering live mode. Verify your settings before ' ...
                    'starting acquisition'];
                msgbox(msg,'Live Mode','warn','modal');  
                acq=defaultLiveAcqSettings;    
                
            case 'Normal'
                acq=defaultNormalAcqSettings;

            otherwise
                warning('Unexpected acqusition mode. What happened?');
        end        
        
    end

% Timer to check on acquisition
acqTimer=timer('Name','iXonAcquisitionWatchTimer','Period',.5,...
    'TimerFcn',@acqTimerFcn,'ExecutionMode','FixedSpacing');

% Timer function checks if acquisition is over and restarts it
    function acqTimerFcn(src,evt)      
        % Camera Status
        [out,outstr]=getCameraStatus;        
        switch outstr
            case 'DRV_IDLE'
                % Grab the images from the camera
                imgs=grabRawImages;       
                
                % Assign images metadata
                mydata=processImages(imgs);

                % Restart Acquisition if desired (auto-stopts)
                if hcAcqRpt.Value
                    startCamera;  
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
                saveData(mydata)     
                
                % Save images to save directory
                if hcauto.Value
                   saveData(mydata,tSaveDir.UserData); 
                end              
                
                % Update live preview if new                
                if ~cAutoUpdate.Value
                    data=mydata;   
                    data=updateImages(data); 
                else
                    % Just update index
                    updateHistoryInd(data);   
                end  
                
            case 'DRV_ACQUIRING'
                % Acquisition is still going.
                
            otherwise
                warning('Acuisition timer has unexpected result');
                stopCamCB;
        end        
    end



%% Image Process Panel

hpADV=uipanel(hF,'units','pixels','backgroundcolor','w',...
    'Position',[0 hpAcq.Position(2)-60 160 80],'title','processing');

% Checkbox for applying point spread function (should this be a separate
% panel?)
hcPSF=uicontrol(hpADV,'style','checkbox','string','sharpen w/ PSF','fontsize',8,...
    'backgroundcolor','w','Position',[5 0 100 20],'callback',@() disp('hi'),...
    'ToolTipString',ttstr,'enable','off');

% Checkbox for new processings
hcMask=uicontrol(hpADV,'style','checkbox','string','objective mask','fontsize',8,...
    'backgroundcolor','w','Position',[5 20 100 20],'callback',@() disp('hi'),...
    'ToolTipString',ttstr,'enable','off');

% Checkbox for enabling 2D gauss fitting
cGaussFilter=uicontrol('style','checkbox','string','gauss filter',...
    'units','pixels','parent',hpADV,'backgroundcolor','w',...
    'value',0);
cGaussFilter.Position=[5 40 80 20];

tblGaussFilter=uitable('parent',hpADV,'units','pixels',...
    'rowname',{},'columnname',{},'Data',.5,'columneditable',[true],...
    'columnwidth',{40},'fontsize',8,'ColumnFormat',{'numeric'});
tblGaussFilter.Position=[80 cGaussFilter.Position(2)-1 45 20];

% Pixels label
uicontrol('parent',hpADV,'units','pixels',...
    'style','text','string','pixels','position',[125 cGaussFilter.Position(2)+1 30 15],...
    'fontsize',8,'backgroundcolor','w');


% process button
hbprocess=uicontrol(hpADV,'style','pushbutton','string','process',...
    'units','pixels','callback',@processCB,'parent',hpADV,'backgroundcolor','w');
hbprocess.Position=[hpADV.Position(3)-45 1 45 15];

    function processCB(~,~)
       disp('reprocessing images...'); 
    end

%% Analysis Panel
hpAnl=uipanel(hF,'units','pixels','backgroundcolor','w','title','analysis');
hpAnl.Position=[0 hpADV.Position(2)-250 160 250];

% Table of ROIs
tblROI=uitable(hpAnl,'units','pixels','ColumnWidth',{30 30 30 30},...
    'ColumnEditable',true(ones(1,4)),'ColumnName',{'X1','X2','Y1','Y2'},...
    'Data',[1 512 1 512],'FontSize',8,...
    'CellEditCallback',@chROI,'RowName',{});
tblROI.Position(3:4)=tblROI.Extent(3:4)+0*[18 0];
tblROI.Position(1:2)=[5 hpAnl.Position(4)-tblROI.Position(4)-20];

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
cdata=imresize(imread(fullfile(mpath,'icons','target.jpg')),[15 15]);
hbSlctROI=uicontrol(hpAnl,'style','pushbutton','Cdata',cdata,'Fontsize',10,...
    'Backgroundcolor','w','Position',[130 tblROI.Position(2)+2 18 18],'Callback',@slctROICB,...
    'ToolTipString',ttstr);


% Callback for selecting an ROI based upon mouse click input.
    function slctROICB(~,~)
        disp(['Selecting display ROI .' ...
            ' Click two points that form the rectangle ROI.']);
        axes(axImg)                 % Select the OD image axis
        [x1,y1]=ginput(1);          % Get a mouse click
        x1=round(x1);y1=round(y1);  % Round to interger        
        p1=plot(x1,y1,'+','color','k','linewidth',1); % Plot it
        
        [x2,y2]=ginput(1);          % Get a mouse click
        x2=round(x2);y2=round(y2);  % Round it        
        p2=plot(x2,y2,'+','color','k','linewidth',1);  % Plot it

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
           warning('Unable to change display ROI.');
        end
        delete(p1);delete(p2);                   % Delete markers

    end




hcBox=uicontrol(hpAnl,'style','checkbox','string','center of mass','fontsize',8,...
    'backgroundcolor','w','Position',[5 60 100 20],'callback',@(~,~) disp('hi'),...
    'ToolTipString',ttstr,'value',1);

hcGauss=uicontrol(hpAnl,'style','checkbox','string','2D gauss','fontsize',8,...
    'backgroundcolor','w','Position',[5 40 100 20],'callback',@(~,~) disp('hi'),...
    'ToolTipString',ttstr,'enable','off');

hcGaussRot=uicontrol(hpAnl,'style','checkbox','string','2D gauss rot','fontsize',8,...
    'backgroundcolor','w','Position',[5 20 100 20],'callback',@(~,~) disp('hi'),...
    'ToolTipString',ttstr,'enable','off');

hcSingleAtoms=uicontrol(hpAnl,'style','checkbox','string','find single atoms','fontsize',8,...
    'backgroundcolor','w','Position',[5 0 100 20],'callback',@(~,~) disp('hi'),...
    'ToolTipString',ttstr,'enable','off');


% Refit button
hbfit=uicontrol(hpAnl,'style','pushbutton','string','analyze',...
    'units','pixels','callback',@cbrefit,'parent',hpAnl,'backgroundcolor','w');
hbfit.Position=[hpAnl.Position(3)-45 1 45 15];

% Callback function for redoing fits button
    function cbrefit(~,~)
        disp('Redoing fits...');
        data=updateImages(data);
    end

%% Display Options Panel

% Initialize panel object
hpDisp=uipanel(hF,'units','pixels','backgroundcolor','w','title','display');
hpDisp.Position=[0 0 160 hpAnl.Position(2)];

%%%%%% Display ROI %%%%%%

% Table for changing display limits
tbl_dispROI=uitable('parent',hpDisp,'units','pixels','RowName',{},'columnname',{},...
    'ColumnEditable',[true true true true],'CellEditCallback',@tbl_dispROICB,...
    'ColumnWidth',{24 24 24 24},'FontSize',8,'Data',[1 size(Z,2) 1 size(Z,1)]);
tbl_dispROI.Position(3:4)=tbl_dispROI.Extent(3:4);
tbl_dispROI.Position(1:2)=[2 hpDisp.Position(4)-tbl_dispROI.Position(4)-20];

% Button for maximizing the display limits
ttstr='Maximize display ROI to full image size.';
cdata=imresize(imread(fullfile(mpath,'icons','fullLim.png')),[15 15]);
hbFullLim=uicontrol(hpDisp,'style','pushbutton','Cdata',cdata,'Fontsize',10,...
    'Backgroundcolor','w','Position',[102 tbl_dispROI.Position(2)+2 18 18],'Callback',@fullDispCB,...
    'ToolTipString',ttstr);

% Button to snap display ROI to the data ROI
ttstr='Snap display ROI to data ROI.';
cdata=imresize(imread(fullfile(mpath,'icons','snapLim.png')),[15 15]);
hbSnapLim=uicontrol(hpDisp,'style','pushbutton','Cdata',cdata,'Fontsize',10,...
    'Backgroundcolor','w','Position',[120 tbl_dispROI.Position(2)+2 18 18],'Callback',@snapDispCB,...
    'ToolTipString',ttstr);

% Button to enable GUI selection of display limits
ttstr='Select the display ROI.';
cdata=imresize(imread(fullfile(mpath,'icons','target.jpg')),[15 15]);
hbSlctLim=uicontrol(hpDisp,'style','pushbutton','Cdata',cdata,'Fontsize',10,...
    'Backgroundcolor','w','Position',[138 tbl_dispROI.Position(2)+2 18 18],'Callback',@slctDispCB,...
    'ToolTipString',ttstr);

% Callback for changing display table ROI
    function tbl_dispROICB(src,evt)
        ROI=src.Data;        % Grab the new ROI     
        % Check that the data is numeric
        if sum(~isnumeric(ROI)) || sum(isinf(ROI)) || sum(isnan(ROI))
            warning('Incorrect data type provided for ROI.');
            src.Data(evt.Indices(2))=evt.PreviousData;
            return;
        end        
        ROI=round(ROI);      % Make sure this ROI are integers   

        % Keep the ROI within image bounds (this is hardcoded and could be
        % changed if we ever implement hardware ROI but want to keep 
        % absolute pixel positions relative to total sensor.)
        if ROI(2)<=ROI(1) || ROI(4)<=ROI(3)
           warning('Bad ROI specification given.');
           ROI(evt.Indices(2))=evt.PreviousData;
        end       
        if ROI(1)<1; ROI(1)=1; end       
        if ROI(3)<1; ROI(3)=1; end   
        if ROI(4)>512; ROI(4)=512;end       
        if ROI(2)>512; ROI(2)=512;end       
        src.Data=ROI;       
        try
            set(axImg,'XLim',ROI(1:2),'YLim',ROI(3:4));
            drawnow;
            resizePlots;
        catch ab
            warning('Unable to change display ROI.');
            src.Data(evt.Indices)=evt.PreviousData;
        end        
    end

% Callback for change display ROI to maximum
    function fullDispCB(~,~)
       ROI=[1 size(Z,2) 1 size(Z,1)];
       tbl_dispROI.Data=ROI;
       tbl_dispROICB(tbl_dispROI);
       resizePlots;
       drawnow;
    end

% Callback for change display to snap to analysis ROI
    function snapDispCB(~,~)
       ROI=[min(tblROI.Data(:,1)) max(tblROI.Data(:,2)) ...
           min(tblROI.Data(:,3)) max(tblROI.Data(:,4))];
       tbl_dispROI.Data=ROI;
       tbl_dispROICB(tbl_dispROI);
       resizePlots;
       drawnow;
    end

% Callback for selecting an ROI based upon mouse click input.
    function slctDispCB(~,~)
        disp(['Selecting display ROI .' ...
            ' Click two points that form the rectangle ROI.']);
        axes(axImg)                 % Select the OD image axis
        [x1,y1]=ginput(1);          % Get a mouse click
        x1=round(x1);y1=round(y1);  % Round to interger        
        p1=plot(x1,y1,'+','color','k','linewidth',1); % Plot it
        
        [x2,y2]=ginput(1);          % Get a mouse click
        x2=round(x2);y2=round(y2);  % Round it        
        p2=plot(x2,y2,'+','color','k','linewidth',1);  % Plot it

        % Create the ROI
        ROI=[min([x1 x2]) max([x1 x2]) min([y1 y2]) max([y1 y2])];

        % Constrain ROI to image
        if ROI(1)<1; ROI(1)=1; end       
        if ROI(3)<1; ROI(3)=1; end   
        if ROI(4)>512; ROI(4)=512; end       
        if ROI(2)>512; ROI(2)=512; end   
        
        % Try to update ROI graphics
        tbl_dispROI.Data=ROI;
        tbl_dispROICB(tbl_dispROI);
        resizePlots;       
        drawnow;        
        delete(p1);delete(p2);                   % Delete markers
    end


%%%%%% Color Limits %%%%%%

% Text label for color limit table on OD image
climtext=uicontrol('parent',hpDisp,'units','pixels','string','color limits : ',...
    'fontsize',8,'backgroundcolor','w','style','text');
climtext.Position(3:4)=climtext.Extent(3:4);
climtext.Position(1:2)=[2 tbl_dispROI.Position(2)-climtext.Position(4)-5];

% Table to adjust color limits on image
climtbl=uitable('parent',hpDisp,'units','pixels','RowName',{},'ColumnName',{},...
    'Data',[0 12000],'ColumnWidth',{40,40},'ColumnEditable',[true true],...
    'CellEditCallback',@climCB);
climtbl.Position(3:4)=climtbl.Extent(3:4);
climtbl.Position(1:2)=[65 tbl_dispROI.Position(2)-climtext.Position(4)-5];

% Callback for changing the color limits table
    function climCB(src,evt)
        try
            axImg.CLim=climtbl.Data;
        catch exception
            warning('Bad OD color limits given. Using old value.');
            src.Data(evt.Indices)=evt.PreviousData;
        end
    end

%%%%%% Plot Options %%%%%%

% Button group for deciding what the X/Y plots show
bgPlot = uibuttongroup(hpDisp,'units','pixels','backgroundcolor','w','BorderType','None',...
    'SelectionChangeFcn',@chPlotCB);  
bgPlot.Position(3:4)=[125 40];
bgPlot.Position(1:2)=[2 climtbl.Position(2)-bgPlot.Position(4)-2];
    
% Radio buttons for cuts vs sum
rbCut=uicontrol(bgPlot,'Style','radiobutton','String','plot cut',...
    'Position',[0 0 120 20],'units','pixels','backgroundcolor','w','Value',1);
rbSum=uicontrol(bgPlot,'Style','radiobutton','String','plot sum',...
    'Position',[0 20 120 20],'units','pixels','backgroundcolor','w','Value',0);

    function chPlotCB(~,~)
%        updatePlots(dstruct); 
        updateImages(data);
    end

%%%%%% Extra Display Stuff %%%%%%

% Checkbox for enabling display of the gaussian reticle
cGaussRet=uicontrol(hpDisp,'style','checkbox','string','show gauss reticle?',...
    'units','pixels','fontsize',8,'backgroundcolor','w','callback',@cGaussRetCB,...
    'enable','off');
cGaussRet.Position=[2 bgPlot.Position(2)-20 125 20];

    function cGaussRetCB(src,~)
       for n=1:size(tblROI.Data,1)
           pGaussRet(n).Visible=src.Value;
       end        
    end


% Checkbox for enabling display of the CoM analysis
cCoMStr=uicontrol(hpDisp,'style','checkbox','string','show CoM result?',...
    'units','pixels','fontsize',8,'backgroundcolor','w','callback',@cCoMCB,...
    'enable','on','value',1);
cCoMStr.Position=[2 cGaussRet.Position(2)-20 125 20];

    function cCoMCB(src,~)        
        set(tCoMAnalysis,'Visible',src.Value);    
    end


% Checkbox for showing/hiding crosshair
cCross=uicontrol(hpDisp,'style','checkbox','string','cross hair?',...
    'units','pixels','fontsize',8,'backgroundcolor','w','callback',@cCrossCB,...
    'enable','on','value',1);
cCross.Position=[2 cCoMStr.Position(2)-20 80 20];

    function cCrossCB(src,~)        
        set(pCrossX,'Visible',src.Value);
        set(pCrossY,'Visible',src.Value);
    end

% Checkbox for dragging cross-hair
cDrag=uicontrol(hpDisp,'style','checkbox','string','can drag?',...
    'units','pixels','fontsize',8,'backgroundcolor','w','callback',@cDragCB,...
    'enable','on','value',1);
cDrag.Position=[cCross.Position(1)+cCross.Position(3)+2 cCross.Position(2) 125 20];

% Callback for dragging the crosshair (matches cut plots with crosshair)
    function cDragCB(src,~)    
        if src.Value            
            pCrossXDrag=draggable(pCrossX);
            pCrossYDrag=draggable(pCrossY);
            pCrossXDrag.on_move_callback=@Xupdate;
            pCrossYDrag.on_move_callback=@Yupdate;
            pCrossYDrag.constraint="horizontal";
            pCrossXDrag.constraint="vertical";

        else
            delete(pCrossXDrag)
            delete(pCrossYDrag)
        end
    end


% Table to adjust cross hai
tblcross=uitable('parent',hpDisp,'units','pixels','RowName',{},'ColumnName',{},...
    'Data',[200 200],'ColumnWidth',{40,40},'ColumnEditable',[true true],...
    'CellEditCallback',@tblcrossCB);
tblcross.Position(3:4)=tblcross.Extent(3:4);
tblcross.Position(1:2)=[20 cCross.Position(2)-tblcross.Position(4)];

    function tblcrossCB(src,evt)
        m=evt.Indices(1); n=evt.Indices(2);
        
        newPos=src.Data(m,n);
        % Check that the data is numeric
        if sum(~isnumeric(newPos)) || sum(isinf(newPos)) || sum(isnan(newPos))
            warning('Incorrect data type provided for cross hair.');
            src.Data(m,n)=evt.PreviousData;
            return;
        end
        
        newPos=round(newPos);      % Make sure the cross hair is an integer

        % Check that limits go from low to high
        if newPos<1 || newPos>512
           warning('Bad cross hair position given.');
           src.Data(m,n)=evt.PreviousData;
           return;
        end
        
        src.Data(m,n)=newPos;
        % Try to update ROI graphics
        try
            % Update the cross hair position
            pCrossY.XData=[1 1]* tblcross.Data(1,1); 
            pCrossX.YData=[1 1]* tblcross.Data(1,2);            
            Zx=data.Z(tblcross.Data(1,1),:);
            set(pX,'XData',data.X,'YData',Zx);
            Zy=data.Z(:,tblcross.Data(1,2));
            set(pY,'XData',Zy,'YData',data.Y);            
        catch
           warning('Unable to change cross hair position.');
        end
    end

%% Tabular Data Results Panel
% Panel for parameters and analysis results.

hpFit=uitabgroup(hF,'units','pixels');
hpFit.Position=[160 0 200 hF.Position(4)-55];

tabs(1)=uitab(hpFit,'Title','acq','units','pixels');
tabs(2)=uitab(hpFit,'Title','param','units','pixels');
tabs(3)=uitab(hpFit,'Title','1','units','pixels');

% Table for acquisition
tbl_acq=uitable(tabs(1),'units','normalized','RowName',{},'fontsize',7,...
    'ColumnName',{},'ColumnWidth',{80 30 69},'columneditable',[false true false],...
    'Position',[0 0 1 1],'celleditcallback',@acqTblCB,'ColumnFormat',{'char','numeric','char'},'enable','off');

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

% % Load data onto tbl_acq (to 
% tbl_acq.Data=[fieldnames(acq), ...
%     struct2cell(acq)];     
% for kk=1:size(tbl_acq.Data,1)    
%     tbl_acq.Data{kk,3}=desc.(tbl_acq.Data{kk,1});   
% end

% Table for run parameters
tbl_params=uitable(tabs(2),'units','normalized','RowName',{},'fontsize',8,...
    'ColumnName',{},'ColumnWidth',{125 50},'columneditable',[false false],...
    'Position',[0 0 1 1]);

% Table for analysis outputs
tbl_analysis(1)=uitable(tabs(3),'units','normalized','RowName',{},'ColumnName',{},...
    'fontsize',8,'ColumnWidth',{60 65 65},'columneditable',false(ones(1,3)),...
    'Position',[0 0 1 1]);


%% Initialize the image panel
hp=uipanel('parent',hF,'units','pixels','backgroundcolor','w',...
    'Position',[400 0 hF.Position(3)-200 hF.Position(4)-130],...
    'bordertype','beveledin');

% Define spacing for images, useful for resizing
l=80;   % Left gap for fitting and data analysis summary

% Resize the X/Y plots and images to maximize area in figure and line up
% properly
    function resizePlots       
        if hp.Position(3)<200 || hp.Position(4)<200
            return;
        end
        
        
        % Resize the image axis     
        axImg.Position=[40 110 hp.Position(3)-200 hp.Position(4)-200];        
        
        % Get the aspect ratio of plot objects
        Rimg=axImg.PlotBoxAspectRatio;Rimg=Rimg(1)/Rimg(2);
        Rax=axImg.Position(3:4);Rax=Rax(1)/Rax(2);
        
        % Size of plot objects (position is weird in axis equal tight);
        if Rax>Rimg
            h1=axImg.Position(4);
            w1=axImg.Position(4)*Rimg;   
            hAxX.Position=[40+(axImg.Position(3)-w1)/2 axImg.Position(2)-l w1 80-15];
            hAxY.Position=[40+(axImg.Position(3)+w1)/2+15 axImg.Position(2) 80 h1];
        else
            w1=axImg.Position(3);
            h1=w1/Rimg;            
            hAxX.Position=[axImg.Position(1) 110+(axImg.Position(4)-h1)/2-l ...
                w1 80-15];
            hAxY.Position=[axImg.Position(1)+axImg.Position(3)+15 ...
                110+(axImg.Position(4)-h1)/2 l h1];            
        end
        
        % Match cut limits with the images limits
        set(hAxX,'XLim',axImg.XLim,'XTick',axImg.XTick);
        set(hAxY,'YLim',axImg.YLim,'YTick',axImg.YTick);
        
        % Move the colorbar
        cBar.Position=[hAxX.Position(1) hAxY.Position(2)+hAxY.Position(4)+23 ...
            hAxX.Position(3) 15]; 
    end

% Initialize image axis
axImg=axes('parent',hp);cla
co=get(gca,'colororder');
hImg=imagesc(data.X,data.Y,data.Z);
set(axImg,'box','on','linewidth',.1,'fontsize',10,'units','pixels',...
    'XAxisLocation','top','colormap',colormap(cmap),...
    'xcolor',co(4,:),'ycolor',co(4,:));
hold on
axImg.Position=[50 150 hp.Position(3)-200 hp.Position(4)-200];
axis equal tight

% Cross Hair Plots
pCrossX=plot([1 512],[512/2 512/2],'-','color',co(5,:),'linewidth',1);
pCrossY=plot([512/2 512/2],[1 512],'-','color',co(5,:),'linewidth',1);

% Make cross hair draggable. This declaration makes this variable global
pCrossXDrag=draggable(pCrossX);
pCrossYDrag=draggable(pCrossY);
pCrossXDrag.on_move_callback=@Xupdate;
pCrossYDrag.on_move_callback=@Yupdate;

% Delete so that it is initially undraggable. 
% delete(pCrossXDrag)
% delete(pCrossYDrag)

% Callback for adjusting the X crosshair
    function Xupdate(g,~)
        % If you drag, you're not plotting a sum
        rbCut.Value=1;
        rbSum.Value=0;
        
        % Get the cross hair poisition
        Ycross=round(g.YData(1));
        Zx=data.Z(Ycross,:);
        set(pX,'XData',data.X,'YData',Zx);
        tblcross.Data(1,2)=Ycross;
        
        drawnow;       
    end

% Callback 
    function Yupdate(g,~)
        % If you drag, you're not plotting a sum
        rbCut.Value=1;
        rbSum.Value=0;
        
        % Get the cross hair poisition
        Xcross=round(g.XData(1));
        Zy=data.Z(:,Xcross);
        set(pY,'YData',data.Y,'XData',Zy);
        tblcross.Data(1,1)=Xcross;

        drawnow;       
    end

% file name string
tImageFile=text(3,3,'FILENAME','units','pixels','fontsize',8,'fontweight','bold',...
    'horizontalalignment','left','verticalalignment','bottom','margin',1,...
    'interpreter','none','backgroundcolor',[1 1 1 .5]);

% box count analysis sytring
tCoMAnalysis=text(.99,0.01,'FILENAME','units','normalized','fontsize',12,'fontweight','bold',...
    'horizontalalignment','right','verticalalignment','bottom','margin',1,...
    'interpreter','latex',...
    'color',co(2,:),'backgroundcolor',[1 1 1 .1]);


% Box for ROI (this will become an array later)
pROI=rectangle('position',[1 1 512 512],'edgecolor',co(1,:),'linewidth',2);
% Reticle for gaussian fit (this will become an array later)
pGaussRet=plot(0,0,'-','linewidth',1,'Visible','off','color',co(1,:));
% Color bar
cBar=colorbar('fontsize',8,'units','pixels','location','northoutside');

axImg.CLim=climtbl.Data;
drawnow;

% X Cut/Sum Axis
hAxX=axes('box','on','linewidth',1,'fontsize',10,...
    'XAxisLocation','Bottom','units','pixels','parent',hp);
hAxX.Position=[axImg.Position(1) axImg.Position(2)-l axImg.Position(3) l];
hold on
% Add X data data and fit plots
pX=plot(data.X,ones(length(data.X),1),'k.-');
pXF=plot(data.X,ones(length(data.X),1),'-','Visible','off','color',co(1,:),'linewidth',2);


% Y Cut/Sum Axis
hAxY=axes('box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'YAxisLocation','Right','YDir','Reverse','parent',hp);
hAxY.Position=[axImg.Position(1)+axImg.Position(3) axImg.Position(2) l axImg.Position(4)];
hold on
% Add Y data data and fit plots
pY=plot(ones(length(data.Y),1),data.Y,'k.-'); 
pYF=plot(data.X,ones(length(data.X),1),'-','Visible','off','color',co(1,:),'linewidth',2);



set(axImg,'XLim',tbl_dispROI.Data(1:2),'YLim',tbl_dispROI.Data(3:4));


%%

function data=updateImages(data)
    % Grab the ROI
    ROI=tblROI.Data;
    data.ROI=ROI;       
    x=ROI(1):ROI(2);
    y=ROI(3):ROI(4); 

    % Create sub image to do center of mass analysis
    Zsub=data.Z(y,x);
    
    % Perform Box Count
    data=boxCount(data);
    
    % Update box count string
    str=[ num2str(max(max(data.Z)),'%.2e') ' max counts ' newline ...
        num2str(data.BoxCount.Nraw,'%.2e') ' counts' newline ...
        '$(X_\mathrm{c},Y_\mathrm{c}) = ' '('  num2str(round(data.BoxCount.Xc,1)) ',' ...
        num2str(round(data.BoxCount.Yc,1)) ')$' newline ...
        '$(\sigma_X,\sigma_Y) = ' '('  num2str(round(data.BoxCount.Xs,1)) ',' ...
        num2str(round(data.BoxCount.Ys,1)) ')$'];    
    
    % Update box count string object
    set(tCoMAnalysis,'String',str);    
    
    % Update X, Y, and Z objects
    set(hImg,'XData',data.X,'YData',data.Y,'CData',data.Z);
    
    % Move cross hair to center of mass
    pCrossX.YData=[1 1]*round(data.BoxCount.Yc);
    pCrossY.XData=[1 1]*round(data.BoxCount.Xc);

    % Update table that trackes cross hair
    tblcross.Data(1,2)=pCrossX.YData(1);
    tblcross.Data(1,1)=pCrossY.XData(1);    
    
    % Update plots if sum
    if rbSum.Value
        Zy=sum(Zsub,2);
        Zx=sum(Zsub,1);          
        set(pX,'XData',x,'YData',Zx);
        set(pY,'XData',Zy,'YData',y);
        drawnow;
    end
    
    % Update plots if cut
    if rbCut.Value
        Zy=data.Z(:,round(data.BoxCount.Xc));
        Zx=data.Z(:,round(data.BoxCount.Yc));
        set(pX,'XData',data.X,'YData',Zx);
        set(pY,'XData',Zy,'YData',data.Y);
        drawnow;
    end
    
    % Update table parametesr
    tbl_params.Data=[fieldnames(data.Params), ...
        struct2cell(data.Params)]; 
    
    updateHistoryInd(data);   
  
end



data=updateImages(data);
drawnow;
SizeChangedFcn
axes(axImg);

%%
function imgs=grabRawImages

    % How many images to grab
    [ret,first,last] = GetNumberNewImages;
    
    numpix=512^2;
    % Grab the data (number just sets buffer size)
    [ret,D] = GetAcquiredData(last*512^2);    

    imgs={};
    imgmats=zeros(512,512,last);

    for j = 1:last % break up into individual images
        ii=double(D((1+(j-1)*numpix):(j*numpix)));
        imgs{j} = reshape(ii,512,512);
        % Subtract off aritifical baseline counts (see manual).
%         imgs{j}=imgs{j}-200;
        imgmats(:,:,j)=imgs{j};
    end 
    
    imgs=imgmats;  
end

function mydata=processImages(imgs)
    mydata=struct;
    % Create the image data structure
    mydata=struct;
    mydata.Date=datevec(now);
    mydata.Name=['iXonUltra_' datestr(mydata.Date,'yyyy-mm-dd_HH-MM-SS')];   

    % Grab the images
    mydata.RawImages=imgs;

    % Add magnification
    mydata.Magnification=mag;
    
    % Add X and Y vectors
    mydata.X=1:size(mydata.RawImages,2);
    mydata.Y=1:size(mydata.RawImages,1);    
    
    % For now always assumed that its PWA,BKGD
    mydata.Z=mydata.RawImages(:,:,2)-mydata.RawImages(:,:,1);
   

    % Grab the sequence parameters
    [mydata.Params,dstr]=grabSequenceParams;        
    mydata.Params.ExecutionDate=dstr;    

    % Append acquisition information
    mydata.CameraInformation=cam_info;
    mydata.AcquisitionInformation=acq;
    mydata.AcquisitionDescription=desc;
end


    function saveData(data,saveDir)
        if nargin==1
           saveDir=historyDir;
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
        % Update image string
        set(tImageFile,'String',data.Name);

        % Upate history list
        filenames=dir([historyDir  filesep '*.mat']);
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
        thistoryInd.String=sprintf('%03d',ind); 
    end
end

%% Camera Functions

% Connect to the Andor iXon Ultra
function out=connectCam
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

% Disconnect from the Andor iXon Ultra
function out=disconnectCam
    disp('Disconnecting from the iXon camera.');
    
    % Shut down cooler
    fprintf('Turning off cooler ... ');
    [ret]=SetCoolerMode(1);     
    disp(error_code(ret));

    % Close the shutter
    fprintf('Closing the shutter ... ');
    [ret]=SetShutter(1,2,0,0);  
    disp(error_code(ret));

    % Shut down the camera
    fprintf('Shutting down camera ... ');
    [ret]=AndorShutDown;        
    disp(error_code(ret))
    
    if isequal(error_code(ret),'DRV_SUCCESS')
        out=1;
    else
        warning('Unable to shut down the iXon camera.');
        out=0;
    end
end

% Set the temperature setpoint
function out=setTemperature(temp)
    fprintf(['Changing temperature set point to ' num2str(temp) ' C ...']);
    ret=SetTemperature(temp);
    disp(error_code(ret))
    
    if isequal(error_code(ret),'DRV_SUCCESS') 
        out=1;
    else
        warning('Unable to change iXon temperature set point.');
        out=0;
    end
end

% Get the temperature
function [out,temp,outstr]=getTemperature
    [ret,temp]=GetTemperature;
    out=1;
    outstr=error_code(ret);
end

% Get the camera status
function [out,outstr]=getCameraStatus
    out=0;
    [ret,outstr]=AndorGetStatus;
    outstr=error_code(outstr);
    
    if isequal(error_code(ret),'DRV_SUCCESS')
        out=1;
    else
        warning('Unable to read iXon status.');
        out=0;
    end
end

% Get Detector
function [out,xpx,ypx]=GetDetectorInfo    
    [ret,xpx,ypx]=GetDetector;
    
    if isequal(error_code(ret),'DRV_SUCCESS')
        out=1;
    else
        warning('Unable to get detector informaton.');
        out=0;
    end
end

% Engage/Disengage TEC
function out=coolCamera(state)
    if state
        fprintf('Engaging TEC to cool sensor ... ');
        ret=CoolerON;
    else
        fprintf('Disengaging TEC to cool sensor ... ');
        ret=CoolerOFF;
    end    
    disp(error_code(ret));
    
    if isequal(error_code(ret),'DRV_SUCCESS')
        out=1;
    else
        warning('Unable to read iXon temperature.');
        out=0;
    end
end

% Set the camera shutter
function out=setCameraShutter(state)
    typ=1; % HIGH TO OPEN
    
    % shutter_mode : 0: Auto, 1:Open ,2:Close    
    if state
        fprintf('Opening shutter ... ');
        shutter_mode=1;
    else
        fprintf('Closing shutter ... ');
        shutter_mode=2;
    end    
    
    ret=SetShutter(typ,shutter_mode,0,0);
    
    disp(error_code(ret));
    
    if isequal(error_code(ret),'DRV_SUCCESS')
        out=1;
    else
        warning('Unable to set iXon shutter.');
        out=0;
    end
end


% Start Acquisition
function out=startCamera
    fprintf('Starting acquisition ... ');
    [ret]=StartAcquisition;
    disp(error_code(ret));
    
    if isequal(error_code(ret),'DRV_SUCCESS')
        out=1;
    else
        warning('Unable to start acquisition.');
        out=0;
    end
end

% Stop Acquisition
function out=stopCamera
    fprintf('Stopping acquisition ... ');
    [ret]=AbortAcquisition;
    disp(error_code(ret));
    
    switch error_code(ret)
        case 'DRV_SUCCESS'
            out=1;
        case 'DRV_IDLE'
            out=1;
            disp('Camera acquisition not running.');
        otherwise
            warning('Error stopping acquisition.');
            out=0;
    end   

end

% Software Trigger
function out=softwareTrigger
    fprintf('Sending software trigger ... ');
    [ret]=SendSoftwareTrigger;
    disp(error_code(ret));
    if isequal(error_code(ret),'DRV_SUCCESS')
        out=1;
    else
        warning('Unable to send software trigger.');
        out=0;
    end
end



%% Analysis Functions

function dstruct=boxCount(dstruct,bgROI)
    fprintf('Performing box count analysis ...');   
    
    if nargin==1
        bgROI=NaN;
    end
  
    
    BoxCount=struct;    
    for k=1:size(dstruct.ROI,1)
        ROI=dstruct.ROI(k,:);
        x=dstruct.X(ROI(1):ROI(2));                 % X vector
        y=dstruct.Y(ROI(3):ROI(4));                 % Y vector
        z=double(dstruct.Z(ROI(3):ROI(4),ROI(1):ROI(2)));
        nbg=0;
        
        if nargin==2
            zbg=double(dstruct.Z(bgROI(3):bgROI(4),bgROI(1):bgROI(2)));
            nbg=sum(sum(zbg))/(size(zbg,1)*size(zbg,2)); % count density
        end 
        
        Nraw=sum(sum(z));
        Nbg=nbg*size(z,1)*size(z,2);  
        
        zNoBg=z-nbg;        
        Ncounts=sum(sum(zNoBg));   
        zY=sum(zNoBg,2)';
        zX=sum(zNoBg,1);
               
        % Calculate center of mass
        Xc=sum(zX.*x)/Ncounts;
        Yc=sum(zY.*y)/Ncounts;
        
        % Calculate central second moment/variance and the standard
        % deviation
        X2=sum(zX.*(x-Xc).^2)/Ncounts; % x variance
        Xs=sqrt(X2); % standard deviation X
        Y2=sum(zY.*(y-Yc).^2)/Ncounts; % x variance
        Ys=sqrt(Y2); % standard deviation Y               

        BoxCount(k).Ncounts=Ncounts;    % Number of counts (w/ bkgd removed)
        BoxCount(k).Nraw=Nraw;          % Raw of number of counts
        BoxCount(k).Nbkgd=Nbg;          % Bakcground number of counts
        BoxCount(k).nbkgd=nbg;          % Background counts/px
        BoxCount(k).bgROI=bgROI;        % ROI for calculating bgkd
        BoxCount(k).Xc=Xc;              % X center of mass
        BoxCount(k).Yc=Yc;              % Y center of mass
        BoxCount(k).Xs=Xs;              % X standard deviation
        BoxCount(k).Ys=Ys;              % Y standard deviation
    end    
    dstruct.BoxCount=BoxCount;
    disp('done');
end


%% HELPER
function s3=getDayDir
    t=now;
    
    d=['Y:\Data'];

    if ~exist(d,'dir')
        s3=pwd;
        return;
    end
    s1=datestr(t,'yyyy');s2=datestr(t,'yyyy.mm');s3=datestr(t,'mm.dd');
    s1=[d filesep s1];s2=[s1 filesep s2];s3=[s2 filesep s3];

    if ~exist(s1,'dir'); mkdir(s1); end
    if ~exist(s2,'dir'); mkdir(s2); end
    if ~exist(s3,'dir'); mkdir(s3); end
end

function [out,dstr]=grabSequenceParams(src)
    if nargin~=1
        src='Y:\_communication\control.txt';
    end
    disp(['Opening information from from ' src]);
    out=struct;
    % Open the control file
    [fid,errmsg] = fopen(src,'rt');
    if ~isempty(errmsg)
       warning('Unable to read control.txt. Aborting association');
       return
    end
    % Read the first six lines (and throw them away)
    fgetl(fid);
    fgetl(fid);
    dstr=fgetl(fid);   % Date line
    dstr=dstr(17:end-1); % get date string (this is a bit risky coding)
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    % Read the parameters (each line looks like "k_cMOT_detuning: 5")
    params = textscan(fid,'%[^:] %*s %s');
    % Close the file
    fclose(fid);
    % Convert the string into a structure
    out=cell2struct(num2cell(str2double(params{2})),params{1});
end

