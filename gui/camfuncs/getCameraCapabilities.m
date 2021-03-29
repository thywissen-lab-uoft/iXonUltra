function cam_skills=getCameraCapabilities
% Author : C Fujiwara
%
% This function interprets the output given by the GetCapabilities function
% outlined by the SDK Pp. 111-129. The function outputs an interger array
% where each bit corresponds to a different setting.  The intepretation of
% each bit is outlined in the SDK.
%
% The primary purpose of this function is to simply reformat the outputs to
% a human readable form using the struct object in MATLAB.

cam_skills=struct;
cam_skills.ulAcqModes=[];
cam_skills.ulReadModes=[];
cam_skills.ulFTReadModes=[];
cam_skills.ulTriggerModes=[];
cam_skills.ulCameraType=[];
cam_skills.ulPixelModes=[];
cam_skills.ulSetFunctions=[];
cam_skills.ulGetFunctions=[];
cam_skills.ulFeatures=[];
cam_skills.ulPCICard=[];
cam_skills.ulEMGainCapability=[];

fprintf('Reading camera capabilities ... ');
[ret,...
   int32_ulAcqModes,...
   int32_ulReadModes,...
   int32_ulFTReadModes,...
   int32_ulTriggerModes,...
   int32_ulCameraType,...
   int32_ulPixelModes,...
   int32_ulSetFunctions,...
   int32_ulGetFunctions,...
   int32_ulFeatures,...
   int32_ulPCICard
   int32_ulEMGainCapability]=GetCapabilities;

% Note that MATLAB is 1 indexed while the codes are 0 indexed

%% Verify good output

if ~isequal(error_code(ret),'DRV_SUCCESS')
   warning('Unable to send software trigger.');
   out=0;
   return;
end


%% AcqModes
ulAcqModes=struct;
ulAcqModes.AC_ACQMODE_SINGLE=bitget(int32_ulAcqModes,1);
ulAcqModes.AC_ACQMODE_VIDEO=bitget(int32_ulAcqModes,2);
ulAcqModes.AC_ACQMODE_ACCUMULATE=bitget(int32_ulAcqModes,3);
ulAcqModes.AC_ACQMODE_KINETIC=bitget(int32_ulAcqModes,4);
ulAcqModes.AC_ACQMODE_FRAMETRANSFER=bitget(int32_ulAcqModes,5);
ulAcqModes.AC_ACQMODE_FASTKINETICS=bitget(int32_ulAcqModes,6);
ulAcqModes.AC_ACQMODE_OVERLAP=bitget(int32_ulAcqModes,7);

%% Read Modes
ulReadModes=struct;
ulReadModes.AC_READMODE_FULLIMAGE=bitget(int32_ulReadModes,1);
ulReadModes.AC_READMODE_SUBIMAGE=bitget(int32_ulReadModes,2);
ulReadModes.AC_READMODE_SINGLETRACK=bitget(int32_ulReadModes,3);
ulReadModes.AC_READMODE_FVB=bitget(int32_ulReadModes,4);
ulReadModes.AC_READMODE_MULTITRACK=bitget(int32_ulReadModes,5);
ulReadModes.AC_READMODE_RANDOMTRACK=bitget(int32_ulReadModes,6);

%% Frame Transfer Read Modes
ulFTReadModes=struct;
ulFTReadModes.AC_READMODE_FULLIMAGE=bitget(int32_ulFTReadModes,1);
ulFTReadModes.AC_READMODE_SUBIMAGE=bitget(int32_ulFTReadModes,2);
ulFTReadModes.AC_READMODE_SINGLETRACK=bitget(int32_ulFTReadModes,3);
ulFTReadModes.AC_READMODE_FVB=bitget(int32_ulFTReadModes,4);
ulFTReadModes.AC_READMODE_MULTITRACK=bitget(int32_ulFTReadModes,5);
ulFTReadModes.AC_READMODE_RANDOMTRACK=bitget(int32_ulFTReadModes,6);

%% Trigger Modes
ulTriggerModes=struct;
ulTriggerModes.AC_TRIGGERMODE_INTERNAL=bitget(int32_ulTriggerModes,1);
ulTriggerModes.AC_TRIGGERMODE_EXTERNAL=bitget(int32_ulTriggerModes,2);
ulTriggerModes.AC_TRIGGERMODE_EXTERNAL_FVB_EM=bitget(int32_ulTriggerModes,3);
ulTriggerModes.AC_TRIGGERMODE_CONTINUOUS=bitget(int32_ulTriggerModes,4);
ulTriggerModes.AC_TRIGGERMODE_EXTERNALSTART=bitget(int32_ulTriggerModes,5);
%     ulTriggerModes.AC_TRIGGERMODE_BULB=bitget(int32_ulTriggerModes,6); depreciated
ulTriggerModes.AC_TRIGGERMODE_EXTERNALEXPOSURE=bitget(int32_ulTriggerModes,6);
ulTriggerModes.AC_TRIGGERMODE_INVERTED=bitget(int32_ulTriggerModes,7);
ulTriggerModes.AC_TRIGGERMODE_EXTERNAL_CHARGESHIFTING=bitget(int32_ulTriggerModes,8);

%% CameraType
camtypes={'Andor PDA','Andor iXon','Andor ICCD','Andor EMCCD','Andor CCD',...
    'Andor iStar','THIRD PARTY','Andor iDus','Andor Newton','Andor Surcam',...
    'Andor USB ICCD','Andor Luca','Reserved','Andor iKon','Andor InGaAs',...
    'Andor iVac','Andor Clara','Andor USB iStar'};
ulCameraType=camtypes{int32_ulCameraType+1};

%% PCI Card
ulPCICard=int32_ulPCICard;


%% EM Gain Capability
ulEMGainCapability=struct;
ulEMGainCapability.AC_EMGAIN_8BIT=bitget(int32_ulEMGainCapability,1);
ulEMGainCapability.AC_EMGAIN_12BIT=bitget(int32_ulEMGainCapability,2);
ulEMGainCapability.AC_EMGAIN_LINEAR12=bitget(int32_ulEMGainCapability,3);
ulEMGainCapability.AC_EMGAIN_REAL12=bitget(int32_ulEMGainCapability,4);
  
        
end
