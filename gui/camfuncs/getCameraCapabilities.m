function cam_skills=getCameraCapabilities
% Author : C Fujiwara
%
% This function interprets the output given by the GetCapabilities function
% outlined by the SDK Pp. 111-129. The function outputs an integer array
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
   int32_ulPCICard,...
   int32_ulEMGainCapability]=GetCapabilities;

nums=[int32_ulAcqModes,int32_ulReadModes,int32_ulFTReadModes,...
    int32_ulTriggerModes,int32_ulCameraType,int32_ulPixelModes,int32_ulSetFunctions,...
    int32_ulGetFunctions,int32_ulFeatures,int32_ulPCICard,int32_ulEMGainCapability];

% Note that MATLAB is 1 indexed while the codes are 0 indexed


%% Verify good output

if ~isequal(error_code(ret),'DRV_SUCCESS')
   warning('Unable to send software trigger.');
   out=0;
   return;
end


%% AcqModes
ulAcqModes=struct;
ulAcqModes.AC_ACQMODE_SINGLE=bitget(int32_ulAcqModes,0+1);
ulAcqModes.AC_ACQMODE_VIDEO=bitget(int32_ulAcqModes,1+1);
ulAcqModes.AC_ACQMODE_ACCUMULATE=bitget(int32_ulAcqModes,2+1);
ulAcqModes.AC_ACQMODE_KINETIC=bitget(int32_ulAcqModes,3+1);
ulAcqModes.AC_ACQMODE_FRAMETRANSFER=bitget(int32_ulAcqModes,4+1);
ulAcqModes.AC_ACQMODE_FASTKINETICS=bitget(int32_ulAcqModes,5+1);
ulAcqModes.AC_ACQMODE_OVERLAP=bitget(int32_ulAcqModes,6+1);

cam_skills.ulAcqModes=ulAcqModes;
%% Read Modes
ulReadModes=struct;
ulReadModes.AC_READMODE_FULLIMAGE=bitget(int32_ulReadModes,0+1);
ulReadModes.AC_READMODE_SUBIMAGE=bitget(int32_ulReadModes,1+1);
ulReadModes.AC_READMODE_SINGLETRACK=bitget(int32_ulReadMode,2+1);
ulReadModes.AC_READMODE_FVB=bitget(int32_ulReadModes,3+1);
ulReadModes.AC_READMODE_MULTITRACK=bitget(int32_ulReadModes,4+1);
ulReadModes.AC_READMODE_RANDOMTRACK=bitget(int32_ulReadModes,5+1);

cam_skills.ulReadModes=ulReadModes;
%% Frame Transfer Read Modes
ulFTReadModes=struct;
ulFTReadModes.AC_READMODE_FULLIMAGE=bitget(int32_ulFTReadModes,0+1);
ulFTReadModes.AC_READMODE_SUBIMAGE=bitget(int32_ulFTReadModes,1+1);
ulFTReadModes.AC_READMODE_SINGLETRACK=bitget(int32_ulFTReadModes,2+1);
ulFTReadModes.AC_READMODE_FVB=bitget(int32_ulFTReadModes,3+1);
ulFTReadModes.AC_READMODE_MULTITRACK=bitget(int32_ulFTReadModes,4+1);
ulFTReadModes.AC_READMODE_RANDOMTRACK=bitget(int32_ulFTReadModes,5+1);

cam_skills.ulFTReadModes=ulFTReadModes;
%% Trigger Modes
ulTriggerModes=struct;
ulTriggerModes.AC_TRIGGERMODE_INTERNAL=bitget(int32_ulTriggerModes,0+1);
ulTriggerModes.AC_TRIGGERMODE_EXTERNAL=bitget(int32_ulTriggerModes,1+1);
ulTriggerModes.AC_TRIGGERMODE_EXTERNAL_FVB_EM=bitget(int32_ulTriggerModes,2+1);
ulTriggerModes.AC_TRIGGERMODE_CONTINUOUS=bitget(int32_ulTriggerModes,3+1);
ulTriggerModes.AC_TRIGGERMODE_EXTERNALSTART=bitget(int32_ulTriggerModes,4+1);
%     ulTriggerModes.AC_TRIGGERMODE_BULB=bitget(int32_ulTriggerModes,5+1); depreciated
ulTriggerModes.AC_TRIGGERMODE_EXTERNALEXPOSURE=bitget(int32_ulTriggerModes,5+1);
ulTriggerModes.AC_TRIGGERMODE_INVERTED=bitget(int32_ulTriggerModes,6+1);
ulTriggerModes.AC_TRIGGERMODE_EXTERNAL_CHARGESHIFTING=bitget(int32_ulTriggerModes,7+1);

cam_skills.ulTriggerModes=ulTriggerModes;
%% CameraType
camtypes={'Andor PDA','Andor iXon','Andor ICCD','Andor EMCCD','Andor CCD',...
    'Andor iStar','THIRD PARTY','Andor iDus','Andor Newton','Andor Surcam',...
    'Andor USB ICCD','Andor Luca','Reserved','Andor iKon','Andor InGaAs',...
    'Andor iVac','Andor Clara','Andor USB iStar','n/a','n/a','n/a','Andor iXon Ultra',};
% i got the additional inex for the ixon ultra by looking at the Ni's group
% github
% int32_ulCameraType

ulCameraType=int32_ulCameraType;
ulCameraType=camtypes{int32_ulCameraType+1};

cam_skills.ulCameraType=ulCameraType;
%% PixelModes
ulPixelModes=struct;
ulPixelModes.AC_PIXELMODE_8BIT=bitget(int32_ulPixelModes,0+1);
ulPixelModes.AC_PIXELMODE_14BIT=bitget(int32_ulPixelModes,1+1);
ulPixelModes.AC_PIXELMODE_16BIT=bitget(int32_ulPixelModes,2+1);
ulPixelModes.AC_PIXELMODE_32BIT=bitget(int32_ulPixelModes,3+1);

pxmodes={'AC_PIXELMODE_MONO','AC_PIXELMODE_RGB','AC_PIXELMODE_CMY'};

ind=bin2dec(strtrim(num2str(bitget(int32_ulPixelModes,17:32))))+1;

cam_skills.ulPixelModes=pxmodes{ind};

%% ulSetFunctions
ulSetFunctions=struct;
ulSetFunctions.AC_SETFUNCTION_VREADOUT=bitget(int32_ulSetFunctions,0+1);
ulSetFunctions.AC_SETFUNCTION_HREADOUT=bitget(int32_ulSetFunctions,1+1);
ulSetFunctions.AC_SETFUNCTION_TEMPERATURE=bitget(int32_ulSetFunctions,2+1);
ulSetFunctions.AC_SETFUNCTION_MCPGAIN=bitget(int32_ulSetFunctions,3+1);
ulSetFunctions.AC_SETFUNCTION_EMCCDGAIN=bitget(int32_ulSetFunctions,4+1);
ulSetFunctions.AC_SETFUNCTION_BASELINECLAMP=bitget(int32_ulSetFunctions,5+1);
ulSetFunctions.AC_SETFUNCTION_VSAMPLITUDE=bitget(int32_ulSetFunctions,6+1);
ulSetFunctions.AC_SETFUNCTION_HIGHCAPACITY=bitget(int32_ulSetFunctions,7+1);

ulSetFunctions.AC_SETFUNCTION_BASELINEOFFSET=bitget(int32_ulSetFunctions,8+1);
ulSetFunctions.AC_SETFUNCTION_PREAMPGAIN=bitget(int32_ulSetFunctions,9+1);
ulSetFunctions.AC_SETFUNCTION_CROPMODE=bitget(int32_ulSetFunctions,10+1);
ulSetFunctions.AC_SETFUNCTION_DMAPARAMETERS=bitget(int32_ulSetFunctions,11+1);
ulSetFunctions.AC_SETFUNCTION_HORIZONTALBIN=bitget(int32_ulSetFunctions,12+1);
ulSetFunctions.AC_SETFUNCTION_MULTITRACKHRANGE=bitget(int32_ulSetFunctions,13+1);
ulSetFunctions.AC_SETFUNCTION_RANDOMTRACKNOGAPS=bitget(int32_ulSetFunctions,14+1);
ulSetFunctions.AC_SETFUNCTION_EMADVANCED=bitget(int32_ulSetFunctions,15+1);

ulSetFunctions.AC_SETFUNCTION_GATEMODE=bitget(int32_ulSetFunctions,16+1);
ulSetFunctions.AC_SETFUNCTION_DDGTIMES=bitget(int32_ulSetFunctions,17+1);
ulSetFunctions.AC_SETFUNCTION_IOC=bitget(int32_ulSetFunctions,18+1);
ulSetFunctions.AC_SETFUNCTION_INTELLIGATE=bitget(int32_ulSetFunctions,19+1);
ulSetFunctions.AC_SETFUNCTION_INSERTION_DELAY=bitget(int32_ulSetFunctions,20+1);
ulSetFunctions.AC_SETFUNCTION_GATESTEP=bitget(int32_ulSetFunctions,21+1);
ulSetFunctions.AC_SETFUNCTION_TRIGGERTERMINATION=bitget(int32_ulSetFunctions,22+1);
ulSetFunctions.AC_SETFUNCTION_EXTENDEDNIR=bitget(int32_ulSetFunctions,23+1);
ulSetFunctions.AC_SETFUNCTION_SPOOLTHREADCOUNT=bitget(int32_ulSetFunctions,24+1);

% In newer version of SDK?
% ulSetFunctions.AC_SETFUNCTION_REGISTERPACK=bitget(int32_ulSetFunctions,25+1);
% ulSetFunctions.AC_SETFUNCTION_PRESCANS=bitget(int32_ulSetFunctions,26+1);
% ulSetFunctions.AC_SETFUNCTION_GATEWIDTHSTEP=bitget(int32_ulSetFunctions,27+1);
% ulSetFunctions.AC_SETFUNCTION_EXTENDED_CROP_MODE=bitget(int32_ulSetFunctions,28+1);

cam_skills.ulSetFunctions=ulSetFunctions;

%% ulGetFunctions
ulGetFunctions=struct;
ulGetFunctions.AC_GETFUNCTION_TEMPERATURE=bitget(int32_ulGetFunctions,0+1);
ulGetFunctions.AC_GETFUNCTION_TEMPERATURERANGE=bitget(int32_ulGetFunctions,2+1);
ulGetFunctions.AC_GETFUNCTION_DETECTORSIZE=bitget(int32_ulGetFunctions,3+1);
ulGetFunctions.AC_GETFUNCTION_MCPGAIN=bitget(int32_ulGetFunctions,4+1);
ulGetFunctions.AC_GETFUNCTION_EMCCDGAIN=bitget(int32_ulGetFunctions,5+1);
ulGetFunctions.AC_GETFUNCTION_GATEMODE=bitget(int32_ulGetFunctions,7+1);
ulGetFunctions.AC_GETFUNCTION_DDGTIMES=bitget(int32_ulGetFunctions,8+1);
ulGetFunctions.AC_GETFUNCTION_IOC=bitget(int32_ulGetFunctions,9+1);

ulGetFunctions.AC_GETFUNCTION_INTELLIGATE=bitget(int32_ulGetFunctions,10+1);
ulGetFunctions.AC_GETFUNCTION_INSERTION_DELAY=bitget(int32_ulGetFunctions,11+1);
ulGetFunctions.AC_GETFUNCTION_PHOSPHORSTATUS=bitget(int32_ulGetFunctions,13+1);
ulGetFunctions.AC_GETFUNCTION_BASELINECLAMP=bitget(int32_ulGetFunctions,15+1);

cam_skills.ulGetFunctions=ulGetFunctions;

%% ulFeatures
ulFeatures=struct;
ulFeatures.AC_FEATURES_POLLING=bitget(int32_ulFeatures,0+1);
ulFeatures.AC_FEATURES_EVENTS=bitget(int32_ulFeatures,1+1);
ulFeatures.AC_FEATURES_SPOOLING=bitget(int32_ulFeatures,2+1);
ulFeatures.AC_FEATURES_SHUTTER=bitget(int32_ulFeatures,3+1);
ulFeatures.AC_FEATURES_SHUTTEREX=bitget(int32_ulFeatures,4+1);
ulFeatures.AC_FEATURES_EXTERNAL_I2C=bitget(int32_ulFeatures,5+1);
ulFeatures.AC_FEATURES_SATURATIONEVENT=bitget(int32_ulFeatures,6+1);

ulFeatures.AC_FEATURES_FANCONTROL=bitget(int32_ulFeatures,7+1);
ulFeatures.AC_FEATURES_MIDFANCONTROL=bitget(int32_ulFeatures,8+1);
ulFeatures.AC_FEATURES_TEMPERATUREDURINGACQUISITION=bitget(int32_ulFeatures,9+1);
ulFeatures.AC_FEATURES_KEEPCLEANCONTROL=bitget(int32_ulFeatures,10+1);
ulFeatures.AC_FEATURES_DDGLITE=bitget(int32_ulFeatures,11+1);
ulFeatures.AC_FEATURES_FTEXTERNALEXPOSURE=bitget(int32_ulFeatures,12+1);
ulFeatures.AC_FEATURES_KINETICEXTERNALEXPOSURE=bitget(int32_ulFeatures,13+1);
ulFeatures.AC_FEATURES_DACCONTROL=bitget(int32_ulFeatures,14+1);

ulFeatures.AC_FEATURES_METADATA=bitget(int32_ulFeatures,15+1);
ulFeatures.AC_FEATURES_IOCONTROL=bitget(int32_ulFeatures,16+1);
ulFeatures.AC_FEATURES_PHOTONCOUNTING=bitget(int32_ulFeatures,17+1);
ulFeatures.AC_FEATURES_COUNTCONVERT=bitget(int32_ulFeatures,18+1);
ulFeatures.AC_FEATURES_DUALMODE=bitget(int32_ulFeatures,19+1);

cam_skills.ulFeatures=ulFeatures;

%%

%% PCI Card
ulPCICard=int32_ulPCICard;

cam_skills.ulPCICard=ulPCICard;
%% EM Gain Capability
ulEMGainCapability=struct;
ulEMGainCapability.AC_EMGAIN_8BIT=bitget(int32_ulEMGainCapability,1);
ulEMGainCapability.AC_EMGAIN_12BIT=bitget(int32_ulEMGainCapability,2);
ulEMGainCapability.AC_EMGAIN_LINEAR12=bitget(int32_ulEMGainCapability,3);
ulEMGainCapability.AC_EMGAIN_REAL12=bitget(int32_ulEMGainCapability,4);
  

cam_skills.ulEMGainCapability=ulEMGainCapability;

%%

disp('%%%%%%%%%%%%%%%%%%%%%%%%%');
fnames=fieldnames(cam_skills);
for kk=1:length(fnames)
    disp(fnames{kk});
    disp(dec2bin(nums(kk),32));
    
   disp(cam_skills.(fnames{kk})); 
   disp('%%%');
end
        
end
