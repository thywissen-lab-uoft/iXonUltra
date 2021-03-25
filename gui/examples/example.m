%% Example Image
% Author : C Fujiwara
%
% This is a minimal example meant to connect to the Andor camera.

% Add the Andor MATLAB drivers to the MATLAB path. You need to have
% installed the MATLAB drivers in order for this to work.
addpath(genpath(fullfile(matlabroot,'toolbox','Andor')));

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath))

% Find available cameras
[ret,ncams]=GetAvailableCameras();
disp(['Detected ' num2str(ncams) ' valid cameras.']);

% Initialize the camera
currDir = pwd;
fileDir = fileparts(mfilename('fullpath'));
cd(fileDir);
fprintf('Connecting to Andor camera ... ');
[ret] = AndorInitialize('');
disp(error_code(ret))
cd(currDir);

% Acquire camera information
cam = getCamInfo;

% Shutdown the camera
[ret]=SetCoolerMode(1);     % Shut down cooler
[ret]=SetShutter(1,2,0,0);  % Close the shutter
[ret]=AndorShutDown;        % Shut down the camera

