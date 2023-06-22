% Disconnect from the Andor iXon Ultra
function out=ixon_disconnectCam
    disp('Disconnecting from the iXon camera.');
    
    % Close the shutter
    fprintf('Closing the shutter ... ');
    [ret]=SetShutter(1,2,0,0);  
    disp(error_code(ret));
    
    % Shut down cooler
    fprintf('Turning off cooler ... ');
    [ret]=SetCoolerMode(1);     
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