
% Start Acquisition
function out=ixon_startCamera
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