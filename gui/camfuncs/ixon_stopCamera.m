% Stop Acquisition
function out=ixon_stopCamera
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