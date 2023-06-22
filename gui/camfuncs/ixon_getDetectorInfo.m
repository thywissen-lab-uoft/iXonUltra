% Get Detector
function [out,xpx,ypx]=ixon_GetDetectorInfo    
    [ret,xpx,ypx]=GetDetector;
    
    if isequal(error_code(ret),'DRV_SUCCESS')
        out=1;
    else
        warning('Unable to get detector informaton.');
        out=0;
    end
end