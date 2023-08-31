% Get the camera status
function [out,outstr]=ixon_getCameraStatus
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