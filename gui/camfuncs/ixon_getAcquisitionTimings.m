
% Get Timing
function out=ixon_getAcquisitionTimings    
    [ret,texp,taccum,tkin] = GetAcquisitionTimings;
    
    if isequal(error_code(ret),'DRV_SUCCESS')
        out = [texp taccum tkin];
    else
        warning('Unable to read iXon timing.');
    end
end