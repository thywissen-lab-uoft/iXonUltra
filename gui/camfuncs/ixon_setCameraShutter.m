% Set the camera shutter
function out=ixon_setCameraShutter(state)
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