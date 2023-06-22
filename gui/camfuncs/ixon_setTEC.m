
% Engage/Disengage TEC
function out=ixon_setTEC(state)
    if state
        fprintf('Engaging TEC to cool sensor ... ');
        ret=CoolerON;
    else
        fprintf('Disengaging TEC to cool sensor ... ');
        ret=CoolerOFF;
    end    
    disp(error_code(ret));    
    if isequal(error_code(ret),'DRV_SUCCESS')
        out=1;
    else
        warning('Unable to set TEC.');
        out=0;
    end
end