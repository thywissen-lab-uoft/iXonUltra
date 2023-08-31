
% Set the temperature setpoint
function out=ixon_setTemperature(temp)
    fprintf(['Changing temperature set point to ' num2str(temp) ' C ...']);
    ret=SetTemperature(temp);
    disp(error_code(ret))    
    if isequal(error_code(ret),'DRV_SUCCESS') 
        out=1;
    else
        warning('Unable to change iXon temperature set point.');
        out=0;
    end
end