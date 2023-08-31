% Get the temperature
function [out,temp,outstr]=ixon_getTemperature
    [ret,temp]=GetTemperature;
    out=1;
    outstr=error_code(ret);
end