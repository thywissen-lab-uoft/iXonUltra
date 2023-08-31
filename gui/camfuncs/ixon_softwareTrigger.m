
% Software Trigger
function out=ixon_softwareTrigger
    fprintf('Sending software trigger ... ');
    [ret]=SendSoftwareTrigger;
    disp(error_code(ret));
    if isequal(error_code(ret),'DRV_SUCCESS')
        out=1;
    else
        warning('Unable to send software trigger.');
        out=0;
    end
end