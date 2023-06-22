function [vals,units,flags]=ixon_grabSequenceParams(src)
    if nargin~=1
        src='Y:\_communication\control2.mat';
    end    
    data=load(src);    
    disp(['Opening information from from ' src]);
    vals=data.vals;
    units=data.units;   
    flags=data.flags;
end