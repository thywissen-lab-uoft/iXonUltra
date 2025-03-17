function [data,focus_data,hF] = ixon_multi_shot_focusing2(data,opts)
% Author : CF Fujiwara
%
% This code compares the focusing properties of subsequent images


qgm_MultiExposures: [NaN 300 1500 1500 NaN 300 1500 1500]
%% Default Settings

if nargin==1; opts=struct;end

% Default ROI
% When choosing an ROI, it is important to NOT allow for clipping induced
% by the image plane iris or any other "hard" features.
if ~isfield(opts,'ROI')
    opts.ROI = [150 250 150 250];
end

% Images to analyze
if ~isfield(opts,'ImageIndeces')
    opts.ImageIndeces = [2 3];
end


end

function piezos = getPiezoValues(data)
    piezos = zeros(length(data),2);
    P=[data.Params];
    i1 = 3;
    i2 = 4;

    for kk=1:length(data)
        vals=[P(kk).qgm_MultiPiezos];
        piezos(kk,1)=vals(i1);
        piezos(kk,2)=vals(i2);
    end
end

