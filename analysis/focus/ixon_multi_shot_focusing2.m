function [data,focus_data,hF] = ixon_multi_shot_focusing2(data,opts)
% Author : CF Fujiwara
%
% This code compares the focusing properties of subsequent images
% 
% It assume that the image stack comprises of three images where image 2
% and image 3 are taken at different focus positions of the microscope.  We
% do not use the 1st image as we allow for atoms leave to exit the trap
% just in case they are not properly cooled. While the analysis compares
% the normalized count profile between the images, it works best if the
% light source (ie. atomic distribution) is identical from shot to shot.

%% Helper functions

% In variable exposure mode, there are a total of eight exposures.
% [Wipe Image1 Image2 Image3 Wipe Image4 Image5 Image6]
%
% As the code is updated, this function may need to change.
function piezos = getPiezoValues(data)
    piezos = zeros(length(data),2);
    P=[data.Params];
    i1 = 3;                     % Image index
    i2 = 4;                     % Image Index
    varName= 'qgm_MultiPiezos'; % This variable stores the piezo values, it is an array
    for kk=1:length(data)
        vals=[P(kk).(varName)]; 
        piezos(kk,1)=vals(i1);
        piezos(kk,2)=vals(i2);
    end
end

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



