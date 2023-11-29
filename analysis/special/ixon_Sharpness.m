function [data] = ixon_Sharpness(data)
% tic
% 
mag = 80;
pxsize = 16;
latt_spacing = 0.532;

site_per_px = pxsize/mag/latt_spacing;

if data(1).ProcessOptions.doScale
    sc = data(1).ProcessOptions.ScaleFactor;
else
    sc = 1;
end

eff_site_per_px = site_per_px/sc;

for kk=1:length(data)
    c1 = find(data(kk).X>=data(kk).ROI(1),1);
    c2 = find(data(kk).X>=data(kk).ROI(2),1);
    r1 = find(data(kk).Y>=data(kk).ROI(3),1);
    r2 = find(data(kk).Y>=data(kk).ROI(4),1);    
    z = data(kk).Z(r1:r2,c1:c2,:);
    znofilter = data(kk).ZNoFilter(r1:r2,c1:c2,:);  
    for nn=1:size(z,3)
        zme1 = z(:,:,nn);
        zme2 = znofilter(:,:,nn);
        % [s, m1] =MLVSharpnessMeasure(zme1);
        % [s, m2] =MLVSharpnessMeasure(zme2);
        zeff = imresize(zme1,eff_site_per_px);
        zeff2 = imresize(zme2,eff_site_per_px);
        
        [s1, m1] =MLVSharpnessMeasure(zeff/sum(zeff,'all'));
        [s2, m2] =MLVSharpnessMeasure(zeff2/sum(zeff2,'all'));
        data(kk).SharpnessScore(nn) = s1;
        data(kk).SharpnessScoreNoFilter(nn) = s2;
    end
end


end

