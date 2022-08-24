function [hF, dark_data] = ixon_analyzeDarkImage(ixondata,opts_darkImage)
dark_data =struct;

c1 = opts_darkImage.AvgLimits;
c2 = opts_darkImage.StdLimits;

fprintf('Performing dark image analysis ...');    

img_dark_avg = zeros(512,512);

img_dark = zeros(512,512,length(ixondata));

for kk=1:length(ixondata)    
    img_dark(:,:,kk) = ixondata(kk).RawImages(:,:,1);           
end


img_dark_avg = mean(img_dark,3);
img_dark_std = std(img_dark,0,3);

img_dark_tot = sum(sum(img_dark,1),2);
img_dark_tot = img_dark_tot(:);


% v1=linspace(c1(1),c1(2),100);
v1=c1(1):1:c1(2);
v2=c2(1):1:c2(2);
% v2=linspace(c2(1),c2(2),100);

hF = figure;
hF.Color='w';
hF.Position = [10 50 600 450];
subplot(221)
imagesc(img_dark_avg);
title('average of dark image');
colorbar
caxis(c1);
axis equal tight

subplot(222)
imagesc(img_dark_std);
title('deviation of dark image');
colorbar
caxis(c2);
axis equal tight

subplot(223)
histogram(img_dark_avg,v1);
title('average of dark image');
xlabel('average of counts');

subplot(224)
histogram(img_dark_std,v2);
title('deviation of dark image');
xlabel('deviation of counts');


hF = figure;
hF.Color='w';
hF.Position = [10 50 600 450];
plot(img_dark_tot);


disp('done');
end

