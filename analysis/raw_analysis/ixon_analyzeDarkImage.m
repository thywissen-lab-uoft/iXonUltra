function [dark_data] = ixon_analyzeDarkImage(ixondata)

fprintf('Performing dark image analysis ...');    

img_dark_avg = zeros(512,512);

img_dark = zeros(512,512,length(ixondata));

for kk=1:length(ixondata)    
    img_dark(:,:,kk) = ixondata(kk).RawImages(:,:,1);           
end


img_dark_avg = mean(img_dark,3);
img_dark_std = std(img_dark,0,3);

c1 = 1500;
c2 = 200;

v1=linspace(0,c1,100);
v2=linspace(0,c2,100);

figure
subplot(221)
imagesc(img_dark_avg);
title('average of dark image');
colorbar
caxis([0 c1]);

subplot(222)
imagesc(img_dark_std);
title('deviation of dark image');
colorbar
caxis([0 c2]);
subplot(223)
histogram(img_dark_avg,v1);
title('average of dark image');
xlabel('average count');
subplot(224)
histogram(img_dark_std,v2);
title('deviation of dark image');
xlabel('\sigma');

% img_dark_counts = sum(sum(img_dark,1),2);
% img_dark_counts=img_dark_counts(:);

% keyboard

end

