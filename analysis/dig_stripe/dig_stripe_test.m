j = 14;

clear datas;
% for n=1:length(Zbin)
    for n=j:j

    disp(n);
    if sum(Zraw(:,:,n),'all')<1e6
        % keyboard
        continue
    end
    Zb=Zbin{n};
    n1 = 1:size(Zb,2);
    n2 = 1:size(Zb,1);
    %%
    
    figure(12)
    clf
    imagesc(Zraw(:,:,n));
    set(gca,'ydir','normal');
    colormap bone
    caxis([20 100]);
    
    [data,hF]=ixon_fitStripe_dig(n1,n2,Zb);
    datas(n)=data;
%     exportgraphics(hF,"parabola.gif","Append",true)

    % waitforbuttonpress
end

%%
% 
% d2 = datas;
% 
% % B=[datas.ModDepth];
% 
% % d2(B<0.4) = [];
% 
% % B2=[d2.ModDepth];
% inds = [B2>0.4];
% 
% % d2(inds) = [];
% % 
% % B3 = [d2.ModDepth];
% % i2 = [B3>0.4];
% 
% figure(20);
% clf
% Y = unwrapPhaseTime(1:length([d2.Phase]),[d2.Phase]);
% Y0 = unwrapPhaseTime(1:length([datas.Phase]),[datas.Phase]);
% 
% subplot(221);
% 
% plot(Y/(2*pi));
% hold on
% % plot(Y0/(2*pi));
% subplot(222);
% cla
% plot(B2.*inds,'o');
% % plot(B2);
% 
% hold on
% % plot(B);

% subplot(224);
% plot(B2<0.4);

