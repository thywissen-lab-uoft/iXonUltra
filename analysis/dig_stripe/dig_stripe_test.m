n = 12;
% n = 5;
clear datas;
for n=1:length(Zbin)
    disp(n);
    if sum(Zraw(:,:,n),'all')<1e6
        % keyboard
        continue
    end

    Zb=Zbin{n};
    n1 = 1:size(Zb,2);
    n2 = 1:size(Zb,1);
    
    Zr = Zraw(:,:,n);
    x = linspace(0,512,size(Zr,2));
    y = linspace(0,512,size(Zr,2));
    % 
    f=figure(20);
    clf
    f.Color='w';

    subplot(121);
    imagesc(x,y,Zr);
    set(gca,'fontsize',12,'ydir','normal')
    colorbar
    caxis([0 200]);
    axis equal tight
    colormap(bone)
    xlabel('px');ylabel('px')

    subplot(122);
    imagesc(n1,n2,Zb);
    set(gca,'fontsize',12,'ydir','normal')
    colorbar
    caxis([0 3000]);
    axis equal tight
    xlabel('site 1')
    ylabel('site 2')
    
    %%
     [data,hF]=ixon_fitStripe_dig(n1,n2,Zb);
    datas(n)=data;

         % exportgraphics(hF,"parabola.gif","Append",true)

    % waitforbuttonpress
end