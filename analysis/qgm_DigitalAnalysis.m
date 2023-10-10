P = [ixondata.Params];
F = [ixondata.Flags];
U = [ixondata.Units];

LatticeDig = [qgmdata.LatticeDig];

digdata = struct;
digdata.FileNames = {ixondata.Name}';
digdata.xVar = ixon_xVar;
digdata.X = [P.(ixon_xVar)];
digdata.Params = P;
digdata.Units = U;
digdata.Flags = F;
digdata.FitType = 'dig';


%% Recenter all digital images

% Get all bounds on the lattice sites indeces
n1i = min(LatticeDig(1).n1);n1f = max(LatticeDig(1).n1);
n2i = min(LatticeDig(1).n2);n2f = max(LatticeDig(1).n2);

% Find the maximu
for n =1:length(LatticeDig)
    n1i = min([LatticeDig(n).n1 n1i]);
    n1f = max([LatticeDig(n).n1 n1f]);
    n2i = min([LatticeDig(n).n2 n2i]);
    n2f = max([LatticeDig(n).n2 n2f]);
end

% Redfine the lattice vectors
n1 = n1i:n1f;
n2 = n2i:n2f;

% Iterate through all images and center the digitzal in the new lattice
% sites
Zdig_all = zeros(length(n2),length(n1),length(LatticeDig));
for n = 1:length(LatticeDig)
    n1i0 = LatticeDig(n).n1(1); 
    n2i0 = LatticeDig(n).n2(1); 
    dI1 = length(LatticeDig(n).n1);
    dI2 = length(LatticeDig(n).n2);    
    i1 = find(n1==n1i0,1);
    i2 = find(n2==n2i0,1);
    Zthis = zeros(length(n2),length(n1));
    Zthis(i2:(i2+dI2-1),i1:(i1+dI1-1)) = LatticeDig(n).Zdig;
    
    Zdig_all(:,:,n) = Zthis;
end

digdata.n1 = n1;
digdata.n2 = n2;
digdata.Z = Zdig_all;

%% Digital ROI
% digROI = [25 155 100 250];

%% Recalculate Momentse and Atom Numbers (just to be sure)

for n = 1:size(digdata.Z,3)  
        Z = digdata.Z(:,:,n);

%     Zc = zeros(size(Z,1),size(Z,2));
%     Zc(digROI(3):digROI(4),digROI(1):digROI(2)) = 1;%  

%     Z = Z.*Zc;
     
    N = sum(sum(Z));
    n1 = digdata.n1;
    n2 = digdata.n2;
    
    xc = sum(sum(Z,1).*n1)./N;
    yc = sum(sum(Z,2).*n2')./N;    
    x2 = sum(sum(Z,1).*n1.^2)./N;
    y2 = sum(sum(Z,2).*n2'.^2)/N; 
    xs = sqrt(x2 - xc^2);
    ys = sqrt(y2 - yc^2);
    
    digdata.Xc(n) = xc;
    digdata.Yc(n) = yc;
    digdata.Xs(n) = xs;
    digdata.Ys(n) = ys;
    digdata.N(n) = N;  
end
%% Save Out
    filename=fullfile(ixon_imgdir,'figures','digdata.mat');
    save(filename,'digdata');
    
    %%
    
    hf2 = figure;
    hf2.Color='w';
    hf2.WindowStyle='docked';
    
    subplot(321)
    plot([digdata.X],[digdata.Xc],'ko');
    ylabel('xc (sites)');
    
    subplot(322)
    plot([digdata.X],[digdata.Yc],'ko');
    ylabel('yc (sites)');

    subplot(323)
    plot([digdata.X],[digdata.N],'ko');
    ylabel('N (atoms)');

    subplot(324)
    plot([digdata.X],[digdata.Xs],'ko');
    ylabel('\sigma_x)');

            subplot(325)
    plot([digdata.X],[digdata.Ys],'ko');
    ylabel('\sigma_y)');

        %%
%{
%% Fit to gaussian distribution

gaussme=@(A,Xc,Yc,s1,s2,xx,yy) A*exp(-( ...
    (1/(2*s1^2))*(xx-Xc).^2 + ...     
     (1/(2*s2^2))*(yy-Yc).^2));        

myfit=fittype(@(A,Xc,Yc,Xs,Ys,xx,yy) gaussme(A,Xc,Yc,Xs,Ys,xx,yy),...
    'independent',{'xx','yy'},'coefficients',{'A','Xc','Yc','Xs','Ys'});
opt=fitoptions(myfit);
opt.StartPoint=[.3 90 100 30 30];

% Make a mesh grid for fitting
[xx,yy]=meshgrid(n1,n2);    

% Copy the data
xx2=xx;yy2=yy;

% Perform the fit
t1=now;
[fout,gof,output]=fit([xx2(:) yy2(:)],Z2(:),myfit,opt);
t2=now;            

% Fit String
fStr=['(' num2str(round(fout.Xc)) ',' ...
    num2str(round(fout.Yc)) ',' num2str(round(fout.Xs)) ',' num2str(round(fout.Ys)) ',' ...
    num2str(fout.A,'%.e')  ')'];  
disp(fStr);

Zfit = feval(fout,xx,yy);
disp(' ');

% Total number of atoms/counts
Ngauss=fout.A*sqrt(2*pi*fout.Xs^2)*sqrt(2*pi*fout.Ys^2);
%% Summary
           
hf = figure(201);
hf.Color='w';
hf.Position(2:4) =[50 640 600];
clf
ax1 = subplot(5,5,[1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19]); 
imagesc(n1,n2,Zdig_all/length(digdata));
colormap parula
xlabel('x site')
ylabel('y site');
axis equal tight
set(gca,'YAxisLocation','left','XAxisLocation','top','YDir','normal',...
    'XDir','normal');

t=uicontrol(hf,'style','text','string',FigLabel,'FontSize',8,...
    'backgroundcolor','w','horizontalalignment','left');
t.Position(3:4) = t.Extent(3:4);
t.Position(1:2) = [1 hf.Position(4) - t.Extent(4)];

Nimg = length(digdata);

% title(['average digitized image N =' num2str(length(digdata))]);
% c = colorbar;
% c.Label.String = 'filling fraction';

strImg = ['$N_{\mathrm{images}} = ' num2str(Nimg) '$'];
text(.98,.98,strImg,'units','normalized','horizontalalignment','right',...
    'verticalalignment','top','color','w','fontsize',12,'interpreter','latex');

strG  = ['Gauss : ' newline ...
    '$n_0=' num2str(fout.A,3) ', \bar{N}=' num2str(round(Ngauss)) '$' newline ...
    '$(\sigma_{x},\sigma_{y}) = (' num2str(round(fout.Xs)) ',' num2str(round(fout.Ys)) ')$'];
text(.02,.02,strG,'units','normalized','horizontalalignment','left',...
    'verticalalignment','bottom','color','w','fontsize',12,'interpreter','latex');
text(.02,.98,strB,'units','normalized','horizontalalignment','left',...
    'verticalalignment','top','color','w','fontsize',12,'interpreter','latex');


% For JHT
t= 600;
mK = 6.64216e-26;
omega = 2*pi*50;
h=6.62607015e-34;
Rc =sqrt(4*t*h/(mK/2)/(omega^2));
Rc_sites = Rc/(532e-9);

xc = 98;
yc  =93;
hold on
tt=linspace(0,2*pi,1000);
plot(xc+Rc_sites*cos(tt),yc+Rc_sites*sin(tt),'r-','linewidth',2);

str = ['$t(2.5 Er) \approx ' num2str(t) '$ Hz' newline ...
    '$\omega \approx 2\pi\cdot' num2str(omega/(2*pi)) '$ Hz' newline ...
    '$\frac{1}{2} m \omega^2R_c^2=4 h t$' newline ...
    '$R_c = ' num2str(round(Rc_sites,1)) 'a_L$'];

text(.99,.01,str,'units','normalized','interpreter','latex',...
    'horizontalalignment','right','verticalalignment','bottom','color','r',...
    'fontsize',16);


% Y plot
subplot(5,5,[5 10 15 20]);
plot(sum(Zdig_all,2)/length(digdata),n2)
ylim([min(n2) max(n2)]);
set(gca,'ydir','normal');
hold on
plot(sum(Zfit,2),n2,'r-');
set(gca,'YaxisLocation','right');
ylabel('site number');
text(2,2,'sum','units','pixels','verticalalignment','bottom');

% X plot
subplot(5,5,[21 22 23 24]);
plot(n1,sum(Zdig_all,1)/length(digdata))
xlim([min(n1) max(n1)]);
hold on
plot(n1,sum(Zfit,1),'r-');
xlabel('site number');
text(2,2,'sum','units','pixels','verticalalignment','bottom');




if ixon_doSave;ixon_saveFigure2(hf,'ixon_dig_avg',saveOpts);end     




%% Animate
filename='animation_qgm.gif';
hF = figure(301);
clf
hImg = imagesc(n1,n2,Zdig_all);
xlabel('x site');
ylabel('y site');
startDelay = 1;
midDelay = 0.2;
endDelay = .5;

for kk=1:length(digdata)   % Iterate over all unique xvalues   


    set(hImg,'CData',digdata(kk).Zdig)  % Image data
    set(gca,'XDir','normal','YDir','normal');
    
    drawnow % update graphcis
    
    
    % Write the image data
    frame = getframe(hF);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);           

    if kk == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',startDelay);
    else
        if kk==29
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',endDelay);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',midDelay);
        end
    end

end
close;

disp('done');


%%
% x = [params.(xVar)];
% n = [digdata.Natoms];
% xc = [digdata.Xc];
% yc = [digdata.Yc];
% xs = [digdata.Xs];
% ys = [digdata.Ys];
% 
% figure(200);
% clf
% 
% subplot(3,2,1);
% histogram(n,7);
% xlabel('atom number');
% 
% subplot(3,2,3);
% histogram(xc,7);
% xlabel('x center (sites)');
% 
% subplot(3,2,4);
% histogram(xs,7);
% xlabel('x size (sites)');
% 
% subplot(3,2,5);
% histogram(yc,7);
% xlabel('y center (sites)');
% 
% subplot(3,2,6);
% histogram(ys,7);
% xlabel('y sizes (sites)');

%}