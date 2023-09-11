digdata = [qgmdata.LatticeDig];
params = [qgmdata.Params];

xVar = 'ExecutionDate';
%%
n1i = min(digdata(1).n1);
n1f = max(digdata(1).n1);

n2i = min(digdata(1).n2);
n2f = max(digdata(1).n2);

for n =1:length(digdata)
    n1i = min([digdata(n).n1 n1i]);
    n1f = max([digdata(n).n1 n1f]);
    n2i = min([digdata(n).n2 n2i]);
    n2f = max([digdata(n).n2 n2f]);
end

n1 = n1i:n1f;
n2 = n2i:n2f;

Zdig_tot = zeros(length(n2),length(n1));
for n = 1:length(digdata)
    n1i0 = digdata(n).n1(1); 
    n2i0 = digdata(n).n2(1); 
    dI1 = length(digdata(n).n1);
    dI2 = length(digdata(n).n2);    
    i1 = find(n1==n1i0,1);
    i2 = find(n2==n2i0,1);
    Zthis = zeros(length(n2),length(n1));
    Zthis(i2:(i2+dI2-1),i1:(i1+dI1-1)) = digdata(n).Zdig;
    Zdig_tot = Zdig_tot +  Zthis;
end
%%


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
    Z2=Zdig_tot/length(digdata);xx2=xx;yy2=yy;
%         
%         
%     % Display initial guess            
%     gStr=[' guess (Xc,Yc,Xs,Ys,A,bg)=(' num2str(round(Xc)) ',' ...
%         num2str(round(Yc)) ',' num2str(round(s1)) ',' num2str(round(s2)) ',' ...
%         num2str(A,'%.e') ',' ')' ];     
%     fprintf([gStr ' ... ']);

    % Perform the fit
    t1=now;
    [fout,gof,output]=fit([xx2(:) yy2(:)],Z2(:),myfit,opt);
    t2=now;            

    % Fit String
    fStr=['(' num2str(round(fout.Xc)) ',' ...
        num2str(round(fout.Yc)) ',' num2str(round(fout.Xs)) ',' num2str(round(fout.Ys)) ',' ...
        num2str(fout.A,'%.e')  ];     
    fprintf([' fit ' fStr]);
                        Zfit = feval(fout,xx,yy);
                        disp(' ');
            %%
            

figure(201)
clf
subplot(5,5,[1 2 3 4 6 7 8 9 11 12 13 14 16 17 18 19]); 
imagesc(n1,n2,Zdig_tot/length(digdata));
colormap parula
xlabel('x site')
ylabel('y site');
title('average digitized image N=29');
c = colorbar;
c.Label.String = 'filling fraction';
N=fout.A*sqrt(2*pi*fout.Xs^2)*sqrt(2*pi*fout.Ys^2);

str  = ['n_0=' num2str(fout.A,3) ', N=' num2str(round(N)) newline ...
    'sx = ' num2str(round(fout.Xs)) ', ' 'sy = ' num2str(round(fout.Ys))];

text(.02,.02,str,'units','normalized','horizontalalignment','left',...
    'verticalalignment','bottom','color','w','fontsize',14);

subplot(5,5,[5 10 15 20]);
i1 = find(n1==901);
plot(sum(Zdig_tot,2)/length(digdata),n2)
ylim([min(n2) max(n2)]);
set(gca,'ydir','reverse');
hold on
plot(sum(Zfit,2),n2,'r-');

subplot(5,5,[21 22 23 24]);
i2 = find(n2==100,1);
plot(n1,sum(Zdig_tot,1)/length(digdata))
xlim([min(n1) max(n1)]);
hold on
plot(n1,sum(Zfit,1),'r-');

%%


%% Animate
filename='animation_qgm.gif';
hF = figure(301);
clf
hImg = imagesc(n1,n2,Zdig_tot);
xlabel('x site');
ylabel('y site');
startDelay = 1;
midDelay = 0.2;
endDelay = .5;

for kk=1:length(digdata)   % Iterate over all unique xvalues
    
    %%%% Update the graphics
    
%     if isequal(xVar,'ExecutionDate')
%             t.String=[xVar ': ' datestr(uxvals(kk),'YYYY-mm-DD_HH-MM-SS')];          % Variable string
%     else
%             t.String=[xVar ': ' num2str(uxvals(kk))];          % Variable string
%     end


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
x = [params.(xVar)];
n = [digdata.Natoms];
xc = [digdata.Xc];
yc = [digdata.Yc];
xs = [digdata.Xs];
ys = [digdata.Ys];

figure(200);
clf

subplot(3,2,1);
histogram(n,7);
xlabel('atom number');

subplot(3,2,3);
histogram(xc,7);
xlabel('x center (sites)');

subplot(3,2,4);
histogram(xs,7);
xlabel('x size (sites)');

subplot(3,2,5);
histogram(yc,7);
xlabel('y center (sites)');

subplot(3,2,6);
histogram(ys,7);
xlabel('y sizes (sites)');
