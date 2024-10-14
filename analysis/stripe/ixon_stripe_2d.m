% %% Sort Data
% 
% params=[ixondata.Params];
% xvals=[params.(ixon_xVar)];
% [xvals,inds]=sort(xvals,'ascend');
% 
% % [xvals,inds]=sort(xvals,'descend');
% 
% ixondata=ixondata(inds);
% 
% %% 2D Stripe Analysis
% % This analzyes the stripes for take from the ixon camera.  This is a 2D
% % fit over the entire cloud.  This fits a 2D gaussian modulated by a sine
% % wave at a particular angle.  This is pariticularly useful to fit the
% % angular dependence of the data. 
% %
% % The input fit parameters are specified in the options structure.
% 
% 
% stripe_2d_opts=struct;
% 
% stripe_2d_opts.xUnit=ixon_unit;
% 
% stripe_2d_opts.ShimFit=0;
% stripe_2d_opts.Theta=[10 190]; % Specify the domain (MUST BE 180 DEGREES)
% 
% stripe_2d_opts.saveAnimation=1;        % save the animation?
% stripe_2d_opts.StartDelay=1;
% stripe_2d_opts.MidDelay=.25;
% stripe_2d_opts.EndDelay=1;
% stripe_2d_opts.CLim = [0 100];
% 
% stripe_2d_opts.ROI = NaN;
% % stripe_2d_opts.ROI = [100 450 25 425];
% 
% stripe_2d_opts.Threshhhold = 20; % Ignore pixels below this threshhhold
% stripe_2d_opts.Threshhhold = NaN;
% 
% stripe_2d_opts.ConstrainWavelength = NaN;
% % stripe_2d_opts.ConstrainWavelength = 72;
% 
% stripe_2d_opts.ConstrainAngle = NaN;
% %  stripe_2d_opts.ConstrainAngle = 88;
% 
% ixondata=flip(ixondata);
% if ixon_doAnalyzeStripes2D
%     [hF_stripe_2d,stripe_data2d]=analyzeStripes2(ixondata,ixon_xVar,stripe_2d_opts);
% 
%     if ixon_doSave
%         ixon_saveFigure(ixondata,hF_stripe_2d,'ixon_field_stripe_2d');        
%     end
% 
%     outdata.stripe_data2d=stripe_data2d;   
% end

%% 2D Stripe Analysis
if ixon_doAnalyzeStripes2D
  
    stripe_2d_opts=struct;
    stripe_2d_opts.xVar = ixon_xVar;
    stripe_2d_opts.xUnit=ixon_unit;
    stripe_2d_opts.ShimFit=0;
    stripe_2d_opts.Theta=[10 190]; % Specify the domain (MUST BE 180 DEGREES)
    
    stripe_2d_opts.Theta=[-70 110]; % Specify the domain (MUST BE 180 DEGREES)
%     stripe_2d_opts.Theta=[85 95]; % Specify the domain (MUST BE 180 DEGREES)
%     stripe_2d_opts.Theta=[45 135]; % Specify the domain (MUST BE 180 DEGREES)
  

    stripe_2d_opts.saveAnimation=1;        % save the animation?
    stripe_2d_opts.StartDelay=1;
    stripe_2d_opts.MidDelay=.1;
    stripe_2d_opts.EndDelay=1;
    stripe_2d_opts.CLim = [0 100];
    
    params=[ixondata.Params];
    xvals=[params.(ixon_xVar)];
    [xvals,inds]=sort(xvals,'ascend');
    
    % Analyze the stripes
    stripes = analyzeStripes3(ixondata,ixon_xVar,stripe_2d_opts);   
end
    %%
if ixon_doAnalyzeStripes2D
stripes_modify = stripes;
Navg = 1;
phi_unwrap = unwrapPhaseTime(xvals,[stripes_modify.phi],Navg);
     for kk=1:length(stripes_modify)
         stripes_modify(kk).phi = phi_unwrap(kk);
        stripes(kk).phi = mod(stripes(kk).phi+pi/2,(2*pi))-pi/2;
     end

    hF_stripe_summary = ixon_stripeSummary(stripes_modify,xvals,stripe_2d_opts);

    if ixon_doSave
        ixon_saveFigure(ixondata,hF_stripe_summary,'ixon_stripe_summary');        
    end
end

%% Stripe focusing
ixon_doFocusStripes = 0;
if ixon_doFocusStripes
    clear focus
    for kk=1:length(ixondata)

        x = ixondata(kk).X;
        y = ixondata(kk).Y;
        z = ixondata(kk).ZNoFilter

        focus(kk) = ixon_focusStripe(x,y,z,stripes(kk));
    end
    
    hFme = figure;
    hFme.Color='w';
    plot(xvals,[focus.y0],'ko');
    xlabel('shim current (mA)')
    ylabel('focus position (px)');  
end

%% Stripe COM
ixon_doStripeCOM = 0;
if ixon_doStripeCOM

    [ixondata,stripes,qgmdata_stripes] = ixon_stripeCOM2(ixondata,stripes);
    
    % Save the output data
    filename=fullfile(ixon_imgdir,'figures','qgmdata_stripes.mat');
    save(filename,'qgmdata_stripes');

end

%% Stability

ixon_doStripeStability = 0;

if ixon_doStripeStability 
   hFme = figure;
   hFme.Color='w';
   hFme.Position=[100 500 1400 300];
   
   P = [ixondata.Params];
   xvals = [P.ExecutionDate];
   
   phi0 = 0.6*2*pi;
   
   phi_err = [stripes.phi_err];
   phi = [stripes.phi];
   p = [ixondata.Params];
   iz = [p.qgm_plane_tilt_dIz];
   
   phi_unwrap = (phi/(2*pi)-round((phi-phi0)/(2*pi)))*2*pi;
   
   ax1=axes('units','normalized');
   ax1.Position = [.05 .15 .7 .8];
      co=get(gca,'colororder');

   yyaxis left
   errorbar(xvals,phi_unwrap/(2*pi),phi_err,'o-','markerfacecolor',co(1,:),...
       'markeredgecolor',co(1,:)*.5,'color',co(1,:),'linewidth',1,...
       'markersize',8);
   datetick('x');
   xlabel('time');
   ylabel('phase (2\pi)');
   
   
   
   hold on
   plot(get(gca,'XLim'),[1 1]*phi0/(2*pi),'k-');
      set(gca,'box','on','linewidth',1,'fontsize',10);
    ylim(phi0/(2*pi)+[-.5 .5]);
   yyaxis right
   plot(xvals,iz*1e3,'o-','markerfacecolor',co(2,:),...
       'markeredgecolor',co(2,:)*.5,'color',co(2,:),'linewidth',1,...
       'markersize',8);
   ylabel('current shift (mA)');
   
  ax2=axes('units','normalized');
    ax2.Position = [ax1.Position(1)+ax1.Position(3)+.05 ...
        ax1.Position(2) .17 ax1.Position(4)];

   h = histogram(phi_unwrap/(2*pi),20);
   bw = h.BinWidth;
   
   pdg = fitdist(phi_unwrap'/(2*pi),'normal');
   hold on
   
   tt=linspace(-5,5,1e3)+phi0/(2*pi);
   plot(tt,normpdf(tt,pdg.mu,pdg.sigma)*bw*numel(phi_unwrap),'r-','linewidth',2);
   
   str = ['$' num2str(round(pdg.mu,3)) '\pm' num2str(round(pdg.sigma,3)) '$'];
   
   text(.01,.99,str,'units','normalized','fontsize',12,'interpreter','latex',...
       'verticalalignment','top');
    
    
    
     set(gca,'yaxislocation','right');
     xlabel('phase (2\pi)');
     ylabel('occurences','fontsize',8);
    xlim(phi0/(2*pi)+[-.5 .5]);
    
   if ixon_doSave
        ixon_saveFigure(ixondata,hFme,'ixon_stripe_stability');        
    end

end

%% Stripe Analysis FFT
ixon_doAnalyzeStripes2D_Focusing=0;
if ixon_doAnalyzeStripes2D_Focusing
    yfoci = zeros(length(ixondata),1);
    thetas = stripe_data2d.Theta(:,1);
    phis =stripe_data2d.Phi(:,1);
    Ls = stripe_data2d.Wavelength(:,1);
            atomcoms=[];

    for kk=1  :length(ixondata)
%     for kk=length(ixondata)  :length(ixondata)
%             for kk=16:16
%             for kk=11:11

        Z = ixondata(kk).ZNoFilter;
        x=ixondata(kk).X;x=x';
        y=ixondata(kk).Y;y=y';    
        
        dX = x(2)-x(1);
        f_max = 1/dX;
        
        [xx,yy]=meshgrid(x,y);
        theta = thetas(kk);
        phi = phis(kk);
        L = Ls(kk);
            
        phaseMap = 2*pi/L*(cosd(theta)*xx+sind(theta)*yy)+phi+pi/2;
        
        pL_N = floor(min(min(phaseMap))/(2*pi))+1;
        pH_N = floor(max(max(phaseMap))/(2*pi))-1;
        sVec = [];
        Nvec = [];
        metric=[];
        myzfs = {};
        ycoms=[];
        atomcoms(end+1) = sum(Z.*yy,'all')/sum(Z,'all');


        nVec = pL_N:pH_N;
        ps={};
        for nn=pL_N:pH_N
            
            ii=[abs(phaseMap-(nn*2*pi))<.1];
            ps{end+1}=polyfit(xx(ii),yy(ii),1);
            
            i1 = (phaseMap>=(nn*2*pi));
            i2 = (phaseMap<=((nn+1)*2*pi));            
            this_map=logical(i1.*i2);
            
            
            Z_onestripe = Z.*this_map;
            Nfft = 2^10;
            
            ycoms(end+1) = sum(Z_onestripe.*yy,'all')/sum(Z_onestripe,'all');

            zf = fft2(Z_onestripe,Nfft,Nfft);
            zf = fftshift(zf);       
            
            f = 1/2*linspace(-f_max,f_max,size(zf,1));    
            zfnorm = abs(zf);
            
            [fxx,fyy]=meshgrid(f,f);
            fmag = sqrt(fxx.^2+fyy.^2);
            iMask = fmag<=.1;
            
            zfnorm(iMask)=0;
            
            zfnorm=zfnorm/sum(sum(zfnorm));                 
            [f1,f2] = meshgrid(f,f);            
            
            v1 = sum(sum(zfnorm.*f1.^2));
            v2 = sum(sum(zfnorm.*f1.^2));            
            s = sqrt(v1+v2);     
            
%             score = brisque(Z_onestripe/sum(sum(Z_onestripe))*1e9);            
%             metric(end+1) = score;
            myzfs{end+1}=zfnorm;
            sVec(end+1)=s;
            Nvec(end+1)=sum(sum(Z_onestripe));
        end
        
        
        hF = figure(917);
        clf
        hF.Color='w';
        co=get(gca,'colororder');

        subplot(131);
        imagesc(x,y,Z);
        set(gca,'ydir','normal');
        caxis([0 100])
        
        t = [0 512];
        hold on
        for i=1:length(nVec)
            p = ps{i};
            plot(t,p(1)*t+p(2),'r-');
            text(10,p(1)*10+p(2)+10,num2str(nVec(i)),'color','w');
        end
        

        subplot(132);
        yyaxis left
        plot(nVec,sVec,'o-','color',co(1,:),'markerfacecolor',co(1,:),...
            'markeredgecolor',co(1,:)*.5,'linewidth',2,'markersize',8);
        set(gca,'YColor',co(1,:));
        ylabel('fft second moment size');
        
        yyaxis right
        plot(nVec,Nvec,'o-','color',co(2,:),'markerfacecolor',co(2,:),...
            'markeredgecolor',co(2,:)*.5,'linewidth',2,'markersize',8);
        ylabel('counts');
        xlabel('stripe index');
        set(gca,'YColor',co(2,:));

        drawnow;
        
         subplot(133);
         plot(ycoms,sVec,'o-');
         xlabel('stripe com (px)');
         ylabel('fft second moment');
         
         [val,ind]=min(sVec);
         hold on
         switch ind
             case 1
                 yfocus = ycoms(ind);
             case length(sVec)
                 yfocus = ycoms(ind);
             otherwise
                 A = [ycoms(ind-1)^2 ycoms(ind-1) 1;
                     ycoms(ind)^2 ycoms(ind) 1;
                     ycoms(ind+1)^2 ycoms(ind+1) 1];
                 b = [sVec(ind-1); sVec(ind); sVec(ind+1)];
                 
                 coeff_vec = A\b;
                 
                 foo = @(y) coeff_vec(1)*y.^2+coeff_vec(2)*y+coeff_vec(3);
                 yvec = linspace(ycoms(ind-1),ycoms(ind+1),100);
                 plot(yvec,foo(yvec),'r-');
                 yfocus = -coeff_vec(2)/(2*coeff_vec(1));
                 
                 
                 
         end
        s=[num2str(round(yfocus))];
        text(.1,.1,s,'units','normalized');
        yfoci(kk)=yfocus;           
        dosave=1;

        if dosave 
             % Write the image data
            frame = getframe(hF);
            im = frame2im(frame);
            [A,map] = rgb2ind(im,256);      
            filename = 'focus_test.gif';
            switch kk
                case 1
                    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',.5);
                case length(ixondata)
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.5);
                otherwise
                    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',.2);
            end
        end

    end
    
    figure
    pp=[ixondata.Params];
    xvals=[pp.qgm_plane_tilt_dIz];
    p1=plot(xvals,yfoci(:,1));
    hold on
    p2=plot(xvals,atomcoms);
    legend([p1,p2],{'focus','atom com'});
    xlabel('shim current (A)');
    ylabel('position (px)');
end