function [hFs,out] = dig_radialAnalysis_average_images(digdata,opts)
    if nargin == 1
        opts = struct;
    end    
 
    X = digdata.X;    
    
    
    if opts.ForceAverage
       Ux = 1;
    else
        Ux = unique(X);
    end
    
    animateFileName = fullfile(opts.saveDir,'radial_animation');            

    
            
        
        
            
        


    
  
    
    %%
for kk = 1:length(Ux)
   if opts.ForceAverage         
        inds = 1:length(X);
   else
        x = Ux(kk);
        inds = find([X==x]);
   end

   nImg = length(inds);                          % Number of images
   Z = mean(digdata.Zdig(:,:,inds),3);          % Average Image for this X variable
   Zr = mean(digdata.Zr(:,inds),2);             % Average Radial data for this X variable
   err = std(digdata.Zr(:,inds),1,2)/sqrt(nImg);  % Standard error at each radius 
   r = mean(digdata.r(:,inds),2);               % Radial vector

   n1 = digdata.n1-mean(digdata.rCenter(:,1));
   n2 = digdata.n2-mean(digdata.rCenter(:,2));

    hFs(kk) = figure;
    hFs(kk).Color='w';
    hFs(kk).Position = [300 100 1200 350];
    W= hFs(kk).Position(3);
    H = hFs(kk).Position(4);
    clf

    if isfield(opts,'FigLabel') && ~isempty(opts.FigLabel)
        tFig=uicontrol('style','text','string',opts.FigLabel,...
            'units','pixels','backgroundcolor',...
            'w','horizontalalignment','left');
        tFig.Position(4)=tFig.Extent(4);
        tFig.Position(3)=tFig.Extent(3);
        tFig.Position(1:2)=[1 1];
    end    

    % 2D Image of Average Image
    ax1=subplot(131);
    ax1.Units='pixels';
    [x0, y0, w, h]=getAxesPos(1,3,1,W,H);
    ax1.Position = [x0 y0 w h];

    imagesc(n1,n2,imboxfilt(Z,opts.BinStep))
    xlabel('position (sites)');
    ylabel('position (sites)');
    
%     xlim([-1 1]*max(digdata.r);

    strBoxCar = ['boxcar avg: ' num2str(opts.BinStep) '\times' num2str(opts.BinStep)];

    strSummary = ['$n=' num2str(nImg) ' $ images' newline ...
        '$N=' num2str(round(mean(digdata.Natoms(inds)))) '\pm' ...        
        num2str(round(std(digdata.Natoms(inds),1))) '$ atoms'];    

    text(1,1,strSummary,'units','pixels','fontsize',10,...
        'horizontalalignment','left','verticalalignment','bottom',...
        'color','r','interpreter','latex');  

    text(.01,.99,strBoxCar,'units','normalized','fontsize',10,...
        'verticalalignment','top','horizontalalignment','left',...
        'color','r');
        ax1.XAxisLocation='Top';
        set(ax1,'FontSize',8);
        
        
        
    if ~opts.ForceAverage        
        str= [digdata.xVar ': ' num2str(x)];
        text(.99,.99,str,'units','normalized','fontsize',8,'color','r',...
            'horizontalalignment','right','verticalalignment','top',...
            'interpreter','none');            
    end
   
           

    colormap bone
    axis equal tight
    cc1=colorbar;
    cc1.Label.String = 'charge occupation';
    title('average image');
    set(gca,'ydir','normal','box','on','linewidth',1);

    if isfield(opts,'nMaxShow')
       caxis([0 opts.nMaxShow]); 
    end

    
    if isfield(opts,'rMaxShow')
        xlim([-1 1]*opts.rMaxShow); 
        ylim([-1 1]*opts.rMaxShow); 

    end

    % Average charge density
    ax2=subplot(132);
    ax2.Units='pixels';
    [x0, y0, w, h]=getAxesPos(1,3,2,W,H);
    ax2.Position = [x0 y0 w h];

    errorbar(r(2:end),Zr(2:end),err(2:end),...
        'ko','markerfacecolor',[.5 .5 .5],...
        'markersize',8,'linewidth',1);
    hold on      
    xlabel('radial distance (sites)');
    ylabel('mean $\{N(r)\}$','interpreter','latex','fontsize',14);    
    title('average radial profile');        
    text(.99,.99,strSummary,'units','normalized','fontsize',12,'horizontalalignment','right',...
        'verticalalignment','top','interpreter','latex');

    if isfield(opts,'nMaxShow')
       ylim([0 opts.nMaxShow]); 
    end

    if isfield(opts,'rMaxShow')
       xlim([0 opts.rMaxShow]); 
    end

    strRadialBin = ['radial bin ' char(916) 'r:' num2str(opts.BinStep)];
    text(.01,.99,strRadialBin,'units','normalized','horizontalalignment','left',...
        'verticalalignment','top','fontsize',10);

    % Plot radial average ndet std
    ax3=subplot(133);
    ax3.Units='pixels';
    [x0, y0, w, h]=getAxesPos(1,3,3,W,H);
    ax3.Position = [x0 y0 w h];

    plot(r(2:end),err(2:end),...
        'ko','markerfacecolor',[.5 .5 .5],...
        'markersize',8,'linewidth',1);
    title('standard error');
    ylabel('std $\{N(r)\}/\sqrt{n}$','interpreter','latex','fontsize',14);
    xlabel('radial distance (sites)');    
    text(.99,.99,strSummary,'units','normalized','fontsize',12,'horizontalalignment','right',...
        'verticalalignment','top','interpreter','latex');

    if isfield(opts,'rMaxShow')
       xlim([0 opts.rMaxShow]); 
    end
    text(.01,.99,strRadialBin,'units','normalized','horizontalalignment','left',...
        'verticalalignment','top','fontsize',10);

    
    if isfield(opts,'doAnimate') && opts.doAnimate
        frame=getframe(hFs(kk));
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256);              
        switch kk
            case 1
                imwrite(A,map,animateFileName,'gif','LoopCount',...
                    Inf,'DelayTime',1);
            case length(Ux)
                imwrite(A,map,animateFileName,'gif','WriteMode',...
                    'append','DelayTime',1);
            otherwise
                imwrite(A,map,animateFileName,'gif','WriteMode',...
                    'append','DelayTime',.2);
         end
    end
end

  
    
    
    
%     %% Assign to output
%     
%     out = struct;
%     out.CroppedImages                  = Zsub;
%     out.AverageImage                   = ZsubBar;
%     out.SiteVector1                    = r(1):r(2);
%     out.SiteVector2                    = r(3):r(4);
%     out.RadialCenter                   = [Xcbar Ycbar];
%     out.RadialVector                   = rVec;
%     out.PotentialVector                = PotentialVector;
%     out.AverageOccupation              = radial_charge_mean;
%     out.AverageOccupationUncertainty   = radial_charge_mean_std/sqrt(nImg);

%% Plotting


end

function [Tics,Average,dev,n]=radial_profile(data,radial_step)
%main axii cpecified:
x=(1:size(data,2))-size(data,2)/2;
y=(1:size(data,1))-size(data,1)/2;
% coordinate grid:
[X,Y]=meshgrid(x,y);
% creating circular layers
Z_integer=round(abs(X+1i*Y)/radial_step)+1;
% % illustrating the principle:
% % figure;imagesc(Z_integer.*data)
% very fast MatLab calculations:
Tics=accumarray(Z_integer(:),abs(X(:)+1i*Y(:)),[],@mean);
Average=accumarray(Z_integer(:),data(:),[],@mean);


% dev=accumarray(Z_integer(:),data(:),[],@std);
dev=accumarray(Z_integer(:),data(:),[],@(x) std(x,1));

n= accumarray(Z_integer(:),data(:),[],@(x) numel(x));

end


function [axX,axY,axWidth,axHeight]=getAxesPos(num_rows,num_cols,index,figW,figH)
yTop=30;
yBot=50;

xLeft=40;
xRight=20;

ySpace=35;
xSpace=70;

rowNumber = floor((index-1)/num_cols)+1;
colNumber = mod(index-1,num_cols)+1;

axHeight = (figH-yTop-yBot-ySpace*(num_rows-1))/num_rows;
axWidth = (figW-xLeft-xRight-xSpace*(num_cols-1))/num_cols;

axX = xLeft + (axWidth+xSpace)*(colNumber-1);
axY=(figH-yTop-axHeight)-(rowNumber-1)*(axHeight+ySpace);

end

