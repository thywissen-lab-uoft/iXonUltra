function data=ixon_fft_multi_shot_focusing(data,opts)
% Author : Cora Fujiwara
%
% This code assume that multiple images are taken of a single plane at
% variable focus positions.  From the different images it looks for the
% most in focus image.
%
% The code functions by placing a notch filter in the momentum domain to
% look for features at a momentum correpsonding to a magnitude at the
% lattice spacing.  When the image is in focus, these are maximized.

%% Default Options

if nargin==1;opts=struct;end

if ~isfield(opts,'kmag')
    opts.kmag = 0.3727;
end

if ~isfield(opts,'kdelta')
    opts.kdelta = 0.01;
end

if ~isfield(opts,'ControlVariable')
    opts.ControlVariable = 'qgm_MultiPiezos';
end

if ~isfield(opts,'doDebug')
    opts.doDebug=1;
end

if ~isfield(opts,'ROI')
    opts.ROI=[260 300 20 500];
end


if ~isfield(opts,'Name')
    opts.Name=[];
end

disp('momentum pixel domain focus analysis');

%% Construct Filters
% Construct the bandpass filter at a radius of the appropriate value
klims = opts.kmag+[-1 1]*opts.kdelta;

fx = data(1).f;
fy = data(1).f;

[fxx,fyy]=meshgrid(fx,fy);

fmat=(fxx.^2+fyy.^2).^(1/2);

LPF = fmat<klims(2);
HPF = fmat>klims(1);
PF = logical(LPF.*HPF);

%% Iterate

for kk=1:length(data)
    X = data(kk).Params.(opts.ControlVariable);
    X(isnan(X))=[];
    
    if data(kk).ProcessOptions.doSubtractBG
        L = length(X)/2;
        X = X(1:L);
    end

    focus = struct;
    focus.ControlVariable = opts.ControlVariable;
    focus.X = X;

    % Calculate FFT on ROI
    ix_1 = find(data(kk).X>=opts.ROI(1),1);
    ix_2 = find(data(kk).X>=opts.ROI(2),1);
    iy_1 = find(data(kk).Y>=opts.ROI(3),1);
    iy_2 = find(data(kk).Y>=opts.ROI(4),1);            
    Xpx = data(kk).X(ix_1:ix_2);
    Ypx = data(kk).Y(iy_1:iy_2);   
    Z = data(kk).Z(iy_1:iy_2,ix_1:ix_2,:);       

    Nfft = 2^10+1;       
    dX = Xpx(2)-Xpx(1);
    f_max = 1/dX;
    f=1/2*linspace(-f_max,f_max,Nfft);
    zfnorm=zeros(Nfft,Nfft,size(Z,3));

    for nn=1:size(Z,3)
        zf = fftshift(fft2(Z(:,:,nn),Nfft,Nfft));
        zfnorm(:,:,nn)=abs(zf);   
    end     


    % Find Scores
    scores=zeros(size(data(kk).ZfNorm,3),1);
    for nn=1:size(zfnorm,3)
        z = zfnorm(:,:,nn);
        zsub2= z.*PF;
        zsub = z(PF);    
        s = std(zsub);    
        zsub2(zsub2<2*s)=0;      
        % scores(nn) = sum(zsub2,'all')/sum(z,'all');
        % zsub(zsub<2*s)=0;
        b = sort(zsub,'descend');
        N= 1000;

        % K means score
        % [IDX, C, SUMD, D]  = kmeans(b(1:N),2);
        % scores(nn) = max(C)/sum(zsub);

        % Just look at the top N pixels
        scores(nn)=sum(b(1:N))/sum(z,'all');        
    end    
    focus.Scores = scores;
    
    [val,ind]=max(scores);
    myfit=fittype('-A*(x-x0).^2+B','independent','x',...
        'coefficients',{'A','B','x0'});
    fitopt = fitoptions(myfit);
    fitopt.Start = [1e-3 1 X(ind)];
    fitopt.Lower = [0 .9 min(X)-.5];
    fitopt.Upper = [inf 1.1 max(X)+.5];
    [fout,gof,~] = fit(X',scores/val,myfit,fitopt);
    focus.Fit = fout;
    focus.FocusCenter = fout.x0;
    focus.Rsquare = gof.rsquare;
    data(kk).KFocusing = focus;
    

    if opts.doDebug    
        if ~isfield(opts,'Parent')
            FigName = 'StripeCircular';
            ff=get(groot,'Children');
        
            fig=[];
            for kk=1:length(ff)
                if isequal(ff(kk).Name, FigName)
                    fig = ff(kk);
                end
            end
            if isempty(fig)
                fig = figure;
            end    
            figure(fig);
            fig.Color='w';
            fig.Name=FigName;
            fig.Position=[5 50 350 300];
            clf(fig);
            figure(fig);
        else
            fig = opts.Parent;
            for mm=1:length(fig.Children)
                delete(fig.Children(1))
            end
        end  


        if ~isempty(data(kk).Name)
            tt=uicontrol('parent',fig,'units','pixels','string',data(kk).Name,...
            'fontsize',8,'horizontalalignment','left','backgroundcolor','w',...
            'style','text');
            tt.Position=[1 1 300 15];
        end
        coP = jet(length(scores));


        [~,indMax]=max(scores);
        cL = [0 max(data(kk).ZNoFilter(:,:,indMax),[],'all')*.5];
        ax_boop=[];
        for nn=1:4
            ax_boop(nn)=subplot(3,2,nn,'parent',fig); 
            imagesc(data(kk).X,data(kk).Y,data(kk).ZNoFilter(:,:,nn),'parent',ax_boop(nn));
            axis(ax_boop(nn),'equal');
            axis(ax_boop(nn),'tight');
            title(ax_boop(nn),...
                ['(' num2str(nn) ') : ' num2str(X(nn)) 'V, score=' num2str(scores(nn),'%.2e')]);
            set(ax_boop(nn),'XLim',256+[-40 40],'YLim',256+[-40 40])
            colorbar(ax_boop(nn))
            caxis(ax_boop(nn),cL);
        end

        ax_piezo=subplot(3,2,5,'parent',fig);
        co2=get(ax_piezo,'colororder');
        for nn=1:length(X)
            plot(nn,X(nn),'ko','markerfacecolor',coP(nn,:),'markersize',8,...
                'linewidth',2,'parent',ax_piezo);
            hold on        
        end

        xlabel(ax_piezo,'image number');
        ylabel(ax_piezo,'piezo (V)')
        set(ax_piezo,'box','on','linewidth',1,'fontsize',10,'Xgrid','on','ygrid','on');


        ax_score=subplot(3,2,6,'parent',fig);
        xR=max(X)-min(X);
        xx=linspace(min(X)-xR*.2,max(X)+xR*.2,100);
        pF=plot(xx,feval(fout,xx)*val,'r-','linewidth',2,'parent',ax_score);
        hold(ax_score,'on');
        for nn=1:length(X)
            plot(X(nn),scores(nn),'ko','markerfacecolor',coP(nn,:),'markersize',8,...
                'linewidth',2,'parent',ax_score);
            hold on
        end

        xlabel(ax_score,'piezo (V)')
        ylabel(ax_score,'focus score (arb)')
        yL = get(ax_score,'YLim');
        % ylim(ax_score,[0 yL(2)]);
        str = ['$V_0 : ' num2str(round(fout.x0,2)) '~\mathrm{V}$' newline '$R^2:' num2str(round(focus.Rsquare,3)) '$'];

        legend(pF,{str},'location','south','interpreter','latex','parent',fig)               
        set(ax_score,'box','on','linewidth',1,'fontsize',10,'xgrid','on','ygrid','on');

    end
%     keyboard
end

end

