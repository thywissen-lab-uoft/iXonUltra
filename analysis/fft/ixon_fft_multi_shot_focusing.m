function data=ixon_fft_multi_shot_focusing(data,opts)

if nargin==1
    opts=struct;
end

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

%% Construct Filters

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

    % Find Scores
    scores=zeros(size(data(kk).ZfNorm,3),1);
    for nn=1:size(data(kk).ZfNorm,3)
        z = data(kk).ZfNorm(:,:,nn);
        zsub = z(PF);    
        N0 = sum(z,'all');
        Nmean = mean(zsub);
        Nstd = std(zsub);
        inds = [z>(Nmean+5*Nstd)];
        inds = inds.*PF;
        inds = logical(inds);
        Nk=sum(z(inds),'all');    
        scores(nn) = Nk/N0;
    end
    focus.Scores = scores;
    
    [val,ind]=max(scores);
    myfit=fittype('-A*(x-x0).^2+B','independent','x',...
        'coefficients',{'A','B','x0'});
    fitopt = fitoptions(myfit);
    fitopt.Start = [1e-3 1 X(ind)];
    fitopt.Lower = [0 .9 min(X)-.5];
    fitopt.Upper = [inf 1.5 max(X)+.5];
    fout = fit(X',scores/val,myfit,fitopt);
    focus.Fit = fout;
    focus.FocusCenter = fout.x0;
    data(kk).KFocusing = focus;

    if opts.doDebug
        FigName = 'FFTFocusing';
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
        fig.Position=[50 50 400 500];
        co2=get(gca,'colororder');

        tt=uicontrol('style','text','string',data(kk).Name','fontsize',8,...
            'horizontalalignment','left','backgroundcolor','w');
        tt.Position=[1 1 300 15];

        subplot(3,1,1,'parent',fig)
        plot(X,'ko','markerfacecolor',[.5 .5 .5],'markersize',8,...
            'linewidth',2);
        xlabel('image number');
        ylabel('piezo (V)')
        set(gca,'box','on','linewidth',1,'fontsize',10);
        grid on

        subplot(3,1,[2 3])
        xx=linspace(min(X)-.1,max(X)+.1,100);
        plot(xx,feval(fout,xx)*val,'r-','linewidth',2)
        hold on
        plot(X,scores,'o','markerfacecolor',co2(1,:),...
            'markeredgecolor',co2(1,:)*.5,'linewidth',2,...
           'markersize',8);
        xlabel('piezo (V)')
        ylabel('momentum peak score (arb.)')
        yL = get(gca,'YLim');
        ylim([0 yL(2)]);
        grid on
        str = ['$V_0 = ' num2str(round(fout.x0,2)) '~\mathrm{V}$'];
        legend({'data',str},'location','south','interpreter','latex')               
        set(gca,'box','on','linewidth',1,'fontsize',10);
    end
end

end

