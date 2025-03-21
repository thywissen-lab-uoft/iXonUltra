function [hF]= ixon_show_multi_shot_focusing(focus,xVar,opts)

if nargin ==1
    xVar = 'ExecutionDate';
end

if nargin < 3
    opts = struct;
end

if ~isfield(opts,'ForceFit')
    opts.ForceFit=0;
end

%% Sort the data by the parameter given

P = [focus.Params];
X = [P.(xVar)];

%%

s1 = [focus.Score1];
s2 = [focus.Score2];
sg = [focus.ScoreGauss];

V1 = [focus.Piezo1];
V2 = [focus.Piezo2];

Vc = [P.objective_piezo];

dV = mean((V2-V1))/2;

dP = V1-V2;

dS = (s1-s2)./sg;

dSdP = dS./dP;


mycorr = [focus.Correlator];

hF = figure;
hF.Color='w';
hF.Position = [50 50 1200 400];
co=get(gca,'colororder');

ax1 = subplot(121);
p1=plot(X,s1,'o-','markerfacecolor',co(4,:),'color',co(4,:)*.5,...
    'markersize',8,'linewidth',1);
hold on
p2=plot(X,s2,'o-','markerfacecolor',co(5,:),'color',co(5,:)*.5,...
    'markersize',8,'linewidth',1);
pg=plot(X,sg,'o-','markerfacecolor',[.5 .5 .5],'color','k',...
    'markersize',8,'linewidth',1);
if isequal(xVar,'ExecutionDate')
    datetick('x');
end

xlabel(xVar,'interpreter','none');
ylabel('$s_i:=\sum\tilde{f_i}(\vec{k})\cdot k$','interpreter','latex')
set(gca,'box','on','linewidth',1,'fontsize',12,'fontname','times')
title('focus score');
hold on
% 
% yyaxis right
% set(gca,'YColor','k')
% plot(X,V1,'.-','color',co(4,:),'linewidth',1);
% plot(X,V2,'.-','color',co(5,:),'linewidth',1);
% plot(X,0.5*(V1+V2),'k--','color','k','linewidth',.5);
% ylabel('piezo (V)')
str1 = ['Image1 (+' num2str(abs(dV)) 'V)'];
str2 = ['Image2 (-' num2str(abs(dV)) 'V)'];

legend([p1 p2 pg],{str1,str2,'Gauss'},'fontsize',8)


ax2 = subplot(122);
yyaxis left
plot(X,dSdP,'o-','markerfacecolor',co(1,:),'color',co(1,:)*.5,...
    'markersize',8,'linewidth',1);
ylabel('$(s_1-s_2)/(V_1-V_2)/s_G$','interpreter','latex');
hold on
yyaxis right
plot(X,mycorr,'o-','markerfacecolor',co(2,:),'color',co(2,:)*.5,...
    'markersize',8,'linewidth',1);
if isequal(xVar,'ExecutionDate')
    datetick('x');
end
xlabel(xVar,'interpreter','none');
set(gca,'box','on','linewidth',1,'fontsize',12,'fontname','times')
title('differential score and correlator');
ylabel('2d correlator (arb.)','interpreter','latex');
hold on
linkaxes([ax1 ax2],'x');
%% Fit it
if isequal(xVar,'objective_piezo') || opts.ForceFit
    x = X(:);
    y1 = dSdP(:);
    y2 = mycorr(:);

    Sg = 0.1;
    amp = 0.5*(max(y1)-min(y1));
    [m1,ind1]=max(y1);
    [m2,ind2]=min(y1);
    Cg = (abs(m1)*x(ind1)+abs(m2)*x(ind2))/(abs(m1)+abs(m2));

    Mg = -2*amp/(x(ind2)-x(ind1));
    Ag = exp(0.5)*amp/Sg;

    % Fit Differential Score
    fit1 = fittype('-A*(x-xc).*exp(-(x-xc).^2/(2*s^2))','independent','x',...
        'coefficients',{'A','s','xc'});
    fitopt1=fitoptions(fit1);
    fitopt1.StartPoint = [Ag, Sg Cg];
    fout1 = fit(x,y1,fit1,fitopt1);

    fit1b = fittype('m*(x-xc)','independent','x',...
        'coefficients',{'m','xc'});
    fitopt1b=fitoptions(fit1b);
    fitopt1b.StartPoint = [Mg Cg];
    indsL = x>x(ind1);
    indsH = x<x(ind2);
    inds = logical(indsL.*indsH);
    fout1b = fit(x(inds),y1(inds),fit1b,fitopt1b);


    axes(ax2);
    yyaxis left
    xF = linspace(min(x),max(x),100);
    pF1=plot(xF,feval(fout1,xF),'--','color',co(1,:),'linewidth',2,...
        'parent',ax2);

    xFb = linspace(x(ind1),x(ind2),5);
    pF2=plot(xFb,feval(fout1b,xFb),'-','color',co(1,:),'linewidth',2,...
        'parent',ax2);

    str1b = ['$m(x-x_C) ~ (' num2str(round(fout1b.m,1)) '/\mathrm{V},' ...
        num2str(round(fout1b.xc,2)) '\mathrm{V})$'];
    legend([pF2],{str1b},'interpreter','latex');

end

