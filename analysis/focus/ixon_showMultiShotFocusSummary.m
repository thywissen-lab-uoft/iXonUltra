function [hF]= ixon_showMultiShotFocusSummary(focus,xVar,opts)

if nargin ==1
    xVar = 'ExecutionDate';
end

if nargin < 3
    opts = struct;
end

%% Sort the data by the parameter given

P = [focus.Params];
X = [P.(xVar)];

%%

Scores = [focus.Scores];
Counts = [focus.Counts];


hF = figure;
hF.Color='w';
hF.Position = [50 50 1200 400];
co=get(gca,'colororder');
clf

% ax1 = subplot(1,2,1,'parent',hF);
legStr={};
% for kk=1:size(Scores,2)
%     ps(kk) =plot(X,Scores(:,kk),'o','markerfacecolor',co(kk,:),'color',co(kk,:)*.5,...
%     'markersize',8,'linewidth',1);
%     hold on
%     legStr{kk}=['image ' num2str(kk)];
% end
% 
% if isequal(xVar,'ExecutionDate')
%     datetick('x');
% end
% 
% xlabel(xVar,'interpreter','none');
% ylabel('$\sum\tilde{f_i}(\vec{k})\cdot k$','interpreter','latex')
% set(gca,'box','on','linewidth',1,'fontsize',12,'fontname','times')
% title('focus score');
% hold on
% 
% legend(ps,legStr,'location','best');

% ax2 = subplot(1,2,2,'parent',hF);
ax2=axes;
yyaxis left
set(gca,'YColor','k');
for kk=1:size(Scores,2)
    ps(kk) =plot(X,Scores(:,kk)./sum(Scores,2),'o','markerfacecolor',co(kk,:),'color',co(kk,:)*.5,...
    'markersize',8,'linewidth',1);
    hold on
        legStr{kk}=['image ' num2str(kk)];

end

if isequal(xVar,'ExecutionDate')
    datetick('x');
end

xlabel(xVar,'interpreter','none');
ylabel('normalized score')
set(gca,'box','on','linewidth',1,'fontsize',12,'fontname','times')
title('focus score');
hold on

yyaxis right
ps(end) =plot(X,sum(Scores,2),'-','color',[.5 .5 .5],...
'markersize',8,'linewidth',1);
hold on
legStr{end}='sum score';
set(gca,'YColor',[.5 .5 .5]);
ylabel('sum score');

text(.01,.01,'$\sum\tilde{f_i}(\vec{k})\cdot k$','interpreter','latex',...
    'units','normalized','horizontalalignment','left','verticalalignment','bottom')

legend(ps,legStr,'location','best');
end
