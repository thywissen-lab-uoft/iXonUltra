function Focusing = bin_multiShotFocusing(n1,n2,Zb,opts)

if ~isfield(opts,'Sites')
    opts.Sites=[50 50];
end

if ~isfield(opts,'XVal')
    opts.XVal = 1:size(Zb,3);
end

%% Get Central Area

X=opts.XVal;
[nn1, nn2]=meshgrid(n1,n2);

n1c = sum(Zb(:,:,1).*nn1,'all')/sum(Zb(:,:,1),'all');
n1c = round(n1c);
n2c = sum(Zb(:,:,1).*nn2,'all')/sum(Zb(:,:,1),'all');
n2c = round(n2c);

n1_a = max([min(n1) n1c-opts.Sites(1)]);
n1_b = min([max(n1) n1c+opts.Sites(1)]);
i1_a = find(n1==n1_a,1);
i1_b = find(n1==n1_b,1);

n2_a = max([min(n2) n2c-opts.Sites(2)]);
n2_b = min([max(n2) n2c+opts.Sites(2)]);
i2_a = find(n2==n2_a,1);
i2_b = find(n2==n2_b,1);

n1_sub = n1_a:n1_b;
n2_sub = n2_a:n2_b;
Zb_sub = Zb(i2_a:i2_b,i1_a:i1_b,:);

ignoreFirst=1;
if ignoreFirst
    X=X(2:end);
    Zb_sub=Zb_sub(:,:,2:end);
end



thresholds=[];
val=[];

Nc = sum(Zb_sub,[1 2]);
Ncbar=mean(Nc);


for jj=1:size(Zb_sub,3)
    data =  Ncbar/Nc(jj).*Zb_sub(:,:,jj);
    [idx,c,sumD,D] = kmeans(data(:),2);

        % Sort by centroid
    [c,inds]=sort(c,'ascend');
    sumD=sumD(inds);

    dd = zeros(numel(c),1);
    thresh = zeros(numel(c),1);
        nCluster=zeros(numel(c),1);

        for bb=1:numel(c)
            ind = inds(bb);
            nThis = sum(idx==ind);
            dd(bb) = sqrt(sumD(bb)/nThis);

            thresh(bb) = round(max(data(idx==ind)));
            nCluster(bb)=nThis;
        end

    thresholds(jj) = thresh(1);
    val(jj) = sum(data>thresh(1),'all');
    
[sharpnessScores(jj), map] =MLVSharpnessMeasure(Zb_sub(:,:,jj));


end

thresh_max=round(max(thresholds));

for jj=1:size(Zb_sub,3)
    Nabove(jj) = sum( Ncbar/Nc(jj).*Zb_sub(:,:,jj)>thresh_max,'all');
end

%%

hF=figure(opts.FigureNumber);
clf
hF.Color='w';
hF.Position=[50 50 900 800];
hF.Name='MultiShotFocusing';
co=get(gca,'colororder');

for jj=1:size(Zb_sub,3)
    subplot(3,3,jj)
    imagesc(n1_sub,n2_sub,Zb_sub(:,:,jj));
    xlabel('n1');
    ylabel('n2');
    title(num2str(X(jj)));
    caxis([0 max(thresholds)]);
    colorbar
    axis equal tight
    colormap(hF,bone)
end


for jj=1:size(Zb_sub,3)
    subplot(3,3,3+jj)
    histogram(Zb_sub(:,:,jj),50);
    xlabel('counts');
    ylabel('occurences');
    title(num2str(X(jj)));
end


subplot(3,3,7);
plot(X,thresholds,'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'linewidth',1,'markersize',10,...
    'linewidth',2);
xlabel('piezo (V)')
ylabel('threshold');
grid on

subplot(3,3,8);
plot(X,Nabove,'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'linewidth',1,'markersize',10,...
    'linewidth',2);
ylabel(['N>' num2str(thresh_max) ]);
xlabel('piezo (V)')
pp=polyfit(X,Nabove,2);
hold on
xx=linspace(min(X),max(X),100);
plot(xx,polyval(pp,xx),'r-','linewidth',1);
[val,ind]=max(polyval(pp,xx));
x0=xx(ind);
text(.05,.05,['peak : ' num2str(round(x0,2)) 'V'],'units','normalized',...
    'horizontalalignment','left','verticalalignment','bottom','fontsize',12)
grid on



subplot(3,3,9);
plot(X,sharpnessScores,'o','markerfacecolor',co(1,:),...
    'markeredgecolor',co(1,:)*.5,'linewidth',1,'markersize',10,...
    'linewidth',2);
ylabel(['sharpness score']);
xlabel('piezo (V)')
pp=polyfit(X,sharpnessScores,2);
hold on
xx=linspace(min(X),max(X),100);
plot(xx,polyval(pp,xx),'r-','linewidth',1);
[val,ind]=max(polyval(pp,xx));
x0=xx(ind);
text(.05,.05,['peak : ' num2str(round(x0,2)) 'V'],'units','normalized',...
    'horizontalalignment','left','verticalalignment','bottom','fontsize',12)
grid on


Focusing = struct;
end

