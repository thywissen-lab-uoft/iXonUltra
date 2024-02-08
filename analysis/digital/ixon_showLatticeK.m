function [out,hF] = ixon_showLatticeK(data,opts)


if nargin == 2 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end
%%
% k1mag = zeros(length(data),length(data(1).LatticeK));
% k1mag_err = zeros(length(data),length(data(1).LatticeK));
k1x = zeros(length(data),length(data(1).LatticeK));
k1y = zeros(length(data),length(data(1).LatticeK));
k1x_err = zeros(length(data),length(data(1).LatticeK));
k1y_err = zeros(length(data),length(data(1).LatticeK));

% k2mag = zeros(length(data),length(data(1).LatticeK));
% k2mag_err = zeros(length(data),length(data(1).LatticeK));
k2x = zeros(length(data),length(data(1).LatticeK));
k2y = zeros(length(data),length(data(1).LatticeK));
k2x_err = zeros(length(data),length(data(1).LatticeK));
k2y_err = zeros(length(data),length(data(1).LatticeK));

for nn = 1 : length(data)
    for kk = 1 :length(data(nn).LatticeK)
        k1x(nn,kk) = data(nn).LatticeK(kk).Fit1.xc;
        k1y(nn,kk) = data(nn).LatticeK(kk).Fit1.yc;        
        c1 = confint(data(nn).LatticeK(kk).Fit1);
        k1x_err(nn,kk) = diff(c1(:,2))*.5;
        k1y_err(nn,kk) = diff(c1(:,3))*.5;

        k2x(nn,kk) = data(nn).LatticeK(kk).Fit2.xc;
        k2y(nn,kk) = data(nn).LatticeK(kk).Fit2.yc;        
        c2 = confint(data(nn).LatticeK(kk).Fit2);
        k2x_err(nn,kk) = diff(c2(:,2))*.5;
        k2y_err(nn,kk) = diff(c2(:,3))*.5;
    end 
end
%%

k1mag = sqrt(k1x.^2+k1y.^2);
k1mag_err = sqrt((k1x./k1mag.*k1x_err).^2+(k1y./k1mag.*k1y_err).^2);
theta1 = atan2(k1y,k1x);
theta1_err = sqrt((k1x.^2.*k1x_err.^2+k1y.^2.*k1y_err.^2)./k1mag.^2);

k1mag = k1mag(:);
k1mag_err = k1mag_err(:);
theta1 = theta1(:);
theta1_err = theta1_err(:);

[k1mag,i1] = rmoutliers(k1mag);
k1mag_err = k1mag_err(~i1);
N_outliers_k1 = sum(i1);

i1mat = reshape(i1,[length(data) length(data(1).LatticeK)]);

[theta1,i2] = rmoutliers(theta1);
theta1_err = theta1_err(~i2);
N_outliers_theta1 = sum(i2);
i2mat = reshape(i2,[length(data) length(data(1).LatticeK)]);

%%
k2mag = sqrt(k2x.^2+k2y.^2);
k2mag_err = sqrt((k2x./k2mag.*k1x_err).^2+(k2y./k2mag.*k1y_err).^2);
theta2 = atan2(k2y,k2x);
theta2_err = sqrt((k2x.^2.*k1x_err.^2+k2y.^2.*k1y_err.^2)./k2mag.^2);


k2mag = k2mag(:);
k2mag_err = k2mag_err(:);
theta2 = theta2(:);
theta2_err = theta2_err(:);

[k2mag,i3] = rmoutliers(k2mag);
k2mag_err = k2mag_err(~i3);
N_outliers_k2 = sum(i3);
i3mat = reshape(i3,[length(data) length(data(1).LatticeK)]);

[theta2,i4] = rmoutliers(theta2);
theta2_err = theta2_err(~i4);
N_outliers_theta2 = sum(i4);
i4mat = reshape(i4,[length(data) length(data(1).LatticeK)]);

%%

hF = figure(5000);
hF.Color='w';
hF.Position = [100 100 600 400];
clf

subplot(221);
histfit(k1mag,50)
hold on
p1=fitdist(k1mag,'normal');
xlabel('$k_1$ magnitude (1/px)','interpreter','latex');
ylabel('occurences','interpreter','latex');
str = ['$' num2str(p1.mu,'%.4f') '\pm' num2str(p1.sigma,'%.1e') '$' newline ...
    '$N_{\mathrm{bad}} = ' num2str(N_outliers_k1) '$'];
text(.02,.98,...
    str,...
    'units','normalized',...
    'verticalalignment','top','horizontalalignment','left',...
    'interpreter','latex','fontsize',10);

subplot(222);
histfit(k2mag,50)
hold on
p2=fitdist(k2mag,'normal');
xlabel('$k_2$ magnitude (1/px)','interpreter','latex');
ylabel('occurences','interpreter','latex');
str = ['$' num2str(p2.mu,'%.4f') '\pm' num2str(p2.sigma,'%.1e') '$' newline ...
    '$N_{\mathrm{bad}} = ' num2str(N_outliers_k2) '$'];
text(.02,.98,...
    str,...
    'units','normalized',...
    'verticalalignment','top','horizontalalignment','left',...
    'interpreter','latex','fontsize',10);


subplot(223);
histfit(theta1/pi*180,50)
hold on
p3=fitdist(theta1/pi*180,'normal');
xlabel('$\theta_1$ (deg.)','interpreter','latex');
ylabel('occurences','interpreter','latex');
str = ['$' num2str(p3.mu,'%.4f') '\pm' num2str(p3.sigma,'%.1e') '$' newline ...
    '$N_{\mathrm{bad}} = ' num2str(N_outliers_theta1) '$'];
text(.02,.98,...
    str,...
    'units','normalized',...
    'verticalalignment','top','horizontalalignment','left',...
    'interpreter','latex','fontsize',10);



subplot(224);
histfit(theta2/pi*180,50)
hold on
p4=fitdist(theta2/pi*180,'normal');
xlabel('$\theta_2$ (deg.)','interpreter','latex');
ylabel('occurences','interpreter','latex');
str = ['$' num2str(p4.mu,'%.4f') '\pm' num2str(p4.sigma,'%.1e') '$' newline ...
    '$N_{\mathrm{bad}} = ' num2str(N_outliers_theta2) '$'];
text(.02,.98,...
    str,...
    'units','normalized',...
    'verticalalignment','top','horizontalalignment','left',...
    'interpreter','latex','fontsize',10);

t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];
    
%% Acrew Output

out = struct;
out.BadLattice = logical(i1mat(:,1)+i2mat(:,1)+i3mat(:,1)+i4mat(:,1));
out.k1x = p1.mu*cosd(p3.mu);
out.k1y = p1.mu*sind(p3.mu);
out.k2x = p2.mu*cosd(p4.mu);
out.k2y = p2.mu*sind(p4.mu);
out.k1 = [out.k1x out.k1y];
out.k2 = [out.k2x out.k2y];

end

