function [hF] = ixon_showLatticeA(data,opts)

if nargin ==1
   opts = struct; 
end


if nargin == 2 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end

%%

a1x = zeros(length(data),length(data(1).LatticePhase));
a1y = zeros(length(data),length(data(1).LatticePhase));
p1 = zeros(length(data),length(data(1).LatticePhase));

a2x = zeros(length(data),length(data(1).LatticePhase));
a2y = zeros(length(data),length(data(1).LatticePhase));
p2  = zeros(length(data),length(data(1).LatticePhase));

for nn = 1 : length(data)
    for kk = 1 :length(data(nn).LatticeK)
        a1x(nn,kk) = data(nn).LatticePhase(kk).a1(1);
        a1y(nn,kk) = data(nn).LatticePhase(kk).a1(2);   
        p1(nn,kk) = data(nn).LatticePhase(kk).p1;   
        
        a2x(nn,kk) = data(nn).LatticePhase(kk).a2(1);
        a2y(nn,kk) = data(nn).LatticePhase(kk).a2(2);   
        p2(nn,kk) = data(nn).LatticePhase(kk).p2;   
    end 
end
%%

a1mag = sqrt(a1x.^2+a1y.^2);
theta1 = atan2(a1y,a1x);

a1mag = a1mag(:);
theta1 = theta1(:);

a2mag = sqrt(a2x.^2+a2y.^2);
theta2 = atan2(a2y,a2x);
a2mag = a2mag(:);
theta2 = theta2(:);
%%

hF = figure(5001);
hF.Color='w';
hF.Position = [100 100 600 400];
clf

hF.Name = 'Lattice Vectors';

subplot(221);
histfit(a1mag,50)
hold on
p1=fitdist(a1mag,'normal');
xlabel('$a_1$ magnitude (px)','interpreter','latex');
ylabel('occurences','interpreter','latex');
str = ['$' num2str(p1.mu,'%.4f') '\pm' num2str(p1.sigma,'%.1e') '$'];
text(.02,.98,...
    str,...
    'units','normalized',...
    'verticalalignment','top','horizontalalignment','left',...
    'interpreter','latex','fontsize',10);

subplot(222);
histfit(a2mag,50)
hold on
p2=fitdist(a2mag,'normal');
xlabel('$a_2$ magnitude (1/px)','interpreter','latex');
ylabel('occurences','interpreter','latex');
str = ['$' num2str(p2.mu,'%.4f') '\pm' num2str(p2.sigma,'%.1e') '$'];
text(.02,.98,...
    str,...
    'units','normalized',...
    'verticalalignment','top','horizontalalignment','left',...
    'interpreter','latex','fontsize',10);
% 

subplot(223);
histfit(theta1/pi*180,50)
hold on
p3=fitdist(theta1/pi*180,'normal');
xlabel('$\theta_1$ (deg.)','interpreter','latex');
ylabel('occurences','interpreter','latex');
str = ['$' num2str(p3.mu,'%.4f') '\pm' num2str(p3.sigma,'%.1e') '$'];
text(.02,.98,...
    str,...
    'units','normalized',...
    'verticalalignment','top','horizontalalignment','left',...
    'interpreter','latex','fontsize',10);
% 


subplot(224);
histfit(theta2/pi*180,50)
hold on
p4=fitdist(theta2/pi*180,'normal');
xlabel('$\theta_2$ (deg.)','interpreter','latex');
ylabel('occurences','interpreter','latex');
str = ['$' num2str(p4.mu,'%.4f') '\pm' num2str(p4.sigma,'%.1e') '$'];
text(.02,.98,...
    str,...
    'units','normalized',...
    'verticalalignment','top','horizontalalignment','left',...
    'interpreter','latex','fontsize',10);
% 
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];
    
%% Acrew Output
% 
% out = struct;
% out.BadLattice = logical(i1mat(:,1)+i2mat(:,1)+i3mat(:,1)+i4mat(:,1));
% out.k1x = p1.mu*cosd(p3.mu);
% out.k1y = p1.mu*sind(p3.mu);
% out.k2x = p2.mu*cosd(p4.mu);
% out.k2y = p2.mu*sind(p4.mu);
% out.k1 = [out.k1x out.k1y];
% out.k2 = [out.k2x out.k2y];

end

