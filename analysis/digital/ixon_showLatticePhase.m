function [hF] = ixon_showLatticePhase(data,opts)

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
a1mag_std = std(a1mag);
theta1 = atan2(a1y,a1x);
theta1_std=std(theta1);

a2mag = sqrt(a2x.^2+a2y.^2);
a2mag_std = std(a2mag);
theta2 = atan2(a2y,a2x);
theta2_std=std(theta1);

%%

hF = figure(5002);
hF.Color='w';
hF.Position = [100 100 600 200];
clf

hF.Name = 'Lattice Phase';

subplot(121);
histogram(p1,50)
hold on
xlabel('$\phi_1$ phase (sites)','interpreter','latex');
ylabel('occurences','interpreter','latex');
xlim([0 1]);
str= ['$a_1 = ' num2str(mean(a1mag),'%.3f') '\pm' num2str(mean(a1mag_std),'%.3f') '$' newline ...
    '$\theta_1 = ' num2str(mean(theta1*180/pi),'%.3f') '\pm' num2str(mean(theta1_std*180/pi),'%.3f') '^\circ$'];
text(.02,.98,...
    str,...
    'units','normalized',...
    'verticalalignment','top','horizontalalignment','left',...
    'interpreter','latex','fontsize',10);


subplot(122);
histogram(p2,50)
hold on
xlabel('$\phi_1$ phase (sites)','interpreter','latex');
ylabel('occurences','interpreter','latex');
xlim([0 1]);
str= ['$a_2 = ' num2str(mean(a2mag),'%.3f') '\pm' num2str(mean(a2mag_std),'%.3f') '$' newline ...
    '$\theta_2 = ' num2str(mean(theta2*180/pi),'%.3f') '^\circ\pm' num2str(mean(theta2_std*180/pi),'%.3f') '^\circ$'];
text(.02,.98,...
    str,...
    'units','normalized',...
    'verticalalignment','top','horizontalalignment','left',...
    'interpreter','latex','fontsize',10);

    
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

