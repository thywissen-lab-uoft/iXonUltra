function [cface,cedge] = ixoncolororder(n)

co=[0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];
co=circshift(co,4,1);
hsl=rgb2hsl(co);
hsl2=rgb2hsl(co);
hsl(:,3)=ones(size(hsl,1),1)*.9;
hsl2(:,2)=ones(size(hsl2,1),1);

coface=hsl2rgb(hsl);
coedge=hsl2rgb(hsl2);

if nargin==1
    cface=coface(mod(n-1,7)+1,:);
    cedge=coedge(mod(n-1,7)+1,:);
else
    cface=coface;
    cedge=coedge;
end

end

