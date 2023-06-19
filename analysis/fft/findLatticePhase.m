function  out=findLatticePhase(x,y,z,k1,k2)

sc = 3;

% Rescale the image to avoid rounding errors when moving pixels
z2 = imresize(z,sc)/(sc^2); % scale amplitude to preserve norm
x2 = linspace(x(1),x(end),size(z2,2));
y2 = linspace(y(1),y(end),size(z2,1));

[X,Y] = meshgrid(x2,y2);

p1=mod(1-(angle(sum(exp(-1i*2*pi*(k1(1).*X+k1(2).*Y)).*z2,'all')))/(2*pi),1);
p2=mod(1-(angle(sum(exp(-1i*2*pi*(k2(1).*X+k2(2).*Y)).*z2,'all')))/(2*pi),1);

Q = [0 -1;1 0];
a1 = Q*k2/(k1'*Q*k2);
a2 = -Q*k1/(-k2'*Q*k1);

out = struct;
out.k1 = k1;
out.k2 = k2;
out.p1 = p1;
out.p2 = p2;
out.a1 = a1;
out.a2 = a2;
    
end

