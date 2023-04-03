
function M = computeMoment(Z,ix,iy)

x=1:size(Z,2);
y=1:size(Z,1);

[xx,yy]=meshgrid(x,y);

M = sum(sum(xx.^(ix).*yy.^(iy).*Z));

end