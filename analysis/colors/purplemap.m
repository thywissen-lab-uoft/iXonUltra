function cm_data = purplemap(m)

% The purple map that we like to use.
c = 0:0.005:1; d = 0.05:0.005:1; 
cm = [(repmat((1-c)',1,3).*repmat([0 0 0],length(c),1) + ...           
            repmat((c)',1,3).*repmat([0.4 0 0.3],length(c),1)); ...
           (repmat((1-d)',1,3).*repmat([0.4 0 0.3],length(d),1) + ...
            repmat((d)',1,3).*repmat([0.8 0.7 0.4],length(d),1)); ...
            (repmat((1-d)',1,3).*repmat([0.8 0.7 0.4],length(d),1) + ...
            repmat((d)',1,3).*repmat([1 1 1],length(d),1))];
                        
            
if nargin < 1
    cm_data = cm;
else
    hsv=rgb2hsv(cm);
    cm_data=interp1(linspace(0,1,size(cm,1)),hsv,linspace(0,1,m));
    cm_data=hsv2rgb(cm_data);
  
end
            
end

