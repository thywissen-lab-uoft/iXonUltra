function [scores,centers]=bin_StripeScore(n1,Zb,seps,threshold,sumdir)
Zb(Zb<threshold(1))=0;

s = sum(Zb,'all');

ii=zeros(length(seps),1);
for kk=1:length(seps)
    ii(kk) = find(n1 == seps(kk),1);
end
scores = zeros(length(seps)-1,1);
centers = zeros(length(seps)-1,1);
for kk=2:(length(seps))
    
    if sumdir == 2
        Zbsub = Zb((ii(kk-1)):ii(kk),:);
    else 
        Zbsub = Zb(:,(ii(kk-1)):ii(kk));
    end
    Nhigh = sum(Zbsub>=threshold(2),'all');
    

    s_sub = sum(Zbsub,'all');
    
    if (s_sub/s)<0.05
       Nhigh = 0; 
    end
    
    scores(kk-1) = Nhigh*(1/s_sub);
    centers(kk-1) = round(0.5*(seps(kk-1) + seps(kk)));
end



end

