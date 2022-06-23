function [dark_data] = ixon_analyzeDarkImage(ixondata,ixon_xVar)

fprintf('Performing dark image analysis ...');    


for kk=1:length(ixondata)
    
    DarkData=struct;    
    for k=1:size(ixondata(kk).ROI,1)
        
        Ncounts = sum(sum(ixondata(kk).RawImages(:,:,1)));      
    end         
    ixondata(kk).BoxCount=DarkData;
end

end

