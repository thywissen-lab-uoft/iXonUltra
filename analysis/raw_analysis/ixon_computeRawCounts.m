function ixondata=ixon_computeRawCounts(ixondata,opts)


for nn=1:length(ixondata)
    for jj=1:size(ixondata(nn).RawImages,3)
        data=ixondata(nn).RawImages(:,:,jj); 
        
        cMax=max(data,[],'all');
        cMin=min(data,[],'all');        
        cMed=median(data,'all');
        cAvg=mean(data,'all');
        cStd=std(data,1,'all');
        
        ixondata(nn).Raw(jj).Maximum=cMax;
        ixondata(nn).Raw(jj).Minimum=cMin;
        ixondata(nn).Raw(jj).Median=cMed;
        ixondata(nn).Raw(jj).Mean=cAvg;
        ixondata(nn).Raw(jj).Std=cStd;
        ixondata(nn).Raw(jj).Total=sum(data,'all');
    end
end

end

