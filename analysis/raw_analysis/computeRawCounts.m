function atomdata=computeRawCounts(atomdata,opts)


for nn=1:length(atomdata)
    for jj=1:size(atomdata(nn).RawImages,3)
        data=atomdata(nn).RawImages(:,:,jj); 
        
        cMax=max(data,[],'all');
        cMin=min(data,[],'all');        
        cMed=median(data,'all');
        cAvg=mean(data,'all');
        cStd=std(data,1,'all');
        
        atomdata(nn).Raw(jj).Maximum=cMax;
        atomdata(nn).Raw(jj).Minimum=cMin;
        atomdata(nn).Raw(jj).Median=cMed;
        atomdata(nn).Raw(jj).Mean=cAvg;
        atomdata(nn).Raw(jj).Std=cStd;
        atomdata(nn).Raw(jj).Total=sum(data,'all');
    end
end

end

