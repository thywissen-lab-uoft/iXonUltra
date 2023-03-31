function ixondata = ixon_matchParamsFlags(ixondata)

pAll = sort(fieldnames(ixondata(1).Params));

if isfield(ixondata(1),'Flags')
    fAll = sort(fieldnames(ixondata(1).Flags));
else
    badFlags=1;
end

badParams = 0;
badFlags  = 0;

% Find Master list of all params and flags
for kk=1:length(ixondata)
    pThis = sort(fieldnames(ixondata(kk).Params));     
    
    fThis = sort(fieldnames(ixondata(kk).Flags));
    
    if ~isequal(pThis,pAll)
        badParams = 1;
        for nn=1:length(pThis)      
            ind=findStr(pAll,pThis{nn});             
            if isempty(ind)
                pAll{end+1}=pThis{nn}; 
            end  
       end
    end 
    if isfield(ixondata(kk),'Flags')

        if ~isequal(fThis,fAll)
            badFlags = 1;
            for nn=1:length(fThis)      
                ind=findStr(fAll,fThis{nn});             
                if isempty(ind)
                    fAll{end+1}=fThis{nn}; 
                end  
           end
        end 
    end
end

if badParams
   warning('Unequal number of parameters detected. Adding NaN'); 
end

if badFlags
   warning('Unequal number of flags detected. Adding NaN'); 
end

% Add blank values from the master list
for kk=1:length(ixondata)
    pThis = sort(fieldnames(ixondata(kk).Params));      
    fThis = sort(fieldnames(ixondata(kk).Flags));
    
    if ~isequal(pAll,pThis)        
       for nn=1:length(pAll)      
            ind=findStr(pThis,pAll{nn});   
            if isempty(ind)   
                ixondata(kk).Params.(pAll{nn}) = NaN ;
                ixondata(kk).Units.(pAll{nn}) = NaN ;
            end       
       end
    end
        if isfield(ixondata(kk),'Flags')

    if ~isequal(fAll,fThis)        
       for nn=1:length(fAll)      
            ind=findStr(fThis,fAll{nn});   
            if isempty(ind)
               ixondata(kk).Flags.(fAll{nn}) = NaN ;
            end       
       end 
    end 
        end
    
end

if badParams
    for kk=1:length(ixondata)
        ixondata(kk).Params = orderfields(ixondata(kk).Params);
        ixondata(kk).Units = orderfields(ixondata(kk).Units);

    end
end

if badFlags
    for kk=1:length(ixondata)
        ixondata(kk).Flags = orderfields(ixondata(kk).Flags);
    end
end

end

function ind=findStr(cellstr,str)

ind=[];
for kk=1:length(cellstr)
   if isequal(str,cellstr{kk})
       ind = kk;
   end
end

end

