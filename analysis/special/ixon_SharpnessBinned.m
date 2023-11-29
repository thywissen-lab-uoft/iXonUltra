function [data] = ixon_SharpnessBinned(data)


% Ignore lattice sites along the edge (LEFT RIGHT BOTTOM TOP
v = [3 4 2 2];

for kk=1:length(data)
    for nn=1:length(data(kk).LatticeBinNoFilter)
        Zsub = data(kk).LatticeBinNoFilter(nn).Zbin;
        Zsub = Zsub(v(3):(end-v(4)),v(1):(end-v(2)));
   

        
        [s1, m1] =MLVSharpnessMeasure(Zsub/sum(Zsub,'all'));
        data(kk).SharpnessScoreBinned(nn) = s1;   
    end
    disp(data(kk).SharpnessScoreBinned )
end



end
