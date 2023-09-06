function s3=ixon_getDayDir(src)    

    % Default root directory
    if nargin ~= 1;src=['X:\Data'];end

    if ~exist(src,'dir')
        s3=pwd;
        return;
    end
    t = now;
    s1=datestr(t,'yyyy');s2=datestr(t,'yyyy.mm');s3=datestr(t,'mm.dd');
    s1=[src filesep s1];s2=[s1 filesep s2];s3=[s2 filesep s3];

    if ~exist(s1,'dir'); mkdir(s1); end
    if ~exist(s2,'dir'); mkdir(s2); end
    if ~exist(s3,'dir'); mkdir(s3); end
end