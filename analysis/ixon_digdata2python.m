function ixon_digdata2python(digdata)

if nargin==0
   [filename,mydir] = uigetfile; 
end

digdata=load(fullfile(mydir,filename));

Z = [digdata.Z];


end

