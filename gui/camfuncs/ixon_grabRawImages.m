function imgs=ixon_grabRawImages
    % How many images to grab
    [ret,first,last] = GetNumberNewImages;    
    numpix=512^2;
    % Grab the data (number just sets buffer size)
    [ret,D] = GetAcquiredData(last*512^2);    
    imgs={};
    imgmats=zeros(512,512,last);
    for j = 1:last % break up into individual images
        ii=double(D((1+(j-1)*numpix):(j*numpix)));
        imgs{j} = reshape(ii,512,512);
        imgmats(:,:,j)=imgs{j};
    end 
    imgs=imgmats;  
end