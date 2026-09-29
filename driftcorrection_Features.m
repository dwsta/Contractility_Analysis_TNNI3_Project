function [xvec,yvec,tvec,repvec,Ushift,Vshift,angles,newIMAGE] = driftcorrection_Features(IMAGES,cfg,reference_frames,vecind)
    
    reference = IMAGES(:,:,reference_frames(1)); %can only have one frame
    
    reference = double(reference); %convert to double for following math steps
    
    % normalize reference and convert back to int (improves feature detection)
    reference = (reference-min(reference(:)))/(max(reference(:)-min(reference(:)))) * 65535; 
    reference = uint16(reference);
    
    %consider replacing with ORB features for open source
    ptsReference = detectSURFFeatures(reference);
    [featuresReference,validPtsReference] = extractFeatures(reference,ptsReference,"Method","SURF");
    
    newIMAGE = uint16(zeros(size(IMAGES)));
    xvec = 1;
    yvec = 1;
    tvec = vecind;
    repvec = 1;
    Ushift = zeros(1,1,size(IMAGES,3));
    Vshift = zeros(1,1,size(IMAGES,3));
    angles = zeros(1,1,size(IMAGES,3));
    
    
    for img = 1:size(IMAGES,3)
        %grab and normalize the image, then convert back to int
        movie = double(IMAGES(:,:,img));
        movie = (movie-min(movie(:)))/(max(movie(:)-min(movie(:)))) * 65535;
        
        movie = uint16(movie);
            
        %detect features
        ptsMovie = detectSURFFeatures(movie);
        [featuresMovie,validPtsMovie] = extractFeatures(movie,ptsMovie,"Method","SURF");
    
        %match reference features to frame features
        indexPairs = matchFeatures(featuresReference,featuresMovie);
        
        matchedReference = validPtsReference(indexPairs(:,1));
        matchedMovie = validPtsMovie(indexPairs(:,2));

        %calculate rigid transform
        [tform, inlierIdx] = estgeotform2d(matchedMovie,matchedReference,'rigid');
        inlierMovie = matchedMovie(inlierIdx,:);
        inlierReference = matchedReference(inlierIdx,:);

        %convert rotation
        invTform = invert(tform);
        Ainv = invTform.A;
        
        ss = Ainv(1,2);
        sc = Ainv(1,1);
        %scaleRecovered = hypot(ss,sc);
        %disp(['Recovered translation: ',num2str(tform.Translation)])
        %disp(['Recovered scale: ', num2str(scaleRecovered)])
        
        % Recover the rotation in which a positive value represents a rotation in
        % the clockwise direction.
        thetaRecovered = atan2d(-ss,sc);
        %disp(['Recovered theta: ', num2str(thetaRecovered)])
        
        %tform.RotationAngle=0;
        
        outputView = imref2d(size(reference));
        %recovered = imwarp(movie,tform,OutputView=outputView);
        recovered = imwarp(IMAGES(:,:,img),tform,OutputView=outputView);
        
        Ushift(img) = -tform.Translation(1);
        Vshift(img) = -tform.Translation(2);
        angles(img) = thetaRecovered;
        newIMAGE(:,:,img) = uint16(recovered);
    
    end

end
